!-------------------------------------------------------------------------------
! $Id$
module extended_atmosphere_module

! Sometimes we wish to compute the radiative processes (eg ozone absorption) at high altitudes
! without undertaking the computational expense of computing CLUBB's equations there.
! To do so, we feed a grid that includes extra levels at high altitudes into the radiation code
! while computing CLUBB's moment equations on the standard grid.
!
! Altitude increases with increasing grid index.
!
! Grid Layout
!---------------------------------------------------------------------------------------
!  Computational  |     Buffer    |    Extended   | Complete                           |
!                 |               |               |                                    |
!                 |               | The           |                                    |
!                 |               | Extended      | The Complete Grid glues together   |
!                 |               | Grid is a     | parts of the Computational         |
!                 | The Buffer    | sounding      | Buffer, and Extended Grids         |
!                 | Grid          | which         | in the following way:              |
!                 | interpolates  | includes      |                                    |
!                 | between the   | grid          |   Extended Grid                    |
!                 | top of the    | heights       |                                    |
!                 | Computational | from the      |        +                           |
!   The           | Grid and      | ground to     |                                    |
!   Computational | the bottom    | above the     |   Buffer Grid                      |
!   Grid is where | of the        | Radiation     |                                    |
!   CLUBB's       | Radiation     | grid. It      |        +                           |
!   dynamical     | Grid          | includes      |                                    |
!   computations  |               | profiles of   |   Computational Grid               |
!   take place    |               | absorbers     |                                    |
!                 |               | (eg ozone)    |                                    |
!                 |               |               |                                    |
! ///////Ground/////////////////////////////////////////////////////////////////////////
!---------------------------------------------------------------------------------------
! Grid Positions and Variables Explained Graphically
!-------------------------------------------------------------------
!  Computational  |     Buffer           |   Extended              |
!                 | (not a separate      | -                       |
!                 |  array)              | extended_atmos_range_size |
!                 |                      | |                       |
!                 |                      | |                       |
!                 |                      | |                       |
!                 |                      | j                       |
!                 | -                    | -                       |
!                 | lin_int_buffer_size  |                         |
!                 | |                    |                         |
!                 | 1                    |                         |
! -               | -                    |                         |
! grid_size       |                      |                         |
! |               |                      |                         |
! |               |                      |                         |
! |               |                      |                         |
! |               |                      |                         |
! 1               |                      |                         |
! -/////Ground//////////////////////////////////////////////////////
!
! References:
!   McClatchey, et al., (1972) _Environmental Research Papers_, No. 411, p.94
!--------------------------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default Scope

  public :: &
    determine_extended_atmos_bounds,&
    finalize_extended_atm

  integer, public :: &
    extended_atmos_dim, & ! Size of Extended Atmosphere
    total_atmos_dim     ! Total Atmosphere Size (grid + buffer + extended atmosphere)
!$omp threadprivate(extended_atmos_dim, total_atmos_dim)

  real( kind = core_rknd ), public, allocatable, dimension(:) :: &
    complete_alt, &     ! Altitude of complete atmosphere in meters
    complete_momentum   ! Altitude of complete momentum grid in meters
!$omp threadprivate(complete_alt, complete_momentum)

  ! Extended Atmosphere variables
  real( kind = core_rknd ), public, allocatable, dimension(:) :: &
    extended_alt, &         ! Altitude, increases with array index    [m]
    extended_T_in_K, &      ! Temperature in degrees Kelvin
    extended_sp_hmdty, &    ! Specific Humidity ( Water Vapor / Density )
    extended_p_in_mb, &     ! Pressure in millibars
    extended_o3l            ! Ozone ( O_3 / Density )
!$omp threadprivate(extended_alt, extended_T_in_K, extended_sp_hmdty)
!$omp threadprivate(extended_p_in_mb, extended_o3l)

  contains

  !------------------------------------------------------------------------------------------------
  subroutine determine_extended_atmos_bounds( grid_size, zt_grid, &
                                            zm_grid, zm_grid_spacing, p_in_Pa_zm, & 
                                            radiation_top, &
                                            extended_atmos_bottom_level, &
                                            extended_atmos_top_level, &
                                            extended_atmos_range_size, &
                                            lin_int_buffer_size )
    ! Description: 
    !   This subroutine determines the bottom and top levels of the
    !   extended atmosphere to use for a given radiation_top and grid. It also
    !   computes the amount of a linear interpolation buffer between the grid and
    !   the extended atmosphere.
    !
    !   The linear interpolation buffer is used to counter cooling spikes caused
    !   by a potentially large difference in temperature between the extended
    !   atmosphere profile and the CLUBB profile.
    !
    ! References:
    !   None
    !
    !----------------------------------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr, & ! Variable(s)
      pascal_per_mb

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! External
    intrinsic :: real, max, int

    ! Input Variable(s)
    integer, intent(in) :: grid_size ! Size of the model grid  [-]

    real( kind = core_rknd ), dimension(grid_size-1), intent(in) :: &
      zt_grid,         & ! Thermodynamic grid        [m]
      zm_grid_spacing    ! Change per zm grid levels [m]

    real( kind = core_rknd ), dimension(grid_size), intent(in) :: &
      zm_grid,         & ! Momentum grid             [m]
      p_in_Pa_zm         ! Pressure                  [Pa]

    real( kind = core_rknd ), intent(in) :: radiation_top ! Maximum height to extended to [m]

    ! Output Variable(s)
    integer, intent(out) :: &
      extended_atmos_bottom_level, & ! Index of lowest point to use for atmosphere extension [-]
      extended_atmos_top_level, &    ! Index of highest point to use for atmosphere extension [-]
      extended_atmos_range_size, &   ! Size of the range between the
      !                              two points bounds for atmosphere extension[-]
      lin_int_buffer_size          ! Size of linear interpolation buffer [-]

    ! Local Variable(s)
    real(kind=core_rknd) :: &
      dz10, dz_model, dz_extension, dz, &
      zm_grid_top, extended_bottom, buffer_size

    integer :: &
      i, & ! Loop index
      j, & ! Array index in extended_alt which is one above the computational domain  []
      k    ! Array index of altitude one less than the top of the radiative domain  []

    ! ---- Begin Code ----

    if ( radiation_top < zm_grid(grid_size) ) then
      write(fstderr,*) "In subroutine determine_extended_atmos_bounds"
      error stop "top of the radiation grid is below the top of the computational grid"
    end if

    ! Determine the bounds to use for the extended atmosphere

    j=1

    ! This code ensures that altitudes monotonically increases with increasing grid level
    do while ( extended_alt(j) < zm_grid(grid_size) .and. j < extended_atmos_dim )
      j = j + 1
    end do

    ! This code ensures that pressure monotonically decreases with increasing grid level
    do while (p_in_Pa_zm(grid_size) < &
            extended_p_in_mb(j) * pascal_per_mb )
      j = j + 1
    end do

    if ( extended_alt(j) < zm_grid(grid_size) ) then
      write(fstderr,*) "In subroutine determine_extended_atmos_bounds"
      error stop "Extended atmosphere is below the top of the computational grid"
    end if

    if ( extended_alt(extended_atmos_dim) < radiation_top ) then
      write(fstderr,*) "In subroutine determine_extended_atmos_bounds"
      write(fstderr,*) "Atmosphere cannot be extended because extension data does ", &
                         "not reach radiation_top"
      error stop
    end if

    if ( p_in_Pa_zm(grid_size) < &
         extended_p_in_mb(j) * pascal_per_mb ) then
      write(fstderr,*) "In subroutine determine_extended_atmos_bounds"
      error stop &
            "pressure at top of computational grid less than pressure at base of radiative grid"
    end if

    k=1

    if ( j <= extended_atmos_dim ) then

      do while( extended_alt(k) < radiation_top .and. k < extended_atmos_dim )
        k= k+1
      end do

      ! It is possible we could be above the specified radiation top, check
      ! and roll back if neccessary
      if( extended_alt(k) > radiation_top ) then
        k = k-1
      end if

    else
      k = j
    end if

    extended_atmos_bottom_level = j
    extended_atmos_top_level = k
    extended_atmos_range_size = k - j + 1

    if ( extended_atmos_range_size < 1 ) then
      write(fstderr,*) "In subroutine determine_extended_atmos_bounds"
      error stop "radiation top below computational grid"
    end if

    ! Get the altitudes for a couple of key points so we can calculate a buffer
    ! size
    zm_grid_top =  zm_grid(grid_size) !Altitude at top of normal grid
    extended_bottom = extended_alt(extended_atmos_bottom_level) !Altitude at bottom of
                                                          !extended atmos
    
    ! Determine the spacing of the lin_int_buffer, it should have no more than
    ! 10 levels.
    dz10 = (extended_bottom - zm_grid_top) / 10._core_rknd
    dz_model = zm_grid_spacing(grid_size-1)
    dz = max( dz10, dz_model )
    ! Calculate the size of the lin_int_buffer
    buffer_size = (extended_bottom - zm_grid_top) / dz
    lin_int_buffer_size = int( buffer_size )

    ! Calculate the dimension of the entire atmosphere
    total_atmos_dim = grid_size + lin_int_buffer_size + extended_atmos_range_size

    ! Build the complete momentum grid
    ! The extended momentum grid contains one level above the
    ! extended thermodynamic grid.
    allocate( complete_momentum(total_atmos_dim + 1) )

    forall ( j=1:grid_size ) ! Loop from the lowest zm level to zm_grid_top
      complete_momentum(j) = zm_grid(j)
    end forall

    ! Interpolate between the top of the computational grid and the bottom
    ! of the extended altitude
    i = 1 ! Tracks the number of interpolation levels used
    do j=grid_size+1, grid_size+lin_int_buffer_size
      complete_momentum(j) = real( zm_grid_top, kind = core_rknd ) + &
        real( dz, kind = core_rknd ) * real( i, kind = core_rknd )
      i = i + 1
    end do

    i = 0 ! Tracks the number of extension levels used
    do j = grid_size + lin_int_buffer_size + 1, total_atmos_dim
      ! Take values from the extended atmosphere    
      complete_momentum(j) = real( extended_alt(extended_atmos_bottom_level + i), kind = core_rknd )
      ! Keep track of where we are in the extended atmosphere
      i = i + 1
      if ( (i-1)+extended_atmos_bottom_level == extended_atmos_dim ) then
        exit
      else
        cycle ! Loop to next point
      end if
    end do

    ! Use a linear extension for the topmost points (generally this is only 1 or 2 points)
    do j = grid_size + lin_int_buffer_size + i, total_atmos_dim + 1
      dz_extension = complete_momentum(j-1) - complete_momentum(j-2)
      complete_momentum(j) = complete_momentum(j-1) + real( dz_extension, kind = core_rknd )
    end do
    
    allocate( complete_alt(total_atmos_dim) )

    ! Build the total atmosphere grid for the zt levels
    forall ( j=1:grid_size-1 )
      complete_alt(j) = zt_grid(j)
    end forall

    forall ( j=grid_size:total_atmos_dim )
      ! Use a linear interpolation above the zt_grid so the points are between the momentum levels
      complete_alt(j) = (complete_momentum(j+1) + complete_momentum(j)) / 2._core_rknd
    end forall

    return
  end subroutine determine_extended_atmos_bounds
  
  !-------------------------------------------------------------------------------------------------
  subroutine finalize_extended_atm( )
    !
    !  Description:
    !    Frees memory used by the extended_atmosphere module.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! External 
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( extended_alt ) ) then
      deallocate( extended_alt )
    end if

    if ( allocated( extended_T_in_K ) ) then
      deallocate( extended_T_in_K )
    end if

    if ( allocated( extended_sp_hmdty ) ) then
      deallocate( extended_sp_hmdty )
    end if

    if ( allocated( extended_p_in_mb ) ) then
      deallocate( extended_p_in_mb )
    end if

    if ( allocated( extended_o3l ) ) then
      deallocate( extended_o3l )
    end if

    if ( allocated( complete_alt ) ) then
      deallocate( complete_alt )
    end if

    if ( allocated( complete_momentum ) ) then
      deallocate( complete_momentum )
    end if

    return
  end subroutine finalize_extended_atm

end module extended_atmosphere_module
