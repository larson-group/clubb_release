!----------------------------------------------------------------------------
! $Id$
module extend_atmosphere_module

  use clubb_precision, only: &
    dp ! double precision

  implicit none

  private ! Default Scope

  public :: &
    load_extend_std_atm, & 
    convert_snd2extend_atm, &
    determine_extend_atmos_bounds,&
    finalize_extend_atm

  ! Size of Extended Atmosphere
  integer, public :: extend_atmos_dim

  ! Total Atmosphere Size (grid + buffer + extended atmosphere)
  integer, public :: total_atmos_dim

  ! Altitude of complete atmosphere in meters
  real, public, target, allocatable, dimension(:) :: complete_alt

  ! Altitude of complete momentum grid in meters
  real, public, target, allocatable, dimension(:) :: complete_momentum

  ! Extended Atmosphere variables

  ! Altitude in meters
  real( kind = dp ), public, target, allocatable, dimension(:) :: extend_alt

  ! Temperature in degrees Kelvin
  real( kind = dp ), public, target, allocatable,dimension(:) :: extend_T_in_K

  ! Specific Humidity ( Water Vapor / Density )
  real( kind = dp ), public, target, allocatable, dimension(:) ::  extend_sp_hmdty

  ! Pressure in millibars
  real( kind = dp ), public, target, allocatable, dimension(:) :: extend_pinmb


  ! Ozone ( O_3 / Density )
  real( kind = dp ), public, target, allocatable, dimension(:) :: extend_o3l

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine convert_snd2extend_atm( iunit, runtype, n_snd_var, p_sfc, zm_init, &
                                     sounding_profiles )
    !
    !  Description: This subroutine converts information retrieved from the
    !  sounding files of a case into a format usable for an extended atmosphere.
    !  The extended atmosphere profile is stored in module variables.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: &
      read_x_profile,  & ! Procedure(s)
      read_one_dim_file

    use input_reader, only: &
      one_dim_read_var ! Type

    use input_interpret, only: read_theta_profile, read_z_profile ! Procedure(s)

    use constants_clubb, only: kappa, p0, fstderr ! Constant(s)

    use input_names, only: z_name, temperature_name, ozone_name, rt_name ! Variable(s)

    use clubb_precision, only: &
      dp ! double precision

    implicit none

    ! External
    intrinsic :: size, trim

    ! Constant Parameters
    integer, parameter :: &
      n_rad_scalars = 1

    ! Input Variable(s)

    integer, intent(in) :: iunit  ! Fortran file unit

    character(len=*), intent(in) ::  runtype

    integer, intent(in) :: n_snd_var   ! Number of variables from sounding [-]

    real, intent(in) :: p_sfc  ! Pressure at the surface [Pa]

    real, intent(in) :: zm_init ! Height at zm(1) [m]

    type(one_dim_read_var), dimension(n_snd_var), intent(in) :: &
      sounding_profiles ! Sounding profile

    ! Local Variables

    real,  dimension(:), allocatable :: alt, theta, p_in_Pa, exner

    type(one_dim_read_var), dimension(n_rad_scalars) :: &
      rad_scalars_sounding_profile ! Ozone Sounding profile

    integer :: i, ivar

    character(len=20) :: alt_type, theta_type

    ! -- Begin Code --

    ! Determine the size of the extended atmosphere buffer
    extend_atmos_dim = size( sounding_profiles(1)%values )
   
    ! Allocate variables
    allocate( extend_alt(extend_atmos_dim) )
    allocate( extend_T_in_K(extend_atmos_dim) )
    allocate( extend_sp_hmdty(extend_atmos_dim) )
    allocate( extend_pinmb(extend_atmos_dim) )
    allocate( extend_o3l(extend_atmos_dim) )

    allocate( alt(extend_atmos_dim) )
    allocate( theta(extend_atmos_dim)  )
    allocate( p_in_Pa(extend_atmos_dim) )
    allocate( exner(extend_atmos_dim) )

    do ivar = 1, n_rad_scalars
      allocate( rad_scalars_sounding_profile(ivar)%values(extend_atmos_dim) )
    end do

    ! Either convert to pressure or from pressure

    call read_z_profile( n_snd_var, extend_atmos_dim, sounding_profiles, p_sfc, zm_init, &
                         alt , p_in_Pa, alt_type )

    extend_alt = dble( alt )

    if ( alt_type == z_name ) then
      write(fstderr,*) "Fatal error in convert_snd2extend_atm."
      stop "Feature not implemented"
    end if

    extend_pinmb = dble( p_in_Pa/ 100. )

    ! Convert to temperature from thlm or theta

    call read_theta_profile( n_snd_var, extend_atmos_dim, sounding_profiles, theta_type, theta  )
    extend_T_in_K = dble( theta )

    if( theta_type /= temperature_name ) then
      exner(1) = ( p_sfc/p0 ) ** kappa
      do i = 2, extend_atmos_dim
        exner(i) = (p_in_Pa(i)/p0) ** kappa
      end do
      extend_T_in_K = extend_T_in_K * dble( exner )
    end if

    ! Convert rtm to specific humidity

    extend_sp_hmdty = dble( read_x_profile( n_snd_var, extend_atmos_dim, rt_name, &
                                            sounding_profiles ) )

    extend_sp_hmdty = extend_sp_hmdty / ( extend_sp_hmdty +1._dp )

    ! Read in radiation scalars sounding (currently it only holds O3)
    call read_one_dim_file( iunit, n_rad_scalars, & ! In
      '../input/case_setups/'//trim( runtype )//'_ozone_sounding.in', & ! In
      rad_scalars_sounding_profile ) ! Out

    ! Set the array holding the values of o3l (ozone)
    extend_o3l = dble( read_x_profile( 1, extend_atmos_dim, ozone_name, &
                                       rad_scalars_sounding_profile ) )

    ! We would add the setting of new radiation scalar arrays like O3 right
    ! after this. New variable names would have to be added to input_names.

    ! Free Memory
    deallocate( alt )
    deallocate( theta )
    deallocate( p_in_Pa )
    deallocate( exner )

    do ivar = 1, n_rad_scalars
      deallocate( rad_scalars_sounding_profile(ivar)%values )
    end do

    return
  end subroutine convert_snd2extend_atm

  !------------------------------------------------------------------------------------------------
  subroutine determine_extend_atmos_bounds( grid_size, zt_grid, &
                                              zm_grid, zm_grid_spacing, & 
                                              radiation_top, &
                                              extend_atmos_bottom_level, &
                                              extend_atmos_top_level, &
                                              extend_atmos_range_size, &
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
    !----------------------------------------------------------------------------------------------


    use constants_clubb, only: &
      fstderr ! Variable(s)

    use clubb_precision, only: &
      dp ! double precision

    implicit none

    ! External
    intrinsic :: real, max, int

    ! Input Variable(s)
    integer, intent(in) :: grid_size ! Size of the model grid  [-]

    real, dimension(grid_size), intent(in) :: &
      zt_grid,       & ! Thermodynamic grid [m]
      zm_grid,       & ! Momentum grid [m]
      zm_grid_spacing  ! Inverse spacing between zm grid levels [m]

    real, intent(in) :: radiation_top ! Maximum height to extend to [m]

    ! Output Variable(s)
    integer, intent(out) :: &
      extend_atmos_bottom_level, & ! Index of lowest point to use for atmosphere extension [-]
      extend_atmos_top_level, &    ! Index of highest point to use for atmosphere extension [-]
      extend_atmos_range_size, &   ! Size of the range between the 
      !                              two points bounds for atmosphere extension[-]
      lin_int_buffer_size          ! Size of linear interpolation buffer [-]

    ! Local Variable(s)
    real(kind=dp) :: &
      dz10, dz_model, dz_extension, dz, &
      zm_grid_top, extend_bottom, buffer_size

    integer :: i, j, k ! Loop indices

    ! ---- Begin Code ----

    ! Determine the bounds to use for the extended atmosphere

    j=1
    do while ( extend_alt(j) < real( zm_grid(grid_size), kind=dp ) .and. j < extend_atmos_dim )
      j= j+1
    end do

    if ( extend_alt(j) < real( zm_grid(grid_size), kind=dp ) ) then
      write(fstderr,*) "In subroutine determine_extend_atmos_bounds"
      stop "Extended atmosphere is below the top of the computational grid"
    end if

    if ( extend_alt(extend_atmos_dim) < real( radiation_top, kind=dp ) ) then
      write(fstderr,*) "In subroutine determine_extend_atmos_bounds"
      write(fstderr,*) "Atmosphere cannot be extended because extension data does ", &
                         "not reach radiation_top"
      stop
    end if

    k=1

    if ( j <= extend_atmos_dim ) then

      do while( extend_alt(k) < real( radiation_top, kind=dp ) .and. k < extend_atmos_dim )
        k= k+1
      end do

      ! It is possible we could be above the specified radiation top, check
      ! and roll back if neccessary
      if( extend_alt(k) > real( radiation_top, kind=dp ) ) then
        k= k-1
      end if

    else
      k = j
    end if


    extend_atmos_bottom_level = j
    extend_atmos_top_level = k
    extend_atmos_range_size = k - j + 1

    if ( extend_atmos_range_size < 1 ) then
      write(fstderr,*) "In subroutine determine_extend_atmos_bounds"
      stop "radiation top below computational grid"
    end if

    ! Get the altitudes for a couple of key points so we can calculate a buffer
    ! size
    zm_grid_top =  real( zm_grid(grid_size), kind=dp ) !Altitude at top of normal grid
    extend_bottom = extend_alt(extend_atmos_bottom_level) !Altitude at bottom of
                                                          !extended atmos
    
    ! Determine the spacing of the lin_int_buffer, it should have no more than
    ! 10 levels.
    dz10 = (extend_bottom - zm_grid_top) / 10._dp
    dz_model = real( zm_grid_spacing(grid_size), kind=dp )
    dz = max( dz10, dz_model )
    ! Calculate the size of the lin_int_buffer
    buffer_size = (extend_bottom - zm_grid_top) / dz
    lin_int_buffer_size = int( buffer_size )

    ! Calculate the dimension of the entire atmosphere
    total_atmos_dim = grid_size + lin_int_buffer_size + extend_atmos_range_size

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
      complete_momentum(j) = real( zm_grid_top ) + real( dz ) * real( i )
      i = i + 1
    end do

    i = 0 ! Tracks the number of extension levels used
    do j=grid_size+lin_int_buffer_size, total_atmos_dim
      ! Take values from the extended atmosphere    
      complete_momentum(j) = real( extend_alt(extend_atmos_bottom_level + i) )
      ! Keep track of where we are in the extended atmosphere
      i = i + 1
      if ( i+extend_atmos_bottom_level == extend_atmos_dim ) then
        exit
      else
        cycle ! Loop to next point
      end if
    end do

    ! Use a linear extension for the topmost points (generally this is only 1 or 2 points)
    do j = i, total_atmos_dim + 1
      dz_extension = real( complete_momentum(j-1) - complete_momentum(j-2), kind=dp )
      complete_momentum(j) = complete_momentum(j-1) + real( dz_extension )
    end do
    
    allocate( complete_alt(total_atmos_dim) )

    ! Build the total atmosphere grid for the zt levels
    forall ( j=1:grid_size )
      complete_alt(j) = zt_grid(j)
    end forall

    forall ( j=grid_size+1:total_atmos_dim )
      ! Use a linear extension above the zt_grid so the points are between the momentum levels
      complete_alt(j) = (complete_momentum(j-1) + complete_momentum(j)) / 2.
    end forall

    return
  end subroutine determine_extend_atmos_bounds

  !-------------------------------------------------------------------------------------------------
  subroutine load_extend_std_atm ( iunit )
    !
    !  Description:
    !    Loads in the U.S. Standard atmosphere data from a file.
    !  References:
    !    McClatchey, et al., (1972) _Environmental Research Papers_,
    !    No. 411, p.94
    !
    !-----------------------------------------------------------------------------------------------

    use input_reader, only: &
      read_x_profile, & ! Procedure(s)
      read_one_dim_file, &
      deallocate_one_dim_vars

    use input_reader, only: &
      one_dim_read_var ! Derived type

    use input_names, only: &
      z_name, &
      ozone_name, &
      temperature_name, &
      press_mb_name, &
      sp_humidity_name

    use clubb_precision, only: &
      dp ! double precision

    implicit none

    ! Constant Parameters
    integer, parameter :: &
      nCol = 5,  &       ! Number of columns in the text file
      std_atmos_dim = 50 ! Dimension of the standard atmosphere table

    character(len=*), parameter :: &
      atm_input_file = "../input/std_atmosphere/atmosphere.in"

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit [-]


    ! Local Variable(s)
    type(one_dim_read_var), dimension(nCol) :: retVars

    ! -- Begin Code --

    extend_atmos_dim = std_atmos_dim

    call read_one_dim_file( iunit, nCol, atm_input_file, retVars )

    ! Allocate and initialize variables for standard atmosphere

    allocate( extend_alt(extend_atmos_dim) )
    allocate( extend_T_in_K(extend_atmos_dim) )
    allocate( extend_sp_hmdty(extend_atmos_dim) )
    allocate( extend_pinmb(extend_atmos_dim) )
    allocate( extend_o3l(extend_atmos_dim) )

    extend_alt = real( read_x_profile( nCol, extend_atmos_dim, z_name, retVars, &
                                 atm_input_file ), kind=dp )

    extend_T_in_K = real( read_x_profile( nCol, extend_atmos_dim, temperature_name, retVars, &
                                    atm_input_file ), kind=dp )

    extend_sp_hmdty = real( read_x_profile( nCol, extend_atmos_dim, sp_humidity_name, retVars, &
                                      atm_input_file ), kind=dp )

    extend_pinmb = real( read_x_profile( nCol, extend_atmos_dim, press_mb_name, retVars, &
                                   atm_input_file ), kind=dp )

    extend_o3l = real( read_x_profile( nCol, extend_atmos_dim, ozone_name, retVars, &
                                 atm_input_file ), kind=dp )

    ! Deallocate memory
    call deallocate_one_dim_vars( nCol, retVars )

    return
  end subroutine load_extend_std_atm
  !-------------------------------------------------------------------------------------------------
  subroutine finalize_extend_atm( )
    !
    !  Description:
    !    Frees memory used by the extend_atmosphere module.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! External 
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( extend_alt ) ) then
      deallocate( extend_alt )
    end if

    if ( allocated( extend_T_in_K ) ) then
      deallocate( extend_T_in_K )
    end if

    if ( allocated( extend_sp_hmdty ) ) then
      deallocate( extend_sp_hmdty )
    end if

    if ( allocated( extend_pinmb ) ) then
      deallocate( extend_pinmb )
    end if

    if ( allocated( extend_o3l ) ) then
      deallocate( extend_o3l )
    end if

    if ( allocated( complete_alt ) ) then
      deallocate( complete_alt )
    end if

    if ( allocated( complete_momentum ) ) then
      deallocate( complete_momentum )
    end if

    return
  end subroutine finalize_extend_atm

end module extend_atmosphere_module
