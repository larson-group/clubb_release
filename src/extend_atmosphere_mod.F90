!----------------------------------------------------------------------------
! $Id$
module extend_atmosphere_mod

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

  ! Flag to signal the use of the U.S. Standard Atmosphere Profile, 1976
  logical, public :: l_use_default_std_atmosphere

  ! Extended Atmosphere variables

  ! Altitude in meters
  double precision, public, target, allocatable, dimension(:) :: extend_alt

  ! Temperature in degrees Kelvin
  double precision, public, target, allocatable,dimension(:) :: extend_T_in_K

  ! Specific Humidity ( Water Vapor / Density )
  double precision, public, target, allocatable, dimension(:) ::  extend_sp_hmdty

  ! Pressure in millibars
  double precision, public, target, allocatable, dimension(:) :: extend_pinmb


  ! Ozone ( O_3 / Density )
  double precision, public, target, allocatable, dimension(:) :: extend_o3l

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine convert_snd2extend_atm( n_snd_var, psfc, zm_init, n_sclr_var, &
                                  sounding_profiles, sclr_sounding_profiles )
    !
    !  Description: This subroutine converts information retrieved from the
    !  sounding files of a case into a format usable for an extended atmosphere.
    !  The extended atmosphere profile is stored in module variables.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var, read_x_profile ! Procedure(s)

    use input_interpret, only: read_theta_profile, read_z_profile ! Procedure(s)

    use constants, only: kappa, p0, fstderr ! Constant(s)

    use input_names, only: z_name, temperature_name, ozone_name, rt_name ! Variable(s)

    implicit none

    ! External
    intrinsic :: size

    ! Input Variable(s)

    integer, intent(in) :: n_snd_var   ! Number of variables from sounding [-]

    integer, intent(in) :: n_sclr_var  ! Number of variables from sclr_sounding [-]

    real, intent(in) :: psfc  ! Pressure at the surface [Pa]

    real, intent(in) :: zm_init ! Height at zm(1) [m]

    type(one_dim_read_var), dimension(n_snd_var), intent(in) :: &
      sounding_profiles ! Sounding profile

    type(one_dim_read_var), dimension(n_sclr_var), intent(in) :: &
      sclr_sounding_profiles ! Sclr Sounding profile

    ! Local Variables

    real,  dimension(:), allocatable :: alt, theta, p_in_Pa, exner

    integer :: i, sclr_sounding_dim

    character(len=20) :: alt_type, theta_type

    ! -- Begin Code --

    ! Determine the size of the extended atmosphere buffer
    extend_atmos_dim = size( sounding_profiles(1)%values )
    sclr_sounding_dim = size( sclr_sounding_profiles(1)%values )

    ! Check for an error condition
    if ( extend_atmos_dim /= sclr_sounding_dim ) then
      write(fstderr,*) "Fatal error in convert_snd2extend_atm."
      write(fstderr,*) "This code assumes the scalar sounding will be the same size"
      write(fstderr,*) "as the main sounding data."
      write(fstderr,*) "standard sounding dimension: ", extend_atmos_dim
      write(fstderr,*) "scalar sounding dimension: ", sclr_sounding_dim
      stop 
    end if

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

    ! Either convert to pressure or from pressure

    call read_z_profile( n_snd_var, extend_atmos_dim, sounding_profiles, psfc, zm_init, &
                         alt , p_in_Pa, alt_type )

    extend_alt = alt

    if ( alt_type == z_name ) then
      write(fstderr,*) "Fatal error in convert_snd2extend_atm."
      stop "Feature not implemented"
    end if

    extend_pinmb = p_in_Pa/ 100.

    ! Convert to temperature from thlm or theta

    call read_theta_profile( n_snd_var, extend_atmos_dim, sounding_profiles, theta_type, theta  )
    extend_T_in_K = theta

    if( theta_type /= temperature_name ) then
      exner(1) = ( psfc/p0 ) ** kappa
      do i = 2, extend_atmos_dim
        exner(i) = (p_in_Pa(i)/p0) ** kappa
      end do
      extend_T_in_K = extend_T_in_K * exner
    end if

    ! Convert rtm to specific humidity

    extend_sp_hmdty = read_x_profile( n_snd_var, extend_atmos_dim, rt_name, &
                      sounding_profiles )

    extend_sp_hmdty = extend_sp_hmdty / ( extend_sp_hmdty +1 )

    ! Read in ozone
    extend_o3l = read_x_profile( n_sclr_var, extend_atmos_dim, ozone_name, &
                                 sclr_sounding_profiles )

    ! Free Memory
    deallocate( alt )
    deallocate( theta )
    deallocate( p_in_Pa )
    deallocate( exner )

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
    !  Description: This subroutine determines the bottom and top levels of the
    !  extended atmosphere to use for a given radiation_top and grid. It also
    !  computes the amount of a linear interpolation buffer between the grid and
    !  the extended atmosphere.
    !
    !  The linear interpolation buffer is used to counter cooling spikes caused
    !  by a potentially large difference in temperature between the extended
    !  atmosphere profile and the CLUBB profile.
    !
    !
    !  References:
    !    None
    !
    !----------------------------------------------------------------------------------------------


    use constants, only: &
    fstderr ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: grid_size ! Size of the model grid  [-]

    real, dimension(grid_size), intent(in) :: zt_grid ! Model grid [m]

    real, dimension(grid_size), intent(in) :: zm_grid ! Model grid [m]

    real, dimension(grid_size), intent(in) :: &
      zm_grid_spacing ! Inverse spacing between grid levels [m]

    real, intent(in) :: radiation_top ! Maximum height to extend to [m]

    ! Output Variable(s)
    integer, intent(out) :: &
      extend_atmos_bottom_level, & ! Index of lowest point to
      !                              use for atmosphere extension [-]
      extend_atmos_top_level, &    ! Index of highest point to 
      !                              use for atmosphere extension [-]
      extend_atmos_range_size, &   ! Size of the range between the 
      !                              two points bounds for atmosphere extension[-]
      lin_int_buffer_size          ! Size of linear interpolation buffer [-]

    ! Local Variable(s)
    integer :: k, j, i
    double precision :: dz10, dz_model, dz_extension, dz
    double precision :: zm_grid_top, extend_bottom, buffer_size

    ! -- Begin Code --

    ! Determine the bounds to use for the extended atmosphere

    j=1
    do while( extend_alt(j) < zm_grid(grid_size) .and. j < extend_atmos_dim )
      j= j+1
    end do

    if(extend_alt(j) < zm_grid(grid_size)) then
      stop "Extended atmosphere is below the top of the computational grid"
    end if

    if( extend_alt(extend_atmos_dim) < radiation_top   ) then
      write(fstderr,*) "Atmosphere cannot be extended because extension data does ", &
                         "not reach radiation_top"
      stop
    end if

    k=1

    if( j <= extend_atmos_dim ) then

      do while( extend_alt(k) < radiation_top .and. k < extend_atmos_dim )
        k= k+1
      end do

      ! It is possible we could be above the specified radiation top, check
      ! and roll back if neccessary
      if( extend_alt(k) > radiation_top) then
        k= k-1
      end if

    else
      k = j
    end if


    extend_atmos_bottom_level = j
    extend_atmos_top_level = k
    extend_atmos_range_size = k - j + 1

    if( extend_atmos_range_size < 1 ) then
      stop "radiation top below computational grid"
    end if

    ! Get the altitudes for a couple of key points so we can calculate a buffer
    ! size
    zm_grid_top =  zm_grid(grid_size) !Altitude at top of normal grid
    extend_bottom = extend_alt(extend_atmos_bottom_level) !Altitude at bottom of
                                                          !extended atmos
    
    ! Determine the spacing of the lin_int_buffer, it should have no more than
    ! 10 levels.
    dz10 = (extend_bottom - zm_grid_top) / 10
    dz_model = zm_grid_spacing(grid_size)
    dz = max(dz10, dz_model)
    ! Calculate the size of the lin_int_buffer
    buffer_size = (extend_bottom - zm_grid_top) / dz
    lin_int_buffer_size = int(buffer_size)

    ! Calculate the dimension of the entire atmosphere
    total_atmos_dim = grid_size + lin_int_buffer_size + extend_atmos_range_size

    ! Build the complete momentum grid
    ! The extended momentum grid contains one level above the
    ! extended thermodynamic grid.
    allocate( complete_momentum(1:total_atmos_dim + 1) )

    i = 0 !The number of extension levels used
    do j=1, total_atmos_dim + 1
      if (j <= grid_size) then
        complete_momentum(j) = zm_grid(j)
      elseif(j > grid_size .and. j <= (grid_size + lin_int_buffer_size)) then 
        !Interpolate between the top of the computational grid and the bottom
        !of the extended altitude
        complete_momentum(j) = real(zm_grid_top + ((extend_bottom - zm_grid_top) / & 
                          (lin_int_buffer_size + 1)) * &
                          (j - grid_size))
      elseif(i < extend_atmos_range_size) then
        !Take values from the extended atmosphere    
        complete_momentum(j) = real(extend_alt(extend_atmos_bottom_level + i))
        ! Keep track of where we are in the extended atmosphere
        i = i + 1
      else
        !We should only get here for the absolute topmost point, which is a
        !linear extension
        dz_extension = complete_momentum(j-1) - complete_momentum(j-2)
        complete_momentum(j) = real(complete_momentum(j-1) + dz_extension)
      endif 
    end do
    
    allocate( complete_alt(1:total_atmos_dim) )

    ! Build the total atmosphere grid
    do j=1, total_atmos_dim
      if (j <= grid_size) then
        complete_alt(j) = zt_grid(j)
      else    
        !Linear extension above zt_grid so points are between momentum levels
        complete_alt(j) = (complete_momentum(j - 1) + complete_momentum(j)) / 2
      endif 
    end do

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
      read_x_profile, &
      one_dim_read_var, &
      read_one_dim_file

    use input_names, only: &
      z_name, &
      ozone_name, &
      temperature_name, &
      press_mb_name, &
      sp_humidity_name

    implicit none

    ! Constant Parameters
    integer, parameter :: nCol = 5

    integer, parameter :: std_atmos_dim = 50

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit [-]


    ! Local Variable(s)
    type(one_dim_read_var), dimension(nCol) :: retVars

    character(len=*), parameter :: atm_input_file = "../input/std_atmosphere/atmosphere.in"

    ! -- Begin Code --

    extend_atmos_dim = std_atmos_dim

    call read_one_dim_file( iunit, nCol, atm_input_file, retVars )

    ! initializing memory

    allocate( extend_alt(1:extend_atmos_dim) )
    allocate( extend_T_in_K(1:extend_atmos_dim) )
    allocate( extend_sp_hmdty(1:extend_atmos_dim) )
    allocate( extend_pinmb(1:extend_atmos_dim) )
    allocate( extend_o3l(1:extend_atmos_dim) )

    extend_alt = read_x_profile( nCol, extend_atmos_dim, z_name, retVars, &
                                 atm_input_file )

    extend_T_in_K = read_x_profile( nCol, extend_atmos_dim, temperature_name, retVars, &
                                    atm_input_file )

    extend_sp_hmdty = read_x_profile( nCol, extend_atmos_dim, sp_humidity_name, retVars, &
                                      atm_input_file )

    extend_pinmb = read_x_profile( nCol, extend_atmos_dim, press_mb_name, retVars, &
                                   atm_input_file )

    extend_o3l = read_x_profile( nCol, extend_atmos_dim, ozone_name, retVars, &
                                 atm_input_file )                           

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

    ! -- Begin Code --

    deallocate( extend_alt )
    deallocate( extend_T_in_K )
    deallocate( extend_sp_hmdty )
    deallocate( extend_pinmb )
    deallocate( extend_o3l )
    deallocate( complete_alt )
    deallocate( complete_momentum )

  end subroutine finalize_extend_atm
end module extend_atmosphere_mod
