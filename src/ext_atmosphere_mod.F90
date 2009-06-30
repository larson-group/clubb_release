!----------------------------------------------------------------------------
! $Id$
module ext_atmosphere_mod

  implicit none

  private ! Default Scope

  public :: &
    ext_atmosphere, &
    load_ext_std_atm, & 
    convert_snd2ext_atm, &
    finalize_ext_atm

  ! Size of Extended Atmosphere
  integer, public :: ext_atmos_dim

  logical, public :: l_use_default_std_atmosphere

  ! Parameters from U.S. Standard Atmosphere, 1976;
  ! Starting at 1Km altitude

  ! Altitude in meters
  double precision, public, target, allocatable, dimension(:) :: ext_alt

  ! Temperature in degrees Kelvin
  double precision, public, target, allocatable,dimension(:) :: ext_T_in_K

  ! Specific Humidity ( Water Vapor / Density )
  double precision, public, target, allocatable, dimension(:) ::  ext_sp_hmdty

  ! Pressure in millibars
  double precision, public, target, allocatable, dimension(:) :: ext_pinmb


  ! Ozone ( O_3 / Density )
  double precision, public, target, allocatable, dimension(:) :: ext_o3l

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine convert_snd2ext_atm( n_snd_var, psfc, zm_init, n_sclr_var, &
                                  sounding_profiles, sclr_sounding_profiles )
   !
   !  Description: This subroutine converts information retrieved from the
   !  sounding files of a case into a format usable for an extended atmosphere.
   !
   !  References:
   !    none
   !   
   !------------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var, read_x_profile

    use input_interpret, only: read_theta_profile, read_z_profile

    use constants, only: kappa, p0

    use input_names, only: z_name, temperature_name, ozone_name, rt_name
  
    implicit none

    ! Input Variables

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

    integer i

    character(len=20) :: alt_type, theta_type

    ext_atmos_dim = size( sounding_profiles(1)%values )

    ! initializing memory
    allocate( ext_alt(1:ext_atmos_dim) )
    allocate( ext_T_in_K(1:ext_atmos_dim) )
    allocate( ext_sp_hmdty(1:ext_atmos_dim) )
    allocate( ext_pinmb(1:ext_atmos_dim) )
    allocate( ext_o3l(1:ext_atmos_dim) )

    allocate( alt(1:ext_atmos_dim) )
    allocate( theta(1:ext_atmos_dim)  )
    allocate( p_in_Pa(1:ext_atmos_dim) )
    allocate( exner(1:ext_atmos_dim) )


    ! Either convert to pressure or from pressure
    
    call read_z_profile( n_snd_var, ext_atmos_dim, sounding_profiles, psfc, zm_init, &
                         alt , p_in_Pa, alt_type )

    ext_alt = alt

    if( alt_type == z_name ) then
      stop "Feature not implemented"
    end if

    ext_pinmb = p_in_Pa/ 100.

    ! Convert to temperature from thlm or theta

    call read_theta_profile( n_snd_var, ext_atmos_dim, sounding_profiles, theta_type, theta  )
    ext_T_in_K = theta

    if( theta_type /= temperature_name ) then
      exner(1) = ( psfc/p0 ) ** kappa
      do i = 2, ext_atmos_dim
        exner(i) = (p_in_Pa(i)/p0) ** kappa
      end do
      ext_T_in_K = ext_T_in_K * exner
    end if

    ! Convert rtm to specific humidity

    ext_sp_hmdty = read_x_profile( n_snd_var, ext_atmos_dim, rt_name, &
                      sounding_profiles )

    ext_sp_hmdty = ext_sp_hmdty/(ext_sp_hmdty +1)

    ! Read in ozone
    ext_o3l = read_x_profile( n_sclr_var, ext_atmos_dim, ozone_name, &
                              sclr_sounding_profiles )

    ! Free Memory
    deallocate( alt )
    deallocate( theta )
    deallocate( p_in_Pa )
    deallocate( exner )

  end subroutine convert_snd2ext_atm

  !-----------------------------------------------------------------------

  subroutine ext_atmosphere( alt, theta, rtm, p_in_Pa )
    !
    !       Description:
    !       Given a specific altitude this subroutine will return interpolated
    !       values for theta and rtm from U.S. Standard Atmosphere data.

    !       References:
    !       McClatchey, et al., (1972) _Environmental Research Papers_,
    !       No. 411, p.94
    !-----------------------------------------------------------------------
    use constants, only:  & 
        p0,  & ! Variable(s) 
        kappa

    use interpolation, only:  & 
        lin_int,  & ! Procedure(s) 
        binary_search

    implicit none

    intrinsic :: real

    ! Input Variable
    real,intent(in) :: alt ! Altitude         [m]

    ! Output Variables
    real,intent(out) ::  & 
    theta,  & ! Potential Temperature         [K]        
    rtm,    & ! Total Water Mixing Ratio      [kg/kg]
    p_in_Pa   ! Pressure                      [Pa]

    ! Local Variables
    real ::  & 
    exner,       & ! Exner function                [-]
    sp_humidity, & ! Specific humidity             [kg/kg]
    tabs0          ! Temperature                   [K]

    ! These variables are used to make the calls to lin_int cleaner
    real, dimension(ext_atmos_dim) :: & 
    T_in_K, & 
    pinmb, & 
    sp_hmdty, & 
    height

    integer :: varindex
    if( l_use_default_std_atmosphere )then

      T_in_K   = real( ext_T_in_K )
      pinmb    = real( ext_pinmb )
      sp_hmdty = real( ext_sp_hmdty )
      height   = real( ext_alt )

      varindex = binary_search( ext_atmos_dim, height, alt )

      if( varindex < 0 ) then
        stop "Cannot find altitude in Standard Atmosphere"
      endif

      ! Compute thlm from Standard Atmosphere

      tabs0 = lin_int( alt, height(varindex), height(varindex-1),  & 
                      T_in_K(varindex), T_in_K(varindex-1) )

      p_in_Pa = 100. *  & 
              lin_int( alt, height(varindex), height(varindex-1), & 
                      pinmb(varindex), pinmb(varindex-1) )

      exner = (p_in_Pa/p0)**kappa

      theta = tabs0/exner

      ! Compute rtm

      sp_humidity = lin_int( alt, height(varindex), height(varindex-1), & 
                    sp_hmdty(varindex), sp_hmdty(varindex-1))

      rtm = sp_humidity/( 1. - sp_humidity)

    end if

    return

  end subroutine ext_atmosphere

  !-------------------------------------------------------------------------------------------------
  subroutine load_ext_std_atm ( iunit )
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
 
    integer, parameter :: nCol = 5

    integer, parameter :: std_atmos_dim = 50

    integer, intent(in) :: iunit

    type(one_dim_read_var), dimension(nCol) :: retVars

    ext_atmos_dim = std_atmos_dim

    call read_one_dim_file( iunit, nCol, &
      "../input/std_atmosphere/atmosphere.in", retVars )

    ! initializing memory

    allocate( ext_alt(1:ext_atmos_dim) )
    allocate( ext_T_in_K(1:ext_atmos_dim) )
    allocate( ext_sp_hmdty(1:ext_atmos_dim) )
    allocate( ext_pinmb(1:ext_atmos_dim) )
    allocate( ext_o3l(1:ext_atmos_dim) )

    ext_alt = read_x_profile( nCol, ext_atmos_dim, z_name, retVars  )

    ext_T_in_K = read_x_profile( nCol, ext_atmos_dim, temperature_name, retVars )

    ext_sp_hmdty = read_x_profile( nCol, ext_atmos_dim, sp_humidity_name, retVars )

    ext_pinmb = read_x_profile( nCol, ext_atmos_dim, press_mb_name, retVars )

    ext_o3l = read_x_profile( nCol, ext_atmos_dim, ozone_name, retVars )

  end subroutine load_ext_std_atm
  !-------------------------------------------------------------------------------------------------
  subroutine finalize_ext_atm( )
    !
    !  Description:
    !    Frees memory used by the ext_atmosphere module.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------
    implicit none
    deallocate( ext_alt )
    deallocate( ext_T_in_K )
    deallocate( ext_sp_hmdty )
    deallocate( ext_pinmb )
    deallocate( ext_o3l )

  end subroutine finalize_ext_atm
end module ext_atmosphere_mod
