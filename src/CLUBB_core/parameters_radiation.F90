!-------------------------------------------------------------------------------
! $Id$
module parameters_radiation

! Description:
!   Parameters for radiation schemes

! References:
!   None
!-------------------------------------------------------------------------------
  implicit none

  character(len=10), public :: & 
    rad_scheme  ! Either BUGSrad, or simplified

  ! For BUGSrad
  double precision, dimension(1), public :: &
    sol_const ! Solar constant

  integer, public :: &
    ext_atmos_buffer

  ! Albedo values (alvdr is used in the simplifed schemes as well)
  double precision, public :: &
    alvdr, &   !Visible direct surface albedo   [-]
    alndr, &   !Near-IR direct surface albedo   [-]
    alvdf, &   !Visible diffuse surface albedo  [-]
    alndf      !Near-IR diffuse surface albedo  [-]


  ! Long-wave constants (simplified radiation)
  real, public :: &
    kappa, & ! A constant (Duynkerke eqn. 5)                   [m^2/kg]
    F0,    & ! Coefficient for cloud top heating (see Stevens) [W/m^2] 
    F1       ! Coefficient for cloud base heating (see Stevens)[W/m^2]

  ! Short-wave constants
  real, public :: &
    eff_drop_radius, & ! Effective droplet radius [m] 
    gc, & ! Asymmetry parameter, "g" in Duynkerke           [-]
    omega ! Single-scattering albedo                        [-] 

  double precision, public :: &
    amu0, & ! Calculated value of cosine of the solar zenith angle
    slr     ! Fraction of daylight

  real, public, dimension(20) :: &
    Fs_list, &          ! List of Fs0 values for simplified radiation
    cos_solar_zen_list  ! List of cosine of the solar zenith angle values

  logical, public :: &
    l_fix_cos_solar_zen

  integer, public :: &
    nparam

  public :: init_radiation

  private ! Default Scope

! OpenMP directives. These cannot be indented.
!$omp threadprivate(rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
!$omp   kappa, F0, F1, eff_drop_radius, gc, omega, ext_atmos_buffer, Fs_list, &
!$omp   cos_solar_zen_list, l_fix_cos_solar_zen, nparam)

  contains

!-------------------------------------------------------------------------------
  subroutine init_radiation( iunit, namelist_file )
! Description:
!   Setup radiation parameters
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    integer, intent(in) :: iunit

    character(len=*), intent(in) :: &
      namelist_file 

    integer :: k

    namelist /radiation_setting/ &
     rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
     kappa, F0, F1, eff_drop_radius, gc, omega, Fs_list, &
     cos_solar_zen_list, ext_atmos_buffer, l_fix_cos_solar_zen, &
     amu0, slr

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist
    rad_scheme = "none"

    ! BUGSrad parameters
    sol_const =  1367.0 ! W/m^2

    alvdf = 0.1 ! Visible diffuse surface albedo       [-]
    alndr = 0.1 ! Near-IR direct surface albedo        [-]
    alndf = 0.1 ! Near-IR diffuse surface albedo       [-]

    ext_atmos_buffer = 10

    ! Variables used by both schemes
    alvdr = 0.1 ! Visible direct surface albedo        [-]

    ! Simplified radiation parameters
    F0    = 100.0  ! Coefficient for cloud top heating (see Stevens) [W/m^2]
    F1    = 20.0   ! Coefficient for cloud base heating (see Stevens)[W/m^2]
    kappa = 119.0  ! A constant (Duynkerke eqn. 5)                   [m^2/kg]
    gc    = 0.86   ! Asymmetry parameter, "g" in Duynkerke           [-]
    omega = 0.9965 ! Single-scattering albedo                        [-]

    slr  = 1.0d0  ! Fraction of daylight
    amu0 = -999.0 ! Calculated value of cosine of the solar zenith angle

    ! Parameters for fixing the value of cosine of the solar zenith angle
    l_fix_cos_solar_zen = .false.

    ! The incident of incoming SW insolation at cloud top the
    ! direction of the incoming beam (not the vertical)   [W/m^2]
    Fs_list(:) = 0.0 

    ! Cosine of the solar zenith angle [-]
    cos_solar_zen_list(:) = -999.0

    eff_drop_radius = 1.e-5 ! Effective droplet radius [m]

    ! Read the namelist values in
    open(unit=iunit, file=namelist_file, status='old',action='read')
    read(iunit, nml=radiation_setting)
    close(unit=iunit)

    do k = 1, size( cos_solar_zen_list )
      if ( cos_solar_zen_list(k) == -999. ) then
        exit
      else
        nparam = k
      end if
    end do

    return
  end subroutine init_radiation

end module parameters_radiation
