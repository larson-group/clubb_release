!-------------------------------------------------------------------------------
! $Id$
module parameters_radiation

! Description:
!   Parameters for radiation schemes

! References:
!   None
!-------------------------------------------------------------------------------
  implicit none

  character(len=20), public :: & 
    rad_scheme  ! Either BUGSrad, simplified, or simplied_bomex

  double precision, dimension(1), public :: &
    sol_const ! Solar constant

  real, public :: &
    radiation_top ! The top of the atmosphere fed into a radiation scheme.
    !               The computational grid should be extended to reach this
    !               altitude.

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
    slr     ! Fraction of daylight

  real, public, dimension(20) :: &
    Fs_values, &           ! List of Fs0 values for simplified radiation
    cos_solar_zen_times, & ! List of cosine of the solar zenith angle times
    cos_solar_zen_values   ! List of cosine of the solar zenith angle values

  logical, public :: &
    l_fix_cos_solar_zen, l_sw_radiation

  logical, public :: &
    l_rad_above_cloud ! Use DYCOMS II RF02 heaviside step function

  integer, public :: &
    nparam

  public :: init_radiation

  private ! Default Scope

! OpenMP directives. These cannot be indented.
!$omp threadprivate(rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
!$omp   kappa, F0, F1, eff_drop_radius, gc, omega, radiation_top, Fs_values, &
!$omp   l_rad_above_cloud, cos_solar_zen_list, l_fix_cos_solar_zen, nparam, &
!$omp   l_sw_raditiation)

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
     kappa, F0, F1, eff_drop_radius, gc, omega, Fs_values, &
     cos_solar_zen_values, cos_solar_zen_times, &
     radiation_top, l_fix_cos_solar_zen, l_sw_radiation, &
     slr, l_rad_above_cloud

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist
    rad_scheme = "none"

    ! BUGSrad parameters
    sol_const =  1367.0 ! W/m^2

    alvdf = 0.1 ! Visible diffuse surface albedo       [-]
    alndr = 0.1 ! Near-IR direct surface albedo        [-]
    alndf = 0.1 ! Near-IR diffuse surface albedo       [-]

    ! 50000m is the top of the U.S. Standard Atmosphere data used
    ! in CLUBB.
    radiation_top = 50000.! [m]

    ! Variables used by both schemes
    alvdr = 0.1 ! Visible direct surface albedo        [-]

    ! Simplified radiation parameters
    F0    = 100.0  ! Coefficient for cloud top heating (see Stevens) [W/m^2]
    F1    = 20.0   ! Coefficient for cloud base heating (see Stevens)[W/m^2]
    kappa = 119.0  ! A constant (Duynkerke eqn. 5)                   [m^2/kg]
    gc    = 0.86   ! Asymmetry parameter, "g" in Duynkerke           [-]
    omega = 0.9965 ! Single-scattering albedo                        [-]

    slr  = 1.0d0  ! Fraction of daylight

    l_rad_above_cloud = .false. ! For the heaviside step function
    l_sw_radiation = .false. ! Set to true to enable shortwave radiation

    ! Parameters for fixing the value of cosine of the solar zenith angle
    l_fix_cos_solar_zen = .false.

    ! The incident of incoming SW insolation at cloud top the
    ! direction of the incoming beam (not the vertical)   [W/m^2]
    Fs_values(:) = 0.0 

    cos_solar_zen_values(:) = -999.0 ! Cosine of the solar zenith angle [-]
    cos_solar_zen_times(:)  = -999.0 ! Simulation times corresponding to above [s]

    eff_drop_radius = 1.e-5 ! Effective droplet radius [m]

    ! Read the namelist values in
    open(unit=iunit, file=namelist_file, status='old',action='read')
    read(iunit, nml=radiation_setting)
    close(unit=iunit)

    do k = 1, size( cos_solar_zen_values )
      if ( cos_solar_zen_values(k) == -999. ) then
        exit
      else
        nparam = k
      end if
    end do

    return
  end subroutine init_radiation

end module parameters_radiation
