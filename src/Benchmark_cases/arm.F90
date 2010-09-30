!----------------------------------------------------------------------
! $Id$
module arm

  !       Description:
  !       Contains subroutines for the GCSS ARM case.
  !----------------------------------------------------------------------

  implicit none

  public :: arm_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine arm_sfclyr( time, z, dn0, thlm_sfc, ubar,  & 
                         wpthlp_sfc, wprtp_sfc, ustar )

  !       Description:
  !       This subroutine computes surface fluxes of horizontal momentum,
  !       heat and moisture according to GCSS ARM specifications
  !----------------------------------------------------------------------

  use constants_clubb, only: Cp, Lv, grav ! Variable(s)

  use stats_precision, only: time_precision ! Variable(s)

  use diag_ustar_module, only: diag_ustar ! Variable(s)

  use surface_flux, only: compute_ht_mostr_flux ! Procedures

  implicit none

  intrinsic :: max, sqrt

  ! Constants
  real, parameter ::  & 
    z0 = 0.035  ! ARM Cu momentum roughness height
  integer, parameter :: &
    ntimes = 7

  ! Input variables
  real(kind=time_precision), intent(in) ::  & 
    time            ! Current time          [s]

  real, intent(in) ::  & 
    z,               & ! Height at zt(2)       [m]
    dn0,             & ! Density at zm(1)      [kg/m^3]
    thlm_sfc,        & ! Theta_l at zt(2)      [K]
    ubar

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar          ! surface friction velocity [m/s]

  ! Local variables
  real ::  & 
    heat_flx, moisture_flx, & 
    heat_flx2, moisture_flx2, & 
    bflx

  !-----------------BEGIN CODE-------------------------

  ! Compute heat and moisture fluxes from ARM data in (W/m2)
   call compute_ht_mostr_flux( time, ntimes, &
                               heat_flx, moisture_flx )

  ! Compute momentum fluxes

  ! Convert heat_flx and moisture_flx to natural units
  heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
  moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

  ! Heat flux in units of (m2/s3) (needed by diag_ustar)
  bflx = grav/thlm_sfc * heat_flx2

  ! Surface winds

  ! Compute ustar
  ustar = diag_ustar( z, bflx, ubar, z0 )

  ! Assign fluxes

  wpthlp_sfc = heat_flx2
  wprtp_sfc  = moisture_flx2

  return
  end subroutine arm_sfclyr

end module arm
