!----------------------------------------------------------------------
! $Id$
module arm

  !       Description:
  !       Contains subroutines for the GCSS ARM case.
  ! 
  !       References:
  !       Brown et al., 2002:  Large-eddy simulation of the diurnal cycle
  !       of shallow cumulus convection over land. Quart. J. Roy. 
  !       Meteor.  Soc., 128, 1075-1093. 
  !       http://www.atmos.washington.edu/~breth/GCSS/Brown_etal_
  !       DiurnalShCu_QJ2002.pdf
  !----------------------------------------------------------------------

  implicit none

  public :: arm_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine arm_sfclyr( ngrdcol, time, z, &
                         thlm_sfc, ubar,  & 
                         wpthlp_sfc, wprtp_sfc, ustar )

  !       Description:
  !       This subroutine computes surface fluxes of horizontal momentum,
  !       heat and moisture according to GCSS ARM specifications
  !   
  !       References:
  !       See module comments.
  !----------------------------------------------------------------------

  use constants_clubb, only: grav ! Variable(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  use diag_ustar_module, only: diag_ustar ! Variable(s)

  use sfc_flux, only: compute_ht_mostr_flux, &
                          convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedures


  implicit none

  intrinsic :: max, sqrt

  ! Constants
  real( kind = core_rknd ), parameter ::  & 
    z0 = 0.035_core_rknd  ! ARM Cu momentum roughness height

  integer, parameter :: &
    ntimes = 7

  ! Input variables
  integer, intent(in) :: &
  ngrdcol

  real(kind=time_precision), intent(in) ::  & 
    time            ! Current time          [s]

  real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
    z,               & ! Height at zt(2)       [m]
    thlm_sfc,        & ! Theta_l at zt(2)      [K]
    ubar

  ! Output variables
  real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar          ! surface friction velocity [m/s]

  ! Local variables
  real( kind = core_rknd ) ::  & 
    heat_flx, moisture_flx, & 
    heat_flx2, moisture_flx2, & 
    bflx, &
    rho_sfc ! Density at zm(1)      [kg/m^3]

  integer :: i

  !-----------------BEGIN CODE-------------------------

  ! Compute heat and moisture fluxes from ARM data in (W/m2)
  call compute_ht_mostr_flux( time, ntimes, &
                              heat_flx, moisture_flx )

  rho_sfc = 1.1_core_rknd

  ! Convert heat_flx and moisture_flx to natural units
  heat_flx2     = convert_sens_ht_to_km_s( heat_flx, rho_sfc )    ! (K m/s)
  moisture_flx2 = convert_latent_ht_to_m_s( moisture_flx, rho_sfc ) ! (m/s)

  ! Compute momentum fluxes
  do i = 1, ngrdcol

    ! Heat flux in units of (m2/s3) (needed by diag_ustar)
    bflx = grav / thlm_sfc(i) * heat_flx2

    ! Surface winds

    ! Compute ustar
    ustar(i) = diag_ustar( z(i), bflx, ubar(i), z0 )

    ! Assign fluxes

    wpthlp_sfc(i) = heat_flx2
    wprtp_sfc(i)  = moisture_flx2

  end do

  return

  end subroutine arm_sfclyr

end module arm
