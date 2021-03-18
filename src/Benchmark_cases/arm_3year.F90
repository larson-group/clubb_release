!----------------------------------------------------------------------
! $Id$
module arm_3year

! Description:
!   Contains subroutine(s) for the ARM 3 Year case.
!----------------------------------------------------------------------


  implicit none

  public :: arm_3year_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine arm_3year_sfclyr( time, z, rho_sfc, & 
                               thlm_sfc, ubar,  & 
                               wpthlp_sfc, wprtp_sfc, ustar )
    ! Description:
    !   This subroutine computes surface fluxes of horizontal momentum,
    !   heat and moisture according to GCSS ARM specifications.
    ! References:
    !   None.
    !----------------------------------------------------------------------

    use constants_clubb, only: grav ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use sfc_flux, only: compute_ht_mostr_flux, &
                            convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedures

    use time_dependent_input, only: time_sfc_given ! Variable(s)

    implicit none

    ! External
    intrinsic :: size

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  &  
      z0    = 0.035_core_rknd   ! ARM Cu mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real( kind = core_rknd ), intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho_sfc,   & ! Density at zm=1     [kg/m^3] 
      ubar,      & ! mean sfc wind speed [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real( kind = core_rknd ) :: bflx, heat_flx, moisture_flx

    
    !-------------BEGIN CODE--------------

    ! Compute heat and moisture fluxes from ARM data in (W/m2)
    call compute_ht_mostr_flux( time, size( time_sfc_given ), &
                                heat_flx, moisture_flx )

    ! Convert W/m^2 into w'thl' w'rt' units
    wpthlp_sfc = convert_sens_ht_to_km_s( heat_flx, rho_sfc )     ! (K m/s)
    wprtp_sfc  = convert_latent_ht_to_m_s( moisture_flx, rho_sfc ) ! (kg m/ kg s)

    ! Compute momentum fluxes using ARM Cu formulae

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )

    return
  end subroutine arm_3year_sfclyr

end module arm_3year
