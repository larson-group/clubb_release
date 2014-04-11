!----------------------------------------------------------------------
! $Id$
module arm_0003

!       Description:
!       Contains subroutines for the March 2000 IOP ARM case.
!
!       References:
!       Xie, S., et al. (2005), Simulations of midlatitude frontal 
!       clouds by single-column and cloud-resolving models during the
!       Atmospheric Radiation Measurement March 2000 cloud intensive
!       operational period, J. Geophys. Res., 110, D15S03, 
!       doi:10.1029/2004JD005119.
!       http://www.agu.org/journals/jd/jd0506/2004JD005119/2004JD005119.pdf
!----------------------------------------------------------------------

  implicit none

  public :: arm_0003_sfclyr

  private ! Default Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_0003_sfclyr( time, z, rho_sfc, & 
                              thlm_sfc, ubar,  & 
                              wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !
    !       References
    !       Xie, S., et al. (2005), Simulations of midlatitude frontal 
    !       clouds by single-column and cloud-resolving models during the
    !       Atmospheric Radiation Measurement March 2000 cloud intensive
    !       operational period, J. Geophys. Res., 110, D15S03, 
    !       doi:10.1029/2004JD005119.
    !       http://www.agu.org/journals/jd/jd0506/2004JD005119/2004JD005119.pdf
    !----------------------------------------------------------------------

    use constants_clubb, only: grav ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use surface_flux, only: compute_ht_mostr_flux, &
                            convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedures

    use time_dependent_input, only: time_sfc_given ! Variable(s)
    
    implicit none

    intrinsic :: max, sqrt

    real( kind = core_rknd ), parameter ::  &  
      z0    = 0.035_core_rknd   ! ARM Cu mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real( kind = core_rknd ), intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho_sfc,      & ! Density at zm=1     [kg/m^3] 
      ubar, &
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real( kind = core_rknd ) :: bflx, heat_flx, moisture_flx

    !----------------------------------------------------------------------

    ! Compute heat and moisture fluxes from ARM data in (W/m2)
    call compute_ht_mostr_flux( time, size( time_sfc_given ), &
                                heat_flx, moisture_flx )

    ! Convert W/m^2 into w'thl' w'rt' units
    wpthlp_sfc = convert_sens_ht_to_km_s( heat_flx, rho_sfc )     ! (K m/s)
    wprtp_sfc  = convert_latent_ht_to_m_s( moisture_flx, rho_sfc )  ! (kg m/ kg s)

    ! Compute momentum fluxes using ARM Cu formulae

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )


    return
  end subroutine arm_0003_sfclyr

!----------------------------------------------------------------------
end module arm_0003
