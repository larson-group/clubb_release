!----------------------------------------------------------------------
! $Id: twp_ice.F90 3363 2009-04-02 21:42:22Z dschanen@uwm.edu $
module twp_ice
  !
  !       Description:
  !       Contains subroutines for the Jan. 2006 TWP_ICE case.
  !----------------------------------------------------------------------

  implicit none

  public :: twp_ice_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine twp_ice_sfclyr( z, sst, exner_sfc, thlm_sfc, & 
                              um_sfc, vm_sfc, rtm, psfc,  & 
                              upwp_sfc, vpwp_sfc, & 
                              wpthlp_sfc, wprtp_sfc, ustar, & 
                              wpsclrp_sfc, wpedsclrp_sfc )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, grav ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use surface_flux, only: compute_ubar, compute_momentum_flux, &
                              compute_wpthlp_sfc, compute_wprtp_sfc
    implicit none

    intrinsic :: max, sqrt

    ! Constants
    real, parameter :: & 
      C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
      C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
      z0      = 0.00015      ! Roughness length, defined by ATEX specification

    real, intent(in) ::  & 
      z,             & ! Height at zt=2      [s] 
      sst,           & ! Sea surface temp    [K]
      exner_sfc,     & ! Exner function at (2) 
      um_sfc,        & ! um at (2)           [m/s]
      vm_sfc,        & ! vm at (2)           [m/s]
      thlm_sfc,      & ! thlm at (2)         [m/s]
      rtm,           & ! rt at (2)           [kg/kg]
      psfc             ! surface pressure    [Pa]

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
      vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    real, intent(out), dimension(sclr_dim) ::  & 
      wpsclrp_sfc      ! Passive scalar surface flux      [units m/s]

    real, intent(out), dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]

    ! Internal variables
    real :: & 
      ubar, & ! This is root (u^2 + v^2), per ATEX and RICO spec.
      Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
      Cq      ! This is C_q_20 scaled to the height of the lowest model level.
    !----------------------------------------------------------------------

    ! Declare the value of ustar.
    ustar = 0.3

    ! Define variable values
    ubar = compute_ubar( um_sfc, vm_sfc )

    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    !Cm   = C_m_20 * ((log(20/z0))/(log(z/z0))) * &
    !       ((log(20/z0))/(log(z/z0)))
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch   = C_h_20 * ((log(20/z0))/(log(z/z0))) * & 
           ((log(20/z0))/(log(z/z0)))
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq   = C_q_20 * ((log(20/z0))/(log(z/z0))) * & 
           ((log(20/z0))/(log(z/z0)))

    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, sst, exner_sfc )
    wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(psfc,sst) )
    !upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
    !vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

    ! Compute momentum fluxes
    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    ! Avoid uninitialized memory
    wpsclrp_sfc(:)   = 0
    wpedsclrp_sfc(:) = 0

    ! Let passive scalars be equal to rt and theta_l for testing
    if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

    if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
    if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

    return
  end subroutine twp_ice_sfclyr
end module twp_ice
