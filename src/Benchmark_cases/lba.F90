!----------------------------------------------------------------------
! $Id$
module lba

  ! Description:
  !   Contains subroutines for the LBA case.
  ! References:
  !   http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: lba_tndcy, lba_sfclyr

  contains

  !----------------------------------------------------------------------
  subroutine lba_tndcy( thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set theta and water tendencies for LBA case.

    !       References:
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    !----------------------------------------------------------------------

    use clubb_api_module, only: gr !  Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)


    implicit none

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing [units vary]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,edsclr_dim) :: & 
      edsclrm_forcing ! Passive eddy-scalar forcing [units vary]

    ! ---- Begin Code ----

    ! Large-scale temperature tendency
    thlm_forcing(:) = 0.0_core_rknd

    ! Large-scale advective moisture tendency
    rtm_forcing(:) = 0.0_core_rknd

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine lba_tndcy

  !----------------------------------------------------------------------
  subroutine lba_sfclyr( time_current, time_initial, z, & 
                         rho_sfc, thlm_sfc, ubar,  & 
                         wpthlp_sfc, wprtp_sfc, ustar )

    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS BOMEX specifications

    !       References:
    !       Grabowski, et al. (2005)
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    !----------------------------------------------------------------------

    use constants_clubb, only: pi, grav, sec_per_hr ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use sfc_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

    implicit none

    intrinsic :: max, sqrt

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  & 
      z0    = 0.035_core_rknd  ! ARM mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time_current, & ! Current time              [s]
      time_initial    ! Start time of model run   [s]

    real( kind = core_rknd ), intent(in) ::  & 
      z,         & ! Height at zt=2      [m] 
      rho_sfc,   & ! Density at zm=1     [kg/m^3] 
      thlm_sfc,  & ! thlm at (2)         [m/s]
      ubar

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real(kind=time_precision) ::  & 
      time      ! Elapsed time of model run    [s]

    real( kind = core_rknd ) :: ft, bflx

    ! Compute heat and moisture fluxes
    ! From Table A.1.
    time = time_current - time_initial
    ft = real( max( 0._core_rknd,  & 
                 cos( 0.5_core_rknd * pi * ( (5.25_core_rknd - &
                 real( time,kind=core_rknd)/sec_per_hr) / 5.25_core_rknd ) ) & 
            ), kind = core_rknd ) ! Known magic number

    wpthlp_sfc =  convert_sens_ht_to_km_s( ( 270._core_rknd * &
           ft**1.5_core_rknd ), rho_sfc ) ! Known magic number
    wprtp_sfc  =  convert_latent_ht_to_m_s( ( 554._core_rknd * &
           ft**1.3_core_rknd ), rho_sfc ) ! Known magic number

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )


    return
  end subroutine lba_sfclyr


end module lba

