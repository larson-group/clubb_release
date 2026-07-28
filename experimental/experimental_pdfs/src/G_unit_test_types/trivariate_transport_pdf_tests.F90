!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module trivariate_transport_pdf_tests

  implicit none

  public :: trivariate_transport_pdf_tests_driver

  private

  contains

  !=============================================================================
  function trivariate_transport_pdf_tests_driver()

    ! Description:
    !   Verifies that the trivariate transport two-Gaussian driver preserves a
    !   supplied full covariance while returning PSD component covariances and
    !   unequal component w-r_t tilts for a positive-covariance state.
    !
    ! References:
    !   See doc/trivariate_transport_pdf_concept.md.
    !---------------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr, & ! Standard error unit
        fstdout    ! Standard output unit

    use clubb_precision, only: &
        core_rknd

    use trivariate_transport_pdf, only: &
        trivariate_transport_pdf_driver, &
        diagnose_transport_implicit_responders

    use bergmann_pdf_module, only: &
        bergmann_psd_guard

    use pdf_closure_module, only: &
        calc_wpxp3_pdf

    implicit none

    !--------------------------- Local Constants ---------------------------
    integer, parameter :: &
      nz = 1,      & ! Number of vertical levels [#]
      ngrdcol = 1    ! Number of grid columns [#]

    real( kind = core_rknd ), parameter :: &
      relative_tolerance = 5.0e-10_core_rknd, & ! Relative reconstruction tolerance [-]
      absolute_tolerance = 1.0e-12_core_rknd, & ! Absolute reconstruction tolerance [units vary]
      tilt_tolerance = 1.0e-7_core_rknd         ! Minimum resolvable component tilt [-]

    !---------------------------- Local Variables --------------------------
    integer :: &
      trivariate_transport_pdf_tests_driver, & ! Test exit code [#]
      total_failures                           ! Number of failed assertions [#]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wm,      & ! Grid-mean vertical velocity [m/s]
      rtm,     & ! Grid-mean total water mixing ratio [kg/kg]
      thlm,    & ! Grid-mean liquid-water potential temperature [K]
      wp2,     & ! Grid w variance [m2/s2]
      rtp2,    & ! Grid r_t variance [(kg/kg)2]
      thlp2,   & ! Grid theta_l variance [K2]
      wprtp,   & ! Grid w-r_t covariance [m/s kg/kg]
      wpthlp,  & ! Grid w-theta_l covariance [m/s K]
      rtpthlp, & ! Grid r_t-theta_l covariance [kg/kg K]
      Skw,     & ! Grid w skewness [-]
      Skrt,    & ! Grid r_t skewness [-]
      Skthl,   & ! Grid theta_l skewness [-]
      w_1,     & ! G1 w mean [m/s]
      w_2,     & ! G2 w mean [m/s]
      rt_1,    & ! G1 r_t mean [kg/kg]
      rt_2,    & ! G2 r_t mean [kg/kg]
      thl_1,   & ! G1 theta_l mean [K]
      thl_2,   & ! G2 theta_l mean [K]
      varnce_w_1,    & ! G1 w variance [m2/s2]
      varnce_w_2,    & ! G2 w variance [m2/s2]
      varnce_rt_1,   & ! G1 r_t variance [(kg/kg)2]
      varnce_rt_2,   & ! G2 r_t variance [(kg/kg)2]
      varnce_thl_1,  & ! G1 theta_l variance [K2]
      varnce_thl_2,  & ! G2 theta_l variance [K2]
      mixt_frac,     & ! G1 mixture weight [-]
      corr_w_rt_1,   & ! G1 w-r_t correlation [-]
      corr_w_rt_2,   & ! G2 w-r_t correlation [-]
      corr_w_thl_1,   & ! G1 w-theta_l correlation [-]
      corr_w_thl_2,   & ! G2 w-theta_l correlation [-]
      corr_rt_thl_1,  & ! G1 r_t-theta_l correlation [-]
      corr_rt_thl_2,  & ! G2 r_t-theta_l correlation [-]
      wp4,            & ! PDF-integrated w fourth moment [m4/s4]
      wprtp2,         & ! PDF-integrated w-r_t variance transport [m/s (kg/kg)2]
      wpthlp2,        & ! PDF-integrated w-theta_l variance transport [m/s K2]
      wprtpthlp,      & ! PDF-integrated r_t-theta_l covariance transport [m/s (kg/kg) K]
      coef_wp4,       & ! Diagnosed w fourth-moment shape factor [-]
      coef_wprtp2,    & ! Diagnosed r_t variance responder [m/s]
      term_wprtp2,    & ! Diagnosed r_t variance residual [m/s (kg/kg)2]
      coef_wpthlp2,   & ! Diagnosed theta_l variance responder [m/s]
      term_wpthlp2,   & ! Diagnosed theta_l variance residual [m/s K2]
      coef_wprtpthlp, & ! Diagnosed r_t-theta_l covariance responder [m/s]
      term_wprtpthlp    ! Diagnosed r_t-theta_l covariance residual [m/s (kg/kg) K]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wprtp3,          & ! PDF-integrated w-r_t scalar-third-moment flux [m/s (kg/kg)^3]
      expected_wprtp3    ! Independently evaluated Gaussian-mixture flux [m/s (kg/kg)^3]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      moisture_tail_gain, & ! G1 moisture-tail tuning [-]
      center_budget,      & ! Center-separation tuning [-]
      w_direction_scale,  & ! w center-direction tuning [-]
      g1_wrt_capture,     & ! G1 w-r_t covariance capture tuning [-]
      plume_structure_strength ! Coherent transport-plume structure tuning [-]

    real( kind = core_rknd ), dimension(3) :: &
      target_mean,        & ! Input grid mean [mixed units]
      component_mean_1,   & ! Reconstructed G1 mean [mixed units]
      component_mean_2,   & ! Reconstructed G2 mean [mixed units]
      reconstructed_mean    ! Mixture-reconstructed mean [mixed units]

    real( kind = core_rknd ), dimension(3,3) :: &
      target_covariance,        & ! Input grid covariance [mixed units squared]
      component_covariance_1,   & ! Reconstructed G1 covariance [mixed units squared]
      component_covariance_2,   & ! Reconstructed G2 covariance [mixed units squared]
      reconstructed_covariance    ! Mixture-reconstructed covariance [mixed units squared]

    real( kind = core_rknd ), dimension(3) :: &
      bergmann_variance ! Bergmann PSD-guard grid variances [mixed units squared]

    real( kind = core_rknd ), dimension(3,3) :: &
      bergmann_covariance,       & ! Bergmann PSD-guard grid covariance [mixed units squared]
      bergmann_covariance_3,     & ! PSD-safe third Gaussian covariance [mixed units squared]
      bergmann_covariance_outer, & ! PSD-safe outer-pair covariance [mixed units squared]
      bergmann_reconstructed       ! Guard-reconstructed grid covariance [mixed units squared]

    real( kind = core_rknd ) :: &
      weight_1, & ! G1 mixture weight [-]
      weight_2, & ! G2 mixture weight [-]
      baseline_thl_separation, & ! G1-G2 theta_l separation before plume structure [K]
      bergmann_delta_effective   ! PSD-safe Bergmann G3 weight [-]

    logical :: &
      bergmann_valid_input ! True when the Bergmann test covariance is realizable

    !----------------------------- Begin Code -----------------------------
    ! This positive-definite state has a positive w-r_t covariance.  Its
    ! nonzero scalar skewnesses activate the transport-center diagnosis.
    wm(1,1) = 0.15_core_rknd
    rtm(1,1) = 0.010_core_rknd
    thlm(1,1) = 300.0_core_rknd
    wp2(1,1) = 1.00_core_rknd
    rtp2(1,1) = 1.00e-6_core_rknd
    thlp2(1,1) = 1.00_core_rknd
    wprtp(1,1) = 4.00e-4_core_rknd
    wpthlp(1,1) = -0.15_core_rknd
    rtpthlp(1,1) = -2.00e-4_core_rknd
    Skw(1,1) = 0.45_core_rknd
    Skrt(1,1) = 1.10_core_rknd
    Skthl(1,1) = -0.45_core_rknd

    moisture_tail_gain(1) = 1.00_core_rknd
    center_budget(1) = 0.72_core_rknd
    w_direction_scale(1) = 0.20_core_rknd
    g1_wrt_capture(1) = 1.00_core_rknd
    plume_structure_strength(1) = 0.00_core_rknd

    call trivariate_transport_pdf_driver( nz, ngrdcol, &
                                           wm, rtm, thlm, & ! In
                                           wp2, rtp2, thlp2, & ! In
                                           wprtp, wpthlp, rtpthlp, & ! In
                                           Skw, Skrt, Skthl, & ! In
                                           moisture_tail_gain, center_budget, & ! In
                                           w_direction_scale, g1_wrt_capture, & ! In
                                           plume_structure_strength, & ! In
                                           w_1, w_2, rt_1, rt_2, & ! Out
                                           thl_1, thl_2, & ! Out
                                           varnce_w_1, varnce_w_2, & ! Out
                                           varnce_rt_1, varnce_rt_2, & ! Out
                                           varnce_thl_1, varnce_thl_2, & ! Out
                                           mixt_frac, & ! Out
                                           corr_w_rt_1, corr_w_rt_2, & ! Out
                                           corr_w_thl_1, corr_w_thl_2, & ! Out
                                           corr_rt_thl_1, corr_rt_thl_2 ) ! Out

    weight_1 = mixt_frac(1,1)
    weight_2 = 1.0_core_rknd - weight_1

    call calc_wpxp3_pdf( nz, ngrdcol, wm, rtm, w_1, w_2, rt_1, rt_2, &
                          varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, &
                          corr_w_rt_1, corr_w_rt_2, mixt_frac, wprtp3 )
    expected_wprtp3 = mixt_frac * ( &
      (w_1-wm) * ((rt_1-rtm)**3 + 3.0_core_rknd*(rt_1-rtm)*varnce_rt_1) &
      + 3.0_core_rknd*corr_w_rt_1*sqrt(varnce_w_1*varnce_rt_1) &
        *((rt_1-rtm)**2 + varnce_rt_1) ) &
      + (1.0_core_rknd-mixt_frac) * ( &
      (w_2-wm) * ((rt_2-rtm)**3 + 3.0_core_rknd*(rt_2-rtm)*varnce_rt_2) &
      + 3.0_core_rknd*corr_w_rt_2*sqrt(varnce_w_2*varnce_rt_2) &
        *((rt_2-rtm)**2 + varnce_rt_2) )
    target_mean = [ wm(1,1), rtm(1,1), thlm(1,1) ]
    target_covariance = reshape( [ wp2(1,1), wprtp(1,1), wpthlp(1,1), &
                                   wprtp(1,1), rtp2(1,1), rtpthlp(1,1), &
                                   wpthlp(1,1), rtpthlp(1,1), thlp2(1,1) ], [ 3, 3 ] )
    component_mean_1 = [ w_1(1,1), rt_1(1,1), thl_1(1,1) ]
    component_mean_2 = [ w_2(1,1), rt_2(1,1), thl_2(1,1) ]
    call build_component_covariance( varnce_w_1(1,1), varnce_rt_1(1,1), & ! In
                                     varnce_thl_1(1,1), corr_w_rt_1(1,1), & ! In
                                     corr_w_thl_1(1,1), corr_rt_thl_1(1,1), & ! In
                                     component_covariance_1 )                  ! Out
    call build_component_covariance( varnce_w_2(1,1), varnce_rt_2(1,1), & ! In
                                     varnce_thl_2(1,1), corr_w_rt_2(1,1), & ! In
                                     corr_w_thl_2(1,1), corr_rt_thl_2(1,1), & ! In
                                     component_covariance_2 )                  ! Out

    reconstructed_mean = weight_1 * component_mean_1 + weight_2 * component_mean_2
    reconstructed_covariance = weight_1 * ( component_covariance_1 &
                                + outer_product_3( component_mean_1 - target_mean ) ) &
                         + weight_2 * ( component_covariance_2 &
                                + outer_product_3( component_mean_2 - target_mean ) )

    total_failures = 0
    if ( weight_1 <= 0.0_core_rknd .or. weight_1 >= 1.0_core_rknd ) then
      call record_failure( "PDF10 returned a mixture weight outside (0,1).", total_failures )
    end if
    if ( .not. all_close( reshape(wprtp3,[1]), reshape(expected_wprtp3,[1]), &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF-integrated <w'rt'^3> does not match Gaussian moment identity.", &
                           total_failures )
    end if
    if ( .not. all_close( reconstructed_mean, target_mean, &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 failed to reconstruct the trivariate mean.", total_failures )
    end if
    if ( .not. all_close( reshape( reconstructed_covariance, [ 9 ] ), &
                          reshape( target_covariance, [ 9 ] ), &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 failed to reconstruct the full trivariate covariance.", &
                           total_failures )
    end if
    if ( .not. covariance_is_psd( component_covariance_1 ) ) then
      call record_failure( "PDF10 returned a non-PSD G1 covariance matrix.", total_failures )
    end if
    if ( .not. covariance_is_psd( component_covariance_2 ) ) then
      call record_failure( "PDF10 returned a non-PSD G2 covariance matrix.", total_failures )
    end if
    if ( abs( corr_w_rt_1(1,1) ) <= tilt_tolerance .or. &
         abs( corr_w_rt_2(1,1) ) <= tilt_tolerance .or. &
         abs( corr_w_rt_1(1,1) - corr_w_rt_2(1,1) ) <= tilt_tolerance ) then
      call record_failure( "PDF10 did not retain unequal nonzero G1/G2 w-r_t tilts.", &
                           total_failures )
    end if

    ! A negative grid w-r_t covariance must keep G1 untilted even when the
    ! diagnosed center vector itself contributes a negative covariance.
    wprtp(1,1) = -1.00e-4_core_rknd
    Skw(1,1) = -0.85_core_rknd
    Skrt(1,1) = 0.85_core_rknd
    Skthl(1,1) = -0.45_core_rknd
    w_direction_scale(1) = 1.20_core_rknd

    call trivariate_transport_pdf_driver( nz, ngrdcol, &
                                           wm, rtm, thlm, & ! In
                                           wp2, rtp2, thlp2, & ! In
                                           wprtp, wpthlp, rtpthlp, & ! In
                                           Skw, Skrt, Skthl, & ! In
                                           moisture_tail_gain, center_budget, & ! In
                                           w_direction_scale, g1_wrt_capture, & ! In
                                           plume_structure_strength, & ! In
                                           w_1, w_2, rt_1, rt_2, & ! Out
                                           thl_1, thl_2, & ! Out
                                           varnce_w_1, varnce_w_2, & ! Out
                                           varnce_rt_1, varnce_rt_2, & ! Out
                                           varnce_thl_1, varnce_thl_2, & ! Out
                                           mixt_frac, & ! Out
                                           corr_w_rt_1, corr_w_rt_2, & ! Out
                                           corr_w_thl_1, corr_w_thl_2, & ! Out
                                           corr_rt_thl_1, corr_rt_thl_2 ) ! Out

    if ( abs( corr_w_rt_1(1,1) ) > tilt_tolerance ) then
      call record_failure( "PDF10 retained a G1 w-r_t tilt for negative grid covariance.", &
                           total_failures )
    end if

    ! A coherent positive-w, positive-w-r_t tail with a negative w-theta_l
    ! relation represents the mature moist/cool transport regime.  The plume
    ! structure branch should rotate G1 onto that signed direction and keep G2
    ! close to the specified well-mixed covariance constraints.
    wprtp(1,1) = 3.50e-4_core_rknd
    wpthlp(1,1) = -0.35_core_rknd
    rtpthlp(1,1) = 1.50e-4_core_rknd
    Skw(1,1) = 0.90_core_rknd
    Skrt(1,1) = -0.70_core_rknd
    Skthl(1,1) = 0.40_core_rknd
    w_direction_scale(1) = 0.20_core_rknd
    plume_structure_strength(1) = 0.00_core_rknd

    call trivariate_transport_pdf_driver( nz, ngrdcol, &
                                           wm, rtm, thlm, & ! In
                                           wp2, rtp2, thlp2, & ! In
                                           wprtp, wpthlp, rtpthlp, & ! In
                                           Skw, Skrt, Skthl, & ! In
                                           moisture_tail_gain, center_budget, & ! In
                                           w_direction_scale, g1_wrt_capture, & ! In
                                           plume_structure_strength, & ! In
                                           w_1, w_2, rt_1, rt_2, & ! Out
                                           thl_1, thl_2, & ! Out
                                           varnce_w_1, varnce_w_2, & ! Out
                                           varnce_rt_1, varnce_rt_2, & ! Out
                                           varnce_thl_1, varnce_thl_2, & ! Out
                                           mixt_frac, & ! Out
                                           corr_w_rt_1, corr_w_rt_2, & ! Out
                                           corr_w_thl_1, corr_w_thl_2, & ! Out
                                           corr_rt_thl_1, corr_rt_thl_2 ) ! Out
    baseline_thl_separation = thl_1(1,1) - thl_2(1,1)

    plume_structure_strength(1) = 1.00_core_rknd
    call trivariate_transport_pdf_driver( nz, ngrdcol, &
                                           wm, rtm, thlm, & ! In
                                           wp2, rtp2, thlp2, & ! In
                                           wprtp, wpthlp, rtpthlp, & ! In
                                           Skw, Skrt, Skthl, & ! In
                                           moisture_tail_gain, center_budget, & ! In
                                           w_direction_scale, g1_wrt_capture, & ! In
                                           plume_structure_strength, & ! In
                                           w_1, w_2, rt_1, rt_2, & ! Out
                                           thl_1, thl_2, & ! Out
                                           varnce_w_1, varnce_w_2, & ! Out
                                           varnce_rt_1, varnce_rt_2, & ! Out
                                           varnce_thl_1, varnce_thl_2, & ! Out
                                           mixt_frac, & ! Out
                                           corr_w_rt_1, corr_w_rt_2, & ! Out
                                           corr_w_thl_1, corr_w_thl_2, & ! Out
                                           corr_rt_thl_1, corr_rt_thl_2 ) ! Out

    if ( baseline_thl_separation <= 0.0_core_rknd .or. rt_1(1,1) <= rt_2(1,1) .or. &
         thl_1(1,1) >= thl_2(1,1) .or. w_1(1,1) <= w_2(1,1) ) then
      call record_failure( "PDF10 plume structure did not select a moist/cool rising G1.", &
                           total_failures )
    end if

    call build_component_covariance( varnce_w_1(1,1), varnce_rt_1(1,1), & ! In
                                     varnce_thl_1(1,1), corr_w_rt_1(1,1), & ! In
                                     corr_w_thl_1(1,1), corr_rt_thl_1(1,1), & ! In
                                     component_covariance_1 )                  ! Out
    call build_component_covariance( varnce_w_2(1,1), varnce_rt_2(1,1), & ! In
                                     varnce_thl_2(1,1), corr_w_rt_2(1,1), & ! In
                                     corr_w_thl_2(1,1), corr_rt_thl_2(1,1), & ! In
                                     component_covariance_2 )                  ! Out
    if ( .not. covariance_is_psd( component_covariance_1 ) .or. &
         .not. covariance_is_psd( component_covariance_2 ) ) then
      call record_failure( "PDF10 plume structure returned a non-PSD component covariance.", &
                           total_failures )
    end if

    ! The semi-implicit responder split must reconstruct every direct PDF
    ! moment at the frozen state.  A deliberately oversized r_t transport
    ! also checks that the speed cap moves only the excess into the residual.
    wp2(1,1) = 4.0_core_rknd
    rtp2(1,1) = 9.0e-6_core_rknd
    thlp2(1,1) = 1.0_core_rknd
    rtpthlp(1,1) = 1.0e-3_core_rknd
    wp4(1,1) = 40.0_core_rknd
    wprtp2(1,1) = 2.0e-5_core_rknd
    wpthlp2(1,1) = -3.0_core_rknd
    wprtpthlp(1,1) = 4.0e-3_core_rknd

    call diagnose_transport_implicit_responders( nz, ngrdcol, &
                                                  wp2, rtp2, thlp2, rtpthlp, & ! In
                                                  wp4, wprtp2, wpthlp2, wprtpthlp, & ! In
                                                  coef_wp4, & ! Out
                                                  coef_wprtp2, term_wprtp2, & ! Out
                                                  coef_wpthlp2, term_wpthlp2, & ! Out
                                                  coef_wprtpthlp, term_wprtpthlp ) ! Out

    if ( .not. all_close( [ coef_wp4(1,1) ], [ 2.5_core_rknd ], &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 responder did not recover the w fourth-moment shape.", &
                           total_failures )
    end if
    if ( .not. all_close( [ coef_wprtp2(1,1) * rtp2(1,1) + term_wprtp2(1,1) ], &
                          [ wprtp2(1,1) ], relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 responder failed to reconstruct w-r_t variance transport.", &
                           total_failures )
    end if
    if ( .not. all_close( [ coef_wpthlp2(1,1) * thlp2(1,1) + term_wpthlp2(1,1) ], &
                          [ wpthlp2(1,1) ], relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 responder failed to reconstruct w-theta_l variance transport.", &
                           total_failures )
    end if
    if ( .not. all_close( [ coef_wprtpthlp(1,1) * rtpthlp(1,1) + term_wprtpthlp(1,1) ], &
                          [ wprtpthlp(1,1) ], relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 responder failed to reconstruct mixed scalar transport.", &
                           total_failures )
    end if

    wprtp2(1,1) = 1.0e-3_core_rknd
    call diagnose_transport_implicit_responders( nz, ngrdcol, &
                                                  wp2, rtp2, thlp2, rtpthlp, & ! In
                                                  wp4, wprtp2, wpthlp2, wprtpthlp, & ! In
                                                  coef_wp4, & ! Out
                                                  coef_wprtp2, term_wprtp2, & ! Out
                                                  coef_wpthlp2, term_wpthlp2, & ! Out
                                                  coef_wprtpthlp, term_wprtpthlp ) ! Out
    if ( .not. all_close( [ coef_wprtp2(1,1) ], [ 8.0_core_rknd ], &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 responder did not cap an excessive r_t variance speed.", &
                           total_failures )
    end if
    if ( .not. all_close( [ coef_wprtp2(1,1) * rtp2(1,1) + term_wprtp2(1,1) ], &
                          [ wprtp2(1,1) ], relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "PDF10 capped responder did not retain its exact residual.", &
                           total_failures )
    end if

    ! The historical independent Bergmann covariance fractions request an
    ! indefinite G3 for this otherwise positive-definite grid covariance.  The
    ! local guard must retain the supplied covariance exactly while returning
    ! positive-semidefinite component and residual covariance matrices.
    bergmann_variance = [ 1.0_core_rknd, 1.0_core_rknd, 1.0_core_rknd ]
    bergmann_covariance = reshape( [ 1.0_core_rknd, 0.75_core_rknd, 0.70_core_rknd, &
                                     0.75_core_rknd, 1.0_core_rknd, 0.80_core_rknd, &
                                     0.70_core_rknd, 0.80_core_rknd, 1.0_core_rknd ], [ 3, 3 ] )
    call bergmann_psd_guard( bergmann_variance, bergmann_covariance, & ! In
                             0.596_core_rknd, 0.437_core_rknd, 0.785_core_rknd, &
                             1.134_core_rknd, 0.778_core_rknd, &
                             bergmann_delta_effective, bergmann_covariance_3, & ! Out
                             bergmann_covariance_outer, bergmann_valid_input )
    bergmann_reconstructed = ( 1.0_core_rknd - bergmann_delta_effective ) &
                               * bergmann_covariance_outer &
                             + bergmann_delta_effective * bergmann_covariance_3
    if ( .not. bergmann_valid_input ) then
      call record_failure( "Bergmann PSD guard rejected a positive-definite grid covariance.", &
                           total_failures )
    end if
    if ( .not. covariance_is_psd( bergmann_covariance_3 ) .or. &
         .not. covariance_is_psd( bergmann_covariance_outer ) ) then
      call record_failure( "Bergmann PSD guard returned an indefinite component covariance.", &
                           total_failures )
    end if
    if ( .not. all_close( reshape( bergmann_reconstructed, [ 9 ] ), &
                          reshape( bergmann_covariance, [ 9 ] ), &
                          relative_tolerance, absolute_tolerance ) ) then
      call record_failure( "Bergmann PSD guard changed the supplied grid covariance.", &
                           total_failures )
    end if
    if ( bergmann_delta_effective > 0.596_core_rknd ) then
      call record_failure( "Bergmann PSD guard increased the requested G3 mixture weight.", &
                           total_failures )
    end if

    if ( total_failures == 0 ) then
      write(fstdout,'(A)') "Trivariate transport PDF unit test: Success!"
      trivariate_transport_pdf_tests_driver = 0
    else
      write(fstderr,'(A,I0,A)') "Trivariate transport PDF unit test: ", total_failures, &
                                " failure(s)."
      trivariate_transport_pdf_tests_driver = 1
    end if

  end function trivariate_transport_pdf_tests_driver

  !=============================================================================
  subroutine build_component_covariance( variance_w, variance_rt, variance_thl, &
                                         corr_w_rt, corr_w_thl, corr_rt_thl, &
                                         covariance )

    ! Description:
    !   Rebuilds one physical component covariance from the variances and
    !   correlations returned by the transport-PDF driver.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), intent(in) :: &
      variance_w,   & ! Component w variance [m2/s2]
      variance_rt,  & ! Component r_t variance [(kg/kg)2]
      variance_thl, & ! Component theta_l variance [K2]
      corr_w_rt,    & ! Component w-r_t correlation [-]
      corr_w_thl,   & ! Component w-theta_l correlation [-]
      corr_rt_thl     ! Component r_t-theta_l correlation [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      covariance ! Component covariance [mixed units squared]

    !----------------------------- Begin Code -----------------------------
    covariance = 0.0_core_rknd
    covariance(1,1) = variance_w
    covariance(2,2) = variance_rt
    covariance(3,3) = variance_thl
    covariance(1,2) = corr_w_rt * sqrt( variance_w * variance_rt )
    covariance(2,1) = covariance(1,2)
    covariance(1,3) = corr_w_thl * sqrt( variance_w * variance_thl )
    covariance(3,1) = covariance(1,3)
    covariance(2,3) = corr_rt_thl * sqrt( variance_rt * variance_thl )
    covariance(3,2) = covariance(2,3)

  end subroutine build_component_covariance

  !=============================================================================
  function outer_product_3( vector ) result( matrix )

    ! Description:
    !   Returns the outer product of one three-component vector.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3), intent(in) :: &
      vector ! Input vector [mixed units]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3,3) :: &
      matrix ! Vector outer product [mixed units squared]

    integer :: &
      row, & ! Matrix row index [#]
      column  ! Matrix column index [#]

    !----------------------------- Begin Code -----------------------------
    do column = 1, 3
      do row = 1, 3
        matrix(row,column) = vector(row) * vector(column)
      end do
    end do

  end function outer_product_3

  !=============================================================================
  function all_close( actual, expected, relative_tolerance, absolute_tolerance ) &
    result( values_match )

    ! Description:
    !   Checks scalar values in an array using combined relative and absolute
    !   tolerances so the small physical r_t moments are checked proportionally.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      actual,   & ! Calculated values [units vary]
      expected    ! Reference values [units vary]

    real( kind = core_rknd ), intent(in) :: &
      relative_tolerance, & ! Relative tolerance [-]
      absolute_tolerance    ! Absolute tolerance [units vary]

    !--------------------------- Output Variables --------------------------
    logical :: &
      values_match ! True when every value is within tolerance [-]

    !----------------------------- Begin Code -----------------------------
    values_match = all( abs( actual - expected ) &
                        <= absolute_tolerance + relative_tolerance * abs( expected ) )

  end function all_close

  !=============================================================================
  function covariance_is_psd( covariance ) result( is_psd )

    ! Description:
    !   Tests all principal minors of a three-by-three covariance after scaling
    !   it to a correlation matrix.  This checks component PSD independently of
    !   the closure's private realizability helper.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance ! Symmetric component covariance [mixed units squared]

    !--------------------------- Output Variables --------------------------
    logical :: &
      is_psd ! True when all covariance principal minors are nonnegative [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), parameter :: &
      tolerance = 1.0e-10_core_rknd ! Principal-minor tolerance [-]

    real( kind = core_rknd ), dimension(3) :: &
      standard_deviation ! Marginal component standard deviations [mixed units]

    real( kind = core_rknd ), dimension(3,3) :: &
      correlation ! Scaled component covariance [-]

    real( kind = core_rknd ) :: &
      determinant ! Determinant of the correlation matrix [-]

    integer :: &
      row, & ! Matrix row index [#]
      column  ! Matrix column index [#]

    !----------------------------- Begin Code -----------------------------
    if ( covariance(1,1) <= tolerance .or. covariance(2,2) <= tolerance .or. &
         covariance(3,3) <= tolerance ) then
      is_psd = .false.
      return
    end if

    standard_deviation = [ sqrt( covariance(1,1) ), sqrt( covariance(2,2) ), &
                           sqrt( covariance(3,3) ) ]
    do column = 1, 3
      do row = 1, 3
        correlation(row,column) = covariance(row,column) / &
                                  ( standard_deviation(row) * standard_deviation(column) )
      end do
    end do
    determinant = correlation(1,1) * ( correlation(2,2) * correlation(3,3) &
                                      - correlation(2,3) * correlation(3,2) ) &
                - correlation(1,2) * ( correlation(2,1) * correlation(3,3) &
                                      - correlation(2,3) * correlation(3,1) ) &
                + correlation(1,3) * ( correlation(2,1) * correlation(3,2) &
                                      - correlation(2,2) * correlation(3,1) )
    is_psd = correlation(1,1) >= -tolerance .and. &
             correlation(2,2) >= -tolerance .and. &
             correlation(3,3) >= -tolerance .and. &
             correlation(1,1) * correlation(2,2) - correlation(1,2)**2 >= -tolerance .and. &
             correlation(1,1) * correlation(3,3) - correlation(1,3)**2 >= -tolerance .and. &
             correlation(2,2) * correlation(3,3) - correlation(2,3)**2 >= -tolerance .and. &
             determinant >= -tolerance

  end function covariance_is_psd

  !=============================================================================
  subroutine record_failure( message, total_failures )

    ! Description:
    !   Records and reports one failed trivariate-transport-PDF assertion.
    !---------------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr

    implicit none

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: &
      message ! Failed assertion description [-]

    !------------------------ Input/Output Variables -----------------------
    integer, intent(inout) :: &
      total_failures ! Accumulated failed assertions [#]

    !----------------------------- Begin Code -----------------------------
    total_failures = total_failures + 1
    write(fstderr,'(A)') trim( message )

  end subroutine record_failure

end module trivariate_transport_pdf_tests
