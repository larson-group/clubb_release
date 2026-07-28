!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module trivariate_transport_pdf

  ! Description:
  !   Moment-driven trivariate two-Gaussian transport PDF for w, r_t, and
  !   theta_l.  The closure diagnoses both component covariance matrices, so
  !   their correlations must be retained by the caller.
  !-------------------------------------------------------------------------------

  implicit none

  public :: trivariate_transport_pdf_driver, &
            diagnose_transport_implicit_responders, &
            cholesky_3x3

  private

  contains

  !=============================================================================
  subroutine trivariate_transport_pdf_driver( nz, ngrdcol,                    &
                                              wm, rtm, thlm,                   &
                                              wp2, rtp2, thlp2,                &
                                              wprtp, wpthlp, rtpthlp,          &
                                              Skw, Skrt, Skthl,                &
                                              moisture_tail_gain,              &
                                              center_budget,                    &
                                              w_direction_scale,               &
                                              g1_wrt_capture,                   &
                                              plume_structure_strength,         &
                                              w_1, w_2, rt_1, rt_2,            &
                                              thl_1, thl_2,                    &
                                              varnce_w_1, varnce_w_2,          &
                                              varnce_rt_1, varnce_rt_2,        &
                                              varnce_thl_1, varnce_thl_2,      &
                                              mixt_frac,                        &
                                              corr_w_rt_1, corr_w_rt_2,        &
                                              corr_w_thl_1, corr_w_thl_2,      &
                                              corr_rt_thl_1, corr_rt_thl_2 )

    ! Description:
    !   Diagnoses the two component means, variances, and all three independent
    !   component correlations from the grid mean, covariance, and marginal
    !   skewnesses.  The algebraic PSD cap limits a requested contrast directly;
    !   no nonlinear solve or iteration is used.
    !
    ! References:
    !   See doc/trivariate_transport_pdf_concept.md.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels [#]
      ngrdcol    ! Number of grid columns [#]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,      & ! Mean vertical velocity [m/s]
      rtm,     & ! Mean total water mixing ratio [kg/kg]
      thlm,    & ! Mean liquid-water potential temperature [K]
      wp2,     & ! Variance of vertical velocity [m2/s2]
      rtp2,    & ! Variance of total water mixing ratio [(kg/kg)2]
      thlp2,   & ! Variance of liquid-water potential temperature [K2]
      wprtp,   & ! Covariance of w and r_t [m/s kg/kg]
      wpthlp,  & ! Covariance of w and theta_l [m/s K]
      rtpthlp, & ! Covariance of r_t and theta_l [kg/kg K]
      Skw,     & ! Skewness of vertical velocity [-]
      Skrt,    & ! Skewness of total water mixing ratio [-]
      Skthl      ! Skewness of liquid-water potential temperature [-]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      moisture_tail_gain, & ! Maps r_t skewness to the G1 weight [-]
      center_budget,      & ! Covariance-metric center-separation budget [-]
      w_direction_scale,  & ! Relative w contribution to center direction [-]
      g1_wrt_capture,     & ! Positive residual w-r_t covariance in G1 [-]
      plume_structure_strength ! Coherent transport-plume structure strength [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      w_1,           & ! Mean w in G1 [m/s]
      w_2,           & ! Mean w in G2 [m/s]
      rt_1,          & ! Mean r_t in G1 [kg/kg]
      rt_2,          & ! Mean r_t in G2 [kg/kg]
      thl_1,         & ! Mean theta_l in G1 [K]
      thl_2,         & ! Mean theta_l in G2 [K]
      varnce_w_1,    & ! Variance of w in G1 [m2/s2]
      varnce_w_2,    & ! Variance of w in G2 [m2/s2]
      varnce_rt_1,   & ! Variance of r_t in G1 [(kg/kg)2]
      varnce_rt_2,   & ! Variance of r_t in G2 [(kg/kg)2]
      varnce_thl_1,  & ! Variance of theta_l in G1 [K2]
      varnce_thl_2,  & ! Variance of theta_l in G2 [K2]
      mixt_frac,     & ! Weight of G1 [-]
      corr_w_rt_1,   & ! Correlation of w and r_t in G1 [-]
      corr_w_rt_2,   & ! Correlation of w and r_t in G2 [-]
      corr_w_thl_1,  & ! Correlation of w and theta_l in G1 [-]
      corr_w_thl_2,  & ! Correlation of w and theta_l in G2 [-]
      corr_rt_thl_1, & ! Correlation of r_t and theta_l in G1 [-]
      corr_rt_thl_2    ! Correlation of r_t and theta_l in G2 [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(3) :: &
      mean,        & ! Grid-mean trivariate state [mixed units]
      variance,    & ! Grid marginal variances [mixed units squared]
      skewness,    & ! Grid marginal skewnesses [-]
      mean_1,      & ! G1 mean trivariate state [mixed units]
      mean_2         ! G2 mean trivariate state [mixed units]

    real( kind = core_rknd ), dimension(3,3) :: &
      covariance,   & ! Grid covariance matrix [mixed units squared]
      covariance_1, & ! G1 covariance matrix [mixed units squared]
      covariance_2    ! G2 covariance matrix [mixed units squared]

    real( kind = core_rknd ) :: &
      weight_1       ! G1 mixture weight [-]

    integer :: &
      i, & ! Grid-column loop index [#]
      k    ! Vertical-level loop index [#]

    !----------------------------- Begin Code -----------------------------
    do k = 1, nz
      do i = 1, ngrdcol
        mean = [ wm(i,k), rtm(i,k), thlm(i,k) ]
        variance = [ wp2(i,k), rtp2(i,k), thlp2(i,k) ]
        skewness = [ Skw(i,k), Skrt(i,k), Skthl(i,k) ]

        covariance = 0.0_core_rknd
        covariance(1,1) = variance(1)
        covariance(2,2) = variance(2)
        covariance(3,3) = variance(3)
        covariance(1,2) = wprtp(i,k)
        covariance(2,1) = wprtp(i,k)
        covariance(1,3) = wpthlp(i,k)
        covariance(3,1) = wpthlp(i,k)
        covariance(2,3) = rtpthlp(i,k)
        covariance(3,2) = rtpthlp(i,k)

        call diagnose_transport_pdf_state( mean, variance, covariance, skewness, & ! In
                                           moisture_tail_gain(i), center_budget(i), & ! In
                                           w_direction_scale(i), g1_wrt_capture(i), & ! In
                                           plume_structure_strength(i),              & ! In
                                           mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                           weight_1 )                                    ! Out

        w_1(i,k) = mean_1(1)
        w_2(i,k) = mean_2(1)
        rt_1(i,k) = mean_1(2)
        rt_2(i,k) = mean_2(2)
        thl_1(i,k) = mean_1(3)
        thl_2(i,k) = mean_2(3)
        varnce_w_1(i,k) = covariance_1(1,1)
        varnce_w_2(i,k) = covariance_2(1,1)
        varnce_rt_1(i,k) = covariance_1(2,2)
        varnce_rt_2(i,k) = covariance_2(2,2)
        varnce_thl_1(i,k) = covariance_1(3,3)
        varnce_thl_2(i,k) = covariance_2(3,3)
        mixt_frac(i,k) = weight_1

        corr_w_rt_1(i,k) = correlation_from_covariance( covariance_1, 1, 2 )
        corr_w_rt_2(i,k) = correlation_from_covariance( covariance_2, 1, 2 )
        corr_w_thl_1(i,k) = correlation_from_covariance( covariance_1, 1, 3 )
        corr_w_thl_2(i,k) = correlation_from_covariance( covariance_2, 1, 3 )
        corr_rt_thl_1(i,k) = correlation_from_covariance( covariance_1, 2, 3 )
        corr_rt_thl_2(i,k) = correlation_from_covariance( covariance_2, 2, 3 )
      end do
    end do

  end subroutine trivariate_transport_pdf_driver

  !=============================================================================
  subroutine diagnose_transport_implicit_responders( nz, ngrdcol, &
                                                      wp2, rtp2, thlp2, rtpthlp, &
                                                      wp4, wprtp2, wpthlp2, wprtpthlp, &
                                                      coef_wp4_implicit, &
                                                      coef_wprtp2_implicit, &
                                                      term_wprtp2_explicit, &
                                                      coef_wpthlp2_implicit, &
                                                      term_wpthlp2_explicit, &
                                                      coef_wprtpthlp_implicit, &
                                                      term_wprtpthlp_explicit )

    ! Description:
    !   Form frozen-geometry, semi-implicit turbulent-advection responders
    !   from the higher-order moments directly integrated from trivariate transport PDF's
    !   diagnosed two-Gaussian state.  The responders preserve each supplied
    !   moment at the coefficient state:
    !
    !   <w'^4> = coef_wp4_implicit * <w'^2>^2,
    !   <w'x'^2> = coef_wpxp2_implicit * <x'^2>
    !                + term_wpxp2_explicit, and
    !   <w'r_t'theta_l'> = coef_wprtpthlp_implicit * <r_t'theta_l'>
    !                       + term_wprtpthlp_explicit.
    !
    !   The coefficient state is frozen during the linear CLUBB advancement.
    !   A bounded responder speed prevents a nearly zero variance/covariance
    !   from producing an arbitrarily large implicit coefficient; the exact
    !   residual remains on the explicit side in that case.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels [#]
      ngrdcol    ! Number of grid columns [#]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wp2,       & ! Grid w variance [m2/s2]
      rtp2,      & ! Grid r_t variance [(kg/kg)2]
      thlp2,     & ! Grid theta_l variance [K2]
      rtpthlp,   & ! Grid r_t-theta_l covariance [(kg/kg) K]
      wp4,       & ! PDF-integrated w fourth moment [m4/s4]
      wprtp2,    & ! PDF-integrated w-r_t variance transport [m/s (kg/kg)2]
      wpthlp2,   & ! PDF-integrated w-theta_l variance transport [m/s K2]
      wprtpthlp    ! PDF-integrated mixed scalar transport [m/s (kg/kg) K]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      coef_wp4_implicit,      & ! <w4> / <w2>2 [-]
      coef_wprtp2_implicit,   & ! Implicit r_t variance responder [m/s]
      term_wprtp2_explicit,   & ! Explicit r_t variance residual [m/s (kg/kg)2]
      coef_wpthlp2_implicit,  & ! Implicit theta_l variance responder [m/s]
      term_wpthlp2_explicit,  & ! Explicit theta_l variance residual [m/s K2]
      coef_wprtpthlp_implicit,& ! Implicit r_t-theta_l responder [m/s]
      term_wprtpthlp_explicit   ! Explicit r_t-theta_l residual [m/s (kg/kg) K]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), parameter :: &
      covariance_relative_floor = 1.0e-6_core_rknd, & ! Relative covariance floor [-]
      responder_speed_factor = 4.0_core_rknd          ! Maximum responder speed / w stdev [-]

    real( kind = core_rknd ) :: &
      speed_cap,        & ! Local maximum responder speed [m/s]
      covariance_floor, & ! Minimum resolved r_t-theta_l covariance [(kg/kg) K]
      raw_coefficient     ! Unbounded mixed-covariance responder [m/s]

    integer :: &
      i, & ! Grid-column loop index [#]
      k    ! Vertical-level loop index [#]

    !----------------------------- Begin Code -----------------------------
    do k = 1, nz
      do i = 1, ngrdcol
        if ( wp2(i,k) > 0.0_core_rknd ) then
          coef_wp4_implicit(i,k) = max( wp4(i,k) / ( wp2(i,k) * wp2(i,k) ), &
                                        0.0_core_rknd )
          speed_cap = responder_speed_factor * sqrt( wp2(i,k) )
        else
          coef_wp4_implicit(i,k) = 0.0_core_rknd
          speed_cap = 0.0_core_rknd
        end if

        call bounded_variance_responder( rtp2(i,k), wprtp2(i,k), speed_cap, &
                                         coef_wprtp2_implicit(i,k), &
                                         term_wprtp2_explicit(i,k) )
        call bounded_variance_responder( thlp2(i,k), wpthlp2(i,k), speed_cap, &
                                         coef_wpthlp2_implicit(i,k), &
                                         term_wpthlp2_explicit(i,k) )

        covariance_floor = covariance_relative_floor * &
                           sqrt( max( rtp2(i,k), 0.0_core_rknd ) * &
                                 max( thlp2(i,k), 0.0_core_rknd ) )
        if ( abs( rtpthlp(i,k) ) > covariance_floor .and. speed_cap > 0.0_core_rknd ) then
          raw_coefficient = wprtpthlp(i,k) / rtpthlp(i,k)
          coef_wprtpthlp_implicit(i,k) = max( -speed_cap, min( speed_cap, raw_coefficient ) )
        else
          coef_wprtpthlp_implicit(i,k) = 0.0_core_rknd
        end if
        term_wprtpthlp_explicit(i,k) = wprtpthlp(i,k) &
          - coef_wprtpthlp_implicit(i,k) * rtpthlp(i,k)
      end do
    end do

  end subroutine diagnose_transport_implicit_responders

  !=============================================================================
  subroutine bounded_variance_responder( variance, transport, speed_cap, &
                                         coefficient, explicit_term )

    ! Description:
    !   Splits a scalar-variance transport moment into a bounded linear
    !   responder and an exact explicit residual at the frozen PDF state.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), intent(in) :: &
      variance,  & ! Scalar variance [units squared]
      transport, & ! Scalar-variance transport [m/s units squared]
      speed_cap    ! Maximum responder speed [m/s]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), intent(out) :: &
      coefficient,  & ! Bounded implicit responder [m/s]
      explicit_term   ! Explicit transport residual [m/s units squared]

    !----------------------------- Begin Code -----------------------------
    if ( variance > 0.0_core_rknd .and. speed_cap > 0.0_core_rknd ) then
      coefficient = max( -speed_cap, min( speed_cap, transport / variance ) )
    else
      coefficient = 0.0_core_rknd
    end if
    explicit_term = transport - coefficient * variance

  end subroutine bounded_variance_responder

  !=============================================================================
  subroutine diagnose_transport_pdf_state( mean, variance, covariance, skewness, &
                                           moisture_tail_gain, center_budget,    &
                                           w_direction_scale, g1_wrt_capture,   &
                                           plume_structure_strength,             &
                                           mean_1, mean_2, covariance_1,        &
                                           covariance_2, weight_1 )

    ! Description:
    !   Diagnoses one physical trivariate state in standardized coordinates.
    !   The component-center direction uses all three skewnesses, while the
    !   covariance contrasts are algebraically capped to keep both components
    !   positive definite.  In a coherent positive-w, positive-w-r_t transport
    !   regime, existing second moments can align G1 with either a warm or cool
    !   plume and request a simple, well-mixed G2 covariance structure.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3), intent(in) :: &
      mean,     & ! Grid mean state [mixed units]
      variance, & ! Grid marginal variances [mixed units squared]
      skewness    ! Grid marginal skewnesses [-]

    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance ! Grid covariance matrix [mixed units squared]

    real( kind = core_rknd ), intent(in) :: &
      moisture_tail_gain, & ! Maps r_t skewness to the G1 weight [-]
      center_budget,      & ! Covariance-metric center-separation budget [-]
      w_direction_scale,  & ! Relative w contribution to center direction [-]
      g1_wrt_capture,     & ! Positive residual w-r_t covariance in G1 [-]
      plume_structure_strength ! Coherent transport-plume structure strength [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3), intent(out) :: &
      mean_1, & ! G1 mean state [mixed units]
      mean_2    ! G2 mean state [mixed units]

    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      covariance_1, & ! G1 covariance matrix [mixed units squared]
      covariance_2    ! G2 covariance matrix [mixed units squared]

    real( kind = core_rknd ), intent(out) :: &
      weight_1 ! G1 mixture weight [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), parameter :: &
      variance_floor = sqrt( tiny( 1.0_core_rknd ) ), & ! Physical variance floor [units squared]
      tilt_tolerance = 1.0e-12_core_rknd, & ! Normalized covariance sign tolerance [-]
      psd_margin = 1.0e-8_core_rknd,   & ! Positive-definite interior margin [-]
      direction_floor = 1.0e-5_core_rknd, & ! Minimum resolved center direction [-]
      plume_gate_scale = 0.25_core_rknd   ! Correlation scale for plume gate [-]

    real( kind = core_rknd ), dimension(3) :: &
      stdev,           & ! Grid marginal standard deviations [mixed units]
      raw_skewness,    & ! Unbounded normalized third moments [-]
      capped_skewness, & ! Soft-capped normalized third moments [-]
      softened_skewness, & ! Bounded direction inputs [-]
      raw_direction,   & ! Unnormalized center direction [-]
      plume_direction, & ! Signed coherent-plume direction [-]
      direction,       & ! Covariance-metric normalized direction [-]
      displacement,    & ! Standardized G1-to-G2 center displacement [-]
      lower_bounds,    & ! Lower component-width-difference limits [-]
      upper_bounds       ! Upper component-width-difference limits [-]

    real( kind = core_rknd ), dimension(3,3) :: &
      correlation,      & ! Standardized grid covariance matrix [-]
      lower,            & ! Cholesky factor of correlation [-]
      residual,         & ! Covariance not used by center separation [-]
      base_contrast,    & ! Negative-tilt allocation contrast [-]
      covariance_1_base, & ! G1 covariance before requested contrast [-]
      covariance_2_base, & ! G2 covariance before requested contrast [-]
      contrast,         & ! Width and positive-tilt contrast [-]
      covariance_1_std, & ! G1 standardized covariance [-]
      covariance_2_std, & ! G2 standardized covariance [-]
      physical_scale      ! Converts normalized covariances to physical units

    real( kind = core_rknd ) :: &
      tail_gain,            & ! Bounded moisture-tail tuning [-]
      bounded_center_budget, & ! Bounded center-budget tuning [-]
      bounded_w_scale,      & ! Bounded w-direction tuning [-]
      bounded_wrt_capture,  & ! Bounded G1 w-r_t capture tuning [-]
      bounded_plume_structure, & ! Bounded coherent-plume tuning [-]
      moisture_signal,      & ! Magnitude of softened r_t skewness [-]
      w_signal,             & ! Magnitude of softened w skewness [-]
      separation_signal,    & ! Center-separation signal [-]
      tail_weight_signal,   & ! Bounded r_t-tail signal [-]
      w_tail_weight,        & ! w-skewness rare-component weight [-]
      plume_gate,           & ! Coherent transport-plume detection [-]
      plume_blend,          & ! Applied plume-structure fraction [-]
      weight_2,             & ! G2 mixture weight [-]
      metric_length,        & ! Covariance-metric length of raw direction [-]
      negative_tilt_scale,  & ! PSD cap for negative w-r_t allocation [-]
      contrast_scale,       & ! PSD cap for requested contrast [-]
      width_difference,     & ! Requested difference between component widths [-]
      width_margin            ! Margin retained in each marginal variance [-]

    integer :: &
      j ! Variable index [#]

    logical :: &
      valid_input ! True when the standardized covariance is positive definite

    !----------------------------- Begin Code -----------------------------
    if ( any( variance <= variance_floor ) ) then
      call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                          mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                          weight_1 )                                  ! Out
      return
    end if

    stdev = sqrt( variance )
    correlation = covariance / outer_product_3( stdev, stdev )
    correlation = 0.5_core_rknd * ( correlation + transpose( correlation ) )
    correlation(1,1) = 1.0_core_rknd
    correlation(2,2) = 1.0_core_rknd
    correlation(3,3) = 1.0_core_rknd

    call cholesky_3x3( correlation, lower, valid_input )
    if ( .not. valid_input ) then
      call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                          mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                          weight_1 )                                  ! Out
      return
    end if

    raw_skewness = skewness
    capped_skewness = min( max( raw_skewness, -12.0_core_rknd ), 12.0_core_rknd )
    softened_skewness = tanh( capped_skewness / 1.25_core_rknd )

    tail_gain = min( max( moisture_tail_gain, 0.0_core_rknd ), 3.0_core_rknd )
    bounded_center_budget = min( max( center_budget, 0.02_core_rknd ), 0.97_core_rknd )
    bounded_w_scale = min( max( w_direction_scale, 0.0_core_rknd ), 1.2_core_rknd )
    bounded_wrt_capture = min( max( g1_wrt_capture, 0.0_core_rknd ), 1.0_core_rknd )
    bounded_plume_structure = min( max( plume_structure_strength, 0.0_core_rknd ), &
                                        1.0_core_rknd )

    moisture_signal = abs( softened_skewness(2) )
    w_signal = abs( softened_skewness(1) )
    tail_weight_signal = tanh( tail_gain * abs( capped_skewness(2) ) / 1.25_core_rknd )
    weight_1 = min( max( 0.5_core_rknd * ( 1.0_core_rknd - tail_weight_signal ), &
                         0.035_core_rknd ), 0.5_core_rknd )

    ! Identify a coherent transport plume from positive w skewness and a
    ! positive w-r_t relation.  The signed w-theta_l relation is deliberately
    ! not gated: fresh warm/moist and mature moist/cool plumes are both valid
    ! members of the same transport family.
    plume_gate = max( softened_skewness(1), 0.0_core_rknd ) &
                 * max( correlation(1,2), 0.0_core_rknd )
    plume_gate = min( plume_gate / plume_gate_scale**2, 1.0_core_rknd )
    plume_blend = bounded_plume_structure * plume_gate

    ! Blend toward the two-point w-skewness mixture weight only inside the
    ! coherent transport-plume regime.  At zero blend, this is exactly the
    ! moisture-tail weight.
    w_tail_weight = 0.5_core_rknd * ( 1.0_core_rknd - abs( capped_skewness(1) ) &
                    / sqrt( 4.0_core_rknd + capped_skewness(1)**2 ) )
    w_tail_weight = min( max( w_tail_weight, 0.035_core_rknd ), 0.5_core_rknd )
    weight_1 = ( 1.0_core_rknd - plume_blend ) * weight_1 + plume_blend * w_tail_weight
    weight_2 = 1.0_core_rknd - weight_1

    raw_direction = [ bounded_w_scale * softened_skewness(1), &
                      softened_skewness(2), softened_skewness(3) ]
    plume_direction = [ 1.0_core_rknd, correlation(1,2), correlation(1,3) ]
    raw_direction = ( 1.0_core_rknd - plume_blend ) * raw_direction &
                    + plume_blend * plume_direction
    direction = 0.0_core_rknd
    if ( sqrt( dot_product( raw_direction, raw_direction ) ) > direction_floor ) then
      call solve_spd_3x3( lower, raw_direction, direction, valid_input )
      if ( .not. valid_input ) then
        call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                            mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                            weight_1 )                                  ! Out
        return
      end if
      metric_length = sqrt( max( dot_product( raw_direction, direction ), psd_margin ) )
      direction = raw_direction / metric_length
    end if

    separation_signal = ( 1.0_core_rknd - plume_blend ) * moisture_signal &
                        + plume_blend * max( moisture_signal, w_signal )
    displacement = bounded_center_budget * separation_signal &
                   / sqrt( weight_1 * weight_2 ) * direction
    residual = correlation - weight_1 * weight_2 &
               * outer_product_3( displacement, displacement )
    residual = 0.5_core_rknd * ( residual + transpose( residual ) )
    call cholesky_3x3( residual, lower, valid_input )
    if ( .not. valid_input ) then
      call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                          mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                          weight_1 )                                  ! Out
      return
    end if

    base_contrast = 0.0_core_rknd
    negative_tilt_scale = 1.0_core_rknd
    ! A negative grid w-r_t covariance is allocated out of G1 whenever that
    ! strict choice remains PSD.  The analytic cap softens only the rare
    ! incompatible cases.
    if ( correlation(1,2) < -tilt_tolerance ) then
      base_contrast(1,2) = -residual(1,2)
      base_contrast(2,1) = -residual(1,2)
      call maximum_contrast_scale( residual, residual, base_contrast, & ! In
                                   weight_1 / weight_2, negative_tilt_scale, valid_input ) ! Out
      if ( .not. valid_input ) then
        call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                            mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                            weight_1 )                                  ! Out
        return
      end if
    end if
    covariance_1_base = residual + negative_tilt_scale * base_contrast
    covariance_2_base = residual - weight_1 / weight_2 * negative_tilt_scale * base_contrast
    covariance_1_base = 0.5_core_rknd * ( covariance_1_base + transpose( covariance_1_base ) )
    covariance_2_base = 0.5_core_rknd * ( covariance_2_base + transpose( covariance_2_base ) )

    contrast = 0.0_core_rknd
    do j = 1, 3
      if ( abs( displacement(j) ) > direction_floor ) then
        width_difference = ( capped_skewness(j) / ( weight_1 * weight_2 ) &
                             - ( weight_2 - weight_1 ) * displacement(j)**3 ) &
                           / ( 3.0_core_rknd * displacement(j) )
        width_margin = min( 1.0e-6_core_rknd, 0.01_core_rknd * residual(j,j) )
        lower_bounds(j) = ( -residual(j,j) + width_margin ) / weight_2
        upper_bounds(j) = ( residual(j,j) - width_margin ) / weight_1
        width_difference = min( max( width_difference, lower_bounds(j) ), upper_bounds(j) )
        contrast(j,j) = weight_2 * width_difference
      end if
    end do

    if ( correlation(1,2) > tilt_tolerance .and. residual(1,2) > tilt_tolerance &
         .and. max( bounded_wrt_capture, plume_blend ) > 0.0_core_rknd ) then
      contrast(1,2) = max( bounded_wrt_capture, plume_blend ) * weight_2 / weight_1 &
                      * residual(1,2)
      contrast(2,1) = contrast(1,2)
    end if

    ! Keep G2 close to a well-mixed background whenever a coherent transport
    ! plume is detected.  The requested contrast removes G2 w-theta_l tilt
    ! toward zero and removes only *positive* G2 r_t-theta_l tilt.  The equal
    ! and opposite weighted G1 contrast preserves the supplied grid covariance.
    ! maximum_contrast_scale below applies one joint PSD cap to these requests
    ! and to the marginal width contrasts.
    contrast(1,3) = plume_blend * weight_2 / weight_1 * covariance_2_base(1,3)
    contrast(3,1) = contrast(1,3)
    contrast(2,3) = plume_blend * weight_2 / weight_1 &
                    * max( covariance_2_base(2,3), 0.0_core_rknd )
    contrast(3,2) = contrast(2,3)

    call maximum_contrast_scale( covariance_1_base, covariance_2_base, contrast, & ! In
                                 weight_1 / weight_2, contrast_scale, valid_input )       ! Out
    if ( .not. valid_input ) then
      call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                          mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                          weight_1 )                                  ! Out
      return
    end if
    covariance_1_std = covariance_1_base + contrast_scale * contrast
    covariance_2_std = covariance_2_base - weight_1 / weight_2 * contrast_scale * contrast
    covariance_1_std = 0.5_core_rknd * ( covariance_1_std + transpose( covariance_1_std ) )
    covariance_2_std = 0.5_core_rknd * ( covariance_2_std + transpose( covariance_2_std ) )
    call cholesky_3x3( covariance_1_std, lower, valid_input )
    if ( valid_input ) call cholesky_3x3( covariance_2_std, lower, valid_input )
    if ( .not. valid_input ) then
      call collapsed_transport_pdf_state( mean, variance, covariance, & ! In
                                          mean_1, mean_2, covariance_1, covariance_2, & ! Out
                                          weight_1 )                                  ! Out
      return
    end if

    physical_scale = outer_product_3( stdev, stdev )
    mean_1 = mean + weight_2 * displacement * stdev
    mean_2 = mean - weight_1 * displacement * stdev
    covariance_1 = covariance_1_std * physical_scale
    covariance_2 = covariance_2_std * physical_scale

  end subroutine diagnose_transport_pdf_state

  !=============================================================================
  subroutine collapsed_transport_pdf_state( mean, variance, covariance, &
                                            mean_1, mean_2, covariance_1, &
                                            covariance_2, weight_1 )

    ! Description:
    !   Returns two identical component copies of a valid input covariance when
    !   the transport diagnosis cannot safely proceed.  A diagonal covariance is
    !   used only when the supplied covariance itself is not positive definite.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3), intent(in) :: &
      mean,     & ! Grid mean state [mixed units]
      variance    ! Grid marginal variances [mixed units squared]

    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance ! Grid covariance matrix [mixed units squared]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3), intent(out) :: &
      mean_1, & ! G1 mean state [mixed units]
      mean_2    ! G2 mean state [mixed units]

    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      covariance_1, & ! G1 covariance matrix [mixed units squared]
      covariance_2    ! G2 covariance matrix [mixed units squared]

    real( kind = core_rknd ), intent(out) :: &
      weight_1 ! G1 mixture weight [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(3,3) :: &
      lower ! Cholesky factor used to check the input covariance [mixed units]

    logical :: &
      valid_input ! True when the supplied covariance is positive definite

    integer :: &
      j ! Variable index [#]

    !----------------------------- Begin Code -----------------------------
    mean_1 = mean
    mean_2 = mean
    covariance_1 = 0.5_core_rknd * ( covariance + transpose( covariance ) )
    call cholesky_3x3( covariance_1, lower, valid_input )
    if ( .not. valid_input ) then
      covariance_1 = 0.0_core_rknd
      do j = 1, 3
        covariance_1(j,j) = max( variance(j), sqrt( tiny( 1.0_core_rknd ) ) )
      end do
    end if
    covariance_2 = covariance_1
    weight_1 = 0.5_core_rknd

  end subroutine collapsed_transport_pdf_state

  !=============================================================================
  subroutine maximum_contrast_scale( covariance_1, covariance_2, contrast, &
                                     covariance_2_ratio, scale, valid_input )

    ! Description:
    !   Finds the direct PSD-safe scale for C1+tA and
    !   C2-(covariance_2_ratio)tA by whitening each component covariance and
    !   limiting t by its most-negative eigenvalue.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance_1, & ! G1 base covariance [-]
      covariance_2, & ! G2 base covariance [-]
      contrast        ! Requested G1 covariance contrast [-]

    real( kind = core_rknd ), intent(in) :: &
      covariance_2_ratio ! G1-to-G2 mixture-weight ratio [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), intent(out) :: &
      scale ! PSD-safe requested contrast scale [-]

    logical, intent(out) :: &
      valid_input ! True when both base covariances are positive definite

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), parameter :: &
      psd_margin = 1.0e-8_core_rknd ! Positive-definite interior margin [-]

    real( kind = core_rknd ), dimension(3,3) :: &
      lower,      & ! Cholesky factor of one base covariance [-]
      direction,  & ! Contrast direction for one component [-]
      whitened      ! Whitened covariance contrast [-]

    real( kind = core_rknd ), dimension(3) :: &
      eigenvalues ! Eigenvalues of whitened contrast [-]

    real( kind = core_rknd ) :: &
      limiting_eigenvalue, & ! Most negative whitening eigenvalue [-]
      limiting_scale         ! Candidate PSD-safe scale [-]

    integer :: &
      component ! Component loop index [#]

    !----------------------------- Begin Code -----------------------------
    if ( maxval( abs( contrast ) ) <= 1.0e-15_core_rknd ) then
      scale = 1.0_core_rknd
      valid_input = .true.
      return
    end if

    scale = huge( scale )
    valid_input = .true.
    do component = 1, 2
      if ( component == 1 ) then
        call cholesky_3x3( covariance_1, lower, valid_input )
        direction = contrast
      else
        call cholesky_3x3( covariance_2, lower, valid_input )
        direction = -covariance_2_ratio * contrast
      end if
      if ( .not. valid_input ) return

      call whiten_symmetric_3x3( lower, direction, whitened )
      call symmetric_eigenvalues_3x3( whitened, eigenvalues )
      limiting_eigenvalue = minval( eigenvalues )
      if ( limiting_eigenvalue < 0.0_core_rknd ) then
        limiting_scale = ( 1.0_core_rknd - psd_margin ) / ( -limiting_eigenvalue )
        scale = min( scale, limiting_scale )
      end if
    end do

    if ( scale >= 1.0_core_rknd ) then
      scale = 1.0_core_rknd
    else
      scale = max( 0.0_core_rknd, 0.995_core_rknd * scale )
    end if

  end subroutine maximum_contrast_scale

  !=============================================================================
  subroutine cholesky_3x3( matrix, lower, valid_input )

    ! Description:
    !   Computes the lower Cholesky factor of a symmetric 3-by-3 matrix and
    !   reports whether the matrix is safely positive definite.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      matrix ! Symmetric matrix to factor [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      lower ! Lower Cholesky factor [-]

    logical, intent(out) :: &
      valid_input ! True when the matrix is positive definite

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), parameter :: &
      psd_margin = 1.0e-12_core_rknd ! Positive pivot floor [-]

    real( kind = core_rknd ) :: &
      pivot ! Current Cholesky pivot [-]

    !----------------------------- Begin Code -----------------------------
    lower = 0.0_core_rknd
    valid_input = .false.

    pivot = matrix(1,1)
    if ( pivot <= psd_margin ) return
    lower(1,1) = sqrt( pivot )

    lower(2,1) = matrix(2,1) / lower(1,1)
    lower(3,1) = matrix(3,1) / lower(1,1)

    pivot = matrix(2,2) - lower(2,1)**2
    if ( pivot <= psd_margin ) return
    lower(2,2) = sqrt( pivot )

    lower(3,2) = ( matrix(3,2) - lower(3,1) * lower(2,1) ) / lower(2,2)

    pivot = matrix(3,3) - lower(3,1)**2 - lower(3,2)**2
    if ( pivot <= psd_margin ) return
    lower(3,3) = sqrt( pivot )

    valid_input = .true.

  end subroutine cholesky_3x3

  !=============================================================================
  subroutine solve_spd_3x3( lower, rhs, solution, valid_input )

    ! Description:
    !   Solves L L^T x = rhs for the Cholesky factor L returned by
    !   cholesky_3x3.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      lower ! Lower Cholesky factor [-]

    real( kind = core_rknd ), dimension(3), intent(in) :: &
      rhs ! Right-hand side [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3), intent(out) :: &
      solution ! SPD-system solution [-]

    logical, intent(out) :: &
      valid_input ! True when the factor has nonzero diagonal entries

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(3) :: &
      intermediate ! Forward-substitution result [-]

    !----------------------------- Begin Code -----------------------------
    valid_input = lower(1,1) > 1.0e-12_core_rknd .and. &
                  lower(2,2) > 1.0e-12_core_rknd .and. &
                  lower(3,3) > 1.0e-12_core_rknd
    if ( .not. valid_input ) then
      solution = 0.0_core_rknd
      return
    end if

    intermediate(1) = rhs(1) / lower(1,1)
    intermediate(2) = ( rhs(2) - lower(2,1) * intermediate(1) ) / lower(2,2)
    intermediate(3) = ( rhs(3) - lower(3,1) * intermediate(1) &
                        - lower(3,2) * intermediate(2) ) / lower(3,3)

    solution(3) = intermediate(3) / lower(3,3)
    solution(2) = ( intermediate(2) - lower(3,2) * solution(3) ) / lower(2,2)
    solution(1) = ( intermediate(1) - lower(2,1) * solution(2) &
                    - lower(3,1) * solution(3) ) / lower(1,1)

  end subroutine solve_spd_3x3

  !=============================================================================
  subroutine whiten_symmetric_3x3( lower, matrix, whitened )

    ! Description:
    !   Forms L^{-1} A L^{-T} for a symmetric matrix A and its lower Cholesky
    !   factor L.  This supplies the eigenvalue test used by the PSD cap.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      lower,  & ! Lower Cholesky factor [-]
      matrix    ! Symmetric matrix to whiten [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      whitened ! Whitened symmetric matrix [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(3,3) :: &
      intermediate ! Intermediate L^{-1} A matrix [-]

    real( kind = core_rknd ), dimension(3) :: &
      rhs,       & ! One row/column used by triangular solve [-]
      solution     ! Triangular-solve result [-]

    integer :: &
      j ! Matrix-column/row index [#]

    !----------------------------- Begin Code -----------------------------
    do j = 1, 3
      rhs = matrix(:,j)
      solution(1) = rhs(1) / lower(1,1)
      solution(2) = ( rhs(2) - lower(2,1) * solution(1) ) / lower(2,2)
      solution(3) = ( rhs(3) - lower(3,1) * solution(1) &
                      - lower(3,2) * solution(2) ) / lower(3,3)
      intermediate(:,j) = solution
    end do

    do j = 1, 3
      rhs = intermediate(j,:)
      solution(1) = rhs(1) / lower(1,1)
      solution(2) = ( rhs(2) - lower(2,1) * solution(1) ) / lower(2,2)
      solution(3) = ( rhs(3) - lower(3,1) * solution(1) &
                      - lower(3,2) * solution(2) ) / lower(3,3)
      whitened(j,:) = solution
    end do
    whitened = 0.5_core_rknd * ( whitened + transpose( whitened ) )

  end subroutine whiten_symmetric_3x3

  !=============================================================================
  subroutine symmetric_eigenvalues_3x3( matrix, eigenvalues )

    ! Description:
    !   Computes the eigenvalues of a real symmetric 3-by-3 matrix with the
    !   closed-form trigonometric formula.  No iterative eigensolver is used.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      matrix ! Real symmetric matrix [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3), intent(out) :: &
      eigenvalues ! Matrix eigenvalues [-]

    !---------------------------- Local Variables --------------------------
    real( kind = core_rknd ), dimension(3,3) :: &
      shifted ! Trace-free, scaled matrix [-]

    real( kind = core_rknd ) :: &
      off_diagonal_norm, & ! Sum of squared off-diagonal entries [-]
      trace_third,       & ! One third of matrix trace [-]
      scale_squared,     & ! Squared eigensystem scaling factor [-]
      scale,             & ! Eigensystem scaling factor [-]
      determinant_half,  & ! Half determinant of the scaled matrix [-]
      angle,             & ! Cubic-eigenvalue trigonometric angle [radians]
      pi                   ! Pi [radians]

    !----------------------------- Begin Code -----------------------------
    off_diagonal_norm = matrix(1,2)**2 + matrix(1,3)**2 + matrix(2,3)**2
    if ( off_diagonal_norm <= 1.0e-30_core_rknd ) then
      eigenvalues = [ matrix(1,1), matrix(2,2), matrix(3,3) ]
      return
    end if

    trace_third = ( matrix(1,1) + matrix(2,2) + matrix(3,3) ) / 3.0_core_rknd
    scale_squared = ( matrix(1,1) - trace_third )**2 &
                  + ( matrix(2,2) - trace_third )**2 &
                  + ( matrix(3,3) - trace_third )**2 &
                  + 2.0_core_rknd * off_diagonal_norm
    scale = sqrt( scale_squared / 6.0_core_rknd )
    if ( scale <= 1.0e-20_core_rknd ) then
      eigenvalues = [ trace_third, trace_third, trace_third ]
      return
    end if

    shifted = matrix / scale
    shifted(1,1) = shifted(1,1) - trace_third / scale
    shifted(2,2) = shifted(2,2) - trace_third / scale
    shifted(3,3) = shifted(3,3) - trace_third / scale
    determinant_half = 0.5_core_rknd * determinant_3x3( shifted )
    determinant_half = min( max( determinant_half, -1.0_core_rknd ), 1.0_core_rknd )
    angle = acos( determinant_half ) / 3.0_core_rknd
    pi = acos( -1.0_core_rknd )

    eigenvalues(1) = trace_third + 2.0_core_rknd * scale * cos( angle )
    eigenvalues(3) = trace_third + 2.0_core_rknd * scale &
                     * cos( angle + 2.0_core_rknd * pi / 3.0_core_rknd )
    eigenvalues(2) = 3.0_core_rknd * trace_third - eigenvalues(1) - eigenvalues(3)

  end subroutine symmetric_eigenvalues_3x3

  !=============================================================================
  function outer_product_3( first, second ) result( product )

    ! Description:
    !   Returns the outer product of two three-element vectors.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3), intent(in) :: &
      first,  & ! First vector [-]
      second    ! Second vector [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(3,3) :: &
      product ! Outer product [-]

    integer :: &
      j ! Matrix-column index [#]

    !----------------------------- Begin Code -----------------------------
    do j = 1, 3
      product(:,j) = first * second(j)
    end do

  end function outer_product_3

  !=============================================================================
  function determinant_3x3( matrix ) result( determinant )

    ! Description:
    !   Returns the determinant of a 3-by-3 matrix by direct expansion.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      matrix ! Matrix whose determinant is needed [-]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ) :: &
      determinant ! Matrix determinant [-]

    !----------------------------- Begin Code -----------------------------
    determinant = matrix(1,1) * ( matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2) ) &
                - matrix(1,2) * ( matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1) ) &
                + matrix(1,3) * ( matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1) )

  end function determinant_3x3

  !=============================================================================
  function correlation_from_covariance( covariance, first, second ) result( correlation )

    ! Description:
    !   Converts one covariance entry to a bounded correlation, returning zero
    !   if either component variance is too small to define a correlation.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      covariance ! Component covariance matrix [mixed units squared]

    integer, intent(in) :: &
      first,  & ! First covariance-variable index [#]
      second    ! Second covariance-variable index [#]

    !--------------------------- Output Variables --------------------------
    real( kind = core_rknd ) :: &
      correlation ! Bounded correlation coefficient [-]

    !----------------------------- Begin Code -----------------------------
    if ( covariance(first,first) > sqrt( tiny( 1.0_core_rknd ) ) .and. &
         covariance(second,second) > sqrt( tiny( 1.0_core_rknd ) ) ) then
      correlation = covariance(first,second) &
                  / ( sqrt( covariance(first,first) ) * sqrt( covariance(second,second) ) )
      correlation = min( max( correlation, -0.999999_core_rknd ), 0.999999_core_rknd )
    else
      correlation = 0.0_core_rknd
    end if

  end function correlation_from_covariance

end module trivariate_transport_pdf
