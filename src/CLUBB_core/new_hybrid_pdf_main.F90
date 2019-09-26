! $Id$
!===============================================================================
module new_hybrid_pdf_main

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! multivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), liquid water potential temperature (thl), and optionally,
  ! the west-east horizontal wind component (u), the south-north horizontal wind
  ! component (v), and passive scalars (sclr).

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: new_hybrid_pdf_driver    ! Procedure(s)

  private :: calc_responder_var, & ! Procedure(s)
             calc_F_w_zeta_w

  private

  contains

  !=============================================================================
  subroutine new_hybrid_pdf_driver( wm, rtm, thlm, um, vm, wp2, rtp2,   & ! In
                                    thlp2, up2, vp2, Skw, wprtp,        & ! In
                                    wpthlp, rtpthlp, upwp, vpwp,        & ! In
                                    sclrm, sclrp2, wpsclrp,             & ! In
                                    Skrt, Skthl, Sku, Skv, Sksclr,      & ! I/O
                                    mu_w_1, mu_w_2,                     & ! Out
                                    mu_rt_1, mu_rt_2,                   & ! Out
                                    mu_thl_1, mu_thl_2,                 & ! Out
                                    mu_u_1, mu_u_2, mu_v_1, mu_v_2,     & ! Out
                                    sigma_w_1_sqd, sigma_w_2_sqd,       & ! Out
                                    sigma_rt_1_sqd, sigma_rt_2_sqd,     & ! Out
                                    sigma_thl_1_sqd, sigma_thl_2_sqd,   & ! Out
                                    sigma_u_1_sqd, sigma_u_2_sqd,       & ! Out
                                    sigma_v_1_sqd, sigma_v_2_sqd,       & ! Out
                                    mu_sclr_1, mu_sclr_2,               & ! Out
                                    sigma_sclr_1_sqd, sigma_sclr_2_sqd, & ! Out
                                    mixt_frac,                          & ! Out
                                    pdf_implicit_coefs_terms,           & ! Out
                                    F_w, min_F_w, max_F_w               ) ! Out
                             

    ! Description:
    ! Calculate the PDF parameters for w (including mixture fraction), rt,
    ! theta-l, and optionally, u, v, and passive scalar variables.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        zero,    & ! Constant(s)
        fstderr

    use new_hybrid_pdf, only: &
        calculate_w_params,         & ! Procedure(s)
        calc_coef_wp4_implicit,     &
        calc_coef_wp2xp_implicit,   &
        calc_coefs_wpxp2_semiimpl,  &
        calc_coefs_wpxpyp_semiimpl

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use model_flags, only: &
        l_explicit_turbulent_adv_wp3,  & ! Variable(s)
        l_explicit_turbulent_adv_wpxp, &
        l_explicit_turbulent_adv_xpyp

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wm,      & ! Mean of w (overall)                 [m/s]
      rtm,     & ! Mean of rt (overall)                [kg/kg]
      thlm,    & ! Mean of thl (overall)               [K]
      wp2,     & ! Variance of w (overall)             [m^2/s^2]
      rtp2,    & ! Variance of rt (overall)            [kg^2/kg^2]
      thlp2,   & ! Variance of thl (overall)           [K^2]
      Skw,     & ! Skewness of w (overall)             [-]
      wprtp,   & ! Covariance of w and rt (overall)    [(m/s)kg/kg]
      wpthlp,  & ! Covariance of w and thl (overall)   [(m/s)K]
      rtpthlp    ! Covariance of rt and thl (overall)  [(kg/kg)K]

    ! Input/Output Variables
    ! These variables are input/output because their values may be clipped.
    ! Otherwise, as long as it is not necessary to clip them, their values
    ! will stay the same.
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Skrt,  & ! Skewness of rt (overall)            [-]
      Skthl    ! Skewness of thl (overall)           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mu_w_1,          & ! Mean of w (1st PDF component)        [m/s]
      mu_w_2,          & ! Mean of w (2nd PDF component)        [m/s]
      mu_rt_1,         & ! Mean of rt (1st PDF component)       [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)       [kg/kg]
      mu_thl_1,        & ! Mean of thl (1st PDF component)      [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)      [K]
      sigma_w_1_sqd,   & ! Variance of w (1st PDF component)    [m^2/s^2]
      sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)    [m^2/s^2]
      sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)   [kg^2/kg^2]
      sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)   [kg^2/kg^2]
      sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)  [K^2]
      sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)  [K^2]
      mixt_frac          ! Mixture fraction                     [-]

    type(implicit_coefs_terms), dimension(gr%nz), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Output only for recording statistics.
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      F_w,     & ! Parameter for the spread of the PDF component means of w  [-]
      min_F_w, & ! Minimum allowable value of parameter F_w                  [-]
      max_F_w    ! Maximum allowable value of parameter F_w                  [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      sigma_w_1, & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2    ! Standard deviation of w (2nd PDF component)      [m/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_sigma_w_1_sqd,   & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
      coef_sigma_w_2_sqd,   & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
      coef_sigma_rt_1_sqd,  & ! sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
      coef_sigma_rt_2_sqd,  & ! sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
      coef_sigma_thl_1_sqd, & ! sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
      coef_sigma_thl_2_sqd    ! sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      zeta_w    ! Parameter for the PDF component variances of w           [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp4_implicit,     & ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]
      coef_wprtp2_implicit,  & ! <w'rt'^2> = coef_wprtp2_implicit*<rt'^2>  [m/s]
      coef_wpthlp2_implicit    ! <w'thl'^2>=coef_wpthlp2_implicit*<thl'^2> [m/s]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      max_corr_w_sclr_sqd    ! Max value of wpsclrp^2 / ( wp2 * sclrp2 )     [-]

    integer :: k, sclr_indx  ! Loop indices


    ! Calculate the maximum value of the square of the correlation of w and a
    ! scalar when scalars are used.
    ! Initialize max_corr_w_sclr_sqd to 0.  It needs to retain this value even
    ! when sclr_dim = 0.
    max_corr_w_sclr_sqd = zero
    if ( sclr_dim > 0 ) then
       do k = 1, gr%nz, 1
          do sclr_indx = 1, sclr_dim, 1
             max_corr_w_sclr_sqd(k) = max( wpsclrp(k,sclr_indx)**2 &
                                           / ( wp2(k) * sclrp2(k,sclr_indx) ), &
                                           max_corr_w_sclr_sqd(k) )
          enddo ! sclr_indx = 1, sclr_dim, 1
       enddo ! k = 1, gr%nz, 1
    endif ! sclr_dim > 0

    ! Calculate the values of PDF tunable parameters F_w and zeta_w.
    ! Vertical velocity, w, will always be the setter variable.
    call calc_F_w_zeta_w( Skw, wprtp, wpthlp, upwp, vpwp, & ! In
                          wp2, rtp2, thlp2, up2, vp2,     & ! In
                          max_corr_w_sclr_sqd,            & ! In
                          F_w, zeta_w, min_F_w, max_F_w   ) ! Out

    ! Calculate the PDF parameters, including mixture fraction, for the
    ! setter variable, w.
    call calculate_w_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                             mu_w_1, mu_w_2, sigma_w_1, & ! Out
                             sigma_w_2, mixt_frac,      & ! Out
                             coef_sigma_w_1_sqd,        & ! Out
                             coef_sigma_w_2_sqd         ) ! Out

    sigma_w_1_sqd = sigma_w_1**2
    sigma_w_2_sqd = sigma_w_2**2

    if ( any( mixt_frac < zero ) ) then
       write(fstderr,*) "Mixture fraction cannot be calculated."
       write(fstderr,*) "The value of F_w must be greater than 0 when " &
                        // "| Skw | > 0."
       stop
    endif

    ! Calculate the PDF parameters for responder variable rt.
    call calc_responder_var( rtm, rtp2, sgn_wprtp, mixt_frac, & ! In
                             coef_spread_DG_means_rt,         & ! In
                             exp_factor_rt,                   & ! In
                             max_Skx2_pos_Skx_sgn_wpxp,       & ! In
                             max_Skx2_neg_Skx_sgn_wpxp,       & ! In
                             Skrt,                            & ! In/Out
                             mu_rt_1, mu_rt_2,                & ! Out
                             sigma_rt_1_sqd, sigma_rt_2_sqd,  & ! Out
                             coef_sigma_rt_1_sqd,             & ! Out
                             coef_sigma_rt_2_sqd,             & ! Out
                             F_rt, min_F_rt, max_F_rt         ) ! Out

    ! Calculate the PDF parameters for responder variable thl.
    call calc_responder_var( thlm, thlp2, sgn_wpthlp, mixt_frac, & ! In
                             coef_spread_DG_means_thl,           & ! In
                             exp_factor_thl,                     & ! In
                             max_Skx2_pos_Skx_sgn_wpxp,          & ! In
                             max_Skx2_neg_Skx_sgn_wpxp,          & ! In
                             Skthl,                              & ! In/Out
                             mu_thl_1, mu_thl_2,                 & ! Out
                             sigma_thl_1_sqd, sigma_thl_2_sqd,   & ! Out
                             coef_sigma_thl_1_sqd,               & ! Out
                             coef_sigma_thl_2_sqd,               & ! Out
                             F_thl, min_F_thl, max_F_thl         ) ! Out


    if ( .not. l_explicit_turbulent_adv_wp3 ) then

       ! Turbulent advection of <w'^3> is being handled implicitly.

       ! <w'^4> = coef_wp4_implicit * <w'^2>^2.
       coef_wp4_implicit &
       = calc_coef_wp4_implicit( mixt_frac, F_w, &
                                 coef_sigma_w_1_sqd, &
                                 coef_sigma_w_2_sqd )

    else ! l_explicit_turbulent_adv_wp3

       ! Turbulent advection of <w'^3> is being handled explicitly.
       coef_wp4_implicit = zero

    endif ! .not. l_explicit_turbulent_adv_wp3

    if ( .not. l_explicit_turbulent_adv_xpyp ) then

       ! Turbulent advection of <rt'^2> and <thl'^2> is being handled
       ! implicitly.  Turbulent advection of <rt'thl'> is being handled
       ! semi-implicitly.

       ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>
       coef_wprtp2_implicit &
       = calc_coef_wpxp2_implicit( wp2, rtp2, wprtp, sgn_wprtp, &
                                   mixt_frac, F_w, F_rt, &
                                   coef_sigma_w_1_sqd, &
                                   coef_sigma_w_2_sqd, &
                                   coef_sigma_rt_1_sqd, &
                                   coef_sigma_rt_2_sqd  )

       ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>
       coef_wpthlp2_implicit &
       = calc_coef_wpxp2_implicit( wp2, thlp2, wpthlp, sgn_wpthlp, &
                                   mixt_frac, F_w, F_thl, &
                                   coef_sigma_w_1_sqd, &
                                   coef_sigma_w_2_sqd, &
                                   coef_sigma_thl_1_sqd, &
                                   coef_sigma_thl_2_sqd  )

       ! <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
       !               + term_wprtpthlp_explicit
       call calc_coefs_wpxpyp_semiimpl( wp2, rtp2, thlp2, wprtp,       & ! In
                                        wpthlp, sgn_wprtp, sgn_wpthlp, & ! In
                                        mixt_frac, F_w, F_rt, F_thl,   & ! In
                                        coef_sigma_w_1_sqd  ,          & ! In
                                        coef_sigma_w_2_sqd,            & ! In
                                        coef_sigma_rt_1_sqd,           & ! In
                                        coef_sigma_rt_2_sqd,           & ! In
                                        coef_sigma_thl_1_sqd,          & ! In
                                        coef_sigma_thl_2_sqd,          & ! In
                                        coef_wprtpthlp_implicit,       & ! Out
                                        term_wprtpthlp_explicit        ) ! Out

    else ! l_explicit_turbulent_adv_xpyp

       ! Turbulent advection of <rt'^2>, <thl'^2>, and <rt'thl'> is being
       ! handled explicitly.
       coef_wprtp2_implicit = zero
       coef_wpthlp2_implicit = zero
       coef_wprtpthlp_implicit = zero
       term_wprtpthlp_explicit = zero

    endif ! .not. l_explicit_turbulent_adv_xpyp

    if ( .not. l_explicit_turbulent_adv_wpxp ) then

       ! Turbulent advection of <w'rt'> and <w'thl'> is being handled
       ! semi-implicitly.

       ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
       call calc_coefs_wp2xp_semiimpl( wp2, rtp2, sgn_wprtp, & ! In
                                       mixt_frac, F_w, F_rt, & ! In
                                       coef_sigma_w_1_sqd,   & ! In
                                       coef_sigma_w_2_sqd,   & ! In
                                       coef_sigma_rt_1_sqd,  & ! In
                                       coef_sigma_rt_2_sqd,  & ! In
                                       coef_wp2rtp_implicit, & ! Out
                                       term_wp2rtp_explicit  ) ! Out

       ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
       call calc_coefs_wp2xp_semiimpl( wp2, thlp2, sgn_wpthlp, & ! In
                                       mixt_frac, F_w, F_thl,  & ! In
                                       coef_sigma_w_1_sqd,     & ! In
                                       coef_sigma_w_2_sqd,     & ! In
                                       coef_sigma_thl_1_sqd,   & ! In
                                       coef_sigma_thl_2_sqd,   & ! In
                                       coef_wp2thlp_implicit,  & ! Out
                                       term_wp2thlp_explicit   ) ! Out

    else ! l_explicit_turbulent_adv_wpxp

       ! Turbulent advection of <w'rt'> and <w'thl'> is being handled
       ! explicitly.
       coef_wp2rtp_implicit = zero
       term_wp2rtp_explicit = zero
       coef_wp2thlp_implicit = zero
       term_wp2thlp_explicit = zero

    endif ! .not. l_explicit_turbulent_adv_wpxp

    ! Pack the implicit coefficients and explicit terms into a single type
    ! variable for output.
    pdf_implicit_coefs_terms%coef_wp4_implicit = coef_wp4_implicit
    pdf_implicit_coefs_terms%coef_wprtp2_implicit = coef_wprtp2_implicit
    pdf_implicit_coefs_terms%coef_wpthlp2_implicit = coef_wpthlp2_implicit
    pdf_implicit_coefs_terms%coef_wp2rtp_implicit = coef_wp2rtp_implicit
    pdf_implicit_coefs_terms%term_wp2rtp_explicit = term_wp2rtp_explicit
    pdf_implicit_coefs_terms%coef_wp2thlp_implicit = coef_wp2thlp_implicit
    pdf_implicit_coefs_terms%term_wp2thlp_explicit = term_wp2thlp_explicit
    pdf_implicit_coefs_terms%coef_wprtpthlp_implicit = coef_wprtpthlp_implicit
    pdf_implicit_coefs_terms%term_wprtpthlp_explicit = term_wprtpthlp_explicit


    return

  end subroutine new_pdf_driver

  !=============================================================================
  elemental subroutine calc_responder_var( xm, xp2, wpxp, wp2, & ! In
                                           mixt_frac, F_w,     & ! In
                                           Skx,                & ! In/Out
                                           mu_x_1, mu_x_2,     & ! Out
                                           sigma_x_1_sqd,      & ! Out
                                           sigma_x_2_sqd,      & ! Out
                                           coef_sigma_x_1_sqd, & ! Out
                                           coef_sigma_x_2_sqd  ) ! Out

    ! Description:
    ! This is the sub-driver for a responder variable.  The limits of the range
    ! of values of Skx that the PDF is able to represent are calculated, and
    ! Skx is clipped when it exceeds the upper or lower limit.  Then, the PDF
    ! parameters for responder variable x are calculated.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three,        & ! Constant(s)
        two,          &
        three_halves, &
        one,          &
        zero

    use new_hybrid_pdf, only: &
        calculate_responder_params    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean of x (overall)                          [units vary]
      xp2,       & ! Variance of x (overall)                  [(units vary)^2]
      wpxp,      & ! Covariance of w and x (overall)        [m/s (units vary)]
      wp2,       & ! Variance of w (overall)                         [m^2/s^2]
      mixt_frac, & ! Mixture fraction                                      [-]
      F_w          ! Parameter for the spread of the PDF comp. means of w  [-]

    ! Input/Output Variable
    real( kind = core_rknd ), intent(inout) :: &
      Skx    ! Skewness of x (overall)              [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)        [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)        [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)    [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)    [(units vary)^2]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>    [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_w_x, & ! Correlation (overall) of w and x                  [-]
      min_Skx,  & ! Minimum Skx that can be represented by the PDF    [-]
      max_Skx     ! Maximum Skx that can be represented by the PDF    [-]


    ! Calculate the limits of Skx.
    if ( F_w > zero ) then

       ! Calculate the overall correlation of w and x.
       if ( wp2 * xp2 > zero ) then

          ! The overall correlation of w and x is:
          !
          ! corr_w_x = <w'x'> / sqrt( <w'^2> * <x'^2> ).
          corr_w_x = wpxp / sqrt( wp2 * xp2 )

       else ! <w'^2> * <x'^2> = 0

          ! When <w'^2> = 0 or <x'^2> = 0, <w'x'> = 0.  The correlation of w
          ! and x is undefined, however, since <w'x'> = 0, Skx = 0.  Setting
          ! corr_w_x = 0 in this scenario will set max_Skx = min_Skx = 0.
          corr_w_x = zero

       endif

       if ( wpxp >= zero ) then

          min_Skx = ( one + mixt_frac ) &
                    / sqrt( mixt_frac * ( one - mixt_frac ) ) &
                    * corr_w_x**3 / F_w**three_halves &
                    - sqrt( mixt_frac / ( one - mixt_frac ) ) &
                      * three * corr_w_x / sqrt( F_w )

          max_Skx = ( mixt_frac - two ) &
                    / sqrt( mixt_frac * ( one - mixt_frac ) ) &
                    * corr_w_x**3 / F_w**three_halves &
                    + sqrt( ( one - mixt_frac ) / mixt_frac ) &
                      * three * corr_w_x / sqrt( F_w )

       else ! <w'x'> < 0

          min_Skx = ( mixt_frac - two ) &
                    / sqrt( mixt_frac * ( one - mixt_frac ) ) &
                    * corr_w_x**3 / F_w**three_halves &
                    + sqrt( ( one - mixt_frac ) / mixt_frac ) &
                      * three * corr_w_x / sqrt( F_w )

          max_Skx = ( one + mixt_frac ) &
                    / sqrt( mixt_frac * ( one - mixt_frac ) ) &
                    * corr_w_x**3 / F_w**three_halves &
                    - sqrt( mixt_frac / ( one - mixt_frac ) ) &
                      * three * corr_w_x / sqrt( F_w )

       endif ! <w'x'> >= 0

    else ! F_w > 0

       ! When F_w = 0, <w'x'> = 0, and Skx = 0.
       min_Skx = zero
       max_Skx = zero

    endif ! F_w > 0

    ! Limit Skx so that min_Skx <= Skx <= max_Skx.
    if ( Skx > max_Skx ) then
       Skx = max_Skx
    elseif ( Skx < min_Skx ) then
       Skx = min_Skx
    endif ! Skx >= max_Skx

    ! Calculate the PDF parameters for responder variable x.
    call calculate_responder_params( xm, xp2, Skx, wpxp,  & ! In
                                     wp2, F_w, mixt_frac, & ! In
                                     mu_x_1, mu_x_2,      & ! Out
                                     sigma_x_1_sqd,       & ! Out
                                     sigma_x_2_sqd,       & ! Out
                                     coef_sigma_x_1_sqd,  & ! Out
                                     coef_sigma_x_2_sqd   ) ! Out


    return

  end subroutine calc_responder_var

  !=============================================================================
  elemental subroutine calc_F_w_zeta_w( Skw, wprtp, wpthlp, upwp, vpwp, & ! In
                                        wp2, rtp2, thlp2, up2, vp2,     & ! In
                                        max_corr_w_sclr_sqd,            & ! In
                                        F_w, zeta_w, min_F_w, max_F_w   ) ! Out

    ! Description:
    ! Calculates the values of F_w and zeta_w for w, which is the setter
    ! variable (which is the variable that sets the mixture fraction).
    !
    ! Based purely on the PDF of w and not considering restrictions caused by
    ! other variables in the multivariate PDF, the value of F_w is calculated
    ! between 0 (min_F_w) and 1 (max_F_w).  However, other variables in the PDF
    ! place restrictions on the minimum value of F_w.  The range of acceptible
    ! values for F_w is given by:
    !
    ! max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all variables x ) <= F_w <= 1.
    !
    ! Additionally, the value of F_w should still increase with an increasing
    ! magnitude of Skw.  The value of F_w can be 0 only when Skw = 0 (as well as
    ! all <w'x'> = 0, as well).  The F_w skewness equation used is:
    !
    ! F_w = max_F_w + ( min_F_w - max_F_w )
    !                 * exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w };
    !
    ! where lambda > 0 and slope_coef_spread_DG_means_w > 0.  As |Skw| goes
    ! toward 0, the value of F_w goes toward min_F_w, and as |Skw| becomes
    ! large, the value of F_w goes toward max_F_w.  When the value of
    ! slope_coef_spread_DG_means_w is small, the value of F_w tends toward
    ! max_F_w, and when the value of slope_coef_spread_DG_means_w is large, the
    ! value of F_w tends toward min_F_w.  When lambda is small, the value of F_w
    ! is less dependent on Skw, and when lambda is large, the value of F_w is
    ! more dependent on Skw.
    !
    ! Mathematically, this equation will always produce a value of F_w that
    ! falls between min_F_w and max_F_w.  However, in order to prevent a value
    ! of F_w from being calculated outside the bounds of min_F_w and max_F_w
    ! owing to numerical underflow or loss of precision, this equation can be
    ! rewritten as:
    !
    ! F_w
    ! = min_F_w * exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w }
    !   + max_F_w * ( 1 - exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w } ).
    !
    ! The value of zeta_w used to adjust the PDF component standard devations:
    !
    ! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
    !              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
    !
    ! where zeta_w > -1.  The sign of zeta_w is used to easily determine if
    ! mixt_frac * sigma_w_1^2 is greater than ( 1 - mixt_frac ) * sigma_w_2^2
    ! (when zeta_w is positive), mixt_frac * sigma_w_1^2 is less than
    ! ( 1 - mixt_frac ) * sigma_w_2^2 (when zeta_w is negative), or
    ! mixt_frac * sigma_w_1^2 is equal to ( 1 - mixt_frac ) * sigma_w_2^2 (when
    ! zeta_w is 0).
    !
    ! In order to allow for a tunable parameter that is the pure ratio of
    ! mixt_frac * sigma_w_1^2 to ( 1 - mixt_frac ) * sigma_w_2^2, zeta_w is
    ! related to the parameter pdf_component_stdev_factor_w, where:
    !
    ! 1 + zeta_w = pdf_component_stdev_factor_w.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use parameters_tunable, only: &
        slope_coef_spread_DG_means_w, & ! Variable(s)
        pdf_component_stdev_factor_w

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw,    & ! Skewness of w (overall)              [-]
      wprtp,  & ! Covariance (overall) of w and rt     [m/s (kg/kg)]
      wpthlp, & ! Covariance (overall) of w and thl    [m/s K]
      upwp,   & ! Covariance (overall) of u and w      [m^2/s^2]
      vpwp,   & ! Covariance (overall) of u and v      [m^2/s^2]
      wp2,    & ! Variance (overall) of w              [m^2/s^2]
      rtp2,   & ! Variance (overall) of rt             [(kg/kg)^2]
      thlp2,  & ! Variance (overall) of thl            [K^2]
      up2,    & ! Variance (overall) of u              [m^2/s^2]
      vp2       ! Variance (overall) of v              [m^2/s^2]

    real( kind = core_rknd ), intent(in) :: &
      max_corr_w_sclr_sqd    ! Max value of wpsclrp^2 / ( wp2 * sclrp2 )     [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      F_w,     & ! Parameter for the spread of the PDF component means of w  [-]
      zeta_w,  & ! Parameter for the PDF component variances of w            [-]
      min_F_w, & ! Minimum allowable value of parameter F_w                  [-]
      max_F_w    ! Maximum allowable value of parameter F_w                  [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_w_rt_sqd,    & ! Value of wprtp^2 / ( wp2 * rtp2 )              [-]
      corr_w_thl_sqd,   & ! Value of wpthlp^2 / ( wp2 * thlp2 )            [-]
      corr_u_w_sqd,     & ! Value of upwp^2 / ( up2 * wp2 )                [-]
      corr_v_w_sqd,     & ! Value of vpwp^2 / ( vp2 * wp2 )                [-]
      max_corr_w_x_sqd    ! Max. val. of wpxp^2/(wp2*xp2) for all vars. x  [-]

    real( kind = core_rknd ) :: &
      exp_Skw_interp_factor    ! Function to interp. between min. and max.   [-]

    real( kind = core_rknd ), parameter :: &
      lambda = 0.5_core_rknd    ! Parameter for Skw dependence    [-]


    ! Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
    ! variables x that are Double Gaussian PDF responder variables.  This
    ! includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
    ! also calculated as part of the PDF, and they are included as well.
    ! Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.

    ! Calculate the squared value of the correlation of w and rt.
    if ( wp2 * rtp2 > zero ) then
       corr_w_rt_sqd = wprtp**2 / ( wp2 * rtp2 )
    else
       corr_w_rt_sqd = zero
    endif

    ! Calculate the squared value of the correlation of w and thl.
    if ( wp2 * thlp2 > zero ) then
       corr_w_thl_sqd = wpthlp**2 / ( wp2 * thlp2 )
    else
       corr_w_thl_sqd = zero
    endif

    ! Calculate the squared value of the correlation of u and w.
    if ( up2 * wp2 > zero ) then
       corr_u_w_sqd = upwp**2 / ( up2 * wp2 )
    else
       corr_u_w_sqd = zero
    endif

    ! Calculate the squared value of the correlation of v and w.
    if ( vp2 * wp2 > zero ) then
       corr_v_w_sqd = vpwp**2 / ( vp2 * wp2 )
    else
       corr_v_w_sqd = zero
    endif

    ! Calculate the maximum value of corr_w_rt^2, corr_w_thl^2, corr_u_w^2, and
    ! corr_v_w^2 in the calculation of max(corr_w_x^2).  Include the maximum
    ! value of corr_w_sclr^2 (for all sclrs) when scalars are found.
    max_corr_w_x_sqd = max( corr_w_rt_sqd, corr_w_thl_sqd, &
                            corr_u_w_sqd, corr_v_w_sqd, max_corr_w_sclr_sqd )

    ! Calculate min_F_w and set max_F_w to 1.
    if ( abs( Skw ) > zero ) then
       min_F_w = max( max_corr_w_x_sqd, 1.0e-3_core_rknd )
    else
       min_F_w = max_corr_w_x_sqd
    endif
    max_F_w = one

    ! F_w must have a value between min_F_w and max_F_w.
    exp_Skw_interp_factor &
    = exp( -abs(Skw)**lambda / slope_coef_spread_DG_means_w )

    F_w = min_F_w * exp_Skw_interp_factor &
          + max_F_w * ( one - exp_Skw_interp_factor )

    ! The value of zeta_w must be greater than -1.
    zeta_w = pdf_component_stdev_factor_w - one


    return

  end subroutine calc_F_w_zeta_w

  !=============================================================================

end module new_hybrid_pdf_main
