! $Id$
!===============================================================================
module new_pdf_main

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! trivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), and liquid water potential temperature (thl).

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: new_pdf_driver    ! Procedure(s)

  private :: calc_responder_var,     & ! Procedure(s)
             calc_F_x_zeta_x_setter, &
             calc_F_x_Skx_fnc

  private

  contains

  !=============================================================================
  subroutine new_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,          & ! In
                             Skw_in, Skrt_in, Skthl_in, wprtp, wpthlp, & ! In
                             mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,         & ! Out
                             mu_thl_1, mu_thl_2, sigma_w_1_sqd,        & ! Out
                             sigma_w_2_sqd, sigma_rt_1_sqd,            & ! Out
                             sigma_rt_2_sqd, sigma_thl_1_sqd,          & ! Out
                             sigma_thl_2_sqd, mixt_frac,               & ! Out
                             F_w, F_rt, F_thl, min_F_w, max_F_w,       & ! Out
                             min_F_rt, max_F_rt, min_F_thl, max_F_thl  ) ! Out
                             

    ! Description:
    ! Selects which variable is used to set the mixture fraction for the PDF
    ! ("the setter") and which variables are handled after that mixture fraction
    ! has been set ("the responders").  Traditionally, w has been used to set
    ! the PDF.  However, here, the variable with the greatest magnitude of
    ! skewness is used to set the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four, & ! Variable(s)
        two,  &
        one,  &
        zero

    use new_pdf, only: &
        calc_setter_var_params,    & ! Procedure(s)
        calc_coef_wp4_implicit,    &
        calc_coef_wpxp2_implicit,  &
        calc_coefs_wp2xp_semiimpl

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wm,       & ! Mean of w (overall)                [m/s]
      rtm,      & ! Mean of rt (overall)               [kg/kg]
      thlm,     & ! Mean of thl (overall)              [K]
      wp2,      & ! Variance of w (overall)            [m^2/s^2]
      rtp2,     & ! Variance of rt (overall)           [kg^2/kg^2]
      thlp2,    & ! Variance of thl (overall)          [K^2]
      Skw_in,   & ! Skewness of w (overall)            [-]
      Skrt_in,  & ! Skewness of rt (overall)           [-]
      Skthl_in, & ! Skewness of thl (overall)          [-]
      wprtp,    & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
      wpthlp      ! Covariance of w and thl (overall)  [(m/s)K]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
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

    ! Output only for recording statistics.
    real( kind = core_rknd ), intent(out) :: &
      F_w,   & ! Parameter for the spread of the PDF component means of w    [-]
      F_rt,  & ! Parameter for the spread of the PDF component means of rt   [-]
      F_thl    ! Parameter for the spread of the PDF component means of thl  [-]

    real( kind = core_rknd ), intent(out) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      sigma_w_1,   & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2,   & ! Standard deviation of w (2nd PDF component)      [m/s]
      sigma_rt_1,  & ! Standard deviation of rt (1st PDF component)     [kg/kg]
      sigma_rt_2,  & ! Standard deviation of rt (2nd PDF component)     [kg/kg]
      sigma_thl_1, & ! Standard deviation of thl (1st PDF component)    [K]
      sigma_thl_2, & ! Standard deviation of thl (2nd PDF component)    [K]
      Skw,         & ! Skewness of w (overall)                          [-]
      Skrt,        & ! Skewness of rt (overall)                         [-]
      Skthl,       & ! Skewness of thl (overall)                        [-]
      sgn_wprtp,   & ! Sign of the covariance of w and rt (overall)     [-]
      sgn_wpthlp     ! Sign of the covariance of w and thl (overall)    [-]

    real( kind = core_rknd ), parameter :: &
      sgn_wp2 = one   ! Sign of the variance of w (overall); always positive [-]

    real( kind = core_rknd ) :: &
      coef_sigma_w_1_sqd,   & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
      coef_sigma_w_2_sqd,   & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
      coef_sigma_rt_1_sqd,  & ! sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
      coef_sigma_rt_2_sqd,  & ! sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
      coef_sigma_thl_1_sqd, & ! sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
      coef_sigma_thl_2_sqd    ! sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]

    real( kind = core_rknd ) :: &
      max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
      max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

    ! The zeta_x parameter is only used when a variable is the setter variable.
    real( kind = core_rknd ) :: &
      zeta_w,   & ! Parameter for the PDF component variances of w           [-]
      zeta_rt,  & ! Parameter for the PDF component variances of rt          [-]
      zeta_thl    ! Parameter for the PDF component variances of thl         [-]

    real ( kind = core_rknd ) :: &
      coef_wp4_implicit,     & ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]
      coef_wprtp2_implicit,  & ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>  [-]
      coef_wpthlp2_implicit    ! <w'thl'^2> = coef_wpthlp2_implicit*<thl'^2> [-]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real ( kind = core_rknd ) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real ( kind = core_rknd ) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    logical, parameter :: &
      l_use_w_setter_var = .false. ! Flag to always use w as the setter variable


    Skw = Skw_in
    Skrt = Skrt_in
    Skthl = Skthl_in

    ! Calculate sgn( <w'rt'> ).
    if ( wprtp >= zero ) then
       sgn_wprtp = one
    else ! wprtp < 0
       sgn_wprtp = -one
    endif ! wprtp >= 0

    ! Calculate sgn( <w'thl'> ).
    if ( wpthlp >= zero ) then
       sgn_wpthlp = one
    else ! wpthlp < 0
       sgn_wpthlp = -one
    endif ! wpthlp >= 0


    ! The variable with the greatest magnitude of skewness will be the setter
    ! variable and the other variables will be responder variables.  When the
    ! l_use_w_setter_var flag is enabled, w will always be the setter variable.
    if ( ( abs( Skw ) >= abs( Skrt ) .and. abs( Skw ) >= abs( Skthl ) ) &
         .or. l_use_w_setter_var ) then

       ! The variable w has the greatest magnitude of skewness or the
       ! l_use_w_setter_var flag is enabled.

       call calc_F_x_zeta_x_setter( Skw, 200.0_core_rknd, & ! In
                                    F_w, zeta_w,          & ! Out
                                    min_F_w, max_F_w      ) ! Out

       ! Calculate the PDF parameters, including mixture fraction, for the
       ! setter variable, w.
       call calc_setter_var_params( wm, wp2, Skw, sgn_wp2,     & ! In
                                    F_w, zeta_w,               & ! In
                                    mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                    sigma_w_2, mixt_frac,      & ! Out
                                    coef_sigma_w_1_sqd,        & ! Out
                                    coef_sigma_w_2_sqd         ) ! Out

       sigma_w_1_sqd = sigma_w_1**2
       sigma_w_2_sqd = sigma_w_2**2

       ! Calculate the upper limit on the magnitude of skewness for responding
       ! variables.
       max_Skx2_pos_Skx_sgn_wpxp = four * ( one - mixt_frac )**2 &
                                   / ( mixt_frac * ( two - mixt_frac ) )

       max_Skx2_neg_Skx_sgn_wpxp = four * mixt_frac**2 / ( one - mixt_frac**2 )

       ! Calculate the PDF parameters for responder variable rt.
       call calc_responder_var( rtm, rtp2, sgn_wprtp, mixt_frac, & ! In
                                200.0_core_rknd,                 & ! In
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
                                200.0_core_rknd,                    & ! In
                                max_Skx2_pos_Skx_sgn_wpxp,          & ! In
                                max_Skx2_neg_Skx_sgn_wpxp,          & ! In
                                Skthl,                              & ! In/Out
                                mu_thl_1, mu_thl_2,                 & ! Out
                                sigma_thl_1_sqd, sigma_thl_2_sqd,   & ! Out
                                coef_sigma_thl_1_sqd,               & ! Out
                                coef_sigma_thl_2_sqd,               & ! Out
                                F_thl, min_F_thl, max_F_thl         ) ! Out

    elseif ( abs( Skrt ) > abs( Skw ) .and. abs( Skrt ) >= abs( Skthl ) ) then

       ! The variable rt has the greatest magnitude of skewness.

       call calc_F_x_zeta_x_setter( Skrt, 200.0_core_rknd,  & ! In
                                    F_rt, zeta_rt,          & ! Out
                                    min_F_rt, max_F_rt      ) ! Out

       ! Calculate the PDF parameters, including mixture fraction, for the
       ! setter variable, rt.
       call calc_setter_var_params( rtm, rtp2, Skrt, sgn_wprtp,   & ! In
                                    F_rt, zeta_rt,                & ! In
                                    mu_rt_1, mu_rt_2, sigma_rt_1, & ! Out
                                    sigma_rt_2, mixt_frac,        & ! Out
                                    coef_sigma_rt_1_sqd,          & ! Out
                                    coef_sigma_rt_2_sqd           ) ! Out

       sigma_rt_1_sqd = sigma_rt_1**2
       sigma_rt_2_sqd = sigma_rt_2**2

       ! Calculate the upper limit on the magnitude of skewness for responding
       ! variables.
       max_Skx2_pos_Skx_sgn_wpxp = four * ( one - mixt_frac )**2 &
                                   / ( mixt_frac * ( two - mixt_frac ) )

       max_Skx2_neg_Skx_sgn_wpxp = four * mixt_frac**2 / ( one - mixt_frac**2 )

       ! Calculate the PDF parameters for responder variable w.
       call calc_responder_var( wm, wp2, sgn_wp2, mixt_frac,  & ! In
                                200.0_core_rknd,              & ! In
                                max_Skx2_pos_Skx_sgn_wpxp,    & ! In
                                max_Skx2_neg_Skx_sgn_wpxp,    & ! In
                                Skw,                          & ! In/Out
                                mu_w_1, mu_w_2,               & ! Out
                                sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                                coef_sigma_w_1_sqd,           & ! Out
                                coef_sigma_w_2_sqd,           & ! Out
                                F_w, min_F_w, max_F_w         ) ! Out

       ! Calculate the PDF parameters for responder variable thl.
       call calc_responder_var( thlm, thlp2, sgn_wpthlp, mixt_frac, & ! In
                                200.0_core_rknd,                    & ! In
                                max_Skx2_pos_Skx_sgn_wpxp,          & ! In
                                max_Skx2_neg_Skx_sgn_wpxp,          & ! In
                                Skthl,                              & ! In/Out
                                mu_thl_1, mu_thl_2,                 & ! Out
                                sigma_thl_1_sqd, sigma_thl_2_sqd,   & ! Out
                                coef_sigma_thl_1_sqd,               & ! Out
                                coef_sigma_thl_2_sqd,               & ! Out
                                F_thl, min_F_thl, max_F_thl         ) ! Out

    else ! abs( Skthl ) > abs( Skw ) .and. abs( Skthl ) > abs( Skrt )

       ! The variable thl has the greatest magnitude of skewness.

       call calc_F_x_zeta_x_setter( Skthl, 200.0_core_rknd,   & ! In
                                    F_thl, zeta_thl,          & ! Out
                                    min_F_thl, max_F_thl      ) ! Out

       ! Calculate the PDF parameters, including mixture fraction, for the
       ! setter variable, thl.
       call calc_setter_var_params( thlm, thlp2, Skthl, sgn_wpthlp,  & ! In
                                    F_thl, zeta_thl,                 & ! In
                                    mu_thl_1, mu_thl_2, sigma_thl_1, & ! Out
                                    sigma_thl_2, mixt_frac,          & ! Out
                                    coef_sigma_thl_1_sqd,            & ! Out
                                    coef_sigma_thl_2_sqd             ) ! Out

       sigma_thl_1_sqd = sigma_thl_1**2
       sigma_thl_2_sqd = sigma_thl_2**2

       ! Calculate the upper limit on the magnitude of skewness for responding
       ! variables.
       max_Skx2_pos_Skx_sgn_wpxp = four * ( one - mixt_frac )**2 &
                                   / ( mixt_frac * ( two - mixt_frac ) )

       max_Skx2_neg_Skx_sgn_wpxp = four * mixt_frac**2 / ( one - mixt_frac**2 )

       ! Calculate the PDF parameters for responder variable w.
       call calc_responder_var( wm, wp2, sgn_wp2, mixt_frac,  & ! In
                                200.0_core_rknd,              & ! In
                                max_Skx2_pos_Skx_sgn_wpxp,    & ! In
                                max_Skx2_neg_Skx_sgn_wpxp,    & ! In
                                Skw,                          & ! In/Out
                                mu_w_1, mu_w_2,               & ! Out
                                sigma_w_1_sqd, sigma_w_2_sqd, & ! Out
                                coef_sigma_w_1_sqd,           & ! Out
                                coef_sigma_w_2_sqd,           & ! Out
                                F_w, min_F_w, max_F_w         ) ! Out

       ! Calculate the PDF parameters for responder variable rt.
       call calc_responder_var( rtm, rtp2, sgn_wprtp, mixt_frac, & ! In
                                200.0_core_rknd,                 & ! In
                                max_Skx2_pos_Skx_sgn_wpxp,       & ! In
                                max_Skx2_neg_Skx_sgn_wpxp,       & ! In
                                Skrt,                            & ! In/Out
                                mu_rt_1, mu_rt_2,                & ! Out
                                sigma_rt_1_sqd, sigma_rt_2_sqd,  & ! Out
                                coef_sigma_rt_1_sqd,             & ! Out
                                coef_sigma_rt_2_sqd,             & ! Out
                                F_rt, min_F_rt, max_F_rt         ) ! Out

    endif ! Find variable with the greatest magnitude of skewness.


    ! <w'^4> = coef_wp4_implicit * <w'^2>^2.
    coef_wp4_implicit &
    = calc_coef_wp4_implicit( mixt_frac, F_w, &
                              coef_sigma_w_1_sqd, &
                              coef_sigma_w_2_sqd )

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


    return

  end subroutine new_pdf_driver

  !=============================================================================
  subroutine calc_responder_var( xm, xp2, sgn_wpxp, mixt_frac, & ! In
                                 coefficient,                  & ! In
                                 max_Skx2_pos_Skx_sgn_wpxp,    & ! In
                                 max_Skx2_neg_Skx_sgn_wpxp,    & ! In
                                 Skx,                          & ! In/Out
                                 mu_x_1, mu_x_2,               & ! Out
                                 sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                 coef_sigma_x_1_sqd,           & ! Out
                                 coef_sigma_x_2_sqd,           & ! Out
                                 F_x, min_F_x, max_F_x         ) ! Out

    ! Description:
    ! This is the sub-driver for a responder variable.  The upper limits of the
    ! magnitude of Skx are calculated, and Skx is clipped when its magnitude
    ! exceeds the upper limits.  The limits of the F_x parameter are calculated,
    ! and the value of F_x is set within those limits.  Then, the PDF parameters
    ! for responder variable x are calculated.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Variable(s)

    use new_pdf, only: &
        calc_limits_F_x_responder, & ! Procedure(s)
        calc_responder_params

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,          & ! Mean of x (overall)                  [units vary]
      xp2,         & ! Variance of x (overall)              [(units vary)^2]
      sgn_wpxp,    & ! Sign of the covariance of w and x    [-]
      mixt_frac,   & ! Mixture fraction                     [-]
      coefficient    ! Coefficient                          [-]

    real( kind = core_rknd ), intent(in) :: &
      max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
      max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

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

    ! Output only for recording statistics.
    real( kind = core_rknd ), intent(out) :: &
      F_x,     & ! Parameter for the spread of the PDF component means of x  [-]
      min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
      max_F_x    ! Maximum allowable value of parameter F_x                  [-]


    ! Calculate the upper limit of the magnitude of Skx.
    if ( Skx * sgn_wpxp >= zero ) then
       if ( Skx**2 >= max_Skx2_pos_Skx_sgn_wpxp ) then
          if ( Skx >= zero ) then
             Skx = sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          else
             Skx = -sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          endif
       endif ! Skx^2 >= max_Skx2_pos_Skx_sgn_wpxp
    else ! Skx * sgn( <w'x'> ) < 0
       if ( Skx**2 >= max_Skx2_neg_Skx_sgn_wpxp ) then
          if ( Skx >= zero ) then
             Skx = sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          else
             Skx = -sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          endif
       endif ! Skx^2 >= max_Skx2_neg_Skx_sgn_wpxp
    endif ! Skx * sgn( <w'x'> ) >= 0

    call calc_limits_F_x_responder( mixt_frac, Skx, sgn_wpxp,  & ! In
                                    max_Skx2_pos_Skx_sgn_wpxp, & ! In
                                    max_Skx2_neg_Skx_sgn_wpxp, & ! In
                                    min_F_x, max_F_x )           ! Out

    ! F_x must have a value between min_F_x and max_F_x.
    F_x = calc_F_x_Skx_fnc( Skx, coefficient, min_F_x, max_F_x )

    call calc_responder_params( xm, xp2, Skx, sgn_wpxp,       & ! In
                                F_x, mixt_frac,               & ! In
                                mu_x_1, mu_x_2,               & ! Out
                                sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                coef_sigma_x_1_sqd,           & ! Out
                                coef_sigma_x_2_sqd            ) ! Out


    return

  end subroutine calc_responder_var

  !=============================================================================
  subroutine calc_F_x_zeta_x_setter( Skx, coefficient, & ! In
                                     F_x, zeta_x,      & ! Out
                                     min_F_x, max_F_x  ) ! Out

    ! Description:
    ! Calculates the values of F_x and zeta_x for the setter variable (which is
    ! the variable that sets the mixture fraction).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx,         & ! Skewness of x (overall)      [-]
      coefficient    ! Coefficient                  [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      F_x,     & ! Parameter for the spread of the PDF component means of x  [-]
      zeta_x,  & ! Parameter for the PDF component variances of x            [-]
      min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
      max_F_x    ! Maximum allowable value of parameter F_x                  [-]


    min_F_x = zero
    max_F_x = one

    F_x = one - exp( -Skx**2 / coefficient )
    zeta_x = zero


    return

  end subroutine calc_F_x_zeta_x_setter

  !=============================================================================
  function calc_F_x_Skx_fnc( Skx, coefficient, min_F_x, max_F_x ) &
  result( F_x )

    ! Description:
    ! Calculates the value of F_x as a function of Skx between min_F_x and
    ! max_F_x.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx,         & ! Skewness of x (overall)                          [-]
      coefficient, & ! Coefficient                                      [-]
      min_F_x,     & ! Minimum allowable value of parameter F_x         [-]
      max_F_x        ! Maximum allowable value of parameter F_x         [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      F_x    ! Parameter for the spread of the PDF component means of x  [-]


    ! F_x must have a value between min_F_x and max_F_x.
    F_x = max_F_x + ( min_F_x - max_F_x ) * exp( -Skx**2 / coefficient )


    return

  end function calc_F_x_Skx_fnc

  !=============================================================================

end module new_pdf_main
