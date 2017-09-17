! $Id$
!===============================================================================
module new_pdf

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! trivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), and liquid water potential temperature (thl).

  ! References:
  ! Griffin and Larson (2018)
  !-------------------------------------------------------------------------

  implicit none

  public :: new_pdf_driver,        & ! Procedure(s)
            calc_mixture_fraction, &
            calc_setter_var_params

  private :: calc_responder_params,     & ! Procedure(s)
             calc_limits_F_x_responder, &
             sort_roots

  private

  contains

  !=============================================================================
  subroutine new_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,          & ! In
                             Skw_in, Skrt_in, Skthl_in, wprtp, wpthlp, & ! In
                             mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,         & ! Out
                             mu_thl_1, mu_thl_2, sigma_w_1_sqd,        & ! Out
                             sigma_w_2_sqd, sigma_rt_1_sqd,            & ! Out
                             sigma_rt_2_sqd, sigma_thl_1_sqd,          & ! Out
                             sigma_thl_2_sqd, mixt_frac                ) ! Out
                             

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

    ! Local Variables
    real( kind = core_rknd ) :: &
      sigma_w_1,  & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2,  & ! Standard deviation of w (2nd PDF component)      [m/s]
      Skw,        & ! Skewness of w (overall)                          [-]
      Skrt,       & ! Skewness of rt (overall)                         [-]
      Skthl,      & ! Skewness of thl (overall)                        [-]
      sgn_wprtp,  & ! Sign of the covariance of w and rt (overall)     [-]
      sgn_wpthlp    ! Sign of the covariance of w and thl (overall)    [-]

    real( kind = core_rknd ), parameter :: &
      sgn_wp2 = one   ! Sign of the variance of w (overall); always positive [-]

    real( kind = core_rknd ) :: &
      max_Skx2_pos_Skx_sgn_wpxp, &
      max_Skx2_neg_Skx_sgn_wpxp

    real( kind = core_rknd ) :: &
      min_F_rt,  & !
      max_F_rt,  & !
      min_F_thl, & !
      max_F_thl    !

    real( kind = core_rknd ) :: &
      F_w,    & ! Parameter for the spread of the PDF component means of w   [-]
      zeta_w, & ! Parameter for the PDF component variances of w             [-]
      F_rt,   & ! Parameter for the spread of the PDF component means of rt  [-]
      F_thl     ! Parameter for the spread of the PDF component means of thl [-]


    Skw = Skw_in
    Skrt = Skrt_in
    Skthl = Skthl_in

    if ( wprtp >= zero ) then
       sgn_wprtp = one
    else ! wprtp < 0
       sgn_wprtp = -one
    endif ! wprtp >= 0

    if ( wpthlp >= zero ) then
       sgn_wpthlp = one
    else ! wpthlp < 0
       sgn_wpthlp = -one
    endif ! wpthlp >= 0

    F_w = one - exp( -Skw**2 / 200.0_core_rknd )
    zeta_w = zero

    call calc_setter_var_params( wm, wp2, Skw, sgn_wp2,     & ! In
                                 F_w, zeta_w,               & ! In
                                 mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                 sigma_w_2, mixt_frac )       ! Out

    sigma_w_1_sqd = sigma_w_1**2
    sigma_w_2_sqd = sigma_w_2**2

    ! Calculate the upper limit on the magnitude of skewness for responding
    ! variables.
    max_Skx2_pos_Skx_sgn_wpxp = four * ( one - mixt_frac )**2 &
                                / ( mixt_frac * ( two - mixt_frac ) )

    max_Skx2_neg_Skx_sgn_wpxp = four * mixt_frac**2 / ( one - mixt_frac**2 )

    ! Calculate the upper limit of the magnitude of Skrt.
    if ( Skrt * sgn_wprtp >= zero ) then
       if ( Skrt**2 >= max_Skx2_pos_Skx_sgn_wpxp ) then
          if ( Skrt >= zero ) then
             Skrt = sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          else
             Skrt = -sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          endif
       endif ! Skrt^2 >= max_Skx2_pos_Skx_sgn_wpxp
    else ! Skrt * sgn( <w'rt'> ) < 0
       if ( Skrt**2 >= max_Skx2_neg_Skx_sgn_wpxp ) then
          if ( Skrt >= zero ) then
             Skrt = sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          else
             Skrt = -sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          endif
       endif ! Skrt^2 >= max_Skx2_neg_Skx_sgn_wpxp
    endif ! Skrt * sgn( <w'rt'> ) >= 0

    call calc_limits_F_x_responder( mixt_frac, Skrt, sgn_wprtp, & ! In
                                    max_Skx2_pos_Skx_sgn_wpxp,  & ! In
                                    max_Skx2_neg_Skx_sgn_wpxp,  & ! In
                                    min_F_rt, max_F_rt )          ! Out

    ! F_rt must have a value between min_F_rt and max_F_rt.
    F_rt = max_F_rt &
           + ( min_F_rt - max_F_rt ) * exp( -Skrt**2 / 200.0_core_rknd )

    call calc_responder_params( rtm, rtp2, Skrt, sgn_wprtp,       & ! In
                                F_rt, mixt_frac,                  & ! In
                                mu_rt_1, mu_rt_2,                 & ! In
                                sigma_rt_1_sqd, sigma_rt_2_sqd )    ! Out

    ! Calculate the upper limit of the magnitude of Skthl.
    if ( Skthl * sgn_wpthlp >= zero ) then
       if ( Skthl**2 >= max_Skx2_pos_Skx_sgn_wpxp ) then
          if ( Skthl >= zero ) then
             Skthl = sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          else
             Skthl = -sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          endif
       endif ! Skthl^2 >= max_Skx2_pos_Skx_sgn_wpxp
    else ! Skthl * sgn( <w'thl'> ) < 0
       if ( Skthl**2 >= max_Skx2_neg_Skx_sgn_wpxp ) then
          if ( Skthl >= zero ) then
             Skthl = sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          else
             Skthl = -sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          endif
       endif ! Skthl^2 >= max_Skx2_neg_Skx_sgn_wpxp
    endif ! Skthl * sgn( <w'thl'> ) >= 0

    call calc_limits_F_x_responder( mixt_frac, Skthl, sgn_wpthlp, & ! In
                                    max_Skx2_pos_Skx_sgn_wpxp,    & ! In
                                    max_Skx2_neg_Skx_sgn_wpxp,    & ! In
                                    min_F_thl, max_F_thl )          ! Out

    ! F_thl must have a value between min_F_thl and max_F_thl.
    F_thl = max_F_thl &
            + ( min_F_thl - max_F_thl ) * exp( -Skthl**2 / 200.0_core_rknd )

    call calc_responder_params( thlm, thlp2, Skthl, sgn_wpthlp,     & ! In
                                F_thl, mixt_frac,                   & ! In
                                mu_thl_1, mu_thl_2,                 & ! In
                                sigma_thl_1_sqd, sigma_thl_2_sqd )    ! Out

   
    return

  end subroutine new_pdf_driver

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =========================================================================
  !
  ! Many times, w has been used as the variable that sets the mixture fraction
  ! for the PDF.  There are five PDF parameters that need to be calculated,
  ! which are mu_w_1 (the mean of w is the 1st PDF component), mu_w_2 (the mean
  ! of w in the 2nd PDF component), sigma_w_1 (the standard deviation of w in
  ! the 1st PDF component), sigma_w_2 (the standard deviation of w in the 2nd
  ! PDF component), and mixt_frac (the mixture fraction, which is the weight of
  ! the 1st PDF component).  In order to solve for these five parameters, five
  ! equations are needed.  These five equations are the equations for <w>,
  ! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally, two
  ! more equations, which involve tunable parameters F_w and zeta_w, and which
  ! are used to help control the spread of the PDF component means and the size
  ! of the PDF component standard deviations compared to each other,
  ! respectively, are used in this equation set.  The five equations are:
  !
  ! <w> = mixt_frac * mu_w_1 + ( 1 - mixt_frac ) * mu_w_2;
  !
  ! <w'^2> = mixt_frac * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 );
  !
  ! <w'^3> = mixt_frac * ( mu_w_1 - <w> )
  !                    * ( ( mu_w_1 - <w> )^2 + 3 * sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> )
  !                              * ( ( mu_w_2 - <w> )^2 + 3 * sigma_w_2^2 );
  !
  ! mu_w_1 - <w> = sqrt(F_w) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <w'^2> );
  !
  ! where 0 <= F_w <= 1; and
  !
  ! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
  !              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
  !
  ! where zeta_w > -1.
  !
  ! Following convention for w, mu_w_1 is defined to be greater than or equal to
  ! mu_w_2 (and is also greater than or equal to <w>, while mu_w_2 is less than
  ! or equal to <w>).  This is relationship is found in the mu_w_1 - <w>
  ! equation above.
  !
  ! The resulting equations for the five PDF parameters are:
  !
  ! mixt_frac
  ! = ( 4 * F_w^3
  !     + 18 * F_w * ( zeta_w + 1 ) * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + 6 * F_w^2 * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + Skw^2
  !     - Skw * sqrt( 4 * F_w^3
  !                   + 12 * F_w^2 * ( 1 - F_w )
  !                   + 36 * F_w * ( zeta_w + 1 ) * ( 1 - F_w )^2
  !                     / ( zeta_w + 2 )^2
  !                   + Skw^2 ) )
  !   / ( 2 * F_w * ( F_w - 3 )^2 + 2 * Skw^2 );
  !
  ! mu_w_1 = <w> + sqrt( F_w * ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> );
  !
  ! mu_w_2 = <w> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_w_1 - <w> );
  !
  ! sigma_w_1
  ! = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
  !         / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> ); and
  !
  ! sigma_w_2
  ! = sqrt( ( mixt_frac * sigma_w_1^2 )
  !         / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) ).
  !
  ! This method works for all values of F_w (where 0 <= F_w <= 1) and zeta_w
  ! (where zeta_w > -1).
  !
  !
  ! Generalized equations for any variable, x, that sets the mixture fraction:
  !
  ! A slight alteration is made to the above equations in order to have any
  ! variable, x, set the mixture fraction.  The same five PDF parameters need to
  ! be calculated for the setting variable, which are mu_x_1 (the mean of x in
  ! the 1st PDF component), mu_x_2 (the mean of x in the 2nd PDF component),
  ! sigma_x_1 (the standard deviation of x in the 1st PDF component), sigma_x_2
  ! (the standard deviation of x in the 2nd PDF component), and mixt_frac (the
  ! mixture fraction).  Again, five equations are needed, and they are the
  ! equations for <x>, <x'^2>, and <x'^3> as found by integrating over the PDF,
  ! as well as the equations that involve tunable parameters F_x and zeta_x.
  ! However, the equation for F_x is multiplied by a new variable,
  ! sgn( <w'x'> ), where <w'x'> is the covariance of w and x, and sgn( <w'x'> )
  ! is given by:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The five equations are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 );
  !
  ! mu_x_1 - <x> = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= F_x <= 1; and
  !
  ! 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
  !              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
  !
  ! where zeta_x > -1.
  !
  ! The only equations that are altered are the equation for mu_x_1 and the
  ! equation for mixt_frac, which now both contain a sgn( <w'x'> ).  The mu_x_2
  ! equation is not altered, but the sign of mu_x_2 - <x> will be the opposite
  ! of the sign of mu_x_1 - <x>.  The resulting equations for the five PDF
  ! parameters are:
  !
  ! mixt_frac
  ! = ( 4 * F_x^3
  !     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + Skx^2
  !     - Skx * sgn( <w'x'> )
  !           * sqrt( 4 * F_x^3
  !                   + 12 * F_x^2 * ( 1 - F_x )
  !                   + 36 * F_x * ( zeta_x + 1 ) * ( 1 - F_x )^2
  !                     / ( zeta_x + 2 )^2
  !                   + Skx^2 ) )
  !   / ( 2 * F_x * ( F_x - 3 )^2 + 2 * Skx^2 );
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
  !
  ! sigma_x_1
  ! = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
  !         / ( ( zeta_x + 2 ) * mixt_frac ) * <x'^2> ); and
  !
  ! sigma_x_2
  ! = sqrt( ( mixt_frac * sigma_x_1^2 )
  !         / ( ( 1 - mixt_frac ) * ( 1 + zeta_x ) ) ).
  !
  ! This method works for all values of F_x (where 0 <= F_x <= 1) and zeta_x
  ! (where zeta_x > -1).
  !
  ! When the generalized form is solved for w (x = w), sgn( <w'^2> ) = 1, and
  ! the equations are unaltered from the equations listed above for w.
  !
  !
  ! Special case:
  !
  ! When Skx = 0 and F_x = 0, the equation for mixt_frac is undefined.  The
  ! equation for mixture fraction in this scenario can be derived by using the
  ! above equation for mixture fraction and then setting Skx = 0.  The resulting
  ! equation becomes:
  !
  ! mixt_frac
  ! = ( 4 * F_x^3
  !     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 ) )
  !   / ( 2 * F_x * ( F_x - 3 )^2 ).
  !
  ! All of the terms in the numerator and denominator contain a F_x, so this
  ! equation can be rewritten as:
  !
  ! mixt_frac
  ! = ( 4 * F_x^2
  !     + 18 * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x * ( 1 - F_x ) / ( zeta_x + 2 ) )
  !   / ( 2 * ( F_x - 3 )^2 ).
  !
  ! Now setting F_x = 0, the equation becomes:
  !
  ! mixt_frac = ( 18 * ( zeta_x + 1 ) / ( zeta_x + 2 ) ) / 18;
  !
  ! which can be rewritten as:
  !
  ! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).
  !
  ! When F_x = 0, Skx must have a value of 0 in order for the PDF to function
  ! correctly.  When F_x = 0, mu_x_1 = mu_x_2.  When the two PDF component means
  ! are equal to each other (and to the overall mean, <x>), the only value of
  ! Skx that can be represented is a value of 0.  In the equation for mixture
  ! fraction, when F_x is set to 0, but | Skx | > 0, mixt_frac will either have
  ! a value of 0 or 1, depending on whether Skx is positive or negative,
  ! respectively.
  !
  ! The value of F_x should be set as a function of Skx.  The value F_x should
  ! go toward 0 as | Skx | (or Skx^2) goes toward 0.  The value of F_x should
  ! go toward 1 as | Skx | (or Skx^2) goes to infinity.
  !
  !
  ! Notes:
  !
  ! When F_x = 0 (which can only happen when Skx = 0), mu_x_1 = mu_x_2, and
  ! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).  When these equations are
  ! substituted into the equations for sigma_x_1 and sigma_x_2, the result is
  ! sigma_x_1 = sigma_x_2 = sqrt( <x'^2> ).  This means that the distribution
  ! becomes a single Gaussian when F_x = 0 (and Skx = 0).  This happens
  ! regardless of the value of zeta_x.
  !
  ! Symmetry
  !
  !
  !
  ! The equations for the PDF component standard deviations can also be
  ! written as:
  !
  ! sigma_x_1 = coef_sigma_x_1 * sqrt( <x'^2> ); and
  !
  ! sigma_x_2 = coef_sigma_x_2 * sqrt( <x'^2> ); where
  !
  ! coef_sigma_x_1 = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
  !                        / ( ( zeta_x + 2 ) * mixt_frac ) ); and
  !
  ! coef_sigma_x_2 = sqrt( ( 1 - F_x )
  !                        / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ) ).
  !
  !=============================================================================
  function calc_mixture_fraction( Skx, F_x, zeta_x, sgn_wpxp ) &
  result( mixt_frac )

    ! Description:
    ! Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        thirty_six, & ! Constant(s)
        eighteen,   &
        twelve,     &
        six,        &
        four,       &
        three,      &
        two,        &
        one,        &
        one_half,   &
        zero,       &
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx,      & ! Skewness of x                                            [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      zeta_x,   & ! Parameter for the PDF component variances of x           [-]
      sgn_wpxp    ! Sign of the covariance of w and x                        [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      mixt_frac    ! Mixture fraction    [-]


    ! Calculate mixture fraction, which is the weight of the 1st PDF component.
    ! The 2nd PDF component has a weight of 1 - mixt_frac.
    if ( F_x > zero ) then

       mixt_frac &
       = ( four * F_x**3 &
           + eighteen * F_x &
             * ( zeta_x + one ) * ( one - F_x ) / ( zeta_x + two ) &
           + six * F_x**2 * ( one - F_x ) / ( zeta_x + two ) &
           + Skx**2 &
           - Skx * sgn_wpxp * sqrt( four * F_x**3 &
                                    + twelve * F_x**2 * ( one - F_x ) &
                                    + thirty_six * F_x &
                                      * ( zeta_x + one ) * ( one - F_x )**2 &
                                      / ( zeta_x + two )**2 &
                                    + Skx**2 ) ) &
         / ( two * F_x * ( F_x - three )**2 + two * Skx**2 )

    else ! F_x = 0

       if ( abs( Skx ) > zero ) then

          write(fstderr,*) "The value of F_x must be greater than 0 when " &
                           // "| Skx | > 0."
          write(fstderr,*) "F_x = ", F_x
          write(fstderr,*) "Skx = ", Skx
          stop

       else ! Skx = 0

          mixt_frac = ( zeta_x + one ) / ( zeta_x + two )

       endif ! | Skx | > 0

    endif ! F_x > 0 


    return

  end function calc_mixture_fraction

  !=============================================================================
  subroutine calc_setter_var_params( xm, xp2, Skx, sgn_wpxp,    & ! In
                                     F_x, zeta_x,               & ! In
                                     mu_x_1, mu_x_2, sigma_x_1, & ! Out
                                     sigma_x_2, mixt_frac )       ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two, & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,       & ! Mean of x (overall)                             [units vary]
      xp2,      & ! Variance of x (overall)                     [(units vary)^2]
      Skx,      & ! Skewness of x                                            [-]
      sgn_wpxp, & ! Sign of the covariance of w and x (overall)              [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      zeta_x      ! Parameter for the PDF component variances of x           [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,    & ! Mean of x (1st PDF component)                  [units vary]
      mu_x_2,    & ! Mean of x (2nd PDF component)                  [units vary]
      sigma_x_1, & ! Standard deviation of x (1st PDF component)    [units vary]
      sigma_x_2, & ! Standard deviation of x (2nd PDF component)    [units vary]
      mixt_frac    ! Mixture fraction                                        [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      coef_sigma_x_1, & !
      coef_sigma_x_2    !


    ! Calculate mixture fraction.
    mixt_frac = calc_mixture_fraction( Skx, F_x, zeta_x, sgn_wpxp )

    mu_x_1 = xm + sqrt( F_x * ( ( one - mixt_frac ) / mixt_frac ) * xp2 ) &
                  * sgn_wpxp

    mu_x_2 = xm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_x_1 - xm )

    sigma_x_1 &
    = sqrt( ( ( zeta_x + one ) * ( one - F_x ) ) &
            / ( ( zeta_x + two ) * mixt_frac ) * xp2 )

    sigma_x_2 &
    = sqrt( ( mixt_frac * sigma_x_1**2 ) &
            / ( ( one - mixt_frac ) * ( one + zeta_x ) ) )

    coef_sigma_x_1 = sqrt( ( zeta_x + one ) * ( one - F_x ) &
                           / ( ( zeta_x + two ) * mixt_frac ) )

    coef_sigma_x_2 = sqrt( ( one - F_x ) &
                           / ( ( zeta_x + two ) * ( one - mixt_frac ) ) )


    return

  end subroutine calc_setter_var_params

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
  ! ======================================================
  !
  ! In order to find equations for the four PDF parameters for each responding
  ! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
  ! for a responding variable here), four equations are needed.  These four
  ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
  ! integrating over the PDF.  Additionally, one more equation, which involves
  ! a tunable parameter F_x, and which is used to help control the spread of the
  ! PDF component means, is used in this equation set.  The four equations are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 ); and
  !
  ! mu_x_1 - <x> = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= F_x <= 1.
  !
  ! The resulting equations for the four PDF parameters are:
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
  !
  ! sigma_x_1^2
  ! = ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !             - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
  !     / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !   * <x'^2>; and
  !
  ! sigma_x_2^2 = ( ( 1 - F_x ) / ( 1 - mixt_frac ) ) * <x'^2>
  !               - ( mixt_frac / ( 1 - mixt_frac ) ) * sigma_x_1^2.
  !
  ! Special case:
  !
  ! When Skx = 0 and F_x = 0, the equations for sigma_x_1^2 and sigma_x_2^2 are
  ! both undefined.  The equations for sigma_x_1^2 and sigma_x_2^2 in this
  ! scenario can be derived by using the above equations for sigma_x_1^2 and
  ! sigma_x_2^2 and then setting Skx = 0.  The resulting equation for
  ! sigma_x_1^2 becomes:
  !
  ! sigma_x_1^2
  ! = ( ( - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
  !     / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !   * <x'^2>.
  !
  ! All of the terms in the numerator and denominator contain a sqrt( F_x ),
  ! so this equation can be rewritten as:
  !
  ! sigma_x_1^2
  ! = ( ( - ( 1 + mixt_frac ) * F_x + 3 * mixt_frac ) / ( 3 * mixt_frac ) )
  !   * <x'^2>.
  !
  ! Now setting F_x = 0, the equation becomes:
  !
  ! sigma_x_1^2 = ( ( 3 * mixt_frac ) / ( 3 * mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_1^2 = <x'^2>.
  !
  ! Substituting the equation for sigma_x_1^2 into the equation for sigma_x_2^2,
  ! and also setting F_x = 0, the equation for sigma_x_2^2 becomes:
  !
  ! sigma_x_2^2 = ( 1 / ( 1 - mixt_frac ) ) * <x'^2>
  !               - ( mixt_frac / ( 1 - mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_2^2
  ! = ( ( 1 / ( 1 - mixt_frac ) ) - ( mixt_frac / ( 1 - mixt_frac ) ) )
  !   * <x'^2>.
  !
  ! This equation becomes:
  !
  ! sigma_x_2^2 = ( ( 1 - mixt_frac ) / ( 1 - mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_2^2 = <x'^2>.
  !
  ! When F_x = 0, Skx must have a value of 0 in order for the PDF to function
  ! correctly.  When F_x = 0, mu_x_1 = mu_x_2.  When the two PDF component means
  ! are equal to each other (and to the overall mean, <x>), the only value of
  ! Skx that can be represented is a value of 0.  The equations that place
  ! limits on F_x for a responding variable (below) calculate the minimum
  ! allowable value of F_x to be greater than 0 when | Skx | > 0.
  !
  ! The value of F_x should be set as a function of Skx.  The value F_x should
  ! go toward 0 as | Skx | (or Skx^2) goes toward 0.  The value of F_x should
  ! go toward 1 as | Skx | (or Skx^2) goes to infinity.  However, the value of
  ! F_x must also be between the minimum and maximum allowable values of F_x for
  ! a responding variable (below).
  !
  ! Limits on F_x:
  !
  ! Since the PDF parameters for this variable need to work with the mixture
  ! fraction that has been provided by the setting variable, the method does
  ! not work for all values of F_x and Skx.  However, the limits of Skx and F_x
  ! can always be calculated.  The limits are based on keeping the values of
  ! sigma_x_1 and sigma_x_2 greater than or equal to 0.  The equation for
  ! keeping the value of sigma_x_1 greater than or equal to 0 is:
  !
  ! - ( 1 + mixt_frac ) * sqrt( F_x )^3 + 3 * mixt_frac * sqrt( F_x )
  ! + sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
  !
  ! The roots of sqrt( F_x ) can be solved by an equation of the form:
  !
  ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
  !
  ! where:
  !
  ! A = - ( 1 + mixt_frac );
  ! B = 0;
  ! C = 3 * mixt_frac; and
  ! D = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
  !
  ! The equation for keeping the value of sigma_x_2 greater than or equal to 0
  ! is:
  !
  ! - ( 2 - mixt_frac ) * sqrt( F_x )^3 + 3 * ( 1 - mixt_frac ) * sqrt( F_x )
  ! - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
  !
  ! The roots of sqrt( F_x ) can be solved by an equation of the form:
  !
  ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
  !
  ! where:
  !
  ! A = - ( 2 - mixt_frac );
  ! B = 0;
  ! C = 3 * ( 1 - mixt_frac ); and
  ! D = - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
  !
  ! After careful analysis of the above equations, the following properties
  ! emerge:
  !
  ! When Skx * sgn( <w'x'> ) >= 0,
  !    Skx^2 < 4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) )
  !    is required; and
  ! when Skx * sgn( <w'x'> ) < 0,
  !    Skx^2 < 4 * mixt_frac^2 / ( 1 - mixt_frac^2 ) is required.
  !
  ! Whenever Skx^2 exceeds these limits, Skx must be limited (preserving its
  ! sign) in order to have any value of F_x that will work in the equation set.
  !
  ! When Skx is found to be within the above limits (or after it has been
  ! limited to fall within its limits), the range of valid values of F_x can be
  ! found according to the following:
  !
  ! When Skx * sgn( <w'x'> ) >= 0:
  !
  !     When  4 * mixt_frac^2 / ( 1 - mixt_frac^2 )  <  Skx^2
  !           <  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the second equation (sigma_x_2
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the second
  !                                equation (sigma_x_2 based) and the only* root
  !                                of the first equation (sigma_x_1 based).
  !
  !     When  Skx^2  <=  4 * mixt_frac^2 / ( 1 - mixt_frac^2 ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the second equation (sigma_x_2
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the second
  !                                equation (sigma_x_2 based) and the largest
  !                                root of the first equation (sigma_x_1 based).
  !
  ! When Skx * sgn( <w'x'> ) < 0:
  !
  !     When  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) )
  !           <  Skx^2  <  4 * mixt_frac^2 / ( 1 - mixt_frac^2 ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the first equation (sigma_x_1
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the first
  !                                equation (sigma_x_1 based) and the only* root
  !                                of the second equation (sigma_x_2 based).
  !
  !     When  Skx^2
  !           <=  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the first equation (sigma_x_1
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the first
  !                                equation (sigma_x_1 based) and the largest
  !                                root of the second equation (sigma_x_2
  !                                based).
  !
  ! Here, "only* root" means the the only root that isn't a complex root.
  !
  ! The value of sqrt( F_x ) is also limited with a minimum of 0 and a maximum
  ! of 1.  The minimum and maximum allowable values of F_x are found by taking
  ! the square of the minimum and maximum allowable values of sqrt( F_x ),
  ! respectively.
  !
  ! Notes:
  !
  ! When F_x = 0 (which can only happen when Skx = 0), mu_x_1 = mu_x_2, and
  ! sigma_x_1 = sigma_x_2 = sqrt( <x'^2> ).  This means that the distribution
  ! becomes a single Gaussian when F_x = 0 (and Skx = 0).
  !
  ! The equations for the PDF component standard deviations can also be
  ! written as:
  !
  ! sigma_x_1^2 = coef_sigma_x_1^2 * <x'^2>; and
  !
  ! sigma_x_2^2 = coef_sigma_x_2^2 * <x'^2>; where
  !
  ! coef_sigma_x_1
  ! = sqrt( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !           / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !         - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
  !         + 1 ); and
  !
  ! coef_sigma_x_2
  ! = sqrt( ( 1 - F_x ) / ( 1 - mixt_frac )
  !         - mixt_frac / ( 1 - mixt_frac )
  !           * ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !                 / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !               - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
  !               + 1 ) ).
  !
  !=============================================================================
  subroutine calc_responder_params( xm, xp2, Skx, sgn_wpxp,         & ! In
                                    F_x, mixt_frac,                 & ! In
                                    mu_x_1, mu_x_2,                 & ! In
                                    sigma_x_1_sqd, sigma_x_2_sqd )    ! Out

    ! Description:
    ! Calculates the PDF component means and the PDF component standard
    ! deviations for a responding variable (a variable that is not used to set
    ! the mixture fraction).

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three, & ! Variable(s)
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,       & ! Mean of x (overall)                             [units vary]
      xp2,      & ! Variance of x (overall)                     [(units vary)^2]
      Skx,      & ! Skewness of x                                            [-]
      sgn_wpxp, & ! Sign of the covariance of w and x (overall)              [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      mixt_frac   ! Mixture fraction                                         [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)      [(units vary)^2]

    ! Local Variables
    real( kind = core_rknd ) :: &
      coef_sigma_x_1, & !
      coef_sigma_x_2    !


    if ( F_x > zero ) then

       mu_x_1 = xm + sqrt( F_x * ( ( one - mixt_frac ) / mixt_frac ) * xp2 ) &
                     * sgn_wpxp

       mu_x_2 = xm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_x_1 - xm )

       sigma_x_1_sqd &
       = ( ( sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp &
                   - ( one + mixt_frac ) * F_x**1.5 &
                   + three * mixt_frac * sqrt( F_x ) ) &
           / ( three * mixt_frac * sqrt( F_x ) ) ) &
         * xp2

       sigma_x_2_sqd = ( ( one - F_x ) / ( one - mixt_frac ) ) * xp2 &
                       - ( mixt_frac / ( one - mixt_frac ) ) * sigma_x_1_sqd

       coef_sigma_x_1 &
       = sqrt( ( sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp &
                 / ( three * mixt_frac * sqrt( F_x ) ) ) &
               - ( one + mixt_frac ) * F_x / ( three * mixt_frac ) &
               + one )

       coef_sigma_x_2 &
       = sqrt( ( one - F_x ) / ( one - mixt_frac ) &
               - mixt_frac / ( one - mixt_frac ) &
                 * ( ( sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp&
                       / ( three * mixt_frac * sqrt( F_x ) ) ) &
                     - ( one + mixt_frac ) * F_x / ( three * mixt_frac ) &
                     + one ) )

    else ! F_x = 0

       mu_x_1 = xm
       mu_x_2 = xm
       sigma_x_1_sqd = xp2
       sigma_x_2_sqd = xp2
       coef_sigma_x_1 = one
       coef_sigma_x_2 = one

    endif ! F_x > 0


    return

  end subroutine calc_responder_params

  !=============================================================================
  subroutine calc_limits_F_x_responder( mixt_frac, Skx, sgn_wpxp,  & ! In
                                        max_Skx2_pos_Skx_sgn_wpxp, & ! In
                                        max_Skx2_neg_Skx_sgn_wpxp, & ! In
                                        min_F, max_F )               ! Out

    ! Description:
    ! Calculates the minimum and maximum allowable values for F_x for a
    ! responding variable.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three, & ! Variable(s)
        two,   &
        one,   &
        zero

    use calc_roots, only: &
        cubic_solve    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, & ! Mixture fraction                 [-]
      Skx,       & ! Skewness of x                    [-]
      sgn_wpxp     ! Sign of covariance of w and x    [-]

    real( kind = core_rknd ), intent(in) :: &
      max_Skx2_pos_Skx_sgn_wpxp, &
      max_Skx2_neg_Skx_sgn_wpxp

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      min_F, & ! Minimum allowable value of F    [-]
      max_F    ! Maximum allowable value of F    [-]

    ! Local Variables
    real( kind = core_rknd ) &
      coef_A_1, & ! Coef. A in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_B_1, & ! Coef. B in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_C_1, & ! Coef. C in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_D_1, & ! Coef. D in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_A_2, & ! Coef. A in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_B_2, & ! Coef. B in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_C_2, & ! Coef. C in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_D_2    ! Coef. D in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]

    complex( kind = core_rknd ), dimension(3) :: &
      sqrt_F_roots_1, & ! Roots of sqrt(F) for sigma_x_1 term    [-]
      sqrt_F_roots_2    ! Roots of sqrt(F) for sigma_x_2 term    [-]

    real( kind = core_rknd ), dimension(3) :: &
      sqrt_F_roots_1_sorted, & ! Sorted roots of sqrt(F) for sigma_x_1 term  [-]
      sqrt_F_roots_2_sorted    ! Sorted roots of sqrt(F) for sigma_x_2 term  [-]

    real( kind = core_rknd ) :: &
      min_sqrt_F, & ! Minimum allowable value of sqrt(F)    [-]
      max_sqrt_F    ! Maximum allowable value of sqrt(F)    [-]

 
    ! Set up the coefficients in the equation for the limit of sqrt(F) based on
    ! the 1st PDF component standard deviation (sigma_x_1) being greater than or
    ! equal to 0.  This equation has the form:
    ! A * sqrt(F)^3 + B * sqrt(F)^2 + C * sqrt(F) + D = 0.
    coef_A_1 = -( one + mixt_frac )
    coef_B_1 = zero
    coef_C_1 = three * mixt_frac
    coef_D_1 = sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp

    ! Solve for the roots (values of sqrt(F)) that satisfy the above equation.
    sqrt_F_roots_1 = cubic_solve( coef_A_1, coef_B_1, coef_C_1, coef_D_1 )

    ! Set up the coefficients in the equation for the limit of sqrt(F) based on
    ! the 2nd PDF component standard deviation (sigma_x_2) being greater than or
    ! equal to 0.  This equation has the form:
    ! A * sqrt(F)^3 + B * sqrt(F)^2 + C * sqrt(F) + D = 0.
    coef_A_2 = -( two - mixt_frac )
    coef_B_2 = zero
    coef_C_2 = three * ( one - mixt_frac )
    coef_D_2 = -sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp

    ! Solve for the roots (values of sqrt(F)) that satisfy the above equation.
    sqrt_F_roots_2 = cubic_solve( coef_A_2, coef_B_2, coef_C_2, coef_D_2 )


    ! Find the minimum and maximum allowable values of sqrt(F) based on Skx
    ! and sgn( <w'x'> ).
    if ( Skx * sgn_wpxp >= zero ) then

       if ( Skx**2 > max_Skx2_neg_Skx_sgn_wpxp ) then

          sqrt_F_roots_2_sorted &
          = sort_roots( real( sqrt_F_roots_2, kind = core_rknd ) )

          min_sqrt_F = sqrt_F_roots_2_sorted(2)
          max_sqrt_F = min( real( sqrt_F_roots_1(1), kind = core_rknd ), &
                            sqrt_F_roots_2_sorted(3) )

       else ! Skx^2 <= max_Skx2_neg_Skx_sgn_wpxp

          sqrt_F_roots_1_sorted &
          = sort_roots( real( sqrt_F_roots_1, kind = core_rknd ) )
          sqrt_F_roots_2_sorted &
          = sort_roots( real( sqrt_F_roots_2, kind = core_rknd ) )

          min_sqrt_F = sqrt_F_roots_2_sorted(2)
          max_sqrt_F = min( sqrt_F_roots_1_sorted(3), sqrt_F_roots_2_sorted(3) )

       endif ! Skx**2 > max_Skx2_neg_Skx_sgn_wpxp

    else ! Skx * sgn( <w'x'> ) < 0 

       if ( Skx**2 > max_Skx2_pos_Skx_sgn_wpxp ) then

          sqrt_F_roots_1_sorted &
          = sort_roots( real( sqrt_F_roots_1, kind = core_rknd ) )

          min_sqrt_F = sqrt_F_roots_1_sorted(2)
          max_sqrt_F = min( real( sqrt_F_roots_2(1), kind = core_rknd ), &
                            sqrt_F_roots_1_sorted(3) )

       else ! Skx^2 <= max_Skx2_pos_Skx_sgn_wpxp

          sqrt_F_roots_1_sorted &
          = sort_roots( real( sqrt_F_roots_1, kind = core_rknd ) )
          sqrt_F_roots_2_sorted &
          = sort_roots( real( sqrt_F_roots_2, kind = core_rknd ) )

          min_sqrt_F = sqrt_F_roots_1_sorted(2)
          max_sqrt_F = min( sqrt_F_roots_1_sorted(3), sqrt_F_roots_2_sorted(3) )

       endif ! Skx**2 > max_Skx2_pos_Skx_sgn_wpxp

    endif ! Skx * sgn( <w'x'> ) >= 0


    ! The minimum and maximum are also limited by 0 and 1, respectively.
    min_sqrt_F = max( min_sqrt_F, zero )
    max_sqrt_F = min( max_sqrt_F, one )

    ! The minimum and maximum allowable values for F are the squares of the
    ! minimum and maximum allowable values for sqrt(F).
    min_F = min_sqrt_F**2
    max_F = max_sqrt_F**2


    return

  end subroutine calc_limits_F_x_responder

  !=============================================================================
  function sort_roots( roots ) &
  result ( roots_sorted )

    ! Description:
    ! Sorts roots from smallest to largest.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = core_rknd ), dimension(3), intent(in) :: &
      roots    ! Roots    [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: &
      roots_sorted    ! Roots sorted from smallest to largest    [-]


    if ( roots(1) <= roots(2) .and. roots(1) <= roots(3) ) then

       ! The value of roots(1) is the smallest root.
       roots_sorted(1) = roots(1)

       if ( roots(2) <= roots(3) ) then

          ! The value of roots(2) is the middle-valued root and the value of
          ! roots(3) is the largest root.
          roots_sorted(2) = roots(2)
          roots_sorted(3) = roots(3)

       else ! roots(3) < roots(2)

          ! The value of roots(3) is the middle-valued root and the value of
          ! roots(2) is the largest root.
          roots_sorted(2) = roots(3)
          roots_sorted(3) = roots(2)

       endif ! roots(2) <= roots(3)

    elseif ( roots(2) < roots(1) .and. roots(2) <= roots(3) ) then

       ! The value of roots(2) is the smallest root.
       roots_sorted(1) = roots(2)

       if ( roots(1) <= roots(3) ) then

          ! The value of roots(1) is the middle-valued root and the value of
          ! roots(3) is the largest root.
          roots_sorted(2) = roots(1)
          roots_sorted(3) = roots(3)

       else ! roots(3) < roots(1)

          ! The value of roots(3) is the middle-valued root and the value of
          ! roots(1) is the largest root.
          roots_sorted(2) = roots(3)
          roots_sorted(3) = roots(1)

       endif ! roots(1) <= roots(3)

    else ! roots(3) < roots(1) .and. roots(3) < roots(2)

       ! The value of roots(3) is the smallest root.
       roots_sorted(1) = roots(3)

       if ( roots(1) <= roots(2) ) then

          ! The value of roots(1) is the middle-valued root and the value of
          ! roots(2) is the largest root.
          roots_sorted(2) = roots(1)
          roots_sorted(3) = roots(2)

       else ! roots(2) < roots(1)

          ! The value of roots(2) is the middle-valued root and the value of
          ! roots(1) is the largest root.
          roots_sorted(2) = roots(2)
          roots_sorted(3) = roots(1)

       endif ! roots(1) <= roots(2)

    endif ! roots(1) <= roots(2) .and. roots(1) <= roots(3)


    return

  end function sort_roots

  !=============================================================================

end module new_pdf
