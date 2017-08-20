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
            calc_mixture_fraction

  private :: calc_setter_var_params, & ! Procedure(s)
             calc_responder_params

  contains

  !=============================================================================
  subroutine new_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,   & ! In
                             Skw, Skrt, Skthl, wprtp, wpthlp,   & ! In
                             mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,  & ! Out
                             mu_thl_1, mu_thl_2, sigma_w_1_sqd, & ! Out
                             sigma_w_2_sqd, sigma_rt_1_sqd,     & ! Out
                             sigma_rt_2_sqd, sigma_thl_1_sqd,   & ! Out
                             sigma_thl_2_sqd, mixt_frac         ) ! Out
                             

    ! Description:
    ! Selects which variable is used to set the mixture fraction for the PDF
    ! ("the setter") and which variables are handled after that mixture fraction
    ! has been set ("the responders").  Traditionally, w has been used to set
    ! the PDF.  However, here, the variable with the greatest magnitude of
    ! skewness is used to set the PDF.

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
      wm,     & ! Mean of w (overall)                [m/s]
      rtm,    & ! Mean of rt (overall)               [kg/kg]
      thlm,   & ! Mean of thl (overall)              [K]
      wp2,    & ! Variance of w (overall)            [m^2/s^2]
      rtp2,   & ! Variance of rt (overall)           [kg^2/kg^2]
      thlp2,  & ! Variance of thl (overall)          [K^2]
      Skw,    & ! Skewness of w (overall)            [-]
      Skrt,   & ! Skewness of rt (overall)           [-]
      Skthl,  & ! Skewness of thl (overall)          [-]
      wprtp,  & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
      wpthlp    ! Covariance of w and thl (overall)  [(m/s)K]

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
      sigma_w_1,  & ! Standard deviation of w (1st PDF component)  [m/s]
      sigma_w_2,  & ! Standard deviation of w (2nd PDF component)  [m/s]
      sgn_wprtp,  & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
      sgn_wpthlp    ! Covariance of w and thl (overall)  [(m/s)K]

    real( kind = core_rknd ) :: &
      F_w,    & !
      zeta_w, & !
      F_rt,   & !
      F_thl     !


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

    F_w = 0.25_core_rknd
    zeta_w = zero
    F_rt = 0.2_core_rknd
    F_thl = 0.2_core_rknd

    call calc_setter_var_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                                 mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                 sigma_w_2, mixt_frac )       ! Out

    sigma_w_1_sqd = sigma_w_1**2
    sigma_w_2_sqd = sigma_w_2**2

    call calc_responder_params( rtm, rtp2, Skrt, sgn_wprtp,       & ! In
                                F_rt, mixt_frac,                  & ! In
                                mu_rt_1, mu_rt_2,                 & ! In
                                sigma_rt_1_sqd, sigma_rt_2_sqd )    ! Out

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
  ! In order to find equations for the five PDF parameters for the setting
  ! variable, which are mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, and mixt_frac,
  ! five equations are needed.  These five equations are the equations for <w>,
  ! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally,
  ! two more equations, which involve tunable parameters F_w and zeta_w, and
  ! which are used to help control the spread of the PDF component means and
  ! the size of the PDF component standard deviations compared to each other,
  ! are used in this equation set.  The five equations are:
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
  ! sigma_w_2 ...
  ! = sqrt( ( mixt_frac * sigma_w_1^2 )
  !         / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) ).
  !
  ! This method works for all values of F_w (where 0 <= F_w <= 1) and zeta_w
  ! (where zeta_w > -1).
  !
  ! Special case:
  !
  ! When Skw = 0 and F_w = 0, the equation for mixt_frac is undefined.  The PDF
  ! component means (mu_w_1 and mu_w_2) are already equal when F_w = 0,
  ! regardless of the value of mixt_frac.  In order to allow for sigma_w_1 to
  ! equal sigma_w_2 when zeta_w = 0 in this special case (allowing the PDF to
  ! return to a single Gaussian), the value of mixt_frac is simply set to 1/2.
  !
  !
  ! The equations for the PDF component standard deviations can also be
  ! written as:
  !
  ! sigma_w_1 = coef_sigma_w_1 * sqrt( <w'^2> ); and
  !
  ! sigma_w_2 = coef_sigma_w_2 * sqrt( <w'^2> ); where
  !
  ! coef_sigma_w_1 = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
  !                        / ( ( zeta_w + 2 ) * mixt_frac ) ); and
  !
  ! coef_sigma_w_2 = sqrt( ( 1 - F_w )
  !                        / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ) ).
  !
  !=============================================================================
  function calc_mixture_fraction( Skw, F_w, zeta_w ) &
  result( mixt_frac )

    ! Description:
    ! Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,     & ! Variable(s)
        three,    &
        two,      &
        one,      &
        one_half, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw,    & !
      F_w,    & !
      zeta_w    !

    ! Return Variable
    real( kind = core_rknd ) :: &
      mixt_frac    ! Mixture fraction    [-]


    ! Calculate mixture fraction, which is the weight of the 1st PDF component.
    ! The 2nd PDF component has weight 1 - mixt_frac.
    if ( abs( Skw ) > zero .or. F_w > zero ) then

       mixt_frac &
       = ( four * F_w**3 &
           + 18.0_core_rknd * F_w &
             * ( zeta_w + one ) * ( one - F_w ) / ( zeta_w + two ) &
           + 6.0_core_rknd * F_w**2 * ( one - F_w ) / ( zeta_w + two ) &
           + Skw**2 &
           - Skw * sqrt( four * F_w**3 &
                         + 12.0_core_rknd * F_w**2 * ( one - F_w ) &
                         + 36.0_core_rknd * F_w &
                           * ( zeta_w + one ) * ( one - F_w )**2 &
                           / ( zeta_w + two )**2 &
                         + Skw**2 ) ) &
         / ( two * F_w * ( F_w - three )**2 + two * Skw**2 )

    else ! Skw = 0 and F_w = 0

       mixt_frac = one_half

    endif ! abs( Skw ) > 0 or F_w > 0 


    return

  end function calc_mixture_fraction

  !=============================================================================
  subroutine calc_setter_var_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                                     mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                     sigma_w_2, mixt_frac )       ! Out

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
      wm,     & !
      wp2,    & !
      Skw,    & !
      F_w,    & !
      zeta_w    !

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_w_1,    & !
      mu_w_2,    & !
      sigma_w_1, & !
      sigma_w_2, & !
      mixt_frac    ! Mixture fraction

    ! Local Variables
    real( kind = core_rknd ) :: &
      coef_sigma_w_1, & !
      coef_sigma_w_2    !


    ! Calculate mixture fraction.
    mixt_frac = calc_mixture_fraction( Skw, F_w, zeta_w )

    mu_w_1 = wm + sqrt( F_w * ( ( one - mixt_frac ) / mixt_frac ) * wp2 )

    mu_w_2 = wm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_w_1 - wm )

    sigma_w_1 &
    = sqrt( ( ( zeta_w + one ) * ( one - F_w ) ) &
            / ( ( zeta_w + two ) * mixt_frac ) * wp2 )

    sigma_w_2 &
    = sqrt( ( mixt_frac * sigma_w_1**2 ) &
            / ( ( one - mixt_frac ) * ( one + zeta_w ) ) )

    coef_sigma_w_1 = sqrt( ( zeta_w + one ) * ( one - F_w ) &
                           / ( ( zeta_w + two ) * mixt_frac ) )

    coef_sigma_w_2 = sqrt( ( one - F_w ) &
                           / ( ( zeta_w + two ) * ( one - mixt_frac ) ) )

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
  ! Since the PDF parmeters for this variable need to work with the mixture
  ! fraction that has been provided by the setting variable, the method does
  ! not work for all values of F_x and Skx.  However, the limits of Skx and F_x
  ! can always be calculated.  For Skx, the magnitude of Skx can be limited.
  ! For F_x, a value can be chosen that is inside the limits.
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
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,       & !
      xp2,      & !
      Skx,      & !
      sgn_wpxp, & !
      F_x,      & !
      mixt_frac   !

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & !
      mu_x_2,        & !
      sigma_x_1_sqd, & !
      sigma_x_2_sqd    !

    ! Local Variables
    real( kind = core_rknd ) :: &
      coef_sigma_x_1, & !
      coef_sigma_x_2    !


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
              * ( ( sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp &
                    / ( three * mixt_frac * sqrt( F_x ) ) ) &
                  - ( one + mixt_frac ) * F_x / ( three * mixt_frac ) &
                  + one ) )


    return

  end subroutine calc_responder_params

  !=============================================================================

end module new_pdf
