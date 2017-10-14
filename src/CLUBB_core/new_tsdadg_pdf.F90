! $Id$
!===============================================================================
module new_tsdadg_pdf

  ! Description:
  ! The new trivariate, skewness-dependent, analytic double Gaussian (TSDADG)
  ! PDF.

  implicit none

  public :: calc_setter_parameters, & ! Procedure(s)
            calc_L_x_Skx_fnc

  private  ! default scope

  contains

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =========================================================================
  !
  ! There are five PDF parameters that need to be calculated for the setting
  ! variable, which are mu_x_1 (the mean of x in the 1st PDF component), mu_x_2
  ! (the mean of x in the 2nd PDF component), sigma_x_1 (the standard deviation
  ! of x in the 1st PDF component), sigma_x_2 (the standard deviation of x in
  ! the 2nd PDF component), and mixt_frac (the mixture fraction).  In order to
  ! solve for these five parameters, five equations are needed.  These five
  ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
  ! integrating over the PDF.  Additionally, two more equations, which involve
  ! tunable parameters L_x_1 and L_x_2, and which are used to help control the
  ! spread of the PDF component means in the 1st PDF component and the 2nd PDF
  ! component, respectively, are used in this equation set.  The five equations
  ! are:
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
  ! mu_x_1 - <x> = L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
  !
  ! mu_x_2 - <x> = -L_x_2
  !                 * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                         / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                 * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= L_x_1 <= 1, 0 <= L_x_2 <= 1, Skx is the skewness of x, such that
  ! Skx = <x'^3> / <x'^2>^(3/2), and sgn( <w'x'> ) is the sign of <w'x'>, such
  ! that:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The resulting equations for the five PDF parameters are:
  !
  ! mu_x_1 = <x> + L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - L_x_2
  !                * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! mixt_frac = 1 / ( 1 + abs( mu_x_1_nrmlized / mu_x_2_nrmlized ) );
  !
  ! sigma_x_1 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     + ( 1 - mixt_frac )
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     - mixt_frac
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); where
  !
  ! mu_x_1_nrmlized = ( mu_x_1 - <x> ) / sqrt( <x'^2> ); and
  !
  ! mu_x_2_nrmlized = ( mu_x_2 - <x> ) / sqrt( <x'^2> ).
  !
  !
  ! Notes:
  !
  ! This method does NOT work for all values of L_x_1 and L_x_2 (where
  ! 0 <= L_x_1 <= 1 and 0 <= L_x_2 <= 1).  Only a subregion of this parameter
  ! space produces valid results.
  !
  ! When both L_x_1 = 0 and L_x_2 = 0, mu_x_1 = mu_x_2 = <x> (which can only
  ! happen when Skx = 0).  In this scenario, the above equations for mixt_frac,
  ! sigma_x_1, and sigma_x_2 are all undefined.  In this special case, the
  ! distribution reduces to a single Gaussian, so the following values are set:
  ! mixt_frac = 1/2, sigma_x_1 = sqrt( <x'^2> ), and sigma_x_2 = sqrt( <x'^2> ).
  !
  !
  ! Tunable parameters:
  !
  ! The parameter L_x_1 controls the 1st PDF component mean while L_x_2 controls
  ! the 2nd PDF component mean.  The equations involving the tunable parameters
  ! L_x_1 and L_x_2 (the mu_x_1 and mu_x_2 equations) are based on the values of
  ! mu_x_1 and mu_x_2 when sigma_x_1 = sigma_x_2 = 0.  In this scenario, the
  ! equation for mixture fraction reduces to:
  !
  ! mixt_frac = (1/2) * ( 1 +/- Skx / sqrt( 4 + Skx^2 ) ).
  !
  ! The +/- is dependent on the sign of ( mu_x_1 - <x> ) vs. ( mu_x_2 - <x> ).
  ! This is dependent on sgn( <w'x'> ), and the mixture fraction equation is
  ! written as:
  !
  ! mixt_frac = (1/2) * ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
  !
  ! Meanwhile, the equation for 1 - mixt_frac is:
  !
  ! 1 - mixt_frac = (1/2) * ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
  !
  ! When sigma_x_1 = sigma_x_2 = 0, the equations for mu_x_1 and mu_x_2 are:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 - mixt_frac ) / mixt_frac ) * sqrt( <x'^2> )
  !                * sgn( <w'x'> ); and
  !
  ! mu_x_2 = <x> - sqrt( mixt_frac / ( 1 - mixt_frac ) ) * sqrt( <x'^2> )
  !                * sgn( <w'x'> ).
  !
  ! Substituting the equations for mixt_frac and 1 - mixt_frac into the
  ! equations for mu_x_1 and mu_x_2 (when sigma_x_1 = sigma_x_2 = 0), the
  ! equations for mu_x_1 and mu_x_2 become:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                      / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
  !
  ! mu_x_2 = <x> - sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                      / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ).
  !
  ! These equations represent the maximum deviation of mu_x_1 and mu_x_2 from
  ! the overall mean, <x>.  The range of parameters of L_x_i is 0 <= L_x_i <= 1.
  ! When L_x_1 = L_x_2 = 0, the value of mu_x_1 = mu_x_2 = <x> (and the
  ! distribution becomes a single Gaussian).  When L_x_i = 1, the value of
  ! mu_x_i - <x> is at its maximum magnitude.
  !
  ! The values of L_x_1 and L_x_2 are also calculated by skewness functions.
  ! Those functions are:
  !
  ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
  ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );
  !
  ! where both l_x_1 and l_x_2 are tunable parameters.
  !
  ! As previously stated, this method does not work for all combinations of
  ! L_x_1 and L_x_2, but rather only for a subregion of parameter space.  This
  ! applies to l_x_1 and l_x_2, as well.  The recommended values are l_x_1 > 2/3
  ! and 0 < l_x_2 < 1.
  !
  !
  ! Equations for PDF component standard deviations:
  !
  ! The equations for the PDF component standard deviations can also be written
  ! as:
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      + ( 1 - mixt_frac )
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ); and
  !
  ! coef_sigma_x_2_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      - mixt_frac
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ).
  !
  ! The above equations can be substituted into an equation for a variable that
  ! has been derived by integrating over the PDF.  Many variables like this are
  ! used in parts of the predictive equation set.  These substitutions allow
  ! some terms to solved implicitly or semi-implicitly in the predictive
  ! equations.
  !
  !
  !=============================================================================
  subroutine calc_setter_parameters( xm, xp2, Skx, sgn_wpxp,    & ! In
                                     big_L_x_1, big_L_x_2,      & ! In
                                     mu_x_1, mu_x_2, sigma_x_1, & ! Out
                                     sigma_x_2, mixt_frac,      & ! Out
                                     coef_sigma_x_1_sqd,        & ! Out
                                     coef_sigma_x_2_sqd         ) ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,  & ! Variable(s)
        three, &
        one,   &
        eps

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean of x (overall)                            [units vary]
      xp2,       & ! Variance of x (overall)                    [(units vary)^2]
      Skx,       & ! Skewness of x                                           [-]
      sgn_wpxp,  & ! Sign of the covariance of w and x (overall)             [-]
      big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
      big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,    & ! Mean of x (1st PDF component)                  [units vary]
      mu_x_2,    & ! Mean of x (2nd PDF component)                  [units vary]
      sigma_x_1, & ! Standard deviation of x (1st PDF component)    [units vary]
      sigma_x_2, & ! Standard deviation of x (2nd PDF component)    [units vary]
      mixt_frac    ! Mixture fraction                                        [-]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      mu_x_1_nrmlized, & ! Normalized mean of x (1st PDF component)          [-]
      mu_x_2_nrmlized    ! Normalized mean of x (2nd PDF component)          [-]

    real( kind = core_rknd ) :: &
      factor_plus,               &
      factor_minus,              &
      sqrt_factor_plus_ov_minus, &
      sqrt_factor_minus_ov_plus


    ! Calculate the factors in the PDF component mean equations.
    factor_plus = one + Skx * sgn_wpxp / sqrt( four + Skx**2 )

    factor_minus = one - Skx * sgn_wpxp / sqrt( four + Skx**2 )

    sqrt_factor_plus_ov_minus = sqrt( factor_plus / factor_minus )
    sqrt_factor_minus_ov_plus = sqrt( factor_minus / factor_plus )

    ! Calculate the normalized mean of x in the 1st PDF component.
    mu_x_1_nrmlized = big_L_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp

    ! Calculate the normalized mean of x in the 2nd PDF component.
    mu_x_2_nrmlized = -big_L_x_2 * sqrt_factor_minus_ov_plus * sgn_wpxp

    ! Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu_x_1_nrmlized * xp2

    ! Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu_x_2_nrmlized * xp2

    ! Calculate the mixture fraction.
    mixt_frac = one / ( one + abs( max( mu_x_1_nrmlized, eps ) &
                                   / max( mu_x_2_nrmlized, eps ) ) )

    ! Calculate the standard deviation of x in the 1st PDF component.
    coef_sigma_x_1_sqd &
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      + ( one - mixt_frac ) &
        * ( Skx / ( three * mixt_frac * max( mu_x_1_nrmlized, eps ) ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_1 = sqrt( coef_sigma_x_1_sqd * xp2 )

    ! Calculate the standard deviation of x in the 2nd PDF component.
    coef_sigma_x_2_sqd & 
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      - mixt_frac &
        * ( Skx / ( three * mixt_frac * max( mu_x_1_nrmlized, eps ) ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_2 = sqrt( coef_sigma_x_2_sqd * xp2 )


    return

  end subroutine calc_setter_parameters

  !=============================================================================
  subroutine calc_L_x_Skx_fnc( Skx, small_l_x_1, small_l_x_2, & ! In
                               big_L_x_1, big_L_x_2           ) ! Out

    ! Description:
    ! Calculates the values of big_L_x_1 and big_L_x_2 as functions of Skx.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx,         & ! Skewness of x (overall)                          [-]
      small_l_x_1, & ! Param. for the spread of the 1st PDF comp. mean of x  [-]
      small_l_x_2    ! Param. for the spread of the 2nd PDF comp. mean of x  [-]

    ! Output Variable
    real( kind = core_rknd ), intent(out) :: &
      big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
      big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]

    ! Local Variable
    real( kind = core_rknd ) :: &
      factor_x


    factor_x = abs( Skx ) / sqrt( 4.0_core_rknd + Skx**2 )

    big_L_x_1 = small_l_x_1 * factor_x
    big_L_x_2 = small_l_x_2 * factor_x


    return

  end subroutine calc_L_x_Skx_fnc

  !=============================================================================

end module new_tsdadg_pdf
