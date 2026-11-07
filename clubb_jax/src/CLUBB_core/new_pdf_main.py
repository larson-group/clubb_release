"""JAX port of new_pdf_main.F90 - the "new" PDF top-level driver helpers.

Description:
The portion of CLUBB's multivariate, two-component PDF that is the
trivariate, two-component normal PDF of vertical velocity (w), total water
mixing ratio (rt), and liquid water potential temperature (thl).

Porting deviations:
- The Fortran nz/ngrdcol loop is expressed by JAX array operations over the
  caller-provided dimensions.
- Fortran intent(out) and intent(inout) arguments are returned in a dictionary;
  the `pdf_implicit_coefs_terms` derived type is updated functionally with
  `.replace(...)`.
- The helper routines are kept private by convention only.
- Fortran constants and flags are imported from the JAX constant modules.
"""
import jax.numpy as jnp


def _spow(x, p):
    """Positive power with a finite gradient at x=0.

    JAX deviation: this is only a numerical helper for the Fortran expression
    `abs(Skx)**lambda`; x is assumed to be non-negative at the call site.
    """
    xp = jnp.where(x > 0.0, x, 1.0)
    pow_x = jnp.where(x > 0.0, xp ** p, 0.0)
    return pow_x


def calc_F_x_zeta_x_setter(Skx, slope_coef_spread_DG_means_x, pdf_component_stdev_factor_x, lambda_):
    """Calculates the values of F_x and zeta_x for the setter variable (which is
    the variable that sets the mixture fraction).

    The value of F_x is calculated between 0 (min_F_x) and 1 (max_F_x).  The
    equation is:

    F_x = max_F_x + ( min_F_x - max_F_x )
                    * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };

    which reduces to:

    F_x = 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };

    where lambda > 0 and slope_coef_spread_DG_means_x > 0.  As |Skx| goes
    toward 0, the value of F_x goes toward 0, and as |Skx| becomes large, the
    value of F_x goes toward 1.  When slope_coef_spread_DG_means_x is small,
    the value of F_x tends toward 1, and when slope_coef_spread_DG_means_x is
    large, the value of F_x tends toward 0.  When lambda is small, the value
    of F_x is less dependent on Skx, and when lambda is large, the value of
    F_x is more dependent on Skx.

    Mathematically, this equation will always produce a value of F_x that
    falls between min_F_x and max_F_x.  However, in order to prevent a value
    of F_x from being calculated outside the bounds of min_F_x and max_F_x
    owing to numerical underflow or loss of precision, this equation can be
    rewritten as:

    F_x
    = min_F_x * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x }
      + max_F_x * ( 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x } ).

    The value of zeta_x used to adjust the PDF component standard devations:

    1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
                 / ( ( 1 - mixt_frac ) * sigma_x_2^2 );

    where zeta_x > -1.  The sign of zeta_x is used to easily determine if
    mixt_frac * sigma_x_1^2 is greater than ( 1 - mixt_frac ) * sigma_x_2^2
    (when zeta_x is positive), mixt_frac * sigma_x_1^2 is less than
    ( 1 - mixt_frac ) * sigma_x_2^2 (when zeta_x is negative), or
    mixt_frac * sigma_x_1^2 is equal to ( 1 - mixt_frac ) * sigma_x_2^2 (when
    zeta_x is 0).

    In order to allow for a tunable parameter that is the pure ratio of
    mixt_frac * sigma_x_1^2 to ( 1 - mixt_frac ) * sigma_x_2^2, zeta_x is
    related to the parameter pdf_component_stdev_factor_x, where:

    1 + zeta_x = pdf_component_stdev_factor_x.

    Returns:
        (F_x, zeta_x, min_F_x, max_F_x)

    Fortran declarations (calc_F_x_zeta_x_setter in src/CLUBB_core/new_pdf_main.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          Skx    ! Skewness of x (overall)              [-]

        real( kind = core_rknd ), intent(in) :: &
          slope_coef_spread_DG_means_x, & ! Slope coef: spread PDF comp. means x [-]
          pdf_component_stdev_factor_x, & ! Param.: PDF comp. standard devs.; x  [-]
          lambda                          ! Param. for Skx dependence            [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          F_x,     & ! Parameter for the spread of the PDF component means of x  [-]
          zeta_x,  & ! Parameter for the PDF component variances of x            [-]
          min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
          max_F_x    ! Maximum allowable value of parameter F_x                  [-]
    """

    absS = jnp.abs(Skx)
    # Set min_F_x to 0 and max_F_x to 1 for the setter variable.
    min_F_x = jnp.where(absS > 0.0, 1.0e-3, 0.0)
    max_F_x = jnp.ones_like(Skx)
    # F_x must have a value between min_F_x and max_F_x.
    e = jnp.exp(-_spow(absS, lambda_) / slope_coef_spread_DG_means_x)
    F_x = min_F_x * e + max_F_x * (1.0 - e)
    # The value of zeta_x must be greater than -1.
    zeta_x = (pdf_component_stdev_factor_x - 1.0) * jnp.ones_like(Skx)
    return F_x, zeta_x, min_F_x, max_F_x


from clubb_jax.src.CLUBB_core.constants_clubb import (  # noqa: E402
    rt_tol,
    thl_tol,
    max_mag_correlation,
)
from clubb_jax.src.CLUBB_core.model_flags import (  # noqa: E402
    l_explicit_turbulent_adv_wp3,
    l_explicit_turbulent_adv_wpxp,
    l_explicit_turbulent_adv_xpyp,
)
from clubb_jax.src.CLUBB_core.new_pdf import (  # noqa: E402
    calc_setter_var_params,
    calc_responder_params,
    calc_limits_F_x_responder,
    calc_coef_wp4_implicit,
    calc_coef_wpxp2_implicit,
    calc_coefs_wp2xp_semiimpl,
    calc_coefs_wpxpyp_semiimpl,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (  # noqa: E402
    islope_coef_spread_DG_means_w,
    ipdf_component_stdev_factor_w,
    icoef_spread_DG_means_rt,
    icoef_spread_DG_means_thl,
)

_LAMBDA_W = 0.5


def calc_F_x_responder(coef_spread_DG_means_x, exp_factor_x, min_F_x, max_F_x):
    """Calculates the value of F_x as a tunable function between min_F_x and
    max_F_x.

    The value of F_x is calculated between min_F_x and max_F_x.  The equation
    is:

    F_x = min_F_x + ( max_F_x - min_F_x )
                    * coef_spread_DG_means_x * exp_factor_x;

    where 0 <= coef_spread_DG_means_x <= 1.  As coef_spread_DG_means_x
    goes toward 0, the value of F_x goes toward min_F_x, and as
    coef_spread_DG_means_x goes toward 1, the value of F_x goes toward
    max_F_x.  The exp_factor_x is a factor of the form 1 - exp{ }.  The range
    of values of exp_factor_x is 0 <= exp_factor_x <= 1.  Here, exp_factor_x
    is used to reduce the value of F_x under special conditions.

    The equation for exp_factor_x depends on which responder variable is being
    solved for.  For rt, the F_rt equation is:

    F_rt = min_F_rt + ( max_F_rt - min_F_rt )
                      * coef_spread_DG_means_rt * exp_factor_rt;

    where exp_factor_rt is given by:

    exp_factor_rt = 1 - exp{ -0.2 * ( adj_corr_rt_thl + 1 )^5 };

    and where the adjusted (overall) correlation of rt and theta-l, denoted
    adj_corr_rt_thl, is given by:

    adj_corr_rt_thl = <rt'thl'> / ( sqrt( <rt'^2> ) * sqrt( <thl'^2> ) )
                      * sgn( <w'rt'> ) * sgn( <w'thl'> ).

    The range of values of the adjusted (overall) correlation of rt and
    theta-l is -1 <= adj_corr_rt_thl <= 1.  The adjusted (overall) correlation
    of rt and theta-l has the same absolute magnitude as the (overall)
    correlation of rt and theta-l, but the sign of the adjusted correlation is
    dependent on the agreement between the signs of <w'rt'> and <w'thl'>.
    When all three covariances are in sign agreement, which is when:

    sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = 1,

    the value of adj_corr_rt_thl is positive.  However, when the three
    covariances aren't consistent in sign, which is when:

    sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = -1,

    the value of adj_corr_rt_thl is negative.

    For theta-l, the F_thl equation is:

    F_thl = min_F_thl + ( max_F_thl - min_F_thl ) * coef_spread_DG_means_thl;

    where exp_factor_thl is set to 1 (1 - exp{-inf} = 1) because reducing
    F_thl through the use of exp_factor_thl is not desired for theta-l.

    The direction of the two PDF component means of x (mu_x_1 and mu_x_2) with
    respect to each other and the overall mean of x (<x>) is given by
    sgn( <w'x'> ).  When sgn( <w'x'> ) = 1, ( mu_x_1 - <x> ) >= 0 and
    ( mu_x_2 - <x> ) <= 0.  When sgn( <w'x'> ) = -1, ( mu_x_1 - <x> ) <= 0 and
    ( mu_x_2 - <x> ) >= 0.  This helps to promote realizability of the PDF
    component correlations (corr_w_x_1 and corr_w_x_2).

    The realizability of the PDF component correlations corr_w_rt_1,
    corr_w_rt_2, corr_w_thl_1, and corr_w_thl_2 is promoted by using
    sgn( <w'rt'> ) and sgn( <w'thl'> ), respectively to choose the direction
    of the PDF component means of rt and theta-l.  However, when the three
    covariances (<w'rt'>, <w'thl'>, and <rt'thl'>) aren't consistent in sign
    (as shown above), it becomes difficult to keep the PDF component
    correlations corr_rt_thl_1 and corr_rt_thl_2 in the realizable range.
    However, in this situation, the above equation for exp_factor_rt brings
    F_rt closer to min_F_rt, which promotes realizability for corr_rt_thl_1
    and corr_rt_thl_2.

    Mathematically, this equation will always produce a value of F_x that
    falls between min_F_x and max_F_x.  However, in order to prevent a value
    of F_x from being calculated outside the bounds of min_F_x and max_F_x
    owing to numerical underflow or loss of precision, this equation can be
    rewritten as:

    F_x = min_F_x * ( 1 - coef_spread_DG_means_x * exp_factor_x )
          + max_F_x * coef_spread_DG_means_x * exp_factor_x.

    Returns:
        F_x

    Fortran declarations (calc_F_x_responder in src/CLUBB_core/new_pdf_main.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz

        real( kind = core_rknd ), intent(in) :: &
          coef_spread_DG_means_x    ! Coef.: spread betw. PDF comp. means of x   [-]

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          exp_factor_x, & ! Factor of the form 1 - exp{}; reduces F_x  [-]
          min_F_x,      & ! Minimum allowable value of parameter F_x   [-]
          max_F_x         ! Maximum allowable value of parameter F_x   [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz) :: &
          F_x    ! Parameter for the spread between the PDF component means of x [-]
    """
    cse = coef_spread_DG_means_x * exp_factor_x
    # F_x must have a value between min_F_x and max_F_x.
    F_x = min_F_x * (1.0 - cse) + max_F_x * cse
    return F_x


def calc_responder_var(xm, xp2, sgn_wpxp, mixt_frac, coef_spread_DG_means_x, exp_factor_x,
                       max_Skx2_pos_Skx_sgn_wpxp, max_Skx2_neg_Skx_sgn_wpxp, Skx):
    """This is the sub-driver for a responder variable.  The upper limits of the
    magnitude of Skx are calculated, and Skx is clipped when its magnitude
    exceeds the upper limits.  The limits of the F_x parameter are calculated,
    and the value of F_x is set within those limits.  Then, the PDF parameters
    for responder variable x are calculated.

    Returns:
        (mu1, mu2, s1sq, s2sq, c1, c2, F_x, min_F, max_F, Skx_c)

    Fortran declarations (calc_responder_var in src/CLUBB_core/new_pdf_main.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          xm,           & ! Mean of x (overall)                       [units vary]
          xp2,          & ! Variance of x (overall)               [(units vary)^2]
          sgn_wpxp,     & ! Sign of the covariance of w and x                  [-]
          mixt_frac,    & ! Mixture fraction                                   [-]
          exp_factor_x    ! Factor of the form 1 - exp{}; reduces F_x          [-]

        real( kind = core_rknd ), intent(in) :: &
          coef_spread_DG_means_x    ! Coef.: spread betw. PDF comp. means of x   [-]

        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
          max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

        !--------------------------- Input/Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(inout) :: &
          Skx    ! Skewness of x (overall)              [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          mu_x_1,        & ! Mean of x (1st PDF component)        [units vary]
          mu_x_2,        & ! Mean of x (2nd PDF component)        [units vary]
          sigma_x_1_sqd, & ! Variance of x (1st PDF component)    [(units vary)^2]
          sigma_x_2_sqd    ! Variance of x (2nd PDF component)    [(units vary)^2]

        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>    [-]
          coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>    [-]

        ! Output only for recording statistics.
        real( kind = core_rknd ), dimension(nz), intent(out) :: &
          F_x,     & ! Param. for the spread betw. the PDF component means of x  [-]
          min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
          max_F_x    ! Maximum allowable value of parameter F_x                  [-]
    """

    # Calculate the upper limit of the magnitude of Skx.
    sk_pos = jnp.where(Skx >= 0.0, jnp.sqrt(0.99 * max_Skx2_pos_Skx_sgn_wpxp),
                       -jnp.sqrt(0.99 * max_Skx2_pos_Skx_sgn_wpxp))
    sk_neg = jnp.where(Skx >= 0.0, jnp.sqrt(0.99 * max_Skx2_neg_Skx_sgn_wpxp),
                       -jnp.sqrt(0.99 * max_Skx2_neg_Skx_sgn_wpxp))
    clip_pos = jnp.where(Skx ** 2 >= max_Skx2_pos_Skx_sgn_wpxp, sk_pos, Skx)
    clip_neg = jnp.where(Skx ** 2 >= max_Skx2_neg_Skx_sgn_wpxp, sk_neg, Skx)
    Skx_c = jnp.where(Skx * sgn_wpxp >= 0.0, clip_pos, clip_neg)

    min_F, max_F = calc_limits_F_x_responder(mixt_frac, Skx_c, sgn_wpxp,
                                             max_Skx2_pos_Skx_sgn_wpxp, max_Skx2_neg_Skx_sgn_wpxp)
    # F_x must have a value between min_F_x and max_F_x.
    F_x = calc_F_x_responder(coef_spread_DG_means_x, exp_factor_x, min_F, max_F)
    mu1, mu2, s1sq, s2sq, c1, c2 = calc_responder_params(xm, xp2, Skx_c, sgn_wpxp, F_x, mixt_frac)
    return mu1, mu2, s1sq, s2sq, c1, c2, F_x, min_F, max_F, Skx_c


def new_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw,
                   wprtp, wpthlp, rtpthlp,
                   clubb_params, Skrt, Skthl,
                   pdf_implicit_coefs_terms):
    """Selects which variable is used to set the mixture fraction for the PDF
    ("the setter") and which variables are handled after that mixture fraction
    has been set ("the responders").  Traditionally, w has been used to set
    the PDF.

    The JAX function returns the Fortran output and input/output values by
    name, including the packed implicit_coefs_terms fields used by the
    turbulent-advection solvers.

    Returns:
        dict with keys Skrt, Skthl, mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1,
          mu_thl_2, sigma_w_1_sqd, sigma_w_2_sqd, sigma_rt_1_sqd,
          sigma_rt_2_sqd, sigma_thl_1_sqd, sigma_thl_2_sqd, mixt_frac,
          pdf_implicit_coefs_terms

    Fortran declarations (new_pdf_driver in src/CLUBB_core/new_pdf_main.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
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

        real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
          clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

        !--------------------------- Input/Output Variables ---------------------------
        ! These variables are input/output because their values may be clipped.
        ! Otherwise, as long as it is not necessary to clip them, their values
        ! will stay the same.
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
          Skrt,  & ! Skewness of rt (overall)            [-]
          Skthl    ! Skewness of thl (overall)           [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
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

        !--------------------------- Input/Output Variables ---------------------------
        type(implicit_coefs_terms), intent(inout) :: &
          pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]
    """

    # ------------------------- Begin Code -------------------------

    # Calculate sgn( <w'rt'> ).
    sgn_wprtp = jnp.where(wprtp >= 0.0, 1.0, -1.0)
    # Calculate sgn( <w'thl'> ).
    sgn_wpthlp = jnp.where(wpthlp >= 0.0, 1.0, -1.0)
    # Sign of the variance of w (overall), which is always positive.
    sgn_wp2 = jnp.ones_like(Skw)

    # Calculate the adjusted (overall) correlation of rt and theta-l, and the
    # value of exp_factor_rt.
    has_var = (rtp2 >= rt_tol ** 2) & (thlp2 >= thl_tol ** 2)
    denom = jnp.sqrt(jnp.where(has_var, rtp2 * thlp2, 1.0))
    adj = jnp.clip(rtpthlp / denom * sgn_wprtp * sgn_wpthlp, -max_mag_correlation, max_mag_correlation)
    exp_factor_rt = jnp.where(has_var, 1.0 - jnp.exp(-0.2 * (adj + 1.0) ** 5), 1.0)
    # The value of F_thl is not reduced by exp_factor_thl.
    exp_factor_thl = jnp.ones_like(Skw)

    slope_w = clubb_params[:, islope_coef_spread_DG_means_w:islope_coef_spread_DG_means_w + 1]
    stdev_w = clubb_params[:, ipdf_component_stdev_factor_w:ipdf_component_stdev_factor_w + 1]
    coef_rt = clubb_params[:, icoef_spread_DG_means_rt:icoef_spread_DG_means_rt + 1]
    coef_thl = clubb_params[:, icoef_spread_DG_means_thl:icoef_spread_DG_means_thl + 1]

    # Vertical velocity, w, will always be the setter variable.
    F_w, zeta_w, min_F_w, max_F_w = calc_F_x_zeta_x_setter(Skw, slope_w, stdev_w, _LAMBDA_W)
    # Calculate the PDF parameters, including mixture fraction, for the
    # setter variable, w.
    (mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac,
     coef_sigma_w_1_sqd, coef_sigma_w_2_sqd) = calc_setter_var_params(
        wm, wp2, Skw, sgn_wp2, F_w, zeta_w)
    sigma_w_1_sqd = sigma_w_1 ** 2
    sigma_w_2_sqd = sigma_w_2 ** 2

    # Calculate the upper limit on the magnitude of skewness for responding
    # variables.
    mf = mixt_frac
    max_Skx2_pos = 4.0 * (1.0 - mf) ** 2 / (mf * (2.0 - mf))
    max_Skx2_neg = 4.0 * mf ** 2 / (1.0 - mf ** 2)

    # Calculate the PDF parameters for responder variable rt.
    (mu_rt_1, mu_rt_2, sigma_rt_1_sqd, sigma_rt_2_sqd,
     coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd,
     F_rt, min_F_rt, max_F_rt, Skrt) = calc_responder_var(
        rtm, rtp2, sgn_wprtp, mf, coef_rt, exp_factor_rt, max_Skx2_pos, max_Skx2_neg, Skrt)
    # Calculate the PDF parameters for responder variable thl.
    (mu_thl_1, mu_thl_2, sigma_thl_1_sqd, sigma_thl_2_sqd,
     coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd,
     F_thl, min_F_thl, max_F_thl, Skthl) = calc_responder_var(
        thlm, thlp2, sgn_wpthlp, mf, coef_thl, exp_factor_thl, max_Skx2_pos, max_Skx2_neg, Skthl)

    zero = jnp.zeros_like(wp2)

    if not l_explicit_turbulent_adv_wp3:
        # Turbulent advection of <w'^3> is being handled implicitly.

        # <w'^4> = coef_wp4_implicit * <w'^2>^2.
        coef_wp4_implicit = calc_coef_wp4_implicit(
            mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)
    else:
        # Turbulent advection of <w'^3> is being handled explicitly.
        coef_wp4_implicit = zero

    if not l_explicit_turbulent_adv_xpyp:
        # Turbulent advection of <rt'^2> and <thl'^2> is being handled
        # implicitly.  Turbulent advection of <rt'thl'> is being handled
        # semi-implicitly.

        # <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>
        coef_wprtp2_implicit = calc_coef_wpxp2_implicit(
            wp2, rtp2, wprtp, sgn_wprtp,
            mixt_frac, F_w, F_rt,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
            coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd)
        # <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>
        coef_wpthlp2_implicit = calc_coef_wpxp2_implicit(
            wp2, thlp2, wpthlp, sgn_wpthlp,
            mixt_frac, F_w, F_thl,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
            coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd)
        # <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
        #               + term_wprtpthlp_explicit
        coef_wprtpthlp_implicit, term_wprtpthlp_explicit = calc_coefs_wpxpyp_semiimpl(
            wp2, rtp2, thlp2, wprtp, wpthlp, sgn_wprtp, sgn_wpthlp,
            mixt_frac, F_w, F_rt, F_thl,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
            coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd,
            coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd)
    else:
        # Turbulent advection of <rt'^2>, <thl'^2>, and <rt'thl'> is being
        # handled explicitly.
        coef_wprtp2_implicit = zero
        coef_wpthlp2_implicit = zero
        coef_wprtpthlp_implicit = zero
        term_wprtpthlp_explicit = zero

    if not l_explicit_turbulent_adv_wpxp:
        # Turbulent advection of <w'rt'> and <w'thl'> is being handled
        # semi-implicitly.

        # <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
        coef_wp2rtp_implicit, term_wp2rtp_explicit = calc_coefs_wp2xp_semiimpl(
            wp2, rtp2, sgn_wprtp,
            mixt_frac, F_w, F_rt,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
            coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd)
        # <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
        coef_wp2thlp_implicit, term_wp2thlp_explicit = calc_coefs_wp2xp_semiimpl(
            wp2, thlp2, sgn_wpthlp,
            mixt_frac, F_w, F_thl,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
            coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd)
    else:
        # Turbulent advection of <w'rt'> and <w'thl'> is being handled
        # explicitly.
        coef_wp2rtp_implicit = zero
        term_wp2rtp_explicit = zero
        coef_wp2thlp_implicit = zero
        term_wp2thlp_explicit = zero

    # Pack the implicit coefficients and explicit terms into a single type
    # variable for output.
    pdf_implicit_coefs_terms = pdf_implicit_coefs_terms.replace(
        coef_wp4_implicit=coef_wp4_implicit,
        coef_wprtp2_implicit=coef_wprtp2_implicit,
        coef_wpthlp2_implicit=coef_wpthlp2_implicit,
        coef_wp2rtp_implicit=coef_wp2rtp_implicit,
        term_wp2rtp_explicit=term_wp2rtp_explicit,
        coef_wp2thlp_implicit=coef_wp2thlp_implicit,
        term_wp2thlp_explicit=term_wp2thlp_explicit,
        coef_wprtpthlp_implicit=coef_wprtpthlp_implicit,
        term_wprtpthlp_explicit=term_wprtpthlp_explicit,
    )

    return {
        "Skrt": Skrt,
        "Skthl": Skthl,
        "mu_w_1": mu_w_1,
        "mu_w_2": mu_w_2,
        "mu_rt_1": mu_rt_1,
        "mu_rt_2": mu_rt_2,
        "mu_thl_1": mu_thl_1,
        "mu_thl_2": mu_thl_2,
        "sigma_w_1_sqd": sigma_w_1_sqd,
        "sigma_w_2_sqd": sigma_w_2_sqd,
        "sigma_rt_1_sqd": sigma_rt_1_sqd,
        "sigma_rt_2_sqd": sigma_rt_2_sqd,
        "sigma_thl_1_sqd": sigma_thl_1_sqd,
        "sigma_thl_2_sqd": sigma_thl_2_sqd,
        "mixt_frac": mixt_frac,
        "pdf_implicit_coefs_terms": pdf_implicit_coefs_terms,
    }
