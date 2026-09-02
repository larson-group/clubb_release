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

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def _spow(x, p):
    """Positive power with a finite gradient at x=0.

    JAX deviation: this is only a numerical helper for the Fortran expression
    `abs(Skx)**lambda`; x is assumed to be non-negative at the call site.
    """
    xp = jnp.where(x > 0.0, x, 1.0)
    return jnp.where(x > 0.0, xp ** p, 0.0)


def calc_F_x_zeta_x_setter(Skx, slope_coef_spread_DG_means_x, pdf_component_stdev_factor_x, lambda_):
    """Calculates the values of F_x and zeta_x for the setter variable (which is
    the variable that sets the mixture fraction).

    ! The value of F_x is calculated between 0 (min_F_x) and 1 (max_F_x).  The
    ! equation is:
    !
    ! F_x = max_F_x + ( min_F_x - max_F_x )
    !                 * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
    !
    ! which reduces to:
    !
    ! F_x = 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
    !
    ! where lambda > 0 and slope_coef_spread_DG_means_x > 0.  As |Skx| goes
    ! toward 0, the value of F_x goes toward 0, and as |Skx| becomes large, the
    ! value of F_x goes toward 1.  When slope_coef_spread_DG_means_x is small,
    ! the value of F_x tends toward 1, and when slope_coef_spread_DG_means_x is
    ! large, the value of F_x tends toward 0.  When lambda is small, the value
    ! of F_x is less dependent on Skx, and when lambda is large, the value of
    ! F_x is more dependent on Skx.
    !
    ! Mathematically, this equation will always produce a value of F_x that
    ! falls between min_F_x and max_F_x.  However, in order to prevent a value
    ! of F_x from being calculated outside the bounds of min_F_x and max_F_x
    ! owing to numerical underflow or loss of precision, this equation can be
    ! rewritten as:
    !
    ! F_x
    ! = min_F_x * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x }
    !   + max_F_x * ( 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x } ).
    !
    ! The value of zeta_x used to adjust the PDF component standard devations:
    !
    ! 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
    !              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
    !
    ! where zeta_x > -1.
    """
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    slope = jnp.asarray(slope_coef_spread_DG_means_x, dtype=jnp.float64)
    stdev_factor = jnp.asarray(pdf_component_stdev_factor_x, dtype=jnp.float64)
    lam = jnp.asarray(lambda_, dtype=jnp.float64)

    absS = jnp.abs(Skx)
    # Set min_F_x to 0 and max_F_x to 1 for the setter variable.
    min_F_x = jnp.where(absS > 0.0, 1.0e-3, 0.0)
    max_F_x = jnp.ones_like(Skx)
    # F_x must have a value between min_F_x and max_F_x.
    e = jnp.exp(-_spow(absS, lam) / slope)
    F_x = min_F_x * e + max_F_x * (1.0 - e)
    # The value of zeta_x must be greater than -1.
    zeta_x = (stdev_factor - 1.0) * jnp.ones_like(Skx)
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

    ! The value of F_x is calculated between min_F_x and max_F_x.  The equation
    ! is:
    !
    ! F_x = min_F_x + ( max_F_x - min_F_x )
    !                 * coef_spread_DG_means_x * exp_factor_x;
    !
    ! where 0 <= coef_spread_DG_means_x <= 1.  As coef_spread_DG_means_x
    ! goes toward 0, the value of F_x goes toward min_F_x, and as
    ! coef_spread_DG_means_x goes toward 1, the value of F_x goes toward
    ! max_F_x.  The exp_factor_x is a factor of the form 1 - exp{ }.  The range
    ! of values of exp_factor_x is 0 <= exp_factor_x <= 1.  Here, exp_factor_x
    ! is used to reduce the value of F_x under special conditions.
    !
    ! F_x = min_F_x * ( 1 - coef_spread_DG_means_x * exp_factor_x )
    !       + max_F_x * coef_spread_DG_means_x * exp_factor_x.
    """
    cse = jnp.asarray(coef_spread_DG_means_x) * jnp.asarray(exp_factor_x)
    # F_x must have a value between min_F_x and max_F_x.
    return jnp.asarray(min_F_x) * (1.0 - cse) + jnp.asarray(max_F_x) * cse


def calc_responder_var(xm, xp2, sgn_wpxp, mixt_frac, coef_spread_DG_means_x, exp_factor_x,
                       max_Skx2_pos_Skx_sgn_wpxp, max_Skx2_neg_Skx_sgn_wpxp, Skx):
    """This is the sub-driver for a responder variable.  The upper limits of the
    magnitude of Skx are calculated, and Skx is clipped when its magnitude
    exceeds the upper limits.  The limits of the F_x parameter are calculated,
    and the value of F_x is set within those limits.  Then, the PDF parameters
    for responder variable x are calculated.
    """
    Skx = jnp.asarray(Skx, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    mp = jnp.asarray(max_Skx2_pos_Skx_sgn_wpxp, dtype=jnp.float64)
    mn = jnp.asarray(max_Skx2_neg_Skx_sgn_wpxp, dtype=jnp.float64)

    # Calculate the upper limit of the magnitude of Skx.
    sk_pos = jnp.where(Skx >= 0.0, jnp.sqrt(0.99 * mp), -jnp.sqrt(0.99 * mp))
    sk_neg = jnp.where(Skx >= 0.0, jnp.sqrt(0.99 * mn), -jnp.sqrt(0.99 * mn))
    clip_pos = jnp.where(Skx ** 2 >= mp, sk_pos, Skx)
    clip_neg = jnp.where(Skx ** 2 >= mn, sk_neg, Skx)
    Skx_c = jnp.where(Skx * sgn >= 0.0, clip_pos, clip_neg)

    min_F, max_F = calc_limits_F_x_responder(mixt_frac, Skx_c, sgn, mp, mn)
    # F_x must have a value between min_F_x and max_F_x.
    F_x = calc_F_x_responder(coef_spread_DG_means_x, exp_factor_x, min_F, max_F)
    mu1, mu2, s1sq, s2sq, c1, c2 = calc_responder_params(xm, xp2, Skx_c, sgn, F_x, mixt_frac)
    return mu1, mu2, s1sq, s2sq, c1, c2, F_x, min_F, max_F, Skx_c


def new_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw,
                   wprtp, wpthlp, rtpthlp,
                   clubb_params, Skrt, Skthl,
                   pdf_implicit_coefs_terms):
    """Selects which variable is used to set the mixture fraction for the PDF
    ("the setter") and which variables are handled after that mixture fraction
    has been set ("the responders").  Traditionally, w has been used to set
    the PDF.

    JAX returns the Fortran output and input/output values by name, including
    the packed implicit_coefs_terms fields used by turbulent-advection solvers.
    """
    wm = jnp.asarray(wm, dtype=jnp.float64); rtm = jnp.asarray(rtm, dtype=jnp.float64); thlm = jnp.asarray(thlm, dtype=jnp.float64)
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); rtp2 = jnp.asarray(rtp2, dtype=jnp.float64); thlp2 = jnp.asarray(thlp2, dtype=jnp.float64)
    Skw = jnp.asarray(Skw, dtype=jnp.float64); rtpthlp = jnp.asarray(rtpthlp, dtype=jnp.float64)
    wprtp = jnp.asarray(wprtp, dtype=jnp.float64); wpthlp = jnp.asarray(wpthlp, dtype=jnp.float64)
    cp = jnp.asarray(clubb_params, dtype=jnp.float64)

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

    slope_w = cp[:, islope_coef_spread_DG_means_w:islope_coef_spread_DG_means_w + 1]
    stdev_w = cp[:, ipdf_component_stdev_factor_w:ipdf_component_stdev_factor_w + 1]
    coef_rt = cp[:, icoef_spread_DG_means_rt:icoef_spread_DG_means_rt + 1]
    coef_thl = cp[:, icoef_spread_DG_means_thl:icoef_spread_DG_means_thl + 1]

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
