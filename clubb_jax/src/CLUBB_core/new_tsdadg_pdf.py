"""JAX port of new_tsdadg_pdf.F90 - the new TSDADG PDF closure helpers.

Description:
The new trivariate, skewness-dependent, analytic double Gaussian (TSDADG)
PDF.

Porting deviations:
- The Fortran nz loop in `tsdadg_pdf_driver` is expressed as JAX vectorized
  three-way selection.
- Fortran intent(out) subroutines return tuples.
- Fortran warnings printed before clipping negative component variances are not
  emitted; the JAX driver returns clipped values functionally.
- The misspelled Fortran name `calc_respnder_parameters` is preserved for API
  parity.

!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
! =========================================================================
!
! There are five PDF parameters that need to be calculated for the setting
! variable, which are mu_x_1 (the mean of x in the 1st PDF component), mu_x_2
! (the mean of x in the 2nd PDF component), sigma_x_1 (the standard deviation
! of x in the 1st PDF component), sigma_x_2 (the standard deviation of x in
! the 2nd PDF component), and mixt_frac (the mixture fraction).
!
! This method does NOT work for all values of L_x_1 and L_x_2 (where
! 0 <= L_x_1 <= 1 and 0 <= L_x_2 <= 1).  Only a subregion of this parameter
! space produces valid results.
!
!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
! ======================================================
!
! In order to find equations for the four PDF parameters for each responding
! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
! for a responding variable here), four equations are needed.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def calc_L_x_Skx_fnc(Skx, sgn_wpxp, small_l_x_1, small_l_x_2):
    """Calculates the values of big_L_x_1 and big_L_x_2 as functions of Skx.

    ! The values of L_x_1 and L_x_2 are calculated by skewness functions.
    ! Those functions are:
    !
    ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ).
    !
    ! When Skx * sgn( <w'x'> ) >= 0:
    ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );
    !
    ! otherwise, when Skx * sgn( <w'x'> ) < 0, switch l_x_1 and l_x_2:
    ! L_x_1 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ).
    """
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    l1 = jnp.asarray(small_l_x_1, dtype=jnp.float64)
    l2 = jnp.asarray(small_l_x_2, dtype=jnp.float64)

    factor = jnp.abs(Skx) / jnp.sqrt(4.0 + Skx ** 2)
    cond = Skx * sgn >= 0.0
    big_L_x_1 = jnp.where(cond, l1, l2) * factor
    big_L_x_2 = jnp.where(cond, l2, l1) * factor
    return big_L_x_1, big_L_x_2


from clubb_jax.src.CLUBB_core.constants_clubb import eps as _EPS  # noqa: E402


# grad-safe sqrt(max(x,0)) — the canonical tracer-toolkit helper.
from clubb_jax.src.CLUBB_core.pdf_utilities import _safe_sqrt as _ssqrt


def calc_setter_parameters(xm, xp2, Skx, sgn_wpxp, big_L_x_1, big_L_x_2):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    L1 = jnp.asarray(big_L_x_1, dtype=jnp.float64); L2 = jnp.asarray(big_L_x_2, dtype=jnp.float64)

    # Calculate the factors in the PDF component mean equations.
    t = Skx * sgn / jnp.sqrt(4.0 + Skx ** 2)
    factor_plus = 1.0 + t
    factor_minus = 1.0 - t
    # Calculate the normalized mean of x in the 1st PDF component.
    mu1n = L1 * jnp.sqrt(factor_plus / factor_minus) * sgn
    # Calculate the normalized mean of x in the 2nd PDF component.
    mu2n = -L2 * jnp.sqrt(factor_minus / factor_plus) * sgn
    sqrt_xp2 = _ssqrt(xp2)
    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu1n * sqrt_xp2
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu2n * sqrt_xp2

    # Calculate the mixture fraction.
    a1 = jnp.abs(mu1n) >= _EPS
    a2 = jnp.abs(mu2n) >= _EPS
    mu2n_safe = jnp.where(a2, mu2n, 1.0)
    mf_both = 1.0 / (1.0 + jnp.abs(mu1n / mu2n_safe))
    mf_1 = 1.0 / (1.0 + jnp.abs(mu1n / _EPS))
    mf_2 = 1.0 / (1.0 + jnp.abs(_EPS / mu2n_safe))
    mixt_frac = jnp.where(a1 & a2, mf_both,
                          jnp.where(a1 & (~a2), mf_1,
                                    jnp.where((~a1) & a2, mf_2, 0.5)))

    # Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    # term in the PDF component variance equations in order to prevent a
    # divide-by-zero error.
    mu1n_thresh = jnp.where(mu1n >= 0.0, jnp.maximum(mu1n, _EPS), jnp.minimum(mu1n, -_EPS))
    common = Skx / (3.0 * mixt_frac * mu1n_thresh) - mu1n ** 2 / 3.0 + mu2n ** 2 / 3.0
    base = 1.0 - mixt_frac * mu1n ** 2 - (1.0 - mixt_frac) * mu2n ** 2
    # Calculate the variance of x in the 1st PDF component.
    coef1 = base + (1.0 - mixt_frac) * common
    # Calculate the variance of x in the 2nd PDF component.
    coef2 = base - mixt_frac * common
    return mu_x_1, mu_x_2, coef1 * xp2, coef2 * xp2, mixt_frac, coef1, coef2


def calc_respnder_parameters(xm, xp2, Skx, sgn_wpxp, mixt_frac, big_L_x_1):
    """new_tsdadg_pdf.F90:calc_respnder_parameters (the misspelling "respnder" is deliberate — it mirrors the
    Fortran subroutine name exactly, like the preserved "derrived" typo in parameters_tunable).
    Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64); L1 = jnp.asarray(big_L_x_1, dtype=jnp.float64)

    # Calculate the factors in the PDF component mean equations.
    t = Skx * sgn / jnp.sqrt(4.0 + Skx ** 2)
    # Calculate the normalized mean of x in the 1st PDF component.
    mu1n = L1 * jnp.sqrt((1.0 + t) / (1.0 - t)) * sgn
    # Calculate the normalized mean of x in the 2nd PDF component.
    mu2n = -(mf / (1.0 - mf)) * mu1n
    sqrt_xp2 = _ssqrt(xp2)
    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu1n * sqrt_xp2
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu2n * sqrt_xp2

    # Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    # term in the PDF component variance equations in order to prevent a
    # divide-by-zero error.
    mu1n_thresh = jnp.where(mu1n >= 0.0, jnp.maximum(mu1n, _EPS), jnp.minimum(mu1n, -_EPS))
    common = Skx / (3.0 * mf * mu1n_thresh) - mu1n ** 2 / 3.0 + mu2n ** 2 / 3.0
    base = 1.0 - mf * mu1n ** 2 - (1.0 - mf) * mu2n ** 2
    # Calculate the variance of x in the 1st PDF component.
    coef1 = base + (1.0 - mf) * common
    # Calculate the variance of x in the 2nd PDF component.
    coef2 = base - mf * common
    return mu_x_1, mu_x_2, coef1 * xp2, coef2 * xp2, coef1, coef2


def tsdadg_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp):
    """Selects which variable is used to set the mixture fraction for the PDF
    ("the setter") and which variables are handled after that mixture fraction
    has been set ("the responders").  Traditionally, w has been used to set
    the PDF.  However, here, the variable with the greatest magnitude of
    skewness is used to set the PDF.
    """
    Skw = jnp.asarray(Skw, dtype=jnp.float64); Skrt = jnp.asarray(Skrt, dtype=jnp.float64)
    Skthl = jnp.asarray(Skthl, dtype=jnp.float64)
    wprtp = jnp.asarray(wprtp, dtype=jnp.float64); wpthlp = jnp.asarray(wpthlp, dtype=jnp.float64)

    # Calculate sgn( <w'rt'> ).
    sgn_wprtp = jnp.where(wprtp >= 0.0, 1.0, -1.0)
    # Calculate sgn( <w'thl'> ).
    sgn_wpthlp = jnp.where(wpthlp >= 0.0, 1.0, -1.0)
    # The sign of the variance of w is always positive.
    sgn_wp2 = jnp.ones_like(Skw)
    l1, l2 = 0.75, 0.5

    bLw1, bLw2 = calc_L_x_Skx_fnc(Skw, sgn_wp2, l1, l2)
    bLrt1, bLrt2 = calc_L_x_Skx_fnc(Skrt, sgn_wprtp, l1, l2)
    bLthl1, bLthl2 = calc_L_x_Skx_fnc(Skthl, sgn_wpthlp, l1, l2)

    aw, art, athl = jnp.abs(Skw), jnp.abs(Skrt), jnp.abs(Skthl)
    # The variable with the greatest magnitude of skewness will be the setter
    # variable and the other variables will be responder variables.
    w_set = (aw >= art) & (aw >= athl)
    rt_set = (art > aw) & (art >= athl)
    thl_set = (~w_set) & (~rt_set)

    # The variable w has the greatest magnitude of skewness.
    sw = calc_setter_parameters(wm, wp2, Skw, sgn_wp2, bLw1, bLw2)
    # The variable rt has the greatest magnitude of skewness.
    srt = calc_setter_parameters(rtm, rtp2, Skrt, sgn_wprtp, bLrt1, bLrt2)
    # The variable thl has the greatest magnitude of skewness.
    sthl = calc_setter_parameters(thlm, thlp2, Skthl, sgn_wpthlp, bLthl1, bLthl2)
    mixt_frac = jnp.where(w_set, sw[4], jnp.where(rt_set, srt[4], sthl[4]))

    rw = calc_respnder_parameters(wm, wp2, Skw, sgn_wp2, mixt_frac, bLw1)
    rrt = calc_respnder_parameters(rtm, rtp2, Skrt, sgn_wprtp, mixt_frac, bLrt1)
    rthl = calc_respnder_parameters(thlm, thlp2, Skthl, sgn_wpthlp, mixt_frac, bLthl1)

    def _sel(mask, setter, responder, i):
        return jnp.where(mask, setter[i], responder[i])

    mu_w_1 = _sel(w_set, sw, rw, 0); mu_w_2 = _sel(w_set, sw, rw, 1)
    sig_w_1 = _sel(w_set, sw, rw, 2); sig_w_2 = _sel(w_set, sw, rw, 3)
    mu_rt_1 = _sel(rt_set, srt, rrt, 0); mu_rt_2 = _sel(rt_set, srt, rrt, 1)
    sig_rt_1 = _sel(rt_set, srt, rrt, 2); sig_rt_2 = _sel(rt_set, srt, rrt, 3)
    mu_thl_1 = _sel(thl_set, sthl, rthl, 0); mu_thl_2 = _sel(thl_set, sthl, rthl, 1)
    sig_thl_1 = _sel(thl_set, sthl, rthl, 2); sig_thl_2 = _sel(thl_set, sthl, rthl, 3)

    # Negative variances are clipped to 0.  The Fortran also writes warnings to
    # stderr before clipping.
    clip = lambda s: jnp.where(s < 0.0, 0.0, s)
    return (mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2,
            clip(sig_w_1), clip(sig_w_2), clip(sig_rt_1), clip(sig_rt_2),
            clip(sig_thl_1), clip(sig_thl_2), mixt_frac)
