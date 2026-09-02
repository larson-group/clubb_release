"""New-hybrid PDF driver - port of new_hybrid_pdf_main.F90.

Description:
The portion of CLUBB's multivariate, two-component PDF that is the
multivariate, two-component normal PDF of vertical velocity (w), total water
mixing ratio (rt), liquid water potential temperature (thl), and optionally,
the west-east horizontal wind component (u), the south-north horizontal wind
component (v), and passive scalars (sclr).

Porting deviations:
- The Fortran ngrdcol/nz/sclr_dim loops are expressed as JAX array operations
  and `vmap` over passive scalars.
- The Fortran `clubb_params` argument is reduced to the already-computed
  tunable inputs needed by this routine.
- Fortran intent(inout) skewness arrays and the implicit-coefficients derived
  type are returned functionally in a dictionary.
- The Fortran error-stop for invalid mixture fraction is not repeated here;
  callers are expected to provide realizable F_w/Skw states.
- Commented-out Fortran interpolation code in `calc_F_w_zeta_w` is documented
  as a deviation; the active ADG1-like gamma form is implemented.
"""
from __future__ import annotations

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_pdf import _ssqrt
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import (
    calculate_w_params,
    calculate_responder_params,
    calculate_coef_wp4_implicit,
    calc_coef_wp2xp_implicit,
    calc_coefs_wpxp2_semiimpl,
    calc_coefs_wpxpyp_semiimpl,
)

# new_hybrid_pdf_main.F90 `use constants_clubb, only: max_mag_correlation_flux`
from clubb_jax.src.CLUBB_core.constants_clubb import (
    max_mag_correlation_flux as _MAX_MAG_CORRELATION_FLUX,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    l_explicit_turbulent_adv_wp3,
    l_explicit_turbulent_adv_wpxp,
    l_explicit_turbulent_adv_xpyp,
)
# The Fortran still declares lambda and exp_Skw_interp_factor for a commented
# interpolation formula.  The active code uses the ADG1-like gamma form below.


def calc_F_w_zeta_w(Skw, wprtp, wpthlp, upwp, vpwp, wp2, rtp2, thlp2, up2, vp2,
                    gamma_Skw_fnc, slope_coef_spread_DG_means_w,
                    pdf_component_stdev_factor_w, max_corr_w_sclr_sqd):
    """Calculates the values of F_w and zeta_w for w, which is the setter
    variable (which is the variable that sets the mixture fraction).

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
    ! all <w'x'> = 0, as well).
    !
    ! Tunable parameter gamma.
    ! When gamma goes to 0, the standard deviations of w in each PDF component
    ! become small, and the spread between the two PDF component means of w
    ! becomes large.  F_w goes to min_F_w.
    ! When gamma goes to 1, the standard deviations of w in each PDF component
    ! become large, and the spread between the two PDF component means of w
    ! becomes small.  F_w goes to max_F_w.
    """
    Skw = jnp.asarray(Skw, dtype=jnp.float64)
    wprtp = jnp.asarray(wprtp, dtype=jnp.float64); wpthlp = jnp.asarray(wpthlp, dtype=jnp.float64)
    upwp = jnp.asarray(upwp, dtype=jnp.float64); vpwp = jnp.asarray(vpwp, dtype=jnp.float64)
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); rtp2 = jnp.asarray(rtp2, dtype=jnp.float64)
    thlp2 = jnp.asarray(thlp2, dtype=jnp.float64); up2 = jnp.asarray(up2, dtype=jnp.float64)
    vp2 = jnp.asarray(vp2, dtype=jnp.float64)
    gamma = jnp.asarray(gamma_Skw_fnc, dtype=jnp.float64)
    pstdev = jnp.asarray(pdf_component_stdev_factor_w, dtype=jnp.float64)
    mcs = jnp.asarray(max_corr_w_sclr_sqd, dtype=jnp.float64)

    def _corr_sqd(cov, va, vb):
        d = va * vb
        return jnp.where(d > 0.0, cov ** 2 / jnp.where(d > 0.0, d, 1.0), 0.0)

    # Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
    # variables x that are Double Gaussian PDF responder variables.  This
    # includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
    # also calculated as part of the PDF, and they are included as well.
    # Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.

    # Calculate the squared value of the correlation of w and rt.
    corr_w_rt_sqd = _corr_sqd(wprtp, wp2, rtp2)
    # Calculate the squared value of the correlation of w and thl.
    corr_w_thl_sqd = _corr_sqd(wpthlp, wp2, thlp2)
    # Calculate the squared value of the correlation of u and w.
    corr_u_w_sqd = _corr_sqd(upwp, up2, wp2)
    # Calculate the squared value of the correlation of v and w.
    corr_v_w_sqd = _corr_sqd(vpwp, vp2, wp2)

    # Calculate the maximum value of corr_w_rt^2, corr_w_thl^2, corr_u_w^2, and
    # corr_v_w^2 in the calculation of max(corr_w_x^2).  Include the maximum
    # value of corr_w_sclr^2 (for all sclrs) when scalars are found.
    max_corr_w_x_sqd = jnp.maximum(
        jnp.maximum(jnp.maximum(corr_w_rt_sqd, corr_w_thl_sqd),
                    jnp.maximum(corr_u_w_sqd, corr_v_w_sqd)), mcs)
    max_corr_w_x_sqd = jnp.minimum(max_corr_w_x_sqd, _MAX_MAG_CORRELATION_FLUX ** 2)

    # Calculate min_F_w and set max_F_w to 1.
    min_F_w = jnp.where(jnp.abs(Skw) > 0.0,
                        jnp.maximum(max_corr_w_x_sqd, 1.0e-3), max_corr_w_x_sqd)
    max_F_w = jnp.ones_like(min_F_w)

    # F_w must have a value between min_F_w and max_F_w.
    # For now, use a formulation similar to what is used for ADG1.
    # Tunable parameter gamma.
    # When gamma goes to 0, the standard deviations of w in each PDF component
    # become small, and the spread between the two PDF component means of w
    # becomes large.  F_w goes to min_F_w.
    # When gamma goes to 1, the standard deviations of w in each PDF component
    # become large, and the spread between the two PDF component means of w
    # becomes small.  F_w goes to max_F_w.
    F_w = max_F_w - gamma * (max_F_w - min_F_w)

    # The value of zeta_w must be greater than -1.
    zeta_w_star = pstdev - 1.0
    # Make the PDF of w symmetric.  In other words, the PDF at a value of
    # positive skewness will look like a mirror image of the PDF at the
    # opposite value of negative skewness.
    zeta_w = jnp.where(Skw >= 0.0, zeta_w_star, -zeta_w_star / (zeta_w_star + 1.0))
    return F_w, zeta_w, min_F_w, max_F_w


def calc_responder_driver(xm, xp2, wpxp, wp2, mixt_frac, F_w, Skx):
    """This is the sub-driver for a responder variable.  The limits of the range
    of values of Skx that the PDF is able to represent are calculated, and
    Skx is clipped when it exceeds the upper or lower limit.  Then, the PDF
    parameters for responder variable x are calculated.
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64); wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64); F_w = jnp.asarray(F_w, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    omf = 1.0 - mf

    # Calculate the limits of Skx.
    F_pos = F_w > 0.0
    F_safe = jnp.where(F_pos, F_w, 1.0)
    wx = wp2 * xp2
    # Calculate the overall correlation of w and x.
    # The overall correlation of w and x is:
    #
    # corr_w_x = <w'x'> / sqrt( <w'^2> * <x'^2> ).
    #
    # When <w'^2> = 0 or <x'^2> = 0, <w'x'> = 0.  The correlation of w
    # and x is undefined, however, since <w'x'> = 0, Skx = 0.  Setting
    # corr_w_x = 0 in this scenario will set max_Skx = min_Skx = 0.
    corr_w_x = jnp.where(wx > 0.0, wpxp / _ssqrt(jnp.where(wx > 0.0, wx, 1.0)), 0.0)

    sqrt_F = _ssqrt(F_safe)
    F_3half = F_safe ** 1.5                                     # F_w**three_halves (gfortran pow)
    base_mfomf = _ssqrt(mf * omf)
    corr3 = corr_w_x * corr_w_x * corr_w_x                      # gfortran expands x**3 to x*x*x
    A = ((1.0 + mf) / base_mfomf * corr3 / F_3half
         - _ssqrt(mf / omf) * 3.0 * corr_w_x / sqrt_F)
    B = ((mf - 2.0) / base_mfomf * corr3 / F_3half
         + _ssqrt(omf / mf) * 3.0 * corr_w_x / sqrt_F)

    wpxp_nonneg = wpxp >= 0.0
    # When F_w = 0, <w'x'> = 0, and Skx = 0.
    min_Skx = jnp.where(F_pos, jnp.where(wpxp_nonneg, A, B), 0.0)
    max_Skx = jnp.where(F_pos, jnp.where(wpxp_nonneg, B, A), 0.0)

    # Limit Skx so that min_Skx <= Skx <= max_Skx.
    Skx_c = jnp.where(Skx > max_Skx, max_Skx, jnp.where(Skx < min_Skx, min_Skx, Skx))

    # Calculate the PDF parameters for responder variable x.
    mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, c1, c2 = calculate_responder_params(
        xm, xp2, Skx_c, wpxp, wp2, F_w, mf)
    return Skx_c, mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, c1, c2


def new_hybrid_pdf_driver(wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2,
                          Skw, wprtp, wpthlp, upwp, vpwp,
                          gamma_Skw_fnc, slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w,
                          Skrt, Skthl, Sku, Skv, pdf_implicit_coefs_terms,
                          sclrm=None, sclrp2=None, wpsclrp=None, Sksclr=None):
    """Calculate the PDF parameters for w (including mixture fraction), rt,
    theta-l, and optionally, u, v, and passive scalar variables.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    slope = jnp.asarray(slope_coef_spread_DG_means_w, dtype=jnp.float64)
    pstdev = jnp.asarray(pdf_component_stdev_factor_w, dtype=jnp.float64)
    if slope.ndim == 1:
        slope = slope[:, None]
    if pstdev.ndim == 1:
        pstdev = pstdev[:, None]

    # Calculate the maximum value of the square of the correlation of w and a
    # scalar when scalars are used.
    # Initialize max_corr_w_sclr_sqd to 0.  It needs to retain this value even
    # when sclr_dim = 0.
    if sclrm is not None and jnp.asarray(sclrm).shape[-1] > 0:
        sclrp2_a = jnp.asarray(sclrp2, dtype=jnp.float64)
        wpsclrp_a = jnp.asarray(wpsclrp, dtype=jnp.float64)
        d = wp2[..., None] * sclrp2_a
        ratio = jnp.where(d > 0.0, wpsclrp_a ** 2 / jnp.where(d > 0.0, d, 1.0), 0.0)
        max_corr_w_sclr_sqd = jnp.max(ratio, axis=-1)
    else:
        max_corr_w_sclr_sqd = jnp.zeros_like(wp2)

    # Calculate the values of PDF tunable parameters F_w and zeta_w.
    # Vertical velocity, w, will always be the setter variable.
    F_w, zeta_w, _min_F_w, _max_F_w = calc_F_w_zeta_w(
        Skw, wprtp, wpthlp, upwp, vpwp, wp2, rtp2, thlp2, up2, vp2,
        gamma_Skw_fnc, slope, pstdev, max_corr_w_sclr_sqd)

    # Calculate the PDF parameters, including mixture fraction, for the
    # setter variable, w.
    (mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac,
     coef_sigma_w_1_sqd, coef_sigma_w_2_sqd) = calculate_w_params(
        wm, wp2, Skw, F_w, zeta_w)

    out = dict(mu_w_1=mu_w_1, mu_w_2=mu_w_2,
               sigma_w_1_sqd=sigma_w_1 ** 2, sigma_w_2_sqd=sigma_w_2 ** 2,
               mixt_frac=mixt_frac, sigma_sqd_w=1.0 - F_w)

    responders = (('rt', rtm, rtp2, wprtp, Skrt), ('thl', thlm, thlp2, wpthlp, Skthl),
                  ('u', um, up2, upwp, Sku), ('v', vm, vp2, vpwp, Skv))
    coef_sigma_1_sqd = {}
    coef_sigma_2_sqd = {}
    for name, xm, xp2, wpxp, Skx in responders:
        # Calculate the PDF parameters for responder variable rt/thl/u/v.
        Skx_c, mu1, mu2, s1, s2, c1, c2 = calc_responder_driver(
            xm, xp2, wpxp, wp2, mixt_frac, F_w, Skx)
        out[f'Sk{name}'] = Skx_c
        out[f'mu_{name}_1'] = mu1; out[f'mu_{name}_2'] = mu2
        out[f'sigma_{name}_1_sqd'] = s1; out[f'sigma_{name}_2_sqd'] = s2
        coef_sigma_1_sqd[name] = c1
        coef_sigma_2_sqd[name] = c2

    zero = jnp.zeros_like(wp2)
    if sclrm is not None and jnp.asarray(sclrm).shape[-1] > 0:
        # Calculate the PDF parameters for responder variable sclr.
        sclrm_a = jnp.asarray(sclrm, dtype=jnp.float64)
        Sksclr_a = jnp.asarray(Sksclr, dtype=jnp.float64)
        (
            out['Sksclr'],
            out['mu_sclr_1'],
            out['mu_sclr_2'],
            out['sigma_sclr_1_sqd'],
            out['sigma_sclr_2_sqd'],
            coef_sigma_sclr_1_sqd,
            coef_sigma_sclr_2_sqd,
        ) = jax.vmap(
            lambda sclrm_s, sclrp2_s, wpsclrp_s, Sksclr_s: calc_responder_driver(
                sclrm_s, sclrp2_s, wpsclrp_s, wp2, mixt_frac, F_w, Sksclr_s,
            ),
            in_axes=(2, 2, 2, 2),
            out_axes=-1,
        )(sclrm_a, sclrp2_a, wpsclrp_a, Sksclr_a)
    else:
        coef_sigma_sclr_1_sqd = None
        coef_sigma_sclr_2_sqd = None

    if not l_explicit_turbulent_adv_wp3:
        # Turbulent advection of <w'^3> is being handled implicitly.

        # <w'^4> = coef_wp4_implicit * <w'^2>^2.
        coef_wp4_implicit = calculate_coef_wp4_implicit(
            mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)
    else:
        # Turbulent advection of <w'^3> is being handled explicitly.
        coef_wp4_implicit = zero

    if not l_explicit_turbulent_adv_wpxp:
        # Turbulent advection of <w'rt'>, <w'thl'>, <u'w'>, <v'w'>, and <w'sclr'>
        # are all being handled implicitly.

        # <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'>
        coef_wp2rtp_implicit = calc_coef_wp2xp_implicit(
            wp2, mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)
        # <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'>;
        # <w'^2 u'> = coef_wp2up_implicit * <u'w'>; and
        # <w'^2 v'> = coef_wp2vp_implicit * <v'w'>;
        # where each coef_wp2xp_implicit is the same as coef_wp2rtp_implicit.
        coef_wp2thlp_implicit = coef_wp2rtp_implicit
        coef_wp2up_implicit = coef_wp2rtp_implicit
        coef_wp2vp_implicit = coef_wp2rtp_implicit
        if coef_sigma_sclr_1_sqd is not None:
            # <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'>;
            # where each coef_wp2xp_implicit is the same as coef_wp2rtp_implicit.
            coef_wp2sclrp_implicit = jnp.broadcast_to(
                coef_wp2rtp_implicit[..., None], coef_sigma_sclr_1_sqd.shape)
            term_wp2sclrp_explicit = jnp.zeros_like(coef_wp2sclrp_implicit)
        else:
            coef_wp2sclrp_implicit = pdf_implicit_coefs_terms.coef_wp2sclrp_implicit
            term_wp2sclrp_explicit = pdf_implicit_coefs_terms.term_wp2sclrp_explicit
    else:
        # Turbulent advection of <w'rt'>, <w'thl'>, <u'w'>, <v'w'>, and <w'sclr'>
        # are all being handled explicitly.
        coef_wp2rtp_implicit = zero
        coef_wp2thlp_implicit = zero
        coef_wp2up_implicit = zero
        coef_wp2vp_implicit = zero
        if coef_sigma_sclr_1_sqd is not None:
            coef_wp2sclrp_implicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            term_wp2sclrp_explicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
        else:
            coef_wp2sclrp_implicit = pdf_implicit_coefs_terms.coef_wp2sclrp_implicit
            term_wp2sclrp_explicit = pdf_implicit_coefs_terms.term_wp2sclrp_explicit

    if not l_explicit_turbulent_adv_xpyp:
        # Turbulent advection of <rt'^2>, <thl'^2>, <rt'thl'>, <u'^2>, <v'^2>,
        # <sclr'^2>, <sclr'rt'>, and <sclr'thl'> are all being handled
        # semi-implicitly.

        # <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
        coef_wprtp2_implicit, term_wprtp2_explicit = calc_coefs_wpxp2_semiimpl(
            wp2, wprtp, mixt_frac, F_w, coef_sigma_1_sqd['rt'], coef_sigma_2_sqd['rt'])
        # <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wprtp2_explicit
        coef_wpthlp2_implicit, term_wpthlp2_explicit = calc_coefs_wpxp2_semiimpl(
            wp2, wpthlp, mixt_frac, F_w, coef_sigma_1_sqd['thl'], coef_sigma_2_sqd['thl'])
        # <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
        #               + term_wprtpthlp_explicit
        coef_wprtpthlp_implicit, term_wprtpthlp_explicit = calc_coefs_wpxpyp_semiimpl(
            wp2, wprtp, wpthlp, mixt_frac, F_w,
            coef_sigma_1_sqd['rt'], coef_sigma_2_sqd['rt'],
            coef_sigma_1_sqd['thl'], coef_sigma_2_sqd['thl'])
        # <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit
        coef_wpup2_implicit, term_wpup2_explicit = calc_coefs_wpxp2_semiimpl(
            wp2, upwp, mixt_frac, F_w, coef_sigma_1_sqd['u'], coef_sigma_2_sqd['u'])
        # <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit
        coef_wpvp2_implicit, term_wpvp2_explicit = calc_coefs_wpxp2_semiimpl(
            wp2, vpwp, mixt_frac, F_w, coef_sigma_1_sqd['v'], coef_sigma_2_sqd['v'])
        if coef_sigma_sclr_1_sqd is not None:
            # <w'sclr'^2> = coef_wpsclrp2_implicit * <sclr'^2>
            #               + term_wpsclrp2_explicit
            coef_wpsclrp2_implicit, term_wpsclrp2_explicit = jax.vmap(
                lambda coef_s1, coef_s2, wpsclrjp: calc_coefs_wpxp2_semiimpl(
                    wp2, wpsclrjp, mixt_frac, F_w, coef_s1, coef_s2,
                ),
                in_axes=(2, 2, 2),
                out_axes=-1,
            )(coef_sigma_sclr_1_sqd, coef_sigma_sclr_2_sqd, wpsclrp_a)
            # <w'rt'sclr'> = coef_wprtpsclrp_implicit * <sclr'rt'>
            #                + term_wprtpsclrp_explicit
            coef_wprtpsclrp_implicit, term_wprtpsclrp_explicit = jax.vmap(
                lambda coef_s1, coef_s2, wpsclrjp: calc_coefs_wpxpyp_semiimpl(
                    wp2, wprtp, wpsclrjp, mixt_frac, F_w,
                    coef_sigma_1_sqd['rt'], coef_sigma_2_sqd['rt'], coef_s1, coef_s2,
                ),
                in_axes=(2, 2, 2),
                out_axes=-1,
            )(coef_sigma_sclr_1_sqd, coef_sigma_sclr_2_sqd, wpsclrp_a)
            # <w'thl'sclr'> = coef_wpthlpsclrp_implicit * <sclr'thl'>
            #                 + term_wpthlpsclrp_explicit
            coef_wpthlpsclrp_implicit, term_wpthlpsclrp_explicit = jax.vmap(
                lambda coef_s1, coef_s2, wpsclrjp: calc_coefs_wpxpyp_semiimpl(
                    wp2, wpthlp, wpsclrjp, mixt_frac, F_w,
                    coef_sigma_1_sqd['thl'], coef_sigma_2_sqd['thl'], coef_s1, coef_s2,
                ),
                in_axes=(2, 2, 2),
                out_axes=-1,
            )(coef_sigma_sclr_1_sqd, coef_sigma_sclr_2_sqd, wpsclrp_a)
        else:
            coef_wpsclrp2_implicit = pdf_implicit_coefs_terms.coef_wpsclrp2_implicit
            term_wpsclrp2_explicit = pdf_implicit_coefs_terms.term_wpsclrp2_explicit
            coef_wprtpsclrp_implicit = pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit
            term_wprtpsclrp_explicit = pdf_implicit_coefs_terms.term_wprtpsclrp_explicit
            coef_wpthlpsclrp_implicit = pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit
            term_wpthlpsclrp_explicit = pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit
    else:
        # Turbulent advection of <rt'^2>, <thl'^2>, <rt'thl'>, <u'^2>, <v'^2>,
        # <sclr'^2>, <sclr'rt'>, and <sclr'thl'> are being handled explicitly.
        coef_wprtp2_implicit = zero
        term_wprtp2_explicit = zero
        coef_wpthlp2_implicit = zero
        term_wpthlp2_explicit = zero
        coef_wprtpthlp_implicit = zero
        term_wprtpthlp_explicit = zero
        coef_wpup2_implicit = zero
        term_wpup2_explicit = zero
        coef_wpvp2_implicit = zero
        term_wpvp2_explicit = zero
        if coef_sigma_sclr_1_sqd is not None:
            coef_wpsclrp2_implicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            term_wpsclrp2_explicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            coef_wprtpsclrp_implicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            term_wprtpsclrp_explicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            coef_wpthlpsclrp_implicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
            term_wpthlpsclrp_explicit = jnp.zeros_like(coef_sigma_sclr_1_sqd)
        else:
            coef_wpsclrp2_implicit = pdf_implicit_coefs_terms.coef_wpsclrp2_implicit
            term_wpsclrp2_explicit = pdf_implicit_coefs_terms.term_wpsclrp2_explicit
            coef_wprtpsclrp_implicit = pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit
            term_wprtpsclrp_explicit = pdf_implicit_coefs_terms.term_wprtpsclrp_explicit
            coef_wpthlpsclrp_implicit = pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit
            term_wpthlpsclrp_explicit = pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit

    # Pack the implicit coefficients and explicit terms into a single type
    # variable for output.
    out['pdf_implicit_coefs_terms'] = pdf_implicit_coefs_terms.replace(
        coef_wp4_implicit=coef_wp4_implicit,
        coef_wp2rtp_implicit=coef_wp2rtp_implicit,
        term_wp2rtp_explicit=zero,
        coef_wp2thlp_implicit=coef_wp2thlp_implicit,
        term_wp2thlp_explicit=zero,
        coef_wp2up_implicit=coef_wp2up_implicit,
        term_wp2up_explicit=zero,
        coef_wp2vp_implicit=coef_wp2vp_implicit,
        term_wp2vp_explicit=zero,
        coef_wprtp2_implicit=coef_wprtp2_implicit,
        term_wprtp2_explicit=term_wprtp2_explicit,
        coef_wpthlp2_implicit=coef_wpthlp2_implicit,
        term_wpthlp2_explicit=term_wpthlp2_explicit,
        coef_wprtpthlp_implicit=coef_wprtpthlp_implicit,
        term_wprtpthlp_explicit=term_wprtpthlp_explicit,
        coef_wpup2_implicit=coef_wpup2_implicit,
        term_wpup2_explicit=term_wpup2_explicit,
        coef_wpvp2_implicit=coef_wpvp2_implicit,
        term_wpvp2_explicit=term_wpvp2_explicit,
        coef_wp2sclrp_implicit=coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit=term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit=coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit=term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit=coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit=term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit=coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit=term_wpthlpsclrp_explicit,
    )

    # The calculation of skewness of rt, thl, u, v, and scalars is hard-wired
    # for use with the ADG1 code, which contains the variable sigma_sqd_w.
    # In order to use an equivalent expression for these skewnesses using the
    # new hybrid PDF (without doing more recoding), set the value of
    # sigma_sqd_w to 1 - F_w.
    return out
