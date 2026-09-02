"""JAX port of Microphys/mixed_moment_PDF_integrals.F90 — the all-mixed-moment hydrometeor-covariance machinery.

FULL port (all 8 functions). The analytic ("evaluated form") moment integrals, from the ground up:
  * univar_N / univar_L                    — univariate normal / lognormal central moments (PDF component)
  * bivar_NL_int_PDF_comp_all_MM           — bivariate normal-lognormal mixed moment (component)
  * bivar_NL_x_hm_all_MM_comp_eq           — per-component <x'^a hm'^b> (4-branch, jnp.where-selected)
  * xp_a_hmpb_integrals_all_MM             — full 2-component <x'^a hm'^b>
  * xphmp_integral_covar / hmxphmyp_integral_covar — streamlined <x'hm'> / <hmx'hmy'> covariances
  * hydrometeor_mixed_moments              — top driver: <rt'hm'>, <thl'hm'>, <w'^2 hm'>, <hmx'hmy'>

Integer exponents are STATIC in the Fortran → the finite sums unroll at trace time; the 4-way σ-vs-tol
dispatches are ported as `jnp.where`, so everything is differentiable w.r.t. the float moment inputs. Pure jnp.

Validated in `tests/test_mixed_moment_pdf_integrals.py` (closed-form binomial/tilting expansions + Monte-Carlo +
branch-selection vs independent NumPy) and `tests/test_hydrometeor_mixed_moments.py` (the driver vs a literal
per-level/per-hydrometeor Fortran-loop transcription).
"""
import math

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.pdf_utilities import (
    compute_mean_binormal, calc_corr_rt_x, calc_corr_thl_x)


def univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp):
    """Central moment of a normally-distributed variable within a PDF component.

    Evaluates  INT(-inf:inf) (x - <x>)^a P_N_i(x) dx
      = SUM(p=0:floor(a/2)) [ a! / ((a-2p)! p!) ] (½ σ_i²)^p (μ_i - <x>)^(a-2p)
    (mixed_moment_PDF_integrals.F90:univar_N_int_PDF_comp_all_MM).

    Args:
        mu_x_i:    Mean of x in the component       [x units].
        sigma_x_i: Std dev of x in the component    [x units].
        x_mean:    Overall mean <x>                 [x units].
        a_exp:     Non-negative integer moment order (static).
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        fac = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        total = total + fac * (0.5 * sigma_x_i ** 2) ** p * (mu_x_i - x_mean) ** (a_exp - 2 * p)
    return total


def univar_L_int_PDF_comp_all_MM(mu_x_i_n, sigma_x_i_n, x_mean, b_exp):
    """Central moment of a lognormally-distributed variable within a PDF component.

    Evaluates  INT(0:inf) (x - <x>)^b P_L_i(x) dx
      = SUM(q=0:b) [ b! / ((b-q)! q!) ] (-<x>)^(b-q) exp( μ_n q + ½ σ_n² q² )
    where (μ_n, σ_n) are the mean/std of ln x in the component
    (mixed_moment_PDF_integrals.F90:univar_L_int_PDF_comp_all_MM).

    Args:
        mu_x_i_n:    Mean of ln x in the component   [ln(x units)].
        sigma_x_i_n: Std dev of ln x in the component [-].
        x_mean:      Overall mean <x>                [x units].
        b_exp:       Non-negative integer moment order (static).
    """
    total = 0.0
    for q in range(b_exp + 1):
        fac = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        total = total + fac * (-x_mean) ** (b_exp - q) * jnp.exp(
            mu_x_i_n * q + 0.5 * sigma_x_i_n ** 2 * q ** 2)
    return total


def bivar_NL_int_PDF_comp_all_MM(mu_x1_i, mu_x2_i_n, sigma_x1_i, sigma_x2_i_n, corr_x1_x2_i_n,
                                 x1_mean, x2_mean, a_exp, b_exp):
    """Bivariate normal-lognormal central mixed moment within a PDF component.

    Evaluates  INT INT (x1-<x1>)^a (x2-<x2>)^b P_NL_i(x1,x2) dx2 dx1, where x1 has a normal marginal and x2 a
    lognormal marginal (in the component), correlated through corr(x1, ln x2):
      = SUM(p=0:floor(a/2)) SUM(q=0:b) [a!/((a-2p)!p!)] [b!/((b-q)!q!)]
          (½ σ_x1²)^p (μ_x1 - <x1> + ρ σ_x1 σ_x2_n q)^(a-2p)
          (-<x2>)^(b-q) exp( μ_x2_n q + ½ σ_x2_n² q² )
    (mixed_moment_PDF_integrals.F90:bivar_NL_int_PDF_comp_all_MM). Equivalent to exponential-tilting the joint
    normal by exp(q·ln x2): the lognormal raw moment factors out and x1's mean shifts by ρ σ_x1 σ_x2_n q.

    Args:
        mu_x1_i:        Mean of x1 in the component              [x1 units].
        mu_x2_i_n:      Mean of ln x2 in the component           [ln(x2 units)].
        sigma_x1_i:     Std dev of x1 in the component           [x1 units].
        sigma_x2_i_n:   Std dev of ln x2 in the component        [-].
        corr_x1_x2_i_n: Correlation of x1 and ln x2              [-].
        x1_mean:        Overall mean <x1>                        [x1 units].
        x2_mean:        Overall mean <x2>                        [x2 units].
        a_exp, b_exp:   Non-negative integer moment orders (static).
    """
    total = 0.0
    for p in range(a_exp // 2 + 1):
        fac_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        for q in range(b_exp + 1):
            fac_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
            shifted = mu_x1_i - x1_mean + corr_x1_x2_i_n * sigma_x1_i * sigma_x2_i_n * q
            total = total + (fac_p * fac_q) * (0.5 * sigma_x1_i ** 2) ** p \
                * shifted ** (a_exp - 2 * p) \
                * (-x2_mean) ** (b_exp - q) \
                * jnp.exp(mu_x2_i_n * q + 0.5 * sigma_x2_i_n ** 2 * q ** 2)
    return total


def bivar_NL_x_hm_all_MM_comp_eq(mu_x_i, mu_hm_i, mu_hm_i_n, sigma_x_i, sigma_hm_i, sigma_hm_i_n,
                                 corr_x_hm_i_n, precip_frac_i, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp):
    """Per-PDF-component contribution to <x'^a hm'^b> (mixed_moment_PDF_integrals.F90:bivar_NL_x_hm_all_MM_comp_eq).

    x is normal-marginal (w/rt/thl/sclr...), hm is a precipitating hydrometeor (lognormal-marginal in-precip).
    The component splits into within-precip (fraction precip_frac_i, where hm follows its lognormal) and
    out-of-precip (1-precip_frac_i, where hm=0 → (hm-<hm>)^b = (-<hm>)^b, x still normal). The Fortran selects
    one of four closed forms depending on whether x and/or hm vary (σ vs tol); here all four are evaluated and
    selected with `jnp.where` (each is finite for finite inputs → safe value+gradient), so the port is
    vectorizable and differentiable while reproducing the exact branch the Fortran takes.

    Args mirror the Fortran; a_exp/b_exp are static non-negative integers.
    """
    one = 1.0
    out_precip = (one - precip_frac_i) * (-hm_mean) ** b_exp     # hm = 0 contribution (common to most branches)

    # (1) both x and hm constant in the component
    both_const = (mu_x_i - x_mean) ** a_exp * (
        precip_frac_i * (mu_hm_i - hm_mean) ** b_exp + out_precip)
    # (2) only x constant
    x_const = (mu_x_i - x_mean) ** a_exp * (
        precip_frac_i * univar_L_int_PDF_comp_all_MM(mu_hm_i_n, sigma_hm_i_n, hm_mean, b_exp) + out_precip)
    # (3) only hm constant
    univar_N_x = univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp)
    hm_const = (precip_frac_i * (mu_hm_i - hm_mean) ** b_exp + out_precip) * univar_N_x
    # (4) both vary
    both_vary = precip_frac_i * bivar_NL_int_PDF_comp_all_MM(
        mu_x_i, mu_hm_i_n, sigma_x_i, sigma_hm_i_n, corr_x_hm_i_n, x_mean, hm_mean, a_exp, b_exp) \
        + out_precip * univar_N_x

    x_is_const = sigma_x_i <= x_tol
    hm_is_const = sigma_hm_i <= hm_tol
    return jnp.where(x_is_const & hm_is_const, both_const,
                     jnp.where(x_is_const, x_const,
                               jnp.where(hm_is_const, hm_const, both_vary)))


def xphmp_integral_covar(mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, sigma_x_1, sigma_x_2, sigma_hm_1, sigma_hm_2,
                         corr_x_hm_1, corr_x_hm_2, mixt_frac, precip_frac_1, precip_frac_2,
                         x_mean, x_tol, hm_tol):
    """Covariance <x'hm'> (mixed_moment_PDF_integrals.F90:xphmp_integral_covar) — the streamlined a=b=1 case.

    x is binormal (w/rt/thl/sclr), hm a precipitating hydrometeor. Within component i (in-precip only, since
    out-of-precip hm=0): E[(x-<x>)hm] = (μ_x_i - <x>) μ_hm_i + corr σ_x_i σ_hm_i. The Fortran's 4-way branch
    (drop the correlation term for whichever component has x or hm constant) decomposes into an independent
    per-component `jnp.where` — differentiable and equivalent (verified in the test). Uses hm's normal-space
    (in-precip) moments directly.
    """
    def term(mu_x, mu_hm, sig_x, sig_hm, corr):
        varies = (sig_x > x_tol) & (sig_hm > hm_tol)
        return (mu_x - x_mean) * mu_hm + jnp.where(varies, corr * sig_x * sig_hm, 0.0)
    return mixt_frac * precip_frac_1 * term(mu_x_1, mu_hm_1, sigma_x_1, sigma_hm_1, corr_x_hm_1) \
        + (1.0 - mixt_frac) * precip_frac_2 * term(mu_x_2, mu_hm_2, sigma_x_2, sigma_hm_2, corr_x_hm_2)


def hmxphmyp_integral_covar(mu_hmx_1, mu_hmx_2, mu_hmy_1, mu_hmy_2, sigma_hmx_1, sigma_hmx_2,
                            sigma_hmy_1, sigma_hmy_2, corr_hmx_hmy_1, corr_hmx_hmy_2, mixt_frac,
                            precip_frac_1, precip_frac_2, hmx_mean, hmy_mean, hmx_tol, hmy_tol):
    """Covariance <hmx'hmy'> of two precipitating hydrometeors
    (mixed_moment_PDF_integrals.F90:hmxphmyp_integral_covar).

    E[hmx·hmy] over the mixture (in-precip only) minus <hmx><hmy>: within component i in-precip,
    E[hmx hmy] = μ_hmx_i μ_hmy_i + corr σ_hmx_i σ_hmy_i. Same per-component `jnp.where` decomposition of the
    Fortran 4-way branch.
    """
    def term(mu_x, mu_y, sig_x, sig_y, corr):
        varies = (sig_x > hmx_tol) & (sig_y > hmy_tol)
        return mu_x * mu_y + jnp.where(varies, corr * sig_x * sig_y, 0.0)
    return mixt_frac * precip_frac_1 * term(mu_hmx_1, mu_hmy_1, sigma_hmx_1, sigma_hmy_1, corr_hmx_hmy_1) \
        + (1.0 - mixt_frac) * precip_frac_2 * term(mu_hmx_2, mu_hmy_2, sigma_hmx_2, sigma_hmy_2, corr_hmx_hmy_2) \
        - hmx_mean * hmy_mean


def xp_a_hmpb_integrals_all_MM(mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, mu_hm_1_n, mu_hm_2_n,
                               sigma_x_1, sigma_x_2, sigma_hm_1, sigma_hm_2, sigma_hm_1_n, sigma_hm_2_n,
                               corr_x_hm_1_n, corr_x_hm_2_n, mixt_frac, precip_frac_1, precip_frac_2,
                               x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp):
    """Any mixed moment <x'^a hm'^b> over the full 2-component PDF
    (mixed_moment_PDF_integrals.F90:xp_a_hmpb_integrals_all_MM). x binormal, hm a precipitating hydrometeor."""
    return mixt_frac * bivar_NL_x_hm_all_MM_comp_eq(
        mu_x_1, mu_hm_1, mu_hm_1_n, sigma_x_1, sigma_hm_1, sigma_hm_1_n,
        corr_x_hm_1_n, precip_frac_1, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp) \
        + (1.0 - mixt_frac) * bivar_NL_x_hm_all_MM_comp_eq(
            mu_x_2, mu_hm_2, mu_hm_2_n, sigma_x_2, sigma_hm_2, sigma_hm_2_n,
            corr_x_hm_2_n, precip_frac_2, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp)


def hydrometeor_mixed_moments(p):
    """Top driver: <rt'hm'>, <thl'hm'>, <w'^2 hm'> and <hmx'hmy'> for all hydrometeors
    (mixed_moment_PDF_integrals.F90:hydrometeor_mixed_moments).

    These appear in the liquid/ice water-loading buoyancy terms of CLUBB's predictive equations. Pure
    orchestration over the already-validated integral functions; vectorized over the `nzt` vertical levels (each
    integral broadcasts over (nzt,) arrays) with static Python loops over the small hydrometeor dimension.

    `p` is a dict of jnp arrays (all in-precip / normal-space moments already unpacked):
      hydromet           (nzt, hm_dim)   overall hydrometeor means
      mu_w_{1,2}, mu_rt_{1,2}, mu_thl_{1,2}, sigma_w_{1,2}, sigma_rt_{1,2}, sigma_thl_{1,2},
      sigma_chi_{1,2}, sigma_eta_{1,2}, mixt_frac, precip_frac_{1,2}, crt_{1,2}, cthl_{1,2}   (nzt,)
      mu_hm_{1,2}, sigma_hm_{1,2}, mu_hm_{1,2}_n, sigma_hm_{1,2}_n,
      corr_chi_hm_{1,2}, corr_eta_hm_{1,2}, corr_w_hm_{1,2}_n   (nzt, hm_dim)
      corr_hmx_hmy_{1,2}  (nzt, hm_dim, hm_dim)
      hydromet_tol        (hm_dim,)
      rt_tol, thl_tol, w_tol   scalars

    Returns dict with rtphmp_zt, thlphmp_zt, wp2hmp  (nzt, hm_dim) and hmxphmyp_zt (nzt, hm_dim, hm_dim)
    (the latter filled on the hmy_idx > hm_idx upper pairs, as in the Fortran). The Fortran's zt->zm stats
    interpolation is left to the caller.
    """
    hm_dim = p['hydromet'].shape[1]
    mixt = p['mixt_frac']
    # Overall means recomputed from PDF params (consistent with the PDF, as in the Fortran).
    rtm = compute_mean_binormal(p['mu_rt_1'], p['mu_rt_2'], mixt)
    thlm = compute_mean_binormal(p['mu_thl_1'], p['mu_thl_2'], mixt)
    wm = compute_mean_binormal(p['mu_w_1'], p['mu_w_2'], mixt)

    rtphmp_cols, thlphmp_cols, wp2hmp_cols = [], [], []
    # hmxphmyp accumulated into a list of (hmy_idx, hm_idx, column) then scattered.
    hmxphmyp = [[None] * hm_dim for _ in range(hm_dim)]

    for hm in range(hm_dim):
        mu_hm_1, mu_hm_2 = p['mu_hm_1'][:, hm], p['mu_hm_2'][:, hm]
        sig_hm_1, sig_hm_2 = p['sigma_hm_1'][:, hm], p['sigma_hm_2'][:, hm]
        hm_tol = p['hydromet_tol'][hm]

        # Correlations of rt/thl with hm from the chi/eta correlations (PDF transform).
        corr_rt_hm_1 = calc_corr_rt_x(p['crt_1'], p['sigma_rt_1'], p['sigma_chi_1'], p['sigma_eta_1'],
                                      p['corr_chi_hm_1'][:, hm], p['corr_eta_hm_1'][:, hm])
        corr_rt_hm_2 = calc_corr_rt_x(p['crt_2'], p['sigma_rt_2'], p['sigma_chi_2'], p['sigma_eta_2'],
                                      p['corr_chi_hm_2'][:, hm], p['corr_eta_hm_2'][:, hm])
        corr_thl_hm_1 = calc_corr_thl_x(p['cthl_1'], p['sigma_thl_1'], p['sigma_chi_1'], p['sigma_eta_1'],
                                        p['corr_chi_hm_1'][:, hm], p['corr_eta_hm_1'][:, hm])
        corr_thl_hm_2 = calc_corr_thl_x(p['cthl_2'], p['sigma_thl_2'], p['sigma_chi_2'], p['sigma_eta_2'],
                                        p['corr_chi_hm_2'][:, hm], p['corr_eta_hm_2'][:, hm])

        rtphmp_cols.append(xphmp_integral_covar(
            p['mu_rt_1'], p['mu_rt_2'], mu_hm_1, mu_hm_2, p['sigma_rt_1'], p['sigma_rt_2'],
            sig_hm_1, sig_hm_2, corr_rt_hm_1, corr_rt_hm_2, mixt, p['precip_frac_1'], p['precip_frac_2'],
            rtm, p['rt_tol'], hm_tol))
        thlphmp_cols.append(xphmp_integral_covar(
            p['mu_thl_1'], p['mu_thl_2'], mu_hm_1, mu_hm_2, p['sigma_thl_1'], p['sigma_thl_2'],
            sig_hm_1, sig_hm_2, corr_thl_hm_1, corr_thl_hm_2, mixt, p['precip_frac_1'], p['precip_frac_2'],
            thlm, p['thl_tol'], hm_tol))
        wp2hmp_cols.append(xp_a_hmpb_integrals_all_MM(
            p['mu_w_1'], p['mu_w_2'], mu_hm_1, mu_hm_2, p['mu_hm_1_n'][:, hm], p['mu_hm_2_n'][:, hm],
            p['sigma_w_1'], p['sigma_w_2'], sig_hm_1, sig_hm_2, p['sigma_hm_1_n'][:, hm], p['sigma_hm_2_n'][:, hm],
            p['corr_w_hm_1_n'][:, hm], p['corr_w_hm_2_n'][:, hm], mixt, p['precip_frac_1'], p['precip_frac_2'],
            wm, p['hydromet'][:, hm], p['w_tol'], hm_tol, 2, 1))

        for hmy in range(hm + 1, hm_dim):
            hmxphmyp[hmy][hm] = hmxphmyp_integral_covar(
                mu_hm_1, mu_hm_2, p['mu_hm_1'][:, hmy], p['mu_hm_2'][:, hmy],
                sig_hm_1, sig_hm_2, p['sigma_hm_1'][:, hmy], p['sigma_hm_2'][:, hmy],
                p['corr_hmx_hmy_1'][:, hm, hmy], p['corr_hmx_hmy_2'][:, hm, hmy],
                mixt, p['precip_frac_1'], p['precip_frac_2'], p['hydromet'][:, hm], p['hydromet'][:, hmy],
                hm_tol, p['hydromet_tol'][hmy])

    nzt = p['hydromet'].shape[0]
    zeros = jnp.zeros(nzt)
    hmxphmyp_zt = jnp.stack([jnp.stack([hmxphmyp[hmy][hm] if hmxphmyp[hmy][hm] is not None else zeros
                                        for hm in range(hm_dim)], axis=1) for hmy in range(hm_dim)], axis=1)
    return {
        'rtphmp_zt': jnp.stack(rtphmp_cols, axis=1),
        'thlphmp_zt': jnp.stack(thlphmp_cols, axis=1),
        'wp2hmp': jnp.stack(wp2hmp_cols, axis=1),
        'hmxphmyp_zt': hmxphmyp_zt,
    }
