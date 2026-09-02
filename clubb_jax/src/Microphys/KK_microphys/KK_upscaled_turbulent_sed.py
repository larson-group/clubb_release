"""JAX port of KK_upscaled_turbulent_sed.F90 — KK sedimentation-velocity covariances.

`KK_sed_vel_covars` computes the covariances of the rain sedimentation velocities with the rain
fields, used in the semi-implicit turbulent-sedimentation term of the rrm/Nrm predictive equations
(advance_microphys). Per KK00:
    <V_rr' r_r'> = 0.012 (1e6) <r_r' R_vr'>,   <V_Nr' N_r'> = 0.007 (1e6) <N_r' R_vr'>,
where R_vr is the rain-drop mean volume radius (= KK_mvr_coef r_r^(1/3) N_r^(-1/3)) and the
covariances <r_r' R_vr'>, <N_r' R_vr'> are bivariate-lognormal integrals.

Each covariance is written `<x' R_vr'> = coefA <x> + termB` so the turbulent-sed term can be handled
SEMI-IMPLICITLY (l_semi_imp_turbulent_sed=.true., the rico default): coefA -> the implicit component
(impc), termB -> the explicit (expc). The partials are the bivariate-lognormal moment with an extra
factor of the differenced variable: vs `bivar_LL_mean` (PDF_integrals_means), the x1 variance term
carries (alpha^2+2 alpha) instead of alpha^2 and the cross term (alpha+1) instead of alpha. The Nr
side reuses the same machinery with r_r<->N_r and the two exponents swapped.

Verified bit-exact vs the rico oracle `rr_KK_mvr_covar_zt` / `Nr_KK_mvr_covar_zt` stats
(tests/test_kk_rico_oracle.py). Exponents alpha=KK_mvr_rr_exp=1/3, beta=KK_mvr_Nr_exp=-1/3. All-jnp,
differentiable.
"""
import jax.numpy as jnp

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_LL_mean_const_all, bivar_LL_mean_const_x1)
# rr_tol from constants_clubb (KK_upscaled_turbulent_sed.F90:1224 `use constants_clubb`); NR_TOL is the
# pi-derived rain-number tolerance, computed full-precision in KK_upscaled_means.
from clubb_jax.src.CLUBB_core.constants_clubb import rr_tol as RR_TOL
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import NR_TOL
# KK_mvr_rr_exp/KK_mvr_Nr_exp from parameters_KK (KK_upscaled_turbulent_sed.F90:349-351 `use parameters_KK`);
# KK_mvr_coef from its KK_tendency_coefs home in KK_microphys_module.
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import (
    KK_mvr_rr_exp as KK_MVR_RR_EXP, KK_mvr_Nr_exp as KK_MVR_NR_EXP)
from clubb_jax.src.Microphys.KK_microphys_module import KK_MVR_COEF

_MICRON_PER_M = 1.0e6


def bivar_LL_covar_partial(mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho, alpha, beta):
    """Both x1, x2 vary. KK_upscaled_turbulent_sed.F90:1540. The (alpha^2+2 alpha) and (alpha+1)
    terms are the `bivar_LL_mean` exponents with the extra x1 factor of the covariance."""
    return jnp.exp(mu_x1_n * alpha + mu_x2_n * beta
                   + 0.5 * sigma_x1_n ** 2 * (alpha ** 2 + 2.0 * alpha)
                   + 0.5 * sigma_x2_n ** 2 * beta ** 2
                   + rho * sigma_x1_n * (alpha + 1.0) * sigma_x2_n * beta)


def bivar_LL_covar_const_x2_partial(mu_x1_n, mu_x2, sigma_x1_n, alpha, beta):
    """x2 constant (sigma_x2 -> 0). KK_upscaled_turbulent_sed.F90:1584.
    = mu_x2^beta exp( mu_x1_n alpha + 1/2 sigma_x1_n^2 (alpha^2 + 2 alpha) )."""
    return mu_x2 ** beta * jnp.exp(mu_x1_n * alpha
                                   + 0.5 * sigma_x1_n ** 2 * (alpha ** 2 + 2.0 * alpha))


def _covar_partial(mu_x1, mu_x2, mu_x1_n, mu_x2_n, sigma_x1, sigma_x2,
                   sigma_x1_n, sigma_x2_n, rho, alpha, beta, x1_tol, x2_tol):
    """Generic partial of <x1' x1^alpha x2^beta> (the coefA building block), 4-way variance
    dispatch (KK_upscaled_turbulent_sed.F90:1209/1355). The x1-constant branch reuses the MEAN's
    const_x1 (the extra x1 factor vanishes when sigma_x1->0); the x2-constant and general branches
    use the covariance partials (which keep the alpha^2+2 alpha / alpha+1 terms)."""
    x1_const = sigma_x1 <= x1_tol
    x2_const = sigma_x2 <= x2_tol
    f_all = bivar_LL_mean_const_all(mu_x1, mu_x2, alpha, beta)
    f_x1 = bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha, beta)
    f_x2 = bivar_LL_covar_const_x2_partial(mu_x1_n, mu_x2, sigma_x1_n, alpha, beta)
    f_gen = bivar_LL_covar_partial(mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho, alpha, beta)
    return jnp.where(x1_const & x2_const, f_all,
                     jnp.where(x1_const, f_x1,
                               jnp.where(x2_const, f_x2, f_gen)))


def bivar_LL_covar_partial_rr(mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, s_rr, s_Nr, s_rr_n, s_Nr_n, rho):
    """<r_r' R_vr'> partial: x1=r_r, x2=N_r, alpha=rr_exp, beta=Nr_exp."""
    a, b = KK_MVR_RR_EXP, KK_MVR_NR_EXP
    mu_x1 = mu_rr if a >= 0.0 else jnp.maximum(mu_rr, RR_TOL)
    mu_x2 = mu_Nr if b >= 0.0 else jnp.maximum(mu_Nr, NR_TOL)
    return _covar_partial(mu_x1, mu_x2, mu_rr_n, mu_Nr_n, s_rr, s_Nr,
                          s_rr_n, s_Nr_n, rho, a, b, RR_TOL, NR_TOL)


def bivar_LL_covar_partial_Nr(mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, s_rr, s_Nr, s_rr_n, s_Nr_n, rho):
    """<N_r' R_vr'> partial: x1=N_r, x2=r_r, alpha=Nr_exp, beta=rr_exp (rr<->Nr swapped)."""
    a, b = KK_MVR_NR_EXP, KK_MVR_RR_EXP          # alpha=Nr exp, beta=rr exp
    mu_x1 = mu_Nr if a >= 0.0 else jnp.maximum(mu_Nr, NR_TOL)
    mu_x2 = mu_rr if b >= 0.0 else jnp.maximum(mu_rr, RR_TOL)
    return _covar_partial(mu_x1, mu_x2, mu_Nr_n, mu_rr_n, s_Nr, s_rr,
                          s_Nr_n, s_rr_n, rho, a, b, NR_TOL, RR_TOL)


def KK_sed_vel_covars(rr_1, rr_2, Nr_1, Nr_2, mvr,
                      mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2,
                      mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
                      sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
                      sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
                      corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac):
    """KK sed-velocity covariances (semi-implicit form). Oracle KK_sed_vel_covars
    (KK_upscaled_turbulent_sed.F90:25), l_semi_imp_turbulent_sed=.true.

    rr_1/rr_2, Nr_1/Nr_2 are the per-component overall means (= precip_frac_i * mu_hm_i, so that
    <r_r> = a rr_1 + (1-a) rr_2). mvr is the rain-drop mean volume radius (KK_mean_vol_rad).

    Returns a dict: rr_KK_mvr_covar, Nr_KK_mvr_covar (the <x' R_vr'> covariances), and the
    Vrrprrp/VNrpNrp implicit (impc) + explicit (expc) components of <V_x' x'>."""
    a = mixt_frac
    rrm = a * rr_1 + (1.0 - a) * rr_2          # within-step-consistent overall means
    Nrm = a * Nr_1 + (1.0 - a) * Nr_2

    pr1 = bivar_LL_covar_partial_rr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n,
                      sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n)
    pr2 = bivar_LL_covar_partial_rr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                      sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n)
    rr_coefA = KK_MVR_COEF * (pr1 + pr2) - mvr
    rr_termB = -KK_MVR_COEF * (a * rr_1 * pr2 + (1.0 - a) * rr_2 * pr1)
    rr_KK_mvr_covar = rr_coefA * rrm + rr_termB

    pn1 = bivar_LL_covar_partial_Nr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n,
                      sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n)
    pn2 = bivar_LL_covar_partial_Nr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                      sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n)
    Nr_coefA = KK_MVR_COEF * (pn1 + pn2) - mvr
    Nr_termB = -KK_MVR_COEF * (a * Nr_1 * pn2 + (1.0 - a) * Nr_2 * pn1)
    Nr_KK_mvr_covar = Nr_coefA * Nrm + Nr_termB

    return {
        'rr_KK_mvr_covar': rr_KK_mvr_covar,
        'Nr_KK_mvr_covar': Nr_KK_mvr_covar,
        'Vrrprrp_impc': -0.012 * _MICRON_PER_M * rr_coefA,
        'Vrrprrp_expc': -0.012 * _MICRON_PER_M * rr_termB,
        'VNrpNrp_impc': -0.007 * _MICRON_PER_M * Nr_coefA,
        'VNrpNrp_expc': -0.007 * _MICRON_PER_M * Nr_termB,
    }
