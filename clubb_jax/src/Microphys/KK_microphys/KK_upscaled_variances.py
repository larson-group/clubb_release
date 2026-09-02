"""JAX port of KK_microphys/KK_upscaled_variances.F90 — variance of the KK rain mean-volume radius.

Mirrors `variance_KK_mvr`: `Var(R_vr) = E[R_vr^2] - E[R_vr]^2` for `R_vr = coef * r_r^(1/3) * N_r^(-1/3)` over the
2-component, in-precip bivariate-lognormal rain PDF. `E[R_vr^2]` is assembled from the already-ported
`bivar_LL_mean_eq` (the `<r_r^a N_r^b>` integral) with the exponents DOUBLED (2*alpha, 2*beta); `E[R_vr]` is the
precomputed `KK_mean_vol_rad` input. Pure jnp → differentiable. Validated against an independent closed-form
lognormal-moment computation + a Monte-Carlo sample (`tests/test_kk_upscaled_variances.py`).
"""
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import bivar_LL_mean_eq

_MVR_RR_EXP, _MVR_NR_EXP = 1.0 / 3.0, -1.0 / 3.0   # parameters_KK.F90: KK_mvr_rr_exp / KK_mvr_Nr_exp


def variance_KK_mvr(mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
                    sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n,
                    sigma_Nr_1_n, sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n,
                    KK_mean_vol_rad, KK_mvr_coef, mixt_frac, precip_frac_1, precip_frac_2):
    """Variance of the KK rain-drop mean volume radius < R_vr'^2 >  [m^2]."""
    a2, b2 = 2.0 * _MVR_RR_EXP, 2.0 * _MVR_NR_EXP
    e_r2 = KK_mvr_coef ** 2 * (
        mixt_frac * precip_frac_1
        * bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1, sigma_Nr_1,
                           sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n, a2, b2)
        + (1.0 - mixt_frac) * precip_frac_2
        * bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2, sigma_Nr_2,
                           sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n, a2, b2))
    return e_r2 - KK_mean_vol_rad ** 2
