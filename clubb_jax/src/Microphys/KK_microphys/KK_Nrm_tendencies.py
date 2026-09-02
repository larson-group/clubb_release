"""JAX port of KK_Nrm_tendencies.F90 — rain-drop-concentration (N_r) KK tendencies.

KK_Nrm_auto_mean        : the N_r autoconversion tendency from the r_r autoconversion rate
                          (each new rain drop has initial radius r_0). KK_Nrm_tendencies.F90:296.
KK_Nrm_evap_local_mean  : the local N_r evaporation tendency from the r_r evaporation rate
                          (KK 2000 Eq. 23; exponent KK_Nrm_evap_nu = 1). KK_Nrm_tendencies.F90:205.

(The upscaled N_r evaporation tendency KK_Nrm_evap_upscaled_mean is a separate larger function,
not yet ported; the local form here is the one used by the over-evaporation limiter in
KK_microphys_adjust.)
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import rho_lw
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import trivar_NLL_mean_eq
# The KK evap exponents / r_0 / KK_Nrm_evap_nu live in their Fortran home parameters_KK.F90
# (KK_Nrm_tendencies.F90:114-118 `use parameters_KK`)
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import (
    KK_evap_Supersat_exp as KK_EVAP_SUPERSAT_EXP, KK_evap_rr_exp as KK_EVAP_RR_EXP,
    KK_evap_Nr_exp as KK_EVAP_NR_EXP, r_0 as _R_0, KK_Nrm_evap_nu as _KK_NRM_EVAP_NU,
)


def KK_Nrm_auto_mean(KK_rrm_auto_tndcy):
    """N_r autoconversion tendency. KK_Nrm_tendencies.F90:296.

    = (dr_r/dt)_auto / ((4/3) pi rho_lw r_0^3) — one new drop of radius r_0 per unit mass."""
    return KK_rrm_auto_tndcy / ((4.0 / 3.0) * jnp.pi * rho_lw * _R_0 ** 3)


def KK_Nrm_evap_local_mean(KK_rrm_evap_tndcy, Nrm, rrm, dt):
    """Local N_r evaporation tendency. KK_Nrm_tendencies.F90:205 (KK 2000 Eq. 23).

    = dt^(nu-1) * (Nrm / rrm^nu) * KK_rrm_evap_tndcy^nu, nu = KK_Nrm_evap_nu = 1,
    which reduces to (Nrm/rrm) * KK_rrm_evap_tndcy."""
    nu = _KK_NRM_EVAP_NU
    rrm_safe = jnp.where(rrm > 0.0, rrm, 1.0)
    return dt ** (nu - 1.0) * (Nrm / rrm_safe ** nu) * KK_rrm_evap_tndcy ** nu


def KK_Nrm_evap_upscaled_mean(mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2,
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
                              sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                              sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n,
                              sigma_Nr_1_n, sigma_Nr_2_n,
                              corr_chi_rr_1_n, corr_chi_rr_2_n,
                              corr_chi_Nr_1_n, corr_chi_Nr_2_n,
                              corr_rr_Nr_1_n, corr_rr_Nr_2_n,
                              KK_evap_coef, mixt_frac, precip_frac_1, precip_frac_2, dt):
    """Upscaled KK N_r evaporation tendency. KK_Nrm_tendencies.F90:16.

    Reuses the trivariate NLL mean (the same integral family as r_r evaporation) with shifted
    exponents alpha = Supersat_exp*nu = 1, beta = (rr_exp-1)*nu = -2/3, gamma = Nr_exp*nu+1 = 5/3,
    and prefactor dt^(nu-1)*KK_evap_coef^nu (= KK_evap_coef since nu=1). Returns (dN_r/dt)_evap."""
    nu = _KK_NRM_EVAP_NU
    a = KK_EVAP_SUPERSAT_EXP * nu
    b = (KK_EVAP_RR_EXP - 1.0) * nu
    g = KK_EVAP_NR_EXP * nu + 1.0
    comp1 = trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n,
                               sigma_chi_1, sigma_rr_1, sigma_Nr_1, sigma_rr_1_n,
                               sigma_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n,
                               corr_rr_Nr_1_n, a, b, g)
    comp2 = trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                               sigma_chi_2, sigma_rr_2, sigma_Nr_2, sigma_rr_2_n,
                               sigma_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n,
                               corr_rr_Nr_2_n, a, b, g)
    return (dt ** (nu - 1.0) * KK_evap_coef ** nu
            * (mixt_frac * precip_frac_1 * comp1
               + (1.0 - mixt_frac) * precip_frac_2 * comp2))
