"""Synthetic PDF preparation for gradient tests; not runtime scheme orchestration."""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import Nc_in_cloud_to_Ncnm
from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_evap_upscaled_mean,
)
from clubb_jax.tests.microphysics_coefficients import kk_evap_coef, kk_auto_coef
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import C_evap as _C_EVAP_DEFAULT  # parameters_KK.F90:48


def _hm_log_moments(mu_hm, sigma_hm):
    """In-precip lognormal moments (mu_n, sigma_n) of a hydrometeor from its linear
    in-precip mean/stdev, with sigma2_on_mu2 = (sigma_hm/mu_hm)^2."""
    mu_safe = jnp.maximum(jnp.abs(mu_hm), 1e-30)
    s2m2 = (sigma_hm / mu_safe) ** 2
    return mean_L2N(mu_safe, s2m2), stdev_L2N(s2m2)


def kk_autoconversion_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac,
                           Nc_in_cloud, cloud_frac_1, cloud_frac_2, rho,
                           const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn, corr_chi_Ncn_n):
    """Mean upscaled-KK rain-water autoconversion tendency <KK_auto> from the PDF state.

    mu_chi_i, sigma_chi_i : chi PDF component means/stdevs (from the CLUBB PDF closure).
    Nc_in_cloud           : in-cloud mean cloud-droplet concentration [num/kg].
    cloud_frac_1/2, mixt_frac : PDF cloud fractions and mixture fraction.
    rho                   : air density [kg/m^3] (for kk_auto_coef).
    const_Ncnp2_on_Ncnm2  : prescribed <Ncn'^2>/<Ncn>^2 (0 => constant N_c).
    const_corr_chi_Ncn    : prescribed LINEAR corr(chi, Ncn) (for the Ncnm inversion).
    corr_chi_Ncn_n        : prescribed NORMAL-space corr(chi, ln Ncn) (rate-function input).
    Returns rrm_auto [(kg/kg)/s]."""
    Ncnm = Nc_in_cloud_to_Ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac,
                               Nc_in_cloud, cloud_frac_1, cloud_frac_2,
                               const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn)
    # N_cn is a single lognormal over the domain: component params are equal.
    sigma_Ncn = jnp.sqrt(const_Ncnp2_on_Ncnm2) * Ncnm
    mu_Ncn_n = mean_L2N(Ncnm, const_Ncnp2_on_Ncnm2)
    sigma_Ncn_n = stdev_L2N(const_Ncnp2_on_Ncnm2)
    coef = kk_auto_coef(rho)
    return KK_auto_upscaled_mean(
        mu_chi_1, mu_chi_2, Ncnm, Ncnm, mu_Ncn_n, mu_Ncn_n,
        sigma_chi_1, sigma_chi_2, sigma_Ncn, sigma_Ncn, sigma_Ncn_n, sigma_Ncn_n,
        corr_chi_Ncn_n, corr_chi_Ncn_n, coef, mixt_frac)


def kk_accretion_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2,
                      mu_rr_1, mu_rr_2, sigma_rr_1, sigma_rr_2,
                      corr_chi_rr_1_n, corr_chi_rr_2_n,
                      mixt_frac, precip_frac_1, precip_frac_2):
    """Mean upscaled-KK rain-water accretion tendency <KK_accr> from the PDF state.

    mu_rr_i, sigma_rr_i : IN-PRECIP r_r component means/stdevs (from calc_comp_mu_sigma_hm).
    corr_chi_rr_i_n     : prescribed NORMAL-space corr(chi, ln r_r). Returns rrm_accr."""
    mu_rr_1_n, sigma_rr_1_n = _hm_log_moments(mu_rr_1, sigma_rr_1)
    mu_rr_2_n, sigma_rr_2_n = _hm_log_moments(mu_rr_2, sigma_rr_2)
    return KK_accr_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
        sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n,
        corr_chi_rr_1_n, corr_chi_rr_2_n, mixt_frac, precip_frac_1, precip_frac_2)


def kk_evaporation_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2,
                        mu_rr_1, mu_rr_2, sigma_rr_1, sigma_rr_2,
                        mu_Nr_1, mu_Nr_2, sigma_Nr_1, sigma_Nr_2,
                        corr_chi_rr_1_n, corr_chi_rr_2_n,
                        corr_chi_Nr_1_n, corr_chi_Nr_2_n,
                        corr_rr_Nr_1_n, corr_rr_Nr_2_n,
                        T_liq, p_in_Pa, C_evap, mixt_frac, precip_frac_1, precip_frac_2,
                        saturation_formula=3):
    """Mean upscaled-KK rain-water evaporation tendency <KK_evap> from the PDF state.

    In-precip r_r and N_r component moments + the 6 prescribed normal-space correlations;
    the thermodynamic coefficient kk_evap_coef(T_liq, p, C_evap). Returns rrm_evap (<0)."""
    mu_rr_1_n, sigma_rr_1_n = _hm_log_moments(mu_rr_1, sigma_rr_1)
    mu_rr_2_n, sigma_rr_2_n = _hm_log_moments(mu_rr_2, sigma_rr_2)
    mu_Nr_1_n, sigma_Nr_1_n = _hm_log_moments(mu_Nr_1, sigma_Nr_1)
    mu_Nr_2_n, sigma_Nr_2_n = _hm_log_moments(mu_Nr_2, sigma_Nr_2)
    coef = kk_evap_coef(T_liq, p_in_Pa, C_evap, saturation_formula)
    return KK_evap_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2,
        mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
        sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n,
        corr_rr_Nr_1_n, corr_rr_Nr_2_n, coef, mixt_frac, precip_frac_1, precip_frac_2)
