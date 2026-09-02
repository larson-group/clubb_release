"""Assembly of the upscaled-KK autoconversion tendency from the PDF state.

This composes the individually-verified pieces into the single function a running model
calls for the rain-water autoconversion source <KK_auto>:

  Nc_in_cloud_to_Ncnm  (Iter117)  -> <Ncn> (= mu_Ncn) from the in-cloud <Nc> and chi PDF
  mean_L2N / stdev_L2N (Iter109)  -> the lognormal N_cn moments mu_Ncn_n, sigma_Ncn_n
  kk_auto_coef(rho)    (Iter108)  -> the autoconversion coefficient
  KK_auto_upscaled_mean(Iter108)  -> the analytic PDF-integrated mean tendency

N_cn is a single lognormal over the domain (mu_Ncn_1 = mu_Ncn_2 = Ncnm, etc.), with its
variance-over-mean^2 prescribed by const_Ncnp2_on_Ncnm2 and its chi correlation prescribed
(normal space) by corr_chi_Ncn_n. When const_Ncnp2_on_Ncnm2 = 0 (constant N_c, e.g. rico),
sigma_Ncn = 0 and the dispatch takes the const_x2 path (correlations irrelevant).

Validated end-to-end against rico's rrm_auto (tests/test_kk_rico_oracle.py extends to this
composed entry point). All-jnp and differentiable.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import Nc_in_cloud_to_Ncnm
from clubb_jax.src.CLUBB_core.pdf_utilities import mean_L2N, stdev_L2N
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_evap_upscaled_mean,
)
from clubb_jax.src.Microphys.KK_microphys_module import kk_evap_coef, kk_auto_coef
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


# KK_microphys_adjust + KK_sedimentation now live in their Fortran-home module
# Microphys/KK_microphys_module.py (mirror-refactor iter 112), imported below.
from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust


# --- prescribed normal-space correlations, derived from the corr_varnce_module default arrays
# (in-cloud == below-cloud for all rate-relevant entries; see corr_varnce_module). ---
from clubb_jax.src.CLUBB_core.corr_varnce_module import kk_prescribed_correlations
_KK_CORRS = kk_prescribed_correlations()
_CORR_CHI_RR_N = _KK_CORRS['corr_chi_rr']    # 0.788
_CORR_CHI_NR_N = _KK_CORRS['corr_chi_Nr']    # 0.675
_CORR_RR_NR_N = _KK_CORRS['corr_rr_Nr']      # 0.821


# KK_sedimentation moved to Microphys/KK_microphys_module.py (mirror-refactor iter 112), imported above.


# hydrometp2_zt (the overall hydrometeor variance, setup_clubb_pdf_params.F90:449) now lives in its
# Fortran-home module CLUBB_core/setup_clubb_pdf_params.py (mirror-refactor iter 192), imported where used.


def compute_kk_microphysics(rrm, Nrm, mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2,
                            mixt_frac, thl_1, thl_2, Nc_in_cloud, cloud_frac_1, cloud_frac_2,
                            precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol,
                            rho, T_liq, p_in_Pa, exner, rcm, dt,
                            rr_ratio=1.25, Nr_ratio=1.25, omicron=0.5, zeta=0.0,
                            C_evap=_C_EVAP_DEFAULT, hm_tol=1.0e-10, l_return_vel_prereqs=False):
    """Full upscaled-KK microphysics step (PDF state + hydromet fields -> state tendencies).

    Composes: hydrometeor in-precip component moments (calc_comp_mu_sigma_hm for r_r and N_r),
    the cloud-nuclei mean (Nc_in_cloud_to_Ncnm), all KK rates (auto/accr/evap mass; auto/evap
    number), and the tendency assembly (KK_microphys_adjust). Prescribed normal-space correlations
    are the corr_varnce_module defaults. N_cn is constant (rico: l_const_Nc_in_cloud).
    Returns (rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc).

    NOTE: validated against rico for the autoconversion-dominated part; the accr/evap contributions
    from the rrm/Nrm FIELDS carry the documented timing confound (calc_comp_mu_sigma_hm) and are only
    fully validatable in a running rico."""
    from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import (
        compute_mean_stdev, norm_transform_mean_stdev, IIPDF_NCN, hydrometp2_zt)
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import (
        KK_Nrm_auto_mean, KK_Nrm_evap_upscaled_mean)
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import KK_evap_upscaled_mean
    from clubb_jax.src.Microphys.KK_microphys_module import kk_evap_coef

    # In-precip r_r and N_r component moments (linear + normal/log space), assembled via the
    # faithful setup_pdf_parameters orchestration. The hydrometeors follow Ncn in pdf order:
    # [chi, eta, w, Ncn, rr, Nr] -> rr at IIPDF_NCN+1, Nr at IIPDF_NCN+2. Ncn/chi are handled
    # separately by the rate entry points here, so dummy Ncn=0 (l_const) is passed.
    hmp2_rr = hydrometp2_zt(rrm, precip_frac, rr_ratio)
    hmp2_Nr = hydrometp2_zt(Nrm, precip_frac, Nr_ratio)
    Nr_tol = hm_tol / ((4.0 / 3.0) * jnp.pi * 1000.0 * (5.0e-3) ** 3)
    iirr, iiNr = IIPDF_NCN + 1, IIPDF_NCN + 2
    mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, hm_1, hm_2, s2m2_1, s2m2_2 = compute_mean_stdev(
        mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, 0.0, 0.0,
        jnp.zeros_like(rrm), 0.0, True,
        [(rrm, hmp2_rr, rr_ratio, hm_tol), (Nrm, hmp2_Nr, Nr_ratio, Nr_tol)],
        thl_1, thl_2, mixt_frac, precip_frac, precip_frac_1, precip_frac_2,
        precip_frac_tol, omicron, zeta)
    mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n = norm_transform_mean_stdev(
        mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, s2m2_1, s2m2_2,
        jnp.zeros_like(rrm), hm_1, hm_2, [hm_tol, Nr_tol], True)
    mu_rr_1, mu_rr_2 = mu_x_1[..., iirr], mu_x_2[..., iirr]
    sigma_rr_1, sigma_rr_2 = sigma_x_1[..., iirr], sigma_x_2[..., iirr]
    mu_Nr_1, mu_Nr_2 = mu_x_1[..., iiNr], mu_x_2[..., iiNr]
    sigma_Nr_1, sigma_Nr_2 = sigma_x_1[..., iiNr], sigma_x_2[..., iiNr]

    one = jnp.ones_like(rrm)
    ccr, ccN, crN = _CORR_CHI_RR_N * one, _CORR_CHI_NR_N * one, _CORR_RR_NR_N * one

    # Rates. Autoconversion via the validated driver (N_cn constant).
    auto = kk_autoconversion_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac,
                                  Nc_in_cloud, cloud_frac_1, cloud_frac_2, rho, 0.0, 0.0, 0.0)
    accr = kk_accretion_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2,
                             mu_rr_1, mu_rr_2, sigma_rr_1, sigma_rr_2, ccr, ccr,
                             mixt_frac, precip_frac_1, precip_frac_2)
    evap = kk_evaporation_mean(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2,
                               mu_rr_1, mu_rr_2, sigma_rr_1, sigma_rr_2,
                               mu_Nr_1, mu_Nr_2, sigma_Nr_1, sigma_Nr_2,
                               ccr, ccr, ccN, ccN, crN, crN, T_liq, p_in_Pa, C_evap,
                               mixt_frac, precip_frac_1, precip_frac_2)
    Nrm_auto = KK_Nrm_auto_mean(auto)
    mrr1n, srr1n = mu_x_1_n[..., iirr], sigma_x_1_n[..., iirr]
    mrr2n, srr2n = mu_x_2_n[..., iirr], sigma_x_2_n[..., iirr]
    mNr1n, sNr1n = mu_x_1_n[..., iiNr], sigma_x_1_n[..., iiNr]
    mNr2n, sNr2n = mu_x_2_n[..., iiNr], sigma_x_2_n[..., iiNr]
    coef = kk_evap_coef(T_liq, p_in_Pa, C_evap)
    Nrm_evap = KK_Nrm_evap_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mrr1n, mrr2n, mNr1n, mNr2n,
        sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
        srr1n, srr2n, sNr1n, sNr2n, ccr, ccr, ccN, ccN, crN, crN, coef,
        mixt_frac, precip_frac_1, precip_frac_2, dt)

    tendencies = KK_microphys_adjust(dt, exner, rcm, rrm, Nrm, evap, auto, accr, Nrm_evap, Nrm_auto)
    if not l_return_vel_prereqs:
        return tendencies

    # Velocity prerequisites for the transport step (advance_one_hydrometeor): the upscaled mean
    # volume radius mvr (KK_upscaled_means.F90:483) + the rr/Nr in-precip component moments (linear +
    # normal space) that KK_sedimentation / KK_sed_vel_covars consume. Computed from the SAME locals
    # the rates use (zero added physics; constant normal-space corr_rr_Nr = _CORR_RR_NR_N).
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import KK_mvr_upscaled_mean
    mvr = KK_mvr_upscaled_mean(mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mrr1n, mrr2n, mNr1n, mNr2n,
                               sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
                               srr1n, srr2n, sNr1n, sNr2n, crN, crN,
                               mixt_frac, precip_frac_1, precip_frac_2)
    prereqs = dict(
        mvr=mvr,
        mu_rr_1=mu_rr_1, mu_rr_2=mu_rr_2, sigma_rr_1=sigma_rr_1, sigma_rr_2=sigma_rr_2,
        mu_Nr_1=mu_Nr_1, mu_Nr_2=mu_Nr_2, sigma_Nr_1=sigma_Nr_1, sigma_Nr_2=sigma_Nr_2,
        mu_rr_1_n=mrr1n, mu_rr_2_n=mrr2n, sigma_rr_1_n=srr1n, sigma_rr_2_n=srr2n,
        mu_Nr_1_n=mNr1n, mu_Nr_2_n=mNr2n, sigma_Nr_1_n=sNr1n, sigma_Nr_2_n=sNr2n,
        corr_rr_Nr_n=crN,
        # Process mean tendencies (auto/accr/evap mass; auto/evap number) — the KK_*_tndcy inputs the
        # second-moment covariance driver (KK_upscaled_covar_driver) needs, plus the process coefficients
        # and Ncnm (= mu_Ncn, constant-Nc so sigma_Ncn=0) it consumes for the auto covariance.
        auto=auto, accr=accr, evap=evap, Nrm_auto=Nrm_auto, Nrm_evap=Nrm_evap,
        coef_auto=kk_auto_coef(rho), coef_accr=67.0 * one, coef_evap=coef,
        Ncnm=Nc_in_cloud_to_Ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac,
                                 Nc_in_cloud, cloud_frac_1, cloud_frac_2, 0.0, 0.0))
    return tendencies, prereqs
