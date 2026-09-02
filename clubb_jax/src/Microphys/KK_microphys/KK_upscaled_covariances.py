import jax.numpy as jnp
"""KK upscaled second-moment microphysics covariances — JAX port of KK_upscaled_covariances.F90.

These give the covariances of the KK process tendencies (auto/accr/evap) with rt/thl/w that become the
SECOND-MOMENT microphysics source terms wprtp_mc/wpthlp_mc/rtp2_mc/thlp2_mc/rtpthlp_mc (the Iter158-localized
gap that keeps KK from being bit-faithful). They compose the validated `trivar_NNL_covar_eq` covariance
integral (PDF_integrals_covar.py) with the already-ported `bivar_NL_mean_eq` mean integral. Differentiable.
"""

# The covar-integral primitives the *_covar_eq routines compose (mirroring the Fortran
# KK_upscaled_covariances USE PDF_integrals_covar relationship; the _eq routines moved here iter 171):
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_covar import (
    quadrivar_NNLL_covar, quadrivar_NNLL_covar_const_all, quadrivar_NNLL_covar_const_x1,
    quadrivar_NNLL_covar_const_x1x2, quadrivar_NNLL_covar_const_x1x3, quadrivar_NNLL_covar_const_x2,
    quadrivar_NNLL_covar_const_x2x3, quadrivar_NNLL_covar_const_x3, quadrivar_NNLL_covar_const_x3x4,
    quadrivar_NNLL_covar_cst_x1x2x3, quadrivar_NNLL_covar_cst_x1x3x4, quadrivar_NNLL_covar_cst_x2x3x4,
    trivar_NNL_covar, trivar_NNL_covar_const_x1, trivar_NNL_covar_const_x1x2,
    trivar_NNL_covar_const_x1x3, trivar_NNL_covar_const_x2, trivar_NNL_covar_const_x2x3,
    trivar_NNL_covar_const_x3)
# KK_upscaled_covariances.F90 `use constants_clubb` (eta tol + parab-cyl input limit + Nc_tol/rr_tol, F90:90-91/1138)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    chi_tol as _CHI_TOL, parab_cyl_max_input as _PARAB_CYL_MAX, Nc_tol as NC_TOL, rr_tol as RR_TOL,
)
# NR_TOL (the pi-derived rain-number tolerance) is computed full-precision in KK_upscaled_means.
from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
    bivar_NL_mean_eq, trivar_NLL_mean_eq, NR_TOL)
# The KK exponents live in their Fortran home parameters_KK.F90 (KK_upscaled_covariances.F90:504-507/1029-1030
# `use parameters_KK`)
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import (
    KK_auto_rc_exp as KK_AUTO_RC_EXP, KK_auto_Nc_exp as KK_AUTO_NC_EXP,
    KK_accr_rc_exp as KK_ACCR_RC_EXP, KK_accr_rr_exp as KK_ACCR_RR_EXP,
    KK_evap_Supersat_exp as KK_EVAP_SUPERSAT_EXP, KK_evap_rr_exp as KK_EVAP_RR_EXP,
    KK_evap_Nr_exp as KK_EVAP_NR_EXP)


def _covar_x_evap_comp(mu_eta, mu_chi, mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, sigma_eta, sigma_chi, sigma_rr, sigma_Nr,
                       sigma_rr_n, sigma_Nr_n, corr_eta_chi, corr_eta_rr_n, corr_eta_Nr_n, corr_chi_rr_n,
                       corr_chi_Nr_n, corr_rr_Nr_n, xm, mu_x, kk_tndcy, kk_coef, eta_tol, c, precip_frac, s):
    """Per-PDF-component Cov(x, KK_evap) for x=r_t (s=+1) / thl (s=-1) — same x'=(1/(2c))(eta'∓chi') form as
    auto/accr but the 4-var subsaturated quadrivar covariance + trivar_NLL (chi^α r_r^β N_r^γ) means."""
    a, bb, gg = KK_EVAP_SUPERSAT_EXP, KK_EVAP_RR_EXP, KK_EVAP_NR_EXP
    quad = quadrivar_NNLL_covar_eq(mu_eta, mu_chi, mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, sigma_eta, sigma_chi, sigma_rr,
                                   sigma_Nr, sigma_rr_n, sigma_Nr_n, corr_eta_chi, corr_eta_rr_n, corr_eta_Nr_n,
                                   corr_chi_rr_n, corr_chi_Nr_n, corr_rr_Nr_n, mu_eta, kk_tndcy, kk_coef,
                                   eta_tol, RR_TOL, NR_TOL, a, bb, gg)
    m_a1 = trivar_NLL_mean_eq(mu_chi, mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, sigma_chi, sigma_rr, sigma_Nr, sigma_rr_n,
                              sigma_Nr_n, corr_chi_rr_n, corr_chi_Nr_n, corr_rr_Nr_n, a + 1.0, bb, gg)
    m_a = trivar_NLL_mean_eq(mu_chi, mu_rr, mu_Nr, mu_rr_n, mu_Nr_n, sigma_chi, sigma_rr, sigma_Nr, sigma_rr_n,
                             sigma_Nr_n, corr_chi_rr_n, corr_chi_Nr_n, corr_rr_Nr_n, a, bb, gg)
    return kk_coef * precip_frac * ((1.0 / (2.0 * c)) * (quad + s * m_a1)
                                    + (mu_x - xm - s * mu_chi / (2.0 * c)) * m_a)


def _covar_x_comp(mu_eta, mu_chi, mu_y, mu_y_n, sigma_eta, sigma_chi, sigma_y, sigma_y_n,
                  corr_eta_chi, corr_eta_y_n, corr_chi_y_n, xm, mu_x, kk_tndcy, kk_coef,
                  eta_tol, y_tol, c, alpha, beta, precip_frac, s):
    """Per-PDF-component contribution to Cov(x, KK_<process>) for x = r_t (s=+1) or thl (s=-1), via the ADG1
    transform x' = (1/(2 c))(eta' ∓ chi'). Auto (y=N_cn, precip_frac=1) and accr (y=r_r, ×precip_frac) share
    this form. KK_upscaled_covariances.F90:covar_{rt,thl}_KK_{auto,accr}."""
    tri = trivar_NNL_covar_eq(mu_eta, mu_chi, mu_y, mu_y_n, sigma_eta, sigma_chi, sigma_y,
                              sigma_y_n, corr_eta_chi, corr_eta_y_n, corr_chi_y_n,
                              mu_eta, kk_tndcy, kk_coef, eta_tol, y_tol, alpha, beta)
    biv_a1 = bivar_NL_mean_eq(mu_chi, mu_y, mu_y_n, sigma_chi, sigma_y, sigma_y_n,
                              corr_chi_y_n, y_tol, alpha + 1.0, beta)
    biv_a = bivar_NL_mean_eq(mu_chi, mu_y, mu_y_n, sigma_chi, sigma_y, sigma_y_n,
                             corr_chi_y_n, y_tol, alpha, beta)
    return kk_coef * precip_frac * ((1.0 / (2.0 * c)) * (tri + s * biv_a1)
                                    + (mu_x - xm - s * mu_chi / (2.0 * c)) * biv_a)


def covar_rt_KK_auto(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                     sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                     sigma_Ncn_1_n, sigma_Ncn_2_n, corr_eta_chi_1, corr_eta_chi_2, corr_eta_Ncn_1_n,
                     corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, rtm, mu_rt_1, mu_rt_2,
                     KK_auto_tndcy, KK_auto_coef, eta_tol, crt1, crt2, mixt_frac):
    """Covariance of r_t with the KK autoconversion tendency. KK_upscaled_covariances.F90:covar_rt_KK_auto."""
    a, b = KK_AUTO_RC_EXP, KK_AUTO_NC_EXP
    c1 = _covar_x_comp(mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_eta_1, sigma_chi_1, sigma_Ncn_1,
                       sigma_Ncn_1_n, corr_eta_chi_1, corr_eta_Ncn_1_n, corr_chi_Ncn_1_n, rtm, mu_rt_1,
                       KK_auto_tndcy, KK_auto_coef, eta_tol, NC_TOL, crt1, a, b, 1.0, 1.0)
    c2 = _covar_x_comp(mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_eta_2, sigma_chi_2, sigma_Ncn_2,
                       sigma_Ncn_2_n, corr_eta_chi_2, corr_eta_Ncn_2_n, corr_chi_Ncn_2_n, rtm, mu_rt_2,
                       KK_auto_tndcy, KK_auto_coef, eta_tol, NC_TOL, crt2, a, b, 1.0, 1.0)
    return mixt_frac * c1 + (1.0 - mixt_frac) * c2


def covar_rt_KK_accr(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                     sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                     sigma_rr_1_n, sigma_rr_2_n, corr_eta_chi_1, corr_eta_chi_2, corr_eta_rr_1_n,
                     corr_eta_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, rtm, mu_rt_1, mu_rt_2,
                     KK_accr_tndcy, KK_accr_coef, eta_tol, crt1, crt2, mixt_frac, precip_frac_1, precip_frac_2):
    """Covariance of r_t with the KK accretion tendency (y=r_r, ×precip_frac). Same composition as auto.
    KK_upscaled_covariances.F90:covar_rt_KK_accr."""
    a, b = KK_ACCR_RC_EXP, KK_ACCR_RR_EXP
    c1 = _covar_x_comp(mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_eta_1, sigma_chi_1, sigma_rr_1,
                       sigma_rr_1_n, corr_eta_chi_1, corr_eta_rr_1_n, corr_chi_rr_1_n, rtm, mu_rt_1,
                       KK_accr_tndcy, KK_accr_coef, eta_tol, RR_TOL, crt1, a, b, precip_frac_1, 1.0)
    c2 = _covar_x_comp(mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_eta_2, sigma_chi_2, sigma_rr_2,
                       sigma_rr_2_n, corr_eta_chi_2, corr_eta_rr_2_n, corr_chi_rr_2_n, rtm, mu_rt_2,
                       KK_accr_tndcy, KK_accr_coef, eta_tol, RR_TOL, crt2, a, b, precip_frac_2, 1.0)
    return mixt_frac * c1 + (1.0 - mixt_frac) * c2


def covar_thl_KK_auto(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                      sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                      sigma_Ncn_1_n, sigma_Ncn_2_n, corr_eta_chi_1, corr_eta_chi_2, corr_eta_Ncn_1_n,
                      corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, thlm, mu_thl_1, mu_thl_2,
                      KK_auto_tndcy, KK_auto_coef, eta_tol, cthl1, cthl2, mixt_frac):
    """Covariance of thl with the KK autoconversion tendency (thl' = (eta'-chi')/(2 c_thl), so s=-1).
    KK_upscaled_covariances.F90:covar_thl_KK_auto."""
    a, b = KK_AUTO_RC_EXP, KK_AUTO_NC_EXP
    c1 = _covar_x_comp(mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_eta_1, sigma_chi_1, sigma_Ncn_1,
                       sigma_Ncn_1_n, corr_eta_chi_1, corr_eta_Ncn_1_n, corr_chi_Ncn_1_n, thlm, mu_thl_1,
                       KK_auto_tndcy, KK_auto_coef, eta_tol, NC_TOL, cthl1, a, b, 1.0, -1.0)
    c2 = _covar_x_comp(mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_eta_2, sigma_chi_2, sigma_Ncn_2,
                       sigma_Ncn_2_n, corr_eta_chi_2, corr_eta_Ncn_2_n, corr_chi_Ncn_2_n, thlm, mu_thl_2,
                       KK_auto_tndcy, KK_auto_coef, eta_tol, NC_TOL, cthl2, a, b, 1.0, -1.0)
    return mixt_frac * c1 + (1.0 - mixt_frac) * c2


def covar_thl_KK_accr(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                      sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                      sigma_rr_1_n, sigma_rr_2_n, corr_eta_chi_1, corr_eta_chi_2, corr_eta_rr_1_n,
                      corr_eta_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, thlm, mu_thl_1, mu_thl_2,
                      KK_accr_tndcy, KK_accr_coef, eta_tol, cthl1, cthl2, mixt_frac, precip_frac_1, precip_frac_2):
    """Covariance of thl with the KK accretion tendency (s=-1, y=r_r, ×precip_frac).
    KK_upscaled_covariances.F90:covar_thl_KK_accr."""
    a, b = KK_ACCR_RC_EXP, KK_ACCR_RR_EXP
    c1 = _covar_x_comp(mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_eta_1, sigma_chi_1, sigma_rr_1,
                       sigma_rr_1_n, corr_eta_chi_1, corr_eta_rr_1_n, corr_chi_rr_1_n, thlm, mu_thl_1,
                       KK_accr_tndcy, KK_accr_coef, eta_tol, RR_TOL, cthl1, a, b, precip_frac_1, -1.0)
    c2 = _covar_x_comp(mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_eta_2, sigma_chi_2, sigma_rr_2,
                       sigma_rr_2_n, corr_eta_chi_2, corr_eta_rr_2_n, corr_chi_rr_2_n, thlm, mu_thl_2,
                       KK_accr_tndcy, KK_accr_coef, eta_tol, RR_TOL, cthl2, a, b, precip_frac_2, -1.0)
    return mixt_frac * c1 + (1.0 - mixt_frac) * c2


def covar_x_KK_auto(mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                    sigma_x_1, sigma_x_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                    sigma_Ncn_1_n, sigma_Ncn_2_n, corr_x_chi_1, corr_x_chi_2, corr_x_Ncn_1_n,
                    corr_x_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, x_mean, KK_auto_tndcy,
                    KK_auto_coef, x_tol, mixt_frac):
    """Covariance of x (=w) with the KK autoconversion tendency — x is a direct PDF variable, so this is the
    plain KK_coef·<mixt blend of trivar_NNL_covar_eq(x, chi, N_cn)>. KK_upscaled_covariances.F90:covar_x_KK_auto."""
    a, b = KK_AUTO_RC_EXP, KK_AUTO_NC_EXP
    t1 = trivar_NNL_covar_eq(mu_x_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_x_1, sigma_chi_1, sigma_Ncn_1,
                             sigma_Ncn_1_n, corr_x_chi_1, corr_x_Ncn_1_n, corr_chi_Ncn_1_n, x_mean,
                             KK_auto_tndcy, KK_auto_coef, x_tol, NC_TOL, a, b)
    t2 = trivar_NNL_covar_eq(mu_x_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_x_2, sigma_chi_2, sigma_Ncn_2,
                             sigma_Ncn_2_n, corr_x_chi_2, corr_x_Ncn_2_n, corr_chi_Ncn_2_n, x_mean,
                             KK_auto_tndcy, KK_auto_coef, x_tol, NC_TOL, a, b)
    return KK_auto_coef * (mixt_frac * t1 + (1.0 - mixt_frac) * t2)


def covar_x_KK_accr(mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                    sigma_x_1, sigma_x_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                    sigma_rr_1_n, sigma_rr_2_n, corr_x_chi_1, corr_x_chi_2, corr_x_rr_1_n,
                    corr_x_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, x_mean, KK_accr_tndcy,
                    KK_accr_coef, x_tol, mixt_frac, precip_frac_1, precip_frac_2):
    """Covariance of x (=w) with the KK accretion tendency. Unlike auto, accretion has the
    out-of-precipitation correction `−(1−precip_frac)·(mu_x−x_mean)·KK_accr_tndcy` per component.
    KK_upscaled_covariances.F90:covar_x_KK_accr."""
    a, b = KK_ACCR_RC_EXP, KK_ACCR_RR_EXP
    t1 = trivar_NNL_covar_eq(mu_x_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_x_1, sigma_chi_1, sigma_rr_1,
                             sigma_rr_1_n, corr_x_chi_1, corr_x_rr_1_n, corr_chi_rr_1_n, x_mean,
                             KK_accr_tndcy, KK_accr_coef, x_tol, RR_TOL, a, b)
    t2 = trivar_NNL_covar_eq(mu_x_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_x_2, sigma_chi_2, sigma_rr_2,
                             sigma_rr_2_n, corr_x_chi_2, corr_x_rr_2_n, corr_chi_rr_2_n, x_mean,
                             KK_accr_tndcy, KK_accr_coef, x_tol, RR_TOL, a, b)
    comp1 = KK_accr_coef * precip_frac_1 * t1 - (1.0 - precip_frac_1) * (mu_x_1 - x_mean) * KK_accr_tndcy
    comp2 = KK_accr_coef * precip_frac_2 * t2 - (1.0 - precip_frac_2) * (mu_x_2 - x_mean) * KK_accr_tndcy
    return mixt_frac * comp1 + (1.0 - mixt_frac) * comp2


def _covar_x_KK_evap(xm, mu_x_1, mu_x_2, c1_, c2_, s, *,
                     mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2,
                     mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, sigma_chi_1,
                     sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n,
                     sigma_Nr_1_n, sigma_Nr_2_n, corr_eta_chi_1, corr_eta_chi_2, corr_eta_rr_1_n, corr_eta_rr_2_n,
                     corr_eta_Nr_1_n, corr_eta_Nr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n,
                     corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_evap_tndcy, KK_evap_coef, eta_tol,
                     mixt_frac, precip_frac_1, precip_frac_2):
    """Shared rt/thl evaporation covariance (s=+1 for r_t, s=-1 for thl)."""
    d1 = _covar_x_evap_comp(mu_eta_1, mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_eta_1, sigma_chi_1,
                            sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_eta_chi_1, corr_eta_rr_1_n,
                            corr_eta_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n, corr_rr_Nr_1_n, xm, mu_x_1,
                            KK_evap_tndcy, KK_evap_coef, eta_tol, c1_, precip_frac_1, s)
    d2 = _covar_x_evap_comp(mu_eta_2, mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_eta_2, sigma_chi_2,
                            sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_eta_chi_2, corr_eta_rr_2_n,
                            corr_eta_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n, corr_rr_Nr_2_n, xm, mu_x_2,
                            KK_evap_tndcy, KK_evap_coef, eta_tol, c2_, precip_frac_2, s)
    return mixt_frac * d1 + (1.0 - mixt_frac) * d2


def covar_rt_KK_evap(rtm, mu_rt_1, mu_rt_2, crt1, crt2, **kw):
    """Cov(r_t, KK evaporation tendency). KK_upscaled_covariances.F90:covar_rt_KK_evap. (Moment kwargs match
    `_covar_x_KK_evap`.)"""
    return _covar_x_KK_evap(rtm, mu_rt_1, mu_rt_2, crt1, crt2, 1.0, **kw)


def covar_thl_KK_evap(thlm, mu_thl_1, mu_thl_2, cthl1, cthl2, **kw):
    """Cov(thl, KK evaporation tendency) (s=-1). KK_upscaled_covariances.F90:covar_thl_KK_evap."""
    return _covar_x_KK_evap(thlm, mu_thl_1, mu_thl_2, cthl1, cthl2, -1.0, **kw)


def covar_x_KK_evap(mu_x_1, mu_x_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n,
                    mu_Nr_1_n, mu_Nr_2_n, sigma_x_1, sigma_x_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                    sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
                    corr_x_chi_1, corr_x_chi_2, corr_x_rr_1_n, corr_x_rr_2_n, corr_x_Nr_1_n, corr_x_Nr_2_n,
                    corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n,
                    corr_rr_Nr_2_n, x_mean, KK_evap_tndcy, KK_evap_coef, x_tol, mixt_frac, precip_frac_1,
                    precip_frac_2):
    """Cov(x=w, KK evaporation tendency) — direct quadrivar + the out-of-precip correction (like accr).
    KK_upscaled_covariances.F90:covar_x_KK_evap."""
    a, b, g = KK_EVAP_SUPERSAT_EXP, KK_EVAP_RR_EXP, KK_EVAP_NR_EXP
    q1 = quadrivar_NNLL_covar_eq(mu_x_1, mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_x_1, sigma_chi_1,
                                 sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_x_chi_1, corr_x_rr_1_n,
                                 corr_x_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n, corr_rr_Nr_1_n, x_mean,
                                 KK_evap_tndcy, KK_evap_coef, x_tol, RR_TOL, NR_TOL, a, b, g)
    q2 = quadrivar_NNLL_covar_eq(mu_x_2, mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_x_2, sigma_chi_2,
                                 sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_x_chi_2, corr_x_rr_2_n,
                                 corr_x_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n, corr_rr_Nr_2_n, x_mean,
                                 KK_evap_tndcy, KK_evap_coef, x_tol, RR_TOL, NR_TOL, a, b, g)
    comp1 = KK_evap_coef * precip_frac_1 * q1 - (1.0 - precip_frac_1) * (mu_x_1 - x_mean) * KK_evap_tndcy
    comp2 = KK_evap_coef * precip_frac_2 * q2 - (1.0 - precip_frac_2) * (mu_x_2 - x_mean) * KK_evap_tndcy
    return mixt_frac * comp1 + (1.0 - mixt_frac) * comp2


# Constants mirror KK_upscaled_covariances.F90 `use constants_clubb, only: Lv, Cp, w_tol, eta_tol` (eta_tol = chi_tol)
from clubb_jax.src.CLUBB_core.constants_clubb import Lv as _LV, Cp as _CP, chi_tol as _ETA_TOL, w_tol as _W_TOL


def KK_upscaled_covar_driver(
        w_mean, rtm, thlm, exner, mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_eta_1, mu_eta_2,
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        mu_Ncn_1_n, mu_Ncn_2_n, sigma_w_1, sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_eta_1, sigma_eta_2,
        sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, sigma_rr_1_n, sigma_rr_2_n,
        sigma_Nr_1_n, sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_w_chi_1, corr_w_chi_2,
        corr_w_rr_1_n, corr_w_rr_2_n, corr_w_Nr_1_n, corr_w_Nr_2_n, corr_w_Ncn_1_n, corr_w_Ncn_2_n,
        corr_chi_eta_1, corr_chi_eta_2, corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n,
        corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, corr_eta_rr_1_n, corr_eta_rr_2_n, corr_eta_Nr_1_n, corr_eta_Nr_2_n,
        corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac, precip_frac_1,
        precip_frac_2, KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy,
        mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2, crt1, crt2, cthl1, cthl2):
    """Assemble the 5 second-moment microphysics tendencies (wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc,
    rtpthlp_mc) from the 9 KK covariances. Faithful port of KK_upscaled_covariances.F90:KK_upscaled_covar_driver.
    Returns (wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc)."""
    # auto covariances (y = N_cn)
    a_w = covar_x_KK_auto(mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                          sigma_w_1, sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                          sigma_Ncn_1_n, sigma_Ncn_2_n, corr_w_chi_1, corr_w_chi_2, corr_w_Ncn_1_n,
                          corr_w_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, w_mean, KK_auto_tndcy,
                          KK_auto_coef, _W_TOL, mixt_frac)
    a_rt = covar_rt_KK_auto(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                            sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                            sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_Ncn_1_n,
                            corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, rtm, mu_rt_1, mu_rt_2,
                            KK_auto_tndcy, KK_auto_coef, _ETA_TOL, crt1, crt2, mixt_frac)
    a_thl = covar_thl_KK_auto(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
                              sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
                              sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_Ncn_1_n,
                              corr_eta_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, thlm, mu_thl_1, mu_thl_2,
                              KK_auto_tndcy, KK_auto_coef, _ETA_TOL, cthl1, cthl2, mixt_frac)
    # accretion covariances (y = r_r, ×precip_frac)
    c_w = covar_x_KK_accr(mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                          sigma_w_1, sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_rr_1_n,
                          sigma_rr_2_n, corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, corr_w_rr_2_n,
                          corr_chi_rr_1_n, corr_chi_rr_2_n, w_mean, KK_accr_tndcy, KK_accr_coef, _W_TOL,
                          mixt_frac, precip_frac_1, precip_frac_2)
    c_rt = covar_rt_KK_accr(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                            sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                            sigma_rr_1_n, sigma_rr_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n,
                            corr_eta_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, rtm, mu_rt_1, mu_rt_2,
                            KK_accr_tndcy, KK_accr_coef, _ETA_TOL, crt1, crt2, mixt_frac, precip_frac_1,
                            precip_frac_2)
    c_thl = covar_thl_KK_accr(mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n,
                              sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                              sigma_rr_1_n, sigma_rr_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n,
                              corr_eta_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, thlm, mu_thl_1, mu_thl_2,
                              KK_accr_tndcy, KK_accr_coef, _ETA_TOL, cthl1, cthl2, mixt_frac, precip_frac_1,
                              precip_frac_2)
    # evaporation covariances (4-var subsaturated, ×precip_frac)
    _evk = dict(mu_eta_1=mu_eta_1, mu_eta_2=mu_eta_2, mu_chi_1=mu_chi_1, mu_chi_2=mu_chi_2, mu_rr_1=mu_rr_1,
                mu_rr_2=mu_rr_2, mu_Nr_1=mu_Nr_1, mu_Nr_2=mu_Nr_2, mu_rr_1_n=mu_rr_1_n, mu_rr_2_n=mu_rr_2_n,
                mu_Nr_1_n=mu_Nr_1_n, mu_Nr_2_n=mu_Nr_2_n, sigma_eta_1=sigma_eta_1, sigma_eta_2=sigma_eta_2,
                sigma_chi_1=sigma_chi_1, sigma_chi_2=sigma_chi_2, sigma_rr_1=sigma_rr_1, sigma_rr_2=sigma_rr_2,
                sigma_Nr_1=sigma_Nr_1, sigma_Nr_2=sigma_Nr_2, sigma_rr_1_n=sigma_rr_1_n, sigma_rr_2_n=sigma_rr_2_n,
                sigma_Nr_1_n=sigma_Nr_1_n, sigma_Nr_2_n=sigma_Nr_2_n, corr_eta_chi_1=corr_chi_eta_1,
                corr_eta_chi_2=corr_chi_eta_2, corr_eta_rr_1_n=corr_eta_rr_1_n, corr_eta_rr_2_n=corr_eta_rr_2_n,
                corr_eta_Nr_1_n=corr_eta_Nr_1_n, corr_eta_Nr_2_n=corr_eta_Nr_2_n, corr_chi_rr_1_n=corr_chi_rr_1_n,
                corr_chi_rr_2_n=corr_chi_rr_2_n, corr_chi_Nr_1_n=corr_chi_Nr_1_n, corr_chi_Nr_2_n=corr_chi_Nr_2_n,
                corr_rr_Nr_1_n=corr_rr_Nr_1_n, corr_rr_Nr_2_n=corr_rr_Nr_2_n, KK_evap_tndcy=KK_evap_tndcy,
                KK_evap_coef=KK_evap_coef, eta_tol=_ETA_TOL, mixt_frac=mixt_frac, precip_frac_1=precip_frac_1,
                precip_frac_2=precip_frac_2)
    e_rt = covar_rt_KK_evap(rtm, mu_rt_1, mu_rt_2, crt1, crt2, **_evk)
    e_thl = covar_thl_KK_evap(thlm, mu_thl_1, mu_thl_2, cthl1, cthl2, **_evk)
    _evkx = {k: v for k, v in _evk.items()
             if k not in ('mu_eta_1', 'mu_eta_2', 'sigma_eta_1', 'sigma_eta_2', 'corr_eta_chi_1',
                          'corr_eta_chi_2', 'corr_eta_rr_1_n', 'corr_eta_rr_2_n', 'corr_eta_Nr_1_n',
                          'corr_eta_Nr_2_n', 'eta_tol')}
    e_w = covar_x_KK_evap(mu_w_1, mu_w_2, sigma_x_1=sigma_w_1, sigma_x_2=sigma_w_2, corr_x_chi_1=corr_w_chi_1,
                          corr_x_chi_2=corr_w_chi_2, corr_x_rr_1_n=corr_w_rr_1_n, corr_x_rr_2_n=corr_w_rr_2_n,
                          corr_x_Nr_1_n=corr_w_Nr_1_n, corr_x_Nr_2_n=corr_w_Nr_2_n, x_mean=w_mean, x_tol=_W_TOL,
                          **_evkx)

    w_tot = a_w + c_w + e_w
    rt_tot = a_rt + c_rt + e_rt
    thl_tot = a_thl + c_thl + e_thl
    L = _LV / (_CP * exner)
    wprtp_mc = -w_tot
    wpthlp_mc = L * w_tot
    rtp2_mc = -2.0 * rt_tot
    thlp2_mc = 2.0 * L * thl_tot
    rtpthlp_mc = L * rt_tot - thl_tot
    return wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc


def quadrivar_NNLL_covar_eq(mu_x_i, mu_chi_i, mu_rr_i, mu_Nr_i, mu_rr_i_n, mu_Nr_i_n,
                            sigma_x_i, sigma_chi_i, sigma_rr_i, sigma_Nr_i, sigma_rr_i_n, sigma_Nr_i_n,
                            corr_x_chi_i, corr_x_rr_i_n, corr_x_Nr_i_n, corr_chi_rr_i_n, corr_chi_Nr_i_n,
                            corr_rr_Nr_i_n, x_mean, mc_tndcy_mean, mc_coef, x_tol, rr_tol, Nr_tol,
                            alpha_exp, beta_exp, gamma_exp):
    """Dispatch wrapper for the quadrivariate covariance Cov_i(x, chi^α r_r^β N_r^γ) — x=w|eta, chi, r_r, N_r —
    used by the KK EVAPORATION covariances. Selects the right of the base + 11 variants by which σ≈0
    (vectorised per-level), exploiting the (x3=r_r,β)↔(x4=N_r,γ) symmetry (an x4-const branch reuses the
    x3-const variant with x3↔x4 args + β↔γ swapped). Faithful port of
    KK_upscaled_covariances.F90:quadrivar_NNLL_covar_eq."""
    mu_x1 = mu_x_i
    mu_x2 = mu_chi_i
    mu_x3 = jnp.where(beta_exp >= 0.0, mu_rr_i, jnp.maximum(mu_rr_i, rr_tol))
    mu_x4 = jnp.where(gamma_exp >= 0.0, mu_Nr_i, jnp.maximum(mu_Nr_i, Nr_tol))
    mu_x3_n, mu_x4_n = mu_rr_i_n, mu_Nr_i_n
    sg1, sg2, sg3, sg4 = sigma_x_i, sigma_chi_i, sigma_rr_i, sigma_Nr_i
    sg3n, sg4n = sigma_rr_i_n, sigma_Nr_i_n
    r12, r13n, r14n = corr_x_chi_i, corr_x_rr_i_n, corr_x_Nr_i_n
    r23n, r24n, r34n = corr_chi_rr_i_n, corr_chi_Nr_i_n, corr_rr_Nr_i_n
    x1m = x_mean
    x2a = mc_tndcy_mean / mc_coef
    x2_tol = _CHI_TOL
    a, b, g = alpha_exp, beta_exp, gamma_exp

    s2g = jnp.where(sg2 > x2_tol, sg2, 1.0)
    s_cc = jnp.where(sg2 > x2_tol, mu_x2 / s2g + r23n * sg3n * b + r24n * sg4n * g, jnp.inf)
    c1 = sg1 <= x_tol
    c2 = (sg2 <= x2_tol) | (jnp.abs(s_cc) > _PARAB_CYL_MAX)
    c3 = sg3 <= rr_tol
    c4 = sg4 <= Nr_tol

    # variants (x3-form and, for the symmetric ones, the x3<->x4 swapped form with b<->g)
    v_all = quadrivar_NNLL_covar_const_all(mu_x1, mu_x2, mu_x3, mu_x4, x1m, x2a, a, b, g)
    v123 = quadrivar_NNLL_covar_cst_x1x2x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sg4n, x1m, x2a, a, b, g)
    v123s = quadrivar_NNLL_covar_cst_x1x2x3(mu_x1, mu_x2, mu_x4, mu_x3_n, sg3n, x1m, x2a, a, g, b)
    v134 = quadrivar_NNLL_covar_cst_x1x3x4(mu_x1, mu_x2, mu_x3, mu_x4, s2g, x1m, x2a, a, b, g)
    v234 = quadrivar_NNLL_covar_cst_x2x3x4(mu_x1, mu_x2, mu_x3, mu_x4, x1m, x2a, a, b, g)
    v12 = quadrivar_NNLL_covar_const_x1x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sg3n, sg4n, r34n, x1m, x2a, a, b, g)
    v13 = quadrivar_NNLL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, mu_x4_n, s2g, sg4n, r24n, x1m, x2a, a, b, g)
    v13s = quadrivar_NNLL_covar_const_x1x3(mu_x1, mu_x2, mu_x4, mu_x3_n, s2g, sg3n, r23n, x1m, x2a, a, g, b)
    v23 = quadrivar_NNLL_covar_const_x2x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sg1, sg4n, r14n, x1m, x2a, a, b, g)
    v23s = quadrivar_NNLL_covar_const_x2x3(mu_x1, mu_x2, mu_x4, mu_x3_n, sg1, sg3n, r13n, x1m, x2a, a, g, b)
    v34 = quadrivar_NNLL_covar_const_x3x4(mu_x1, mu_x2, mu_x3, mu_x4, sg1, s2g, r12, x1m, x2a, a, b, g)
    v_x1 = quadrivar_NNLL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, mu_x4_n, s2g, sg3n, sg4n, r23n, r24n, r34n, x1m, x2a, a, b, g)
    v_x2 = quadrivar_NNLL_covar_const_x2(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sg1, sg3n, sg4n, r13n, r14n, r34n, x1m, x2a, a, b, g)
    v_x3 = quadrivar_NNLL_covar_const_x3(mu_x1, mu_x2, mu_x3, mu_x4_n, sg1, s2g, sg4n, r12, r14n, r24n, x1m, x2a, a, b, g)
    v_x3s = quadrivar_NNLL_covar_const_x3(mu_x1, mu_x2, mu_x4, mu_x3_n, sg1, s2g, sg3n, r12, r13n, r23n, x1m, x2a, a, g, b)
    v_base = quadrivar_NNLL_covar(mu_x1, mu_x2, mu_x3_n, mu_x4_n, sg1, s2g, sg3n, sg4n,
                                  r12, r13n, r14n, r23n, r24n, r34n, x1m, x2a, a, b, g)

    out = v_base
    out = jnp.where(c4, v_x3s, out)
    out = jnp.where(c3, v_x3, out)
    out = jnp.where(c2, v_x2, out)
    out = jnp.where(c1, v_x1, out)
    out = jnp.where(c3 & c4, v34, out)
    out = jnp.where(c2 & c4, v23s, out)
    out = jnp.where(c2 & c3, v23, out)
    out = jnp.where(c1 & c4, v13s, out)
    out = jnp.where(c1 & c3, v13, out)
    out = jnp.where(c1 & c2, v12, out)
    out = jnp.where(c2 & c3 & c4, v234, out)
    out = jnp.where(c1 & c3 & c4, v134, out)
    out = jnp.where(c1 & c2 & c4, v123s, out)
    out = jnp.where(c1 & c2 & c3, v123, out)
    out = jnp.where(c1 & c2 & c3 & c4, v_all, out)
    return out

def trivar_NNL_covar_eq(mu_x_i, mu_chi_i, mu_y_i, mu_y_i_n, sigma_x_i, sigma_chi_i,
                        sigma_y_i, sigma_y_i_n, corr_x_chi_i, corr_x_y_i_n, corr_chi_y_i_n,
                        x_mean, mc_tndcy_mean, mc_coef, x_tol, y_tol, alpha_exp, beta_exp):
    """Dispatch wrapper for the trivariate covariance Cov_i(x, chi^α y^β) — x=w|eta, chi, y=Ncn|rr.
    Maps the component moments to the integral inputs and selects the right form by which σ≈0
    (vectorised: per-level masks via jnp.where). Faithful port of
    KK_upscaled_covariances.F90:trivar_NNL_covar_eq."""
    mu_x1 = mu_x_i
    mu_x2 = mu_chi_i
    mu_x3 = jnp.where(beta_exp >= 0.0, mu_y_i, jnp.maximum(mu_y_i, y_tol))
    mu_x3_n = mu_y_i_n
    sigma_x1, sigma_x2, sigma_x3, sigma_x3_n = sigma_x_i, sigma_chi_i, sigma_y_i, sigma_y_i_n
    rho_x1x2, rho_x1x3_n, rho_x2x3_n = corr_x_chi_i, corr_x_y_i_n, corr_chi_y_i_n
    x1m = x_mean
    x2a = mc_tndcy_mean / mc_coef
    x2_tol = _CHI_TOL

    # Guard the σ_x2 denominator for the forms that use mu_x2/σ_x2 (only SELECTED when σ_x2>tol).
    s2g = jnp.where(sigma_x2 > x2_tol, sigma_x2, 1.0)
    s_c = jnp.where(sigma_x2 > x2_tol, mu_x2 / s2g + rho_x2x3_n * sigma_x3_n * beta_exp, jnp.inf)
    c1 = sigma_x1 <= x_tol
    c2 = (sigma_x2 <= x2_tol) | (jnp.abs(s_c) > _PARAB_CYL_MAX)
    c3 = sigma_x3 <= y_tol

    v_x1x2 = trivar_NNL_covar_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n, x1m, x2a, alpha_exp, beta_exp)
    v_x1x3 = trivar_NNL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, s2g, x1m, x2a, alpha_exp, beta_exp)
    v_x2x3 = trivar_NNL_covar_const_x2x3(mu_x1, mu_x2, mu_x3, x1m, x2a, alpha_exp, beta_exp)
    v_x1 = trivar_NNL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, s2g, sigma_x3_n, rho_x2x3_n, x1m, x2a, alpha_exp, beta_exp)
    v_x2 = trivar_NNL_covar_const_x2(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n, rho_x1x3_n, x1m, x2a, alpha_exp, beta_exp)
    v_x3 = trivar_NNL_covar_const_x3(mu_x1, mu_x2, mu_x3, sigma_x1, s2g, rho_x1x2, x1m, x2a, alpha_exp, beta_exp)
    v_base = trivar_NNL_covar(mu_x1, mu_x2, mu_x3_n, sigma_x1, s2g, sigma_x3_n,
                              rho_x1x2, rho_x1x3_n, rho_x2x3_n, x1m, x2a, alpha_exp, beta_exp)

    out = v_base
    out = jnp.where(c3, v_x3, out)
    out = jnp.where(c2, v_x2, out)
    out = jnp.where(c1, v_x1, out)
    out = jnp.where(c2 & c3, v_x2x3, out)
    out = jnp.where(c1 & c3, v_x1x3, out)
    out = jnp.where(c1 & c2, v_x1x2, out)
    out = jnp.where(c1 & c2 & c3, v_x2x3, out)   # const_all == const_x2x3
    return out
