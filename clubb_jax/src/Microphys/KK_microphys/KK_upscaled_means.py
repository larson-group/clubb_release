"""JAX port of KK_upscaled_means.F90 — upscaled (analytic-PDF) KK rate means.

This iteration ports the AUTOCONVERSION mean:
  KK_auto_upscaled_mean  (KK_upscaled_means.F90:289) and its per-component dispatch
  bivar_NL_mean_eq       (KK_upscaled_means.F90:864).

The upscaled KK scheme (l_local_kk=.false., used by every precipitating CLUBB case)
integrates the local Khairoutdinov-Kogan rates analytically over CLUBB's two-component
chi / N_cn / r_r PDF. The autoconversion mean is

  < KK_auto > = KK_auto_coef * ( mixt_frac * < chi^alpha N_cn^beta >_1
                                 + (1-mixt_frac) * < chi^alpha N_cn^beta >_2 ),

with the KK exponents alpha = KK_auto_rc_exp = 2.47, beta = KK_auto_Nc_exp = -1.79, and
each component mean given by bivar_NL_mean_eq (one of four reduced/general NL forms,
selected by which component variances vanish). See PDF_integrals_means.py.

Constants (parameters_KK.F90, constants_clubb.F90):
  KK_auto_rc_exp = 2.47, KK_auto_Nc_exp = -1.79
  chi_tol = max(1e-8, eps) = 1e-8 ;  Nc_tol = 1e2 ;  parab_cyl_max_input = 49
  KK_auto_coef = 1350 * (rho / cm3_per_m3)^KK_auto_Nc_exp   (cm3_per_m3 = 1e6)
"""
import jax.numpy as jnp

# chi_tol/rho_lw/mvr_rain_max/Nc_tol/rr_tol/parab_cyl_max_input all mirror constants_clubb
# (KK_upscaled_means.F90 `use constants_clubb`).
from clubb_jax.src.CLUBB_core.constants_clubb import (
    chi_tol as CHI_TOL, rho_lw as _RHO_LW, mvr_rain_max as _MVR_RAIN_MAX,
    Nc_tol as NC_TOL, rr_tol as RR_TOL, parab_cyl_max_input as PARAB_CYL_MAX_INPUT,
)

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_NL_mean,
    bivar_NL_mean_const_x1,
    bivar_NL_mean_const_x2,
    bivar_NL_mean_const_all,
    trivar_NLL_mean,
    trivar_NLL_mean_const_x1,
    trivar_NLL_mean_const_x2,
    trivar_NLL_mean_const_x1x2,
    trivar_NLL_mean_const_x2x3,
    trivar_NLL_mean_const_all,
    bivar_LL_mean,
    bivar_LL_mean_const_x1,
    bivar_LL_mean_const_all,
)

# --- parameters --- the KK exponents live in their Fortran home parameters_KK.F90 (`use parameters_KK`) ----
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import (
    KK_auto_rc_exp as KK_AUTO_RC_EXP, KK_auto_Nc_exp as KK_AUTO_NC_EXP,
    KK_accr_rc_exp as KK_ACCR_RC_EXP, KK_accr_rr_exp as KK_ACCR_RR_EXP,
    KK_evap_Supersat_exp as KK_EVAP_SUPERSAT_EXP, KK_evap_rr_exp as KK_EVAP_RR_EXP,
    KK_evap_Nr_exp as KK_EVAP_NR_EXP, KK_mvr_rr_exp as KK_MVR_RR_EXP, KK_mvr_Nr_exp as KK_MVR_NR_EXP,
)
# KK_ACCR_COEF/KK_MVR_COEF live in their Fortran home KK_microphys_module (KK_tendency_coefs, F90:1185/1188).
from clubb_jax.src.Microphys.KK_microphys_module import KK_ACCR_COEF, KK_MVR_COEF
# Nr_tol = rr_tol / ((4/3) pi rho_lw mvr_rain_max^3). constants_clubb.F90:306  (rho_lw/mvr_rain_max imported above)
NR_TOL = RR_TOL / ((4.0 / 3.0) * jnp.pi * _RHO_LW * _MVR_RAIN_MAX ** 3)  # ~1.9099e-7
# kk_auto_coef lives in its Fortran-home Microphys/KK_microphys_module.py (computed inline at
# KK_microphys_module.F90:1182), not here — KK_upscaled_means.F90 only takes KK_auto_coef as an input arg.


def bivar_NL_mean_eq(mu_chi_i, mu_y_i, mu_y_i_n, sigma_chi_i,
                     sigma_y_i, sigma_y_i_n, corr_chi_y_i_n,
                     y_tol, alpha_exp, beta_exp):
    """Per-component < chi^alpha y^beta >_i with the 4-way variance dispatch.

    KK_upscaled_means.F90:864. y = N_cn (auto) or r_r (accr). Selects:
      both vary                -> bivar_NL_mean
      sigma_x1 ~ 0 or |s_c|>49 -> bivar_NL_mean_const_x1
      sigma_x2 ~ 0             -> bivar_NL_mean_const_x2
      both ~ 0                 -> bivar_NL_mean_const_all
    Implemented branch-free (all forms computed, selected with jnp.where) so it
    vmaps/jits over a grid; divisions in the varying form use a clamped sigma_x1
    so the unused (const) branch never poisons values/gradients with inf/nan.
    """
    mu_x1 = mu_chi_i
    mu_x2 = jnp.where(beta_exp >= 0.0, mu_y_i, jnp.maximum(mu_y_i, y_tol))
    mu_x2_n = mu_y_i_n
    sigma_x1 = sigma_chi_i
    sigma_x2 = sigma_y_i
    sigma_x2_n = sigma_y_i_n
    rho = corr_chi_y_i_n

    # s_c for the parabolic-cylinder input; clamp sigma_x1 away from 0 in the
    # division so the const-x1 branch (selected when sigma_x1 <= CHI_TOL) is clean.
    sig1_safe = jnp.maximum(sigma_x1, CHI_TOL)
    s_c = (mu_x1 / sig1_safe) + rho * sigma_x2_n * beta_exp

    x1_const = (sigma_x1 <= CHI_TOL) | (jnp.abs(s_c) > PARAB_CYL_MAX_INPUT)
    x2_const = sigma_x2 <= y_tol

    f_all = bivar_NL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp)
    f_x1 = bivar_NL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp)
    f_x2 = bivar_NL_mean_const_x2(mu_x1, mu_x2, sig1_safe, alpha_exp, beta_exp)
    f_gen = bivar_NL_mean(mu_x1, mu_x2_n, sig1_safe, sigma_x2_n, rho, alpha_exp, beta_exp)

    # Mirror the Fortran if/elseif/elseif/else ladder.
    return jnp.where(x1_const & x2_const, f_all,
                     jnp.where(x1_const, f_x1,
                               jnp.where(x2_const, f_x2, f_gen)))


def KK_auto_upscaled_mean(mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2,
                          mu_Ncn_1_n, mu_Ncn_2_n, sigma_chi_1, sigma_chi_2,
                          sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n,
                          corr_chi_Ncn_1_n, corr_chi_Ncn_2_n,
                          KK_auto_coef_val, mixt_frac):
    """Mean upscaled KK autoconversion tendency. KK_upscaled_means.F90:289.

    KK_auto_coef_val is KK_auto_coef (a function of rho; see kk_auto_coef)."""
    alpha = KK_AUTO_RC_EXP
    beta = KK_AUTO_NC_EXP
    comp1 = bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1,
                             sigma_Ncn_1, sigma_Ncn_1_n, corr_chi_Ncn_1_n,
                             NC_TOL, alpha, beta)
    comp2 = bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2,
                             sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n,
                             NC_TOL, alpha, beta)
    return KK_auto_coef_val * (mixt_frac * comp1 + (1.0 - mixt_frac) * comp2)


def KK_accr_upscaled_mean(mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2,
                          mu_rr_1_n, mu_rr_2_n, sigma_chi_1, sigma_chi_2,
                          sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n,
                          corr_chi_rr_1_n, corr_chi_rr_2_n,
                          mixt_frac, precip_frac_1, precip_frac_2):
    """Mean upscaled KK rain-water accretion tendency. KK_upscaled_means.F90:386.

    Same NL bivariate dispatch as autoconversion but with y = r_r (in-precip,
    tolerance rr_tol), exponents alpha=beta=KK_accr_rc/rr_exp=1.15, a CONSTANT
    coefficient KK_accr_coef=67, and a per-component precip_frac factor (r_r exists
    only in precipitation)."""
    alpha = KK_ACCR_RC_EXP
    beta = KK_ACCR_RR_EXP
    comp1 = bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1,
                             sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n,
                             RR_TOL, alpha, beta)
    comp2 = bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2,
                             sigma_rr_2, sigma_rr_2_n, corr_chi_rr_2_n,
                             RR_TOL, alpha, beta)
    return KK_ACCR_COEF * (mixt_frac * precip_frac_1 * comp1
                           + (1.0 - mixt_frac) * precip_frac_2 * comp2)


def trivar_NLL_mean_eq(mu_chi_i, mu_rr_i, mu_Nr_i, mu_rr_i_n, mu_Nr_i_n,
                       sigma_chi_i, sigma_rr_i, sigma_Nr_i, sigma_rr_i_n,
                       sigma_Nr_i_n, corr_chi_rr_i_n, corr_chi_Nr_i_n,
                       corr_rr_Nr_i_n, alpha_exp, beta_exp, gamma_exp):
    """Per-component < chi^alpha r_r^beta N_r^gamma >_i with the 8-way variance dispatch.

    KK_upscaled_means.F90:586. x1=chi, x2=r_r, x3=N_r. Selects one of 6 trivariate
    forms (const_x1x2 and const_x2 are reused with swapped (x2,x3) args for the
    chi+Nr and Nr-only cases). Branch-free via jnp.where; divisions use a clamped
    sigma_x1 so the const-x1 branches never poison values/gradients."""
    mu_x1 = mu_chi_i
    mu_x2 = jnp.where(beta_exp >= 0.0, mu_rr_i, jnp.maximum(mu_rr_i, RR_TOL))
    mu_x3 = jnp.where(gamma_exp >= 0.0, mu_Nr_i, jnp.maximum(mu_Nr_i, NR_TOL))
    mu_x2_n, mu_x3_n = mu_rr_i_n, mu_Nr_i_n
    sigma_x1, sigma_x2, sigma_x3 = sigma_chi_i, sigma_rr_i, sigma_Nr_i
    sigma_x2_n, sigma_x3_n = sigma_rr_i_n, sigma_Nr_i_n
    rho12, rho13, rho23 = corr_chi_rr_i_n, corr_chi_Nr_i_n, corr_rr_Nr_i_n

    sig1_safe = jnp.maximum(sigma_x1, CHI_TOL)
    s_cc = ((mu_x1 / sig1_safe)
            + rho12 * sigma_x2_n * beta_exp + rho13 * sigma_x3_n * gamma_exp)

    x1_const = (sigma_x1 <= CHI_TOL) | (jnp.abs(s_cc) > PARAB_CYL_MAX_INPUT)
    x2_const = sigma_x2 <= RR_TOL
    x3_const = sigma_x3 <= NR_TOL

    f_all = trivar_NLL_mean_const_all(mu_x1, mu_x2, mu_x3, alpha_exp, beta_exp, gamma_exp)
    f_x1x2 = trivar_NLL_mean_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n,
                                        alpha_exp, beta_exp, gamma_exp)
    f_x1x3 = trivar_NLL_mean_const_x1x2(mu_x1, mu_x3, mu_x2_n, sigma_x2_n,
                                        alpha_exp, gamma_exp, beta_exp)
    f_x2x3 = trivar_NLL_mean_const_x2x3(mu_x1, mu_x2, mu_x3, sig1_safe,
                                        alpha_exp, beta_exp, gamma_exp)
    f_x1 = trivar_NLL_mean_const_x1(mu_x1, mu_x2_n, mu_x3_n, sigma_x2_n, sigma_x3_n,
                                    rho23, alpha_exp, beta_exp, gamma_exp)
    f_x2 = trivar_NLL_mean_const_x2(mu_x1, mu_x2, mu_x3_n, sig1_safe, sigma_x3_n,
                                    rho13, alpha_exp, beta_exp, gamma_exp)
    f_x3 = trivar_NLL_mean_const_x2(mu_x1, mu_x3, mu_x2_n, sig1_safe, sigma_x2_n,
                                    rho12, alpha_exp, gamma_exp, beta_exp)
    f_gen = trivar_NLL_mean(mu_x1, mu_x2_n, mu_x3_n, sig1_safe, sigma_x2_n, sigma_x3_n,
                            rho12, rho13, rho23, alpha_exp, beta_exp, gamma_exp)

    # Mirror the Fortran if/elseif ladder (top-priority condition first).
    return jnp.where(x1_const & x2_const & x3_const, f_all,
            jnp.where(x1_const & x2_const, f_x1x2,
             jnp.where(x1_const & x3_const, f_x1x3,
              jnp.where(x2_const & x3_const, f_x2x3,
               jnp.where(x1_const, f_x1,
                jnp.where(x2_const, f_x2,
                 jnp.where(x3_const, f_x3, f_gen)))))))


def KK_evap_upscaled_mean(mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2,
                          mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
                          sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
                          sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n,
                          sigma_Nr_1_n, sigma_Nr_2_n,
                          corr_chi_rr_1_n, corr_chi_rr_2_n,
                          corr_chi_Nr_1_n, corr_chi_Nr_2_n,
                          corr_rr_Nr_1_n, corr_rr_Nr_2_n,
                          KK_evap_coef_val, mixt_frac, precip_frac_1, precip_frac_2):
    """Mean upscaled KK rain-water evaporation tendency. KK_upscaled_means.F90:167.

    < KK_evap > = KK_evap_coef * ( mixt_frac*precip_frac_1*<chi^1 r_r^(1/3) N_r^(2/3)>_1
                                   + (1-mixt_frac)*precip_frac_2*<...>_2 ), integrated over
    chi<0 (subsaturated). KK_evap_coef_val (=3 C_evap G_T_p(T,p)) is passed in (its G_T_p
    thermodynamic factor is a separate dependency, like KK_accr/auto coefficients)."""
    a, b, g = KK_EVAP_SUPERSAT_EXP, KK_EVAP_RR_EXP, KK_EVAP_NR_EXP
    comp1 = trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n,
                               sigma_chi_1, sigma_rr_1, sigma_Nr_1, sigma_rr_1_n,
                               sigma_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n,
                               corr_rr_Nr_1_n, a, b, g)
    comp2 = trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                               sigma_chi_2, sigma_rr_2, sigma_Nr_2, sigma_rr_2_n,
                               sigma_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n,
                               corr_rr_Nr_2_n, a, b, g)
    return KK_evap_coef_val * (mixt_frac * precip_frac_1 * comp1
                               + (1.0 - mixt_frac) * precip_frac_2 * comp2)


def bivar_LL_mean_eq(mu_rr_i, mu_Nr_i, mu_rr_i_n, mu_Nr_i_n, sigma_rr_i,
                     sigma_Nr_i, sigma_rr_i_n, sigma_Nr_i_n, corr_rr_Nr_i_n,
                     alpha_exp, beta_exp):
    """Per-component < r_r^alpha N_r^beta >_i with the 4-way variance dispatch.

    KK_upscaled_means.F90:1064. x1=r_r, x2=N_r (both lognormal). When sigma_x2 is the
    constant one, const_x1 is reused with (x2,x1) args swapped. Branch-free via jnp.where."""
    mu_x1 = jnp.where(alpha_exp >= 0.0, mu_rr_i, jnp.maximum(mu_rr_i, RR_TOL))
    mu_x2 = jnp.where(beta_exp >= 0.0, mu_Nr_i, jnp.maximum(mu_Nr_i, NR_TOL))
    mu_x1_n, mu_x2_n = mu_rr_i_n, mu_Nr_i_n
    sigma_x1, sigma_x2 = sigma_rr_i, sigma_Nr_i
    sigma_x1_n, sigma_x2_n = sigma_rr_i_n, sigma_Nr_i_n
    rho = corr_rr_Nr_i_n

    x1_const = sigma_x1 <= RR_TOL
    x2_const = sigma_x2 <= NR_TOL

    f_all = bivar_LL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp)
    f_x1 = bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp)
    f_x2 = bivar_LL_mean_const_x1(mu_x2, mu_x1_n, sigma_x1_n, beta_exp, alpha_exp)
    f_gen = bivar_LL_mean(mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho, alpha_exp, beta_exp)

    return jnp.where(x1_const & x2_const, f_all,
                     jnp.where(x1_const, f_x1,
                               jnp.where(x2_const, f_x2, f_gen)))


def KK_mvr_upscaled_mean(mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n,
                         mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, sigma_rr_2, sigma_Nr_1,
                         sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
                         sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n,
                         mixt_frac, precip_frac_1, precip_frac_2):
    """Mean upscaled KK rain-drop mean volume radius. KK_upscaled_means.F90:483.

    < KK_mvr > = KK_mvr_coef * ( mixt_frac*precip_frac_1*<r_r^(1/3) N_r^(-1/3)>_1 + ... ),
    a bivariate-lognormal mean (both r_r, N_r lognormal). KK_mvr_coef=((4/3)pi rho_lw)^(-1/3)."""
    a, b = KK_MVR_RR_EXP, KK_MVR_NR_EXP
    comp1 = bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                             sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n, a, b)
    comp2 = bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                             sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n, a, b)
    return KK_MVR_COEF * (mixt_frac * precip_frac_1 * comp1
                          + (1.0 - mixt_frac) * precip_frac_2 * comp2)
