// Native KK_upscaled_variances numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline real_t variance_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_mean_vol_rad,
                              real_t KK_mvr_coef, real_t mixt_frac, real_t precip_frac_1,
                              real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the variance of KK mean volume radius of rain
  // drops (R_vr), which can be written as < R_vr'^2 >.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t variance_KK_mvr;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the covariance of r_r and KK mean volume radius of rain drops.
  variance_KK_mvr = (((KK_mvr_coef * KK_mvr_coef) *
                      (((mixt_frac * precip_frac_1) *
                        bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                                         sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n,
                                         (two * alpha_exp), (two * beta_exp), kk)) +
                       (((one - mixt_frac) * precip_frac_2) *
                        bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                                         sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n,
                                         (two * alpha_exp), (two * beta_exp), kk)))) -
                     (KK_mean_vol_rad * KK_mean_vol_rad));

  return variance_KK_mvr;
}
} // namespace clubb_kokkos::kk_engine
