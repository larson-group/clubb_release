// Native KK_local_means numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline real_t KK_evap_local_mean(real_t mean_s, real_t rrm, real_t Nrm, real_t KK_evap_coef,
                                 const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (22) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_evap_local_mean;
  real_t alpha_exp;
  real_t beta_exp;
  real_t gamma_exp;

  // Values of the KK exponents.
  alpha_exp = KK_evap_Supersat_exp;
  beta_exp  = KK_evap_rr_exp;
  gamma_exp = KK_evap_Nr_exp;

  // Calculate the local KK evaporation tendency.
  if ((mean_s <= zero)) {

    KK_evap_local_mean =
        (((KK_evap_coef * Kokkos::pow(mean_s, alpha_exp)) * Kokkos::pow(rrm, beta_exp)) *
         Kokkos::pow(Nrm, gamma_exp));

  } else {

    // The air is cloudy, and rain does not evaporate.
    KK_evap_local_mean = zero;
  }

  return KK_evap_local_mean;
}
inline real_t KK_auto_local_mean(real_t mean_s, real_t Ncm, real_t KK_auto_coef,
                                 const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (29) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_auto_local_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_auto_rc_exp;
  beta_exp  = KK_auto_Nc_exp;

  // Calculate the local KK autoconversion tendency.
  if ((mean_s > zero)) {

    KK_auto_local_mean =
        ((KK_auto_coef * Kokkos::pow(mean_s, alpha_exp)) * Kokkos::pow(Ncm, beta_exp));

  } else {

    // The air is clear, and rain cannot be produced.
    KK_auto_local_mean = zero;
  }

  return KK_auto_local_mean;
}
inline real_t KK_accr_local_mean(real_t mean_s, real_t rrm, real_t KK_accr_coef,
                                 const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (33) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_accr_local_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_accr_rc_exp;
  beta_exp  = KK_accr_rr_exp;

  // Calculate the local KK accretion tendency.
  if ((mean_s > zero)) {

    KK_accr_local_mean =
        ((KK_accr_coef * Kokkos::pow(mean_s, alpha_exp)) * Kokkos::pow(rrm, beta_exp));

  } else {

    // The air is clear, and rain cannot be produced.
    KK_accr_local_mean = zero;
  }

  return KK_accr_local_mean;
}
inline real_t KK_mvr_local_mean(real_t rrm, real_t Nrm, real_t KK_mvr_coef,
                                const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (3) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
  // Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
  // Mon. Wea. Rev., 128, 229--243.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_mvr_local_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the KK mean volume radius of rain drops
  KK_mvr_local_mean = ((KK_mvr_coef * Kokkos::pow(rrm, alpha_exp)) *
                       Kokkos::pow(((Nrm) > (Nr_tol) ? (Nrm) : (Nr_tol)), beta_exp));

  return KK_mvr_local_mean;
}
} // namespace clubb_kokkos::kk_engine
