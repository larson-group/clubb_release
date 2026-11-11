// Constants and case-owned adjustable parameters for KK microphysics.
#pragma once
#include "CLUBB_core/constants_clubb.hpp"
namespace clubb_kokkos {
inline constexpr real_t KK_evap_Supersat_exp = one;
inline constexpr real_t KK_evap_rr_exp       = one_third;
inline constexpr real_t KK_evap_Nr_exp       = two_thirds;
inline constexpr real_t KK_auto_rc_exp       = real_t(2.47);
inline constexpr real_t KK_auto_Nc_exp       = real_t(-1.79);
inline constexpr real_t KK_accr_rc_exp       = real_t(1.15);
inline constexpr real_t KK_accr_rr_exp       = real_t(1.15);
inline constexpr real_t KK_mvr_rr_exp        = one_third;
inline constexpr real_t KK_mvr_Nr_exp        = -one_third;
inline constexpr real_t KK_Nrm_evap_nu       = one;
struct parameters_KK {
  real_t r_0    = real_t(25.0e-6);
  real_t C_evap = real_t(0.86);
  // Native Parabolic accuracy selector; case owned, not a namelist option.
  bool l_high_accuracy_parab_cyl_fnc = false;
};
} // namespace clubb_kokkos
