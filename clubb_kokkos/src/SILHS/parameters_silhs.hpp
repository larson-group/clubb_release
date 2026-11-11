// Case-owned SILHS sampling parameters and Native configuration flags.
// The three flag APIs set defaults, assemble the flag value and print it.
// print_silhs_config_flags_api writes Native rows to a caller-owned C++ stream.
// Prescribed probabilities and decorrelation settings belong to one caller;
// immutable allocation selectors and the inverse-CDF cutoff are shared constants.
#pragma once
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <iosfwd>
#include "clubb_real.hpp"

namespace clubb_kokkos {
inline constexpr int eight_cluster_allocation_opt = 1;
inline constexpr int four_cluster_allocation_opt  = 2;
inline constexpr int two_cluster_cp_nocp_opt      = 3;

struct eight_cluster_presc_probs_type {
  real_t cloud_precip_comp1     = real_t(0.15);
  real_t cloud_precip_comp2     = real_t(0.15);
  real_t nocloud_precip_comp1   = real_t(0.15);
  real_t nocloud_precip_comp2   = real_t(0.15);
  real_t cloud_noprecip_comp1   = real_t(0.15);
  real_t cloud_noprecip_comp2   = real_t(0.15);
  real_t nocloud_noprecip_comp1 = real_t(0.05);
  real_t nocloud_noprecip_comp2 = real_t(0.05);
};

// BEGIN GENERATED SILHS CONFIG FLAGS TYPE
// Mirrors src/SILHS/parameters_silhs.F90; arrays use Fortran-order LayoutLeft.
struct silhs_config_flags_type {
  int cluster_allocation_strategy{};
  std::uint8_t l_lh_importance_sampling{};
  std::uint8_t l_Lscale_vert_avg{};
  std::uint8_t l_lh_straight_mc{};
  std::uint8_t l_lh_clustered_sampling{};
  std::uint8_t l_rcm_in_cloud_k_lh_start{};
  std::uint8_t l_random_k_lh_start{};
  std::uint8_t l_max_overlap_in_cloud{};
  std::uint8_t l_lh_instant_var_covar_src{};
  std::uint8_t l_lh_limit_weights{};
  std::uint8_t l_lh_var_frac{};
  std::uint8_t l_lh_normalize_weights{};
  std::uint8_t l_lh_deterministic_test{};
};

// END GENERATED SILHS CONFIG FLAGS TYPE

struct parameters_silhs_state {
  eight_cluster_presc_probs_type eight_cluster_presc_probs;
  real_t importance_prob_thresh = real_t(1.e-8); // Minimum category PDF probability [-].
  real_t vert_decorr_coef       = real_t(0.03); // Empirical vertical decorrelation coefficient [-].
};

// Inverse-CDF inputs lie in [single_prec_thresh,1-single_prec_thresh], since the
// algorithm only guarantees single-precision accuracy even with double inputs.
inline constexpr real_t single_prec_thresh = real_t(3.e-8);

void set_default_silhs_config_flags_api(int& cluster_allocation_strategy, // Out
                                        bool& l_lh_importance_sampling,   // Out
                                        bool& l_Lscale_vert_avg,          // Out
                                        bool& l_lh_straight_mc,           // Out
                                        bool& l_lh_clustered_sampling,    // Out
                                        bool& l_rcm_in_cloud_k_lh_start,  // Out
                                        bool& l_random_k_lh_start,        // Out
                                        bool& l_max_overlap_in_cloud,     // Out
                                        bool& l_lh_instant_var_covar_src, // Out
                                        bool& l_lh_limit_weights,         // Out
                                        bool& l_lh_var_frac,              // Out
                                        bool& l_lh_normalize_weights,     // Out
                                        bool& l_lh_deterministic_test     // Out
);                                                                        // Out

void initialize_silhs_config_flags_type_api(int cluster_allocation_strategy,            // In
                                            bool l_lh_importance_sampling,              // In
                                            bool l_Lscale_vert_avg,                     // In
                                            bool l_lh_straight_mc,                      // In
                                            bool l_lh_clustered_sampling,               // In
                                            bool l_rcm_in_cloud_k_lh_start,             // In
                                            bool l_random_k_lh_start,                   // In
                                            bool l_max_overlap_in_cloud,                // In
                                            bool l_lh_instant_var_covar_src,            // In
                                            bool l_lh_limit_weights,                    // In
                                            bool l_lh_var_frac,                         // In
                                            bool l_lh_normalize_weights,                // In
                                            bool l_lh_deterministic_test,               // In
                                            silhs_config_flags_type& silhs_config_flags // Out
);                                                                                      // Out

// iunit is the C++ output stream replacing the Native logical-unit selector.
void print_silhs_config_flags_api(std::ostream& iunit,
                                  const silhs_config_flags_type& silhs_config_flags);

} // namespace clubb_kokkos
