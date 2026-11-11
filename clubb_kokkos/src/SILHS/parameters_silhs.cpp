// Defaults, assembly and printing for caller-owned SILHS flags.
// Public functions mirror the three same-named Native configuration APIs.
// print_silhs_config_flags_api writes Native rows to a caller-owned C++ stream.
#include "SILHS/parameters_silhs.hpp"
#include <ostream>
#include <sstream>

namespace clubb_kokkos {

//-----------------------------------------------------------------------
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
) {
  // Set all SILHS flags to their Native default settings.
  // cluster_allocation_strategy: Two clusters, one containing all categories with either cloud or
  //   precip, and the other containing categories with neither
  // l_lh_importance_sampling: Limit noise by performing importance sampling
  // l_Lscale_vert_avg: Calculate Lscale_vert_avg in generate_silhs_sample
  // l_lh_straight_mc: Use true Monte Carlo sampling with no Latin hypercube sampling and no
  //   importance sampling
  // l_lh_clustered_sampling: Use the "new" SILHS importance sampling scheme with prescribed
  //   probabilities
  // l_rcm_in_cloud_k_lh_start: Determine k_lh_start based on maximum within-cloud rcm
  // l_random_k_lh_start: Place k_lh_start at a random grid level between maximum rcm and maximum
  //   rcm_in_cloud
  // l_max_overlap_in_cloud: Assume maximum vertical overlap when grid-box rcm exceeds cloud
  //   threshold
  // l_lh_instant_var_covar_src: Produces "instantaneous" variance-covariance microphysical source
  //   terms, ignoring discretization effects
  // l_lh_limit_weights: Limit SILHS sample point weights for stability
  // l_lh_var_frac: Prescribe variance fractions
  // l_lh_normalize_weights: Scale sample point weights to sum to num_samples (the "ratio estimate")
  // l_lh_deterministic_test: Use deterministic sampling inputs for repeatable testing
  cluster_allocation_strategy = two_cluster_cp_nocp_opt;
  l_lh_importance_sampling    = true;
  l_Lscale_vert_avg           = false;
  l_lh_straight_mc            = false;
  l_lh_clustered_sampling     = true;
  l_rcm_in_cloud_k_lh_start   = true;
  l_random_k_lh_start         = false;
  l_max_overlap_in_cloud      = true;
  l_lh_instant_var_covar_src  = true;
  l_lh_limit_weights          = true;
  l_lh_var_frac               = false;
  l_lh_normalize_weights      = true;
  l_lh_deterministic_test     = false;
}

//-----------------------------------------------------------------------
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
) {
  // Assemble the SILHS flags value from all caller-supplied settings.
  // cluster_allocation_strategy: Two clusters, one containing all categories with either cloud or
  //   precip, and the other containing categories with neither
  // l_lh_importance_sampling: Limit noise by performing importance sampling
  // l_Lscale_vert_avg: Calculate Lscale_vert_avg in generate_silhs_sample
  // l_lh_straight_mc: Use true Monte Carlo sampling with no Latin hypercube sampling and no
  //   importance sampling
  // l_lh_clustered_sampling: Use the "new" SILHS importance sampling scheme with prescribed
  //   probabilities
  // l_rcm_in_cloud_k_lh_start: Determine k_lh_start based on maximum within-cloud rcm
  // l_random_k_lh_start: Place k_lh_start at a random grid level between maximum rcm and maximum
  //   rcm_in_cloud
  // l_max_overlap_in_cloud: Assume maximum vertical overlap when grid-box rcm exceeds cloud
  //   threshold
  // l_lh_instant_var_covar_src: Produces "instantaneous" variance-covariance microphysical source
  //   terms, ignoring discretization effects
  // l_lh_limit_weights: Limit SILHS sample point weights for stability
  // l_lh_var_frac: Prescribe variance fractions
  // l_lh_normalize_weights: Scale sample point weights to sum to num_samples (the "ratio estimate")
  // l_lh_deterministic_test: Use deterministic sampling inputs for repeatable testing
  // silhs_config_flags: Derived type holding all configurable SILHS flags
  silhs_config_flags.cluster_allocation_strategy = cluster_allocation_strategy;
  silhs_config_flags.l_lh_importance_sampling    = l_lh_importance_sampling;
  silhs_config_flags.l_Lscale_vert_avg           = l_Lscale_vert_avg;
  silhs_config_flags.l_lh_straight_mc            = l_lh_straight_mc;
  silhs_config_flags.l_lh_clustered_sampling     = l_lh_clustered_sampling;
  silhs_config_flags.l_rcm_in_cloud_k_lh_start   = l_rcm_in_cloud_k_lh_start;
  silhs_config_flags.l_random_k_lh_start         = l_random_k_lh_start;
  silhs_config_flags.l_max_overlap_in_cloud      = l_max_overlap_in_cloud;
  silhs_config_flags.l_lh_instant_var_covar_src  = l_lh_instant_var_covar_src;
  silhs_config_flags.l_lh_limit_weights          = l_lh_limit_weights;
  silhs_config_flags.l_lh_var_frac               = l_lh_var_frac;
  silhs_config_flags.l_lh_normalize_weights      = l_lh_normalize_weights;
  silhs_config_flags.l_lh_deterministic_test     = l_lh_deterministic_test;
}

//-----------------------------------------------------------------------
void print_silhs_config_flags_api(std::ostream& iunit,
                                  const silhs_config_flags_type& silhs_config_flags) {
  // Description:
  //   Print the thirteen Native SILHS flag rows to the caller-owned output stream.
  // References: None.
  // iunit: Output stream replacing the Native logical-unit selector.
  // silhs_config_flags: Caller-owned configurable SILHS flags [-].
  // Render decimal integers and T/F logicals independently of caller stream flags.
  std::ostringstream output;
  output << "cluster_allocation_strategy = " << silhs_config_flags.cluster_allocation_strategy
         << '\n';
  output << "l_lh_importance_sampling = "
         << (silhs_config_flags.l_lh_importance_sampling ? "T" : "F") << '\n';
  output << "l_Lscale_vert_avg = " << (silhs_config_flags.l_Lscale_vert_avg ? "T" : "F") << '\n';
  output << "l_lh_straight_mc = " << (silhs_config_flags.l_lh_straight_mc ? "T" : "F") << '\n';
  output << "l_lh_clustered_sampling = " << (silhs_config_flags.l_lh_clustered_sampling ? "T" : "F")
         << '\n';
  output << "l_rcm_in_cloud_k_lh_start = "
         << (silhs_config_flags.l_rcm_in_cloud_k_lh_start ? "T" : "F") << '\n';
  output << "l_random_k_lh_start = " << (silhs_config_flags.l_random_k_lh_start ? "T" : "F") << '\n';
  output << "l_max_overlap_in_cloud = " << (silhs_config_flags.l_max_overlap_in_cloud ? "T" : "F")
         << '\n';
  output << "l_lh_instant_var_covar_src = "
         << (silhs_config_flags.l_lh_instant_var_covar_src ? "T" : "F") << '\n';
  output << "l_lh_limit_weights = " << (silhs_config_flags.l_lh_limit_weights ? "T" : "F") << '\n';
  output << "l_lh_var_frac = " << (silhs_config_flags.l_lh_var_frac ? "T" : "F") << '\n';
  output << "l_lh_normalize_weights = " << (silhs_config_flags.l_lh_normalize_weights ? "T" : "F")
         << '\n';
  output << "l_lh_deterministic_test = " << (silhs_config_flags.l_lh_deterministic_test ? "T" : "F")
         << '\n';
  const auto text = output.str();
  iunit.write(text.data(), static_cast<std::streamsize>(text.size()));
}

} // namespace clubb_kokkos
