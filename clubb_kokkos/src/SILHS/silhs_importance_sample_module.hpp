// SILHS cloud/precipitation/mixture categories and importance sampling.
// Public functions define/compute category probabilities, pick/scale sample
// categories and generate ordered stratified cloud-weighted/importance profiles.
#pragma once
#include "scratch_workspace.hpp"
#include <Kokkos_Core.hpp>
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "SILHS/parameters_silhs.hpp"
#include "SILHS/mt95.hpp"

namespace clubb_kokkos {

inline constexpr int num_importance_categories = 8;
struct importance_category_type {
  bool l_in_cloud;
  bool l_in_precip;
  bool l_in_component_1;
};

KOKKOS_INLINE_FUNCTION
Kokkos::Array<importance_category_type, num_importance_categories> define_importance_categories() {
  // Create the eight cloud/precipitation/mixture-component categories, retaining
  // Native statements and order. No RNG draw or mutable module state is involved.
  Kokkos::Array<importance_category_type, num_importance_categories> importance_categories;
  importance_categories[0].l_in_cloud       = true;
  importance_categories[0].l_in_precip      = true;
  importance_categories[0].l_in_component_1 = true;
  importance_categories[1].l_in_cloud       = true;
  importance_categories[1].l_in_precip      = true;
  importance_categories[1].l_in_component_1 = false;
  importance_categories[2].l_in_cloud       = false;
  importance_categories[2].l_in_precip      = true;
  importance_categories[2].l_in_component_1 = true;
  importance_categories[3].l_in_cloud       = false;
  importance_categories[3].l_in_precip      = true;
  importance_categories[3].l_in_component_1 = false;
  importance_categories[4].l_in_cloud       = true;
  importance_categories[4].l_in_precip      = false;
  importance_categories[4].l_in_component_1 = true;
  importance_categories[5].l_in_cloud       = true;
  importance_categories[5].l_in_precip      = false;
  importance_categories[5].l_in_component_1 = false;
  importance_categories[6].l_in_cloud       = false;
  importance_categories[6].l_in_precip      = false;
  importance_categories[6].l_in_component_1 = true;
  importance_categories[7].l_in_cloud       = false;
  importance_categories[7].l_in_precip      = false;
  importance_categories[7].l_in_component_1 = false;
  return importance_categories;
}

using SilhsCategories    = Kokkos::Array<importance_category_type, 8>;
using SilhsProbabilities = Kokkos::Array<real_t, 8>;
using SilhsHostReal1D    = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using SilhsHostReal2D    = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using SilhsHostInt1D     = Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
using SilhsHostInt2D     = Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace>;

// Ordered importance profiles share the case RNG with LH permutation work.
// Host profile Views are an explicit execution boundary, not Fortran forwarding.
void importance_sampling_driver(ScratchWorkspace& workspace,                // InOut
                                int num_samples,                            // In
                                real_t cloud_frac_1, real_t cloud_frac_2,   // In
                                real_t mixt_frac,                           // In
                                real_t precip_frac_1, real_t precip_frac_2, // In
                                int cluster_allocation_strategy, bool l_lh_clustered_sampling, // In
                                bool l_lh_limit_weights, bool l_lh_var_frac,
                                bool l_lh_normalize_weights, // In
                                SilhsHostReal1D X_u_chi_one_lev,
                                SilhsHostReal1D X_u_dp1_one_lev,         // InOut
                                SilhsHostReal1D X_u_dp2_one_lev,         // InOut
                                SilhsHostReal1D lh_sample_point_weights, // Out
                                mt95& rng,                               // Case state
                                const parameters_silhs_state& silhs      // Case state
);

SilhsProbabilities compute_category_real_probs(const SilhsCategories& importance_categories, // In
                                               real_t cloud_frac_1, real_t cloud_frac_2,     // In
                                               real_t mixt_frac,                             // In
                                               real_t precip_frac_1, real_t precip_frac_2    // In
);

// Returned profiles use the caller's active workspace scope; category slots are 0..7.
SilhsHostInt1D pick_sample_categories(ScratchWorkspace& workspace, int num_samples,
                                      const SilhsProbabilities& category_prescribed_probs, // In
                                      SilhsHostReal1D rand_vect                            // In
);

void cloud_weighted_sampling_driver(ScratchWorkspace& workspace, int num_samples,
                                    SilhsHostInt1D p_matrix_chi,
                                    SilhsHostInt1D p_matrix_dp1, // In
                                    real_t cloud_frac_1, real_t cloud_frac_2,
                                    real_t mixt_frac,                                 // In
                                    SilhsHostReal1D X_u_chi, SilhsHostReal1D X_u_dp1, // InOut
                                    SilhsHostReal1D lh_sample_point_weights,          // Out
                                    mt95& rng                                         // Case state
);

SilhsHostReal1D generate_strat_uniform_variate(ScratchWorkspace& workspace, // InOut
                                               int num_samples,             // In
                                               mt95& rng                    // Case state
);

SilhsHostInt1D determine_sample_categories(ScratchWorkspace& workspace, int num_samples,
                                           int pdf_dim,
                                           const hm_metadata_type& hm_metadata, // In
                                           SilhsHostReal2D X_nl_one_lev,        // In
                                           SilhsHostInt1D X_mixt_comp_one_lev,
                                           const SilhsCategories& importance_categories // In
);

} // namespace clubb_kokkos
