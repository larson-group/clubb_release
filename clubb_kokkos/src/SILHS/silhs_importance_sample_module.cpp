// SILHS importance-sampling profiles, ordered by the Native source.
// RNG-coupled profile control runs on host-owned Kokkos Views; no Native
// scientific routine is forwarded. PDF transforms/reductions use execution Views.
#include "scratch_workspace.hpp"
#include "CLUBB_core/error_code.hpp"
#include "SILHS/silhs_importance_sample_module.hpp"
#include "SILHS/generate_uniform_sample_module.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace clubb_kokkos {
namespace {
SilhsProbabilities
compute_category_sample_weights(const SilhsProbabilities& category_real_probs,      // In
                                const SilhsProbabilities& category_prescribed_probs // In
);
void limit_category_weights(const SilhsProbabilities& category_real_probs, // In
                            SilhsProbabilities& category_prescribed_probs  // InOut
);
void scale_sample_to_category(const importance_category_type& category, real_t cloud_frac_1,
                              real_t cloud_frac_2,                              // In
                              real_t mixt_frac,                                 // In
                              real_t precip_frac_1, real_t precip_frac_2,       // In
                              real_t& X_u_chi, real_t& X_u_dp1, real_t& X_u_dp2 // InOut
);
SilhsProbabilities two_cluster_cp_nocp(ScratchWorkspace& workspace,                   // InOut
                                       const SilhsCategories& importance_categories,  // In
                                       const SilhsProbabilities& category_real_probs, // In
                                       bool l_lh_var_frac,                            // In
                                       const parameters_silhs_state& silhs            // Case state
);
SilhsProbabilities eight_cluster_allocation(ScratchWorkspace& workspace,                   // InOut
                                            const SilhsCategories& importance_categories,  // In
                                            const SilhsProbabilities& category_real_probs, // In
                                            bool l_lh_var_frac,                            // In
                                            const parameters_silhs_state& silhs // Case state
);
SilhsProbabilities four_cluster_no_precip(ScratchWorkspace& workspace,                   // InOut
                                          const SilhsCategories& importance_categories,  // In
                                          const SilhsProbabilities& category_real_probs, // In
                                          bool l_lh_var_frac,                            // In
                                          const parameters_silhs_state& silhs // Case state
);
SilhsProbabilities compute_clust_category_probs(ScratchWorkspace& workspace,
                                                const SilhsProbabilities& category_real_probs,
                                                int num_clusters,                  // In
                                                int max_num_categories_in_cluster, // In
                                                SilhsHostInt1D num_categories_in_cluster,
                                                SilhsHostInt2D cluster_categories,  // In
                                                SilhsHostReal1D cluster_fractions,  // In
                                                bool l_lh_var_frac,                 // In
                                                const parameters_silhs_state& silhs // Case state
);
SilhsProbabilities clust_cat_probs_frm_var_fracs(ScratchWorkspace& workspace,
                                                 const SilhsProbabilities& category_real_probs,
                                                 int num_clusters,                  // In
                                                 int max_num_categories_in_cluster, // In
                                                 SilhsHostInt1D num_categories_in_cluster,
                                                 SilhsHostInt2D cluster_categories,         // In
                                                 SilhsHostReal1D cluster_variance_fractions // In
);
SilhsProbabilities clust_cat_probs_frm_presc_prb(ScratchWorkspace& workspace,
                                                 const SilhsProbabilities& category_real_probs,
                                                 int num_clusters,                  // In
                                                 int max_num_categories_in_cluster, // In
                                                 SilhsHostInt1D num_categories_in_cluster,
                                                 SilhsHostInt2D cluster_categories,        // In
                                                 SilhsHostReal1D cluster_prescribed_probs, // In
                                                 const parameters_silhs_state& silhs // Case state
);
SilhsProbabilities cloud_importance_sampling(const SilhsCategories& importance_categories,  // In
                                             const SilhsProbabilities& category_real_probs, // In
                                             real_t cloud_frac_1, real_t cloud_frac_2,      // In
                                             real_t mixt_frac                               // In
);
void importance_sampling_assertions(
    int num_samples, const SilhsCategories& importance_categories,                      // In
    const SilhsProbabilities& category_real_probs,                                      // In
    const SilhsProbabilities& category_prescribed_probs,                                // In
    const SilhsProbabilities& category_sample_weights, SilhsHostReal1D X_u_chi_one_lev, // In
    SilhsHostReal1D X_u_dp1_one_lev, SilhsHostReal1D X_u_dp2_one_lev,                   // In
    SilhsHostReal1D lh_sample_point_weights, SilhsHostInt1D int_sample_category,        // In
    real_t cloud_frac_1, real_t cloud_frac_2,                                           // In
    real_t mixt_frac,                                                                   // In
    real_t precip_frac_1, real_t precip_frac_2,                                         // In
    bool l_lh_normalize_weights,                                                        // In
    bool& l_error                                                                       // Out
);
void choose_X_u_scaled(bool l_cloudy_sample,                             // In
                       int p_matrix_element, int num_samples,            // In
                       real_t cloud_frac_1, real_t cloud_frac_2,         // In
                       real_t mixt_frac, real_t mixt_rand_element,       // In
                       real_t& X_u_dp1_element, real_t& X_u_chi_element, // Out
                       mt95& rng                                         // Case state
);
} // namespace

//-----------------------------------------------------------------------
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
) {
  // Description:
  //   Applies importance sampling to a single vertical level !
  // References:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:importance_sampling
  //   clubb:ticket:736 !
  // num_samples: Number of SILHS sample points
  // cluster_allocation_strategy: Strategy for distributing sample points
  // l_lh_clustered_sampling: Use prescribed probability sampling with clusters (SILHS)
  // l_lh_limit_weights: Ensure weights stay under a given value
  // l_lh_var_frac: Prescribe variance fractions
  // l_lh_normalize_weights: Normalize weights to sum to num_samples
  // X_u_chi_one_lev: These uniform variates are scaled by the importance
  // X_u_dp1_one_lev: sampling process.
  auto workspace_scope             = workspace.scope();
  bool l_error                     = false;
  const auto importance_categories = define_importance_categories();
  const auto category_real_probs   = compute_category_real_probs(importance_categories,       // In
                                                                 cloud_frac_1, cloud_frac_2,  // In
                                                                 mixt_frac,                   // In
                                                                 precip_frac_1, precip_frac_2 // In
  );
  SilhsProbabilities category_prescribed_probs;
  if (l_lh_clustered_sampling) {
    switch (cluster_allocation_strategy) {
    case eight_cluster_allocation_opt:
      category_prescribed_probs =
          eight_cluster_allocation(workspace, importance_categories, category_real_probs, // In
                                   l_lh_var_frac,                                         // In
                                   silhs // Case state
          );
      break;
    case four_cluster_allocation_opt:
      category_prescribed_probs =
          four_cluster_no_precip(workspace, importance_categories, category_real_probs, // In
                                 l_lh_var_frac,                                         // In
                                 silhs // Case state
          );
      break;
    case two_cluster_cp_nocp_opt:
      category_prescribed_probs =
          two_cluster_cp_nocp(workspace, importance_categories, category_real_probs, // In
                              l_lh_var_frac,                                         // In
                              silhs                                                  // Case state
          );
      break;
    default:
      std::cerr << "Unsupported allocation strategy:" << cluster_allocation_strategy << '\n';
      throw std::runtime_error("Fatal error in importance_sampling_driver");
    }
    if (l_lh_limit_weights)
      limit_category_weights(category_real_probs, category_prescribed_probs // In / InOut
      );
  } else {
    // Place half of all sample points in cloud.
    category_prescribed_probs =
        cloud_importance_sampling(importance_categories, category_real_probs, // In
                                  cloud_frac_1, cloud_frac_2,                 // In
                                  mixt_frac                                   // In
        );
  }
  const auto category_sample_weights =
      compute_category_sample_weights(category_real_probs, category_prescribed_probs // In
      );
  const auto rand_vect = generate_strat_uniform_variate(workspace, num_samples, // In
                                                        rng                     // Case state
  );
  const auto int_sample_category =
      pick_sample_categories(workspace, num_samples, category_prescribed_probs, // In
                             rand_vect                                          // In
      );
  for (int sample = 0; sample < num_samples; ++sample) {
    scale_sample_to_category(
        importance_categories[int_sample_category(sample)], cloud_frac_1, cloud_frac_2, // In
        mixt_frac,                                                                      // In
        precip_frac_1, precip_frac_2,                                                   // In
        X_u_chi_one_lev(sample), X_u_dp1_one_lev(sample), X_u_dp2_one_lev(sample)       // InOut
    );
    lh_sample_point_weights(sample) = category_sample_weights[int_sample_category(sample)];
  }
  if (l_lh_normalize_weights) {
    real_t weights_sum = zero;
    for (int sample = 0; sample < num_samples; ++sample)
      weights_sum = weights_sum + lh_sample_point_weights(sample);
    for (int sample = 0; sample < num_samples; ++sample)
      lh_sample_point_weights(sample) =
          lh_sample_point_weights(sample) * real_t(num_samples) / weights_sum;
  }
  if (clubb_at_least_debug_level_api(2)) {
    importance_sampling_assertions(
        num_samples, importance_categories, category_real_probs,                        // In
        category_prescribed_probs, category_sample_weights, X_u_chi_one_lev,            // In
        X_u_dp1_one_lev, X_u_dp2_one_lev, lh_sample_point_weights, int_sample_category, // In
        cloud_frac_1, cloud_frac_2,                                                     // In
        mixt_frac,                                                                      // In
        precip_frac_1, precip_frac_2,                                                   // In
        l_lh_normalize_weights,                                                         // In
        l_error                                                                         // Out
    );
    if (l_error)
      throw std::runtime_error("Fatal error in importance_sampling_driver");
  }
}

// define_importance_categories is device-inline in the public header.

//-----------------------------------------------------------------------
SilhsProbabilities compute_category_real_probs(const SilhsCategories& importance_categories, // In
                                               real_t cloud_frac_1, real_t cloud_frac_2,     // In
                                               real_t mixt_frac,                             // In
                                               real_t precip_frac_1, real_t precip_frac_2    // In
) {
  // Description:
  //   Computes the real PDF probability associated with each importance sampling
  //   category.
  //   For example, if a category is in cloud, out of precipitation, and in mixture
  //   component two, then without importance sampling, the probability that a
  //   point will appear in that category is:
  //   P(cloud,noprecip,comp2) = cloud_frac_2 * (1-precip_frac_2) * (1-mixt_frac)
  // References:
  //   None
  // importance_categories: A list of importance categories
  SilhsProbabilities category_real_probs;
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    real_t cloud_frac_i, precip_frac_i, component_factor;
    if (importance_categories[icategory].l_in_component_1) {
      cloud_frac_i     = cloud_frac_1;
      precip_frac_i    = precip_frac_1;
      component_factor = mixt_frac;
    } else {
      cloud_frac_i     = cloud_frac_2;
      precip_frac_i    = precip_frac_2;
      component_factor = one - mixt_frac;
    }
    real_t cloud_factor, precip_factor;
    if (importance_categories[icategory].l_in_cloud)
      cloud_factor = cloud_frac_i;
    else
      cloud_factor = one - cloud_frac_i;
    if (importance_categories[icategory].l_in_precip)
      precip_factor = precip_frac_i;
    else
      precip_factor = one - precip_frac_i;
    category_real_probs[icategory] = component_factor * cloud_factor * precip_factor;
  }
  return category_real_probs;
}

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities
compute_category_sample_weights(const SilhsProbabilities& category_real_probs,      // In
                                const SilhsProbabilities& category_prescribed_probs // In
) {
  // Description:
  //   Compute the sample point weights for a sample point in each category based
  //   on the PDF probability and the modified probability from importance
  //   sampling
  // References:
  //   None
  // category_real_probs: The actual PDF probability of each category
  // category_prescribed_probs: The modified probability of each category due to importance sampling
  SilhsProbabilities category_sample_weights;
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    if (std::abs(category_prescribed_probs[icategory]) < std::numeric_limits<real_t>::epsilon())
      category_sample_weights[icategory] = real_t(1);
    else
      category_sample_weights[icategory] =
          category_real_probs[icategory] / category_prescribed_probs[icategory];
  }
  return category_sample_weights;
}
} // namespace

namespace {

//-----------------------------------------------------------------------
void limit_category_weights(const SilhsProbabilities& category_real_probs, // In
                            SilhsProbabilities& category_prescribed_probs  // InOut
) {
  // Description:
  //   Modifies the category prescribed probabilities such that an
  //   invariant maximum weight is achieved. The weight for each category
  //   j is equal to category_real_probs(j) / category_prescribed_probs(j).
  //   In this subroutine, category_prescribed_probs(j) is increased in
  //   each necessary category such that no category has a weight that
  //   exceeds the maximum.
  // References:
  //   none
  // category_real_probs: The PDF probability for each importance category
  // category_prescribed_probs: Prescribed probability for each category; these will be modified to
  //   ensure the minimum weight.
  constexpr real_t max_weight                    = real_t(2);
  constexpr real_t real_prob_thresh_for_transfer = real_t(1.e-8);
  SilhsProbabilities min_presc_probs, min_presc_prob_diff;
  real_t total_diff_under = zero;
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    if (category_real_probs[icategory] >= real_prob_thresh_for_transfer)
      min_presc_probs[icategory] = category_real_probs[icategory] / max_weight;
    else
      min_presc_probs[icategory] = category_real_probs[icategory];
    min_presc_prob_diff[icategory] =
        category_prescribed_probs[icategory] - min_presc_probs[icategory];
    if (min_presc_prob_diff[icategory] < zero)
      total_diff_under = total_diff_under + std::abs(min_presc_prob_diff[icategory]);
  }
  if (total_diff_under > zero) {
    real_t total_diff_over = zero;
    for (int icategory = 0; icategory < num_importance_categories; ++icategory)
      if (min_presc_prob_diff[icategory] > zero)
        total_diff_over = total_diff_over + min_presc_prob_diff[icategory];
    if (total_diff_under > total_diff_over) {
      std::cerr << "The sample point weights could not be limited to the maximum value.\n";
      throw std::runtime_error("Fatal error in limit_category_weights");
    }
    for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
      if (min_presc_prob_diff[icategory] > zero)
        category_prescribed_probs[icategory] =
            category_prescribed_probs[icategory] -
            (min_presc_prob_diff[icategory] / total_diff_over) * total_diff_under;
      else if (min_presc_prob_diff[icategory] < zero)
        category_prescribed_probs[icategory] = min_presc_probs[icategory];
    }
  }
  if (clubb_at_least_debug_level_api(2)) {
    for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
      if (category_prescribed_probs[icategory] > zero &&
          category_real_probs[icategory] >= real_prob_thresh_for_transfer) {
        const real_t weight = category_real_probs[icategory] / category_prescribed_probs[icategory];
        if (weight > max_weight) {
          std::cerr << "In limit_category_weights, a weight was not limited.\n"
                    << "category_real_prob = " << category_real_probs[icategory] << '\n'
                    << "category_prescribed_prob = " << category_prescribed_probs[icategory] << '\n'
                    << "weight = " << weight << '\n'
                    << "min_presc_probs = " << min_presc_probs[icategory] << '\n'
                    << "min_presc_prob_diff = " << min_presc_prob_diff[icategory] << '\n';
          throw std::runtime_error("Fatal error in limit_category_weights");
        }
      }
    }
  }
}
} // namespace

//-----------------------------------------------------------------------
SilhsHostInt1D pick_sample_categories(ScratchWorkspace& workspace, int num_samples,
                                      const SilhsProbabilities& category_prescribed_probs, // In
                                      SilhsHostReal1D rand_vect                            // In
) {
  // Description:
  //   Picks a category for each sample point, based on the given probabilities,
  //   such that the distribution of categories of the sample points
  //   approximates the probabilities that are given.
  // References:
  //   None
  // num_samples: Number of sample points to be picked
  // category_prescribed_probs: Prescribed probability for each category
  // rand_vect: A sample of num_samples values from the uniform distribution in the range (0,1).
  //   This will be used to pick the categories.
  // Returned category slots are 0..7; -1 marks an unselected sample before validation.
  SilhsProbabilities category_cumulative_probs;
  category_cumulative_probs[0] = zero;
  for (int category = 1; category < num_importance_categories; ++category)
    category_cumulative_probs[category] =
        category_cumulative_probs[(category - 1)] + category_prescribed_probs[(category - 1)];
  SilhsHostInt1D int_sample_category = workspace.view<SilhsHostInt1D>(
      Kokkos::WithoutInitializing, "int_sample_category", num_samples);
  for (int sample = 0; sample < num_samples; ++sample) {
    int_sample_category(sample) = -1;
    for (int category = 0; category < num_importance_categories; ++category) {
      if (category < (num_importance_categories - 1)) {
        if (rand_vect(sample) >= category_cumulative_probs[category] &&
            rand_vect(sample) < category_cumulative_probs[(category + 1)]) {
          int_sample_category(sample) = category;
          break;
        }
      } else if (category == (num_importance_categories - 1)) {
        // Native ticket805 includes exactly one in the last category.
        if (rand_vect(sample) >= category_cumulative_probs[category] && rand_vect(sample) <= one)
          int_sample_category(sample) = category;
      }
    }
    if (int_sample_category(sample) == -1) {
      if (clubb_at_least_debug_level_api(0)) {
        std::cerr << "Invalid rand_vect number in pick_sample_categories\n"
                  << "rand_vect(sample) = " << rand_vect(sample) << '\n';
        throw std::runtime_error("Fatal error");
      } else
        int_sample_category(sample) = 0;
    }
  }
  return int_sample_category;
}

namespace {

//-----------------------------------------------------------------------
void scale_sample_to_category(const importance_category_type& category, real_t cloud_frac_1,
                              real_t cloud_frac_2,                              // In
                              real_t mixt_frac,                                 // In
                              real_t precip_frac_1, real_t precip_frac_2,       // In
                              real_t& X_u_chi, real_t& X_u_dp1, real_t& X_u_dp2 // InOut
) {
  // Description:
  //   Scale and transpose a sample point to reside in the specified category
  // References:
  //   None
  // category: Scale the sample point to reside in this category
  // X_u_chi: Uniform sample of extend cloud water mixing ratio
  // X_u_dp1: Uniform sample of the d+1 variate
  // X_u_dp2: Uniform sample of the d+2 variate
  real_t cloud_frac_i, precip_frac_i;
  // Scale dp1 into component1 (0,mixt_frac), or component2 (mixt_frac,1).
  if (category.l_in_component_1) {
    X_u_dp1       = X_u_dp1 * mixt_frac;
    cloud_frac_i  = cloud_frac_1;
    precip_frac_i = precip_frac_1;
  } else {
    X_u_dp1       = X_u_dp1 * (one - mixt_frac) + mixt_frac;
    cloud_frac_i  = cloud_frac_2;
    precip_frac_i = precip_frac_2;
  }
  // Scale dp2 into or out of precipitation.
  if (category.l_in_precip)
    X_u_dp2 = X_u_dp2 * precip_frac_i;
  else
    X_u_dp2 = X_u_dp2 * (one - precip_frac_i) + precip_frac_i;
  // Scale chi into cloud (1-cloud_frac_i,1) or clear air (0,1-cloud_frac_i).
  if (category.l_in_cloud)
    X_u_chi = X_u_chi * cloud_frac_i + (one - cloud_frac_i);
  else
    X_u_chi = X_u_chi * (one - cloud_frac_i);
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities two_cluster_cp_nocp(ScratchWorkspace& workspace,                   // InOut
                                       const SilhsCategories& importance_categories,  // In
                                       const SilhsProbabilities& category_real_probs, // In
                                       bool l_lh_var_frac,                            // In
                                       const parameters_silhs_state& silhs            // Case state
) {
  // Description:
  //   Clusters importance categories into two clusters: categories that contain
  //   either cloud or precipitation (or both), and clusters that contain
  //   neither.
  // References:
  //   clubb:ticket:740
  // importance_categories: A list of importance categories
  // category_real_probs: The real probability for each category
  // l_lh_var_frac: Prescribe variance fractions
  auto workspace_scope       = workspace.scope();
  constexpr int num_clusters = 2, max_num_categories_in_cluster = 6;
  constexpr int iicld_or_precip = 0, iinocld_precip = 1;
  constexpr real_t prob_cld_or_precip = real_t(1), prob_nocld_precip = real_t(0);
  SilhsHostInt1D num_categories_in_cluster =
      workspace_scope.view<SilhsHostInt1D>("num_categories_in_cluster", num_clusters);
  SilhsHostInt2D cluster_categories =
      workspace_scope.view<SilhsHostInt2D>(Kokkos::WithoutInitializing, "cluster_categories",
                                           num_clusters, max_num_categories_in_cluster);
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    if (importance_categories[icategory].l_in_cloud ||
        importance_categories[icategory].l_in_precip) {
      ++num_categories_in_cluster(iicld_or_precip);
      cluster_categories(iicld_or_precip, num_categories_in_cluster(iicld_or_precip) - 1) =
          icategory;
    } else {
      ++num_categories_in_cluster(iinocld_precip);
      cluster_categories(iinocld_precip, num_categories_in_cluster(iinocld_precip) - 1) = icategory;
    }
  }
  if (clubb_at_least_debug_level_api(2))
    if (num_categories_in_cluster(iinocld_precip) != 2 ||
        num_categories_in_cluster(iicld_or_precip) != 6)
      throw std::runtime_error("Invalid categories in two_cluster_cp_nocp");
  SilhsHostReal1D cluster_prescribed_probs = workspace_scope.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "cluster_prescribed_probs", num_clusters);
  cluster_prescribed_probs(iicld_or_precip) = prob_cld_or_precip;
  cluster_prescribed_probs(iinocld_precip)  = prob_nocld_precip;
  return compute_clust_category_probs(
      workspace, category_real_probs, num_clusters, max_num_categories_in_cluster, // In
      num_categories_in_cluster, cluster_categories, cluster_prescribed_probs,     // In
      l_lh_var_frac,                                                               // In
      silhs                                                                        // Case state
  );
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities eight_cluster_allocation(ScratchWorkspace& workspace,                   // InOut
                                            const SilhsCategories& importance_categories,  // In
                                            const SilhsProbabilities& category_real_probs, // In
                                            bool l_lh_var_frac,                            // In
                                            const parameters_silhs_state& silhs // Case state
) {
  // Description:
  //   Clusters importance categories such that each of the eight importance
  //   categories has its own cluster. Effectively, there are no clusters.
  // References:
  //   clubb:ticket:752
  // importance_categories: A list of importance categories
  // category_real_probs: The real probability for each category
  // l_lh_var_frac: Prescribe variance fractions
  auto workspace_scope       = workspace.scope();
  constexpr int num_clusters = 8, max_num_categories_in_cluster = 1;
  SilhsHostInt1D num_categories_in_cluster = workspace_scope.view<SilhsHostInt1D>(
      Kokkos::WithoutInitializing, "num_categories_in_cluster", num_clusters);
  Kokkos::deep_copy(num_categories_in_cluster, 1);
  SilhsHostInt2D cluster_categories =
      workspace_scope.view<SilhsHostInt2D>(Kokkos::WithoutInitializing, "cluster_categories",
                                           num_clusters, max_num_categories_in_cluster);
  SilhsHostReal1D cluster_prescribed_probs = workspace_scope.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "cluster_prescribed_probs", num_clusters);
  const auto& eight_cluster_presc_probs = silhs.eight_cluster_presc_probs;
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    cluster_categories(icategory, 0) = icategory;
    const bool l_in_cloud            = importance_categories[icategory].l_in_cloud;
    const bool l_in_component_1      = importance_categories[icategory].l_in_component_1;
    const bool l_in_precip           = importance_categories[icategory].l_in_precip;
    if (l_in_cloud && l_in_precip && l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.cloud_precip_comp1;
    else if (l_in_cloud && l_in_precip && !l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.cloud_precip_comp2;
    else if (!l_in_cloud && l_in_precip && l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.nocloud_precip_comp1;
    else if (!l_in_cloud && l_in_precip && !l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.nocloud_precip_comp2;
    else if (l_in_cloud && !l_in_precip && l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.cloud_noprecip_comp1;
    else if (l_in_cloud && !l_in_precip && !l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.cloud_noprecip_comp2;
    else if (!l_in_cloud && !l_in_precip && l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.nocloud_noprecip_comp1;
    else if (!l_in_cloud && !l_in_precip && !l_in_component_1)
      cluster_prescribed_probs(icategory) = eight_cluster_presc_probs.nocloud_noprecip_comp2;
    else
      throw std::runtime_error("Invalid category in eight_cluster_allocation");
  }
  return compute_clust_category_probs(
      workspace, category_real_probs, num_clusters, max_num_categories_in_cluster, // In
      num_categories_in_cluster, cluster_categories, cluster_prescribed_probs,     // In
      l_lh_var_frac,                                                               // In
      silhs                                                                        // Case state
  );
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities four_cluster_no_precip(ScratchWorkspace& workspace,                   // InOut
                                          const SilhsCategories& importance_categories,  // In
                                          const SilhsProbabilities& category_real_probs, // In
                                          bool l_lh_var_frac,                            // In
                                          const parameters_silhs_state& silhs // Case state
) {
  // Description:
  //   Clusters categories into four clusters for the four combinations of
  //   cloud/no cloud and comp 1/comp 2. Precip fraction is effectively
  //   ignored
  // References:
  //   clubb:ticket:752
  // importance_categories: A list of importance categories
  // category_real_probs: The real probability for each category
  // l_lh_var_frac: Prescribe variance fractions
  auto workspace_scope       = workspace.scope();
  constexpr int num_clusters = 4, max_num_categories_in_cluster = 2;
  constexpr int iicld_comp1 = 0, iicld_comp2 = 1, iincld_comp1 = 2, iincld_comp2 = 3;
  constexpr real_t cloud_comp1 = real_t(.30), cloud_comp2 = real_t(.30),
                   nocloud_comp1 = real_t(.20), nocloud_comp2 = real_t(.20);
  SilhsHostInt1D num_categories_in_cluster =
      workspace_scope.view<SilhsHostInt1D>("num_categories_in_cluster", num_clusters);
  SilhsHostInt2D cluster_categories =
      workspace_scope.view<SilhsHostInt2D>(Kokkos::WithoutInitializing, "cluster_categories",
                                           num_clusters, max_num_categories_in_cluster);
  SilhsHostReal1D cluster_prescribed_probs = workspace_scope.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "cluster_prescribed_probs", num_clusters);
  for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
    const bool l_in_cloud       = importance_categories[icategory].l_in_cloud;
    const bool l_in_component_1 = importance_categories[icategory].l_in_component_1;
    if (l_in_cloud && l_in_component_1) {
      cluster_prescribed_probs(iicld_comp1) = cloud_comp1;
      ++num_categories_in_cluster(iicld_comp1);
      cluster_categories(iicld_comp1, num_categories_in_cluster(iicld_comp1) - 1) = icategory;
    } else if (l_in_cloud && !l_in_component_1) {
      cluster_prescribed_probs(iicld_comp2) = cloud_comp2;
      ++num_categories_in_cluster(iicld_comp2);
      cluster_categories(iicld_comp2, num_categories_in_cluster(iicld_comp2) - 1) = icategory;
    } else if (!l_in_cloud && l_in_component_1) {
      cluster_prescribed_probs(iincld_comp1) = nocloud_comp1;
      ++num_categories_in_cluster(iincld_comp1);
      cluster_categories(iincld_comp1, num_categories_in_cluster(iincld_comp1) - 1) = icategory;
    } else if (!l_in_cloud && !l_in_component_1) {
      cluster_prescribed_probs(iincld_comp2) = nocloud_comp2;
      ++num_categories_in_cluster(iincld_comp2);
      cluster_categories(iincld_comp2, num_categories_in_cluster(iincld_comp2) - 1) = icategory;
    } else
      throw std::runtime_error("Invalid category in four_cluster_no_precip");
  }
  if (clubb_at_least_debug_level_api(2)) {
    bool invalid = false;
    for (int i = 0; i < num_clusters; ++i)
      invalid = invalid || num_categories_in_cluster(i) != 2;
    if (invalid) {
      std::cerr << "Not all clusters have two categories\n";
      throw std::runtime_error("Fatal error in four_cluster_no_precip");
    }
  }
  return compute_clust_category_probs(
      workspace, category_real_probs, num_clusters, max_num_categories_in_cluster, // In
      num_categories_in_cluster, cluster_categories, cluster_prescribed_probs,     // In
      l_lh_var_frac,                                                               // In
      silhs                                                                        // Case state
  );
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities compute_clust_category_probs(ScratchWorkspace& workspace,
                                                const SilhsProbabilities& category_real_probs,
                                                int num_clusters,                  // In
                                                int max_num_categories_in_cluster, // In
                                                SilhsHostInt1D num_categories_in_cluster,
                                                SilhsHostInt2D cluster_categories,  // In
                                                SilhsHostReal1D cluster_fractions,  // In
                                                bool l_lh_var_frac,                 // In
                                                const parameters_silhs_state& silhs // Case state
) {
  // Description:
  //   Calls clust_cat_probs_frm_presc_prb or clust_cat_probs_frm_var_fracs!
  // References:
  //   clubb:ticket:740
  // category_real_probs: The real probability for each category
  // num_clusters: The number of clusters to sample from
  // max_num_categories_in_cluster: The max number of categories in each cluster
  // num_categories_in_cluster: The number of categories in each cluster
  // cluster_categories: An integer matrix containing indices corresponding to the members of the
  //   clusters
  // cluster_fractions: Prescribed fraction of some sort for each cluster
  // l_lh_var_frac: Prescribe variance fractions
  if (l_lh_var_frac)
    return clust_cat_probs_frm_var_fracs(
        workspace, category_real_probs, num_clusters, max_num_categories_in_cluster, // In
        num_categories_in_cluster, cluster_categories, cluster_fractions             // In
    );
  return clust_cat_probs_frm_presc_prb(
      workspace, category_real_probs, num_clusters, max_num_categories_in_cluster, // In
      num_categories_in_cluster, cluster_categories, cluster_fractions,            // In
      silhs                                                                        // Case state
  );
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities clust_cat_probs_frm_var_fracs(ScratchWorkspace& workspace,
                                                 const SilhsProbabilities& category_real_probs,
                                                 int num_clusters,                  // In
                                                 int max_num_categories_in_cluster, // In
                                                 SilhsHostInt1D num_categories_in_cluster,
                                                 SilhsHostInt2D cluster_categories,         // In
                                                 SilhsHostReal1D cluster_variance_fractions // In
) {
  // Description:
  //   This is a generalized algorithm that takes as input a set of "clusters"
  //   of the importance categories and a variance fraction      for each
  //   cluster, and computes the prescribed probabilities for each category.
  // References:
  //   clubb:ticket:740
  // category_real_probs: The real probability for each category
  // num_clusters: The number of clusters to sample from
  // max_num_categories_in_cluster: The max number of categories in each cluster
  // num_categories_in_cluster: The number of categories in each cluster
  // cluster_categories: An integer matrix containing indices corresponding to the members of the
  //   clusters
  // cluster_variance_fractions: Prescribed variance fraction for each cluster
  auto workspace_scope = workspace.scope();
  SilhsHostReal1D cluster_real_probs =
      workspace_scope.view<SilhsHostReal1D>("cluster_real_probs", num_clusters);
  SilhsHostReal1D cluster_prescribed_probs = workspace_scope.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "cluster_prescribed_probs", num_clusters);
  for (int icluster = 0; icluster < num_clusters; ++icluster)
    for (int icategory = 0; icategory < num_categories_in_cluster(icluster); ++icategory)
      cluster_real_probs(icluster) = cluster_real_probs(icluster) +
                                     category_real_probs[cluster_categories(icluster, icategory)];
  real_t pdf_prob_var_frac_prod_sum = zero;
  for (int icluster = 0; icluster < num_clusters; ++icluster)
    pdf_prob_var_frac_prod_sum =
        pdf_prob_var_frac_prod_sum +
        cluster_real_probs(icluster) * cluster_variance_fractions(icluster);
  if (std::abs(pdf_prob_var_frac_prod_sum) < eps) {
    for (int icluster = 0; icluster < num_clusters; ++icluster)
      cluster_prescribed_probs(icluster) = cluster_real_probs(icluster);
  } else {
    for (int icluster = 0; icluster < num_clusters; ++icluster)
      cluster_prescribed_probs(icluster) = cluster_real_probs(icluster) *
                                           cluster_variance_fractions(icluster) /
                                           pdf_prob_var_frac_prod_sum;
  }
  SilhsProbabilities category_prescribed_probs;
  for (int icluster = 0; icluster < num_clusters; ++icluster) {
    for (int icategory = 0; icategory < num_categories_in_cluster(icluster); ++icategory) {
      const int cat_idx = cluster_categories(icluster, icategory);
      if (std::abs(cluster_real_probs(icluster)) < eps)
        category_prescribed_probs[cat_idx] = zero;
      else
        category_prescribed_probs[cat_idx] =
            (category_real_probs[cat_idx] / cluster_real_probs(icluster)) *
            cluster_prescribed_probs(icluster);
    }
  }
  return category_prescribed_probs;
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities clust_cat_probs_frm_presc_prb(ScratchWorkspace& workspace,
                                                 const SilhsProbabilities& category_real_probs,
                                                 int num_clusters,                  // In
                                                 int max_num_categories_in_cluster, // In
                                                 SilhsHostInt1D num_categories_in_cluster,
                                                 SilhsHostInt2D cluster_categories,        // In
                                                 SilhsHostReal1D cluster_prescribed_probs, // In
                                                 const parameters_silhs_state& silhs // Case state
) {
  // Description:
  //   This is a generalized algorithm that takes as input a set of "clusters"
  //   of the importance categories and a prescribed probability for each
  //   cluster, and computes the prescribed probabilities for each category
  //   such that the sum of the prescribed probabilities of every category
  //   within a cluster is equal to the prescribed probability of the cluster.
  // References:
  //   clubb:ticket:752
  // category_real_probs: The real probability for each category
  // num_clusters: The number of clusters to sample from
  // max_num_categories_in_cluster: The max number of categories in each cluster
  // num_categories_in_cluster: The number of categories in each cluster
  // cluster_categories: An integer matrix containing indices corresponding to the members of the
  //   clusters
  // cluster_prescribed_probs: Prescribed probability sum for each cluster
  auto workspace_scope = workspace.scope();
  SilhsHostReal1D cluster_real_probs =
      workspace_scope.view<SilhsHostReal1D>("cluster_real_probs", num_clusters);
  SilhsHostReal1D cluster_prescribed_probs_mod = workspace_scope.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "cluster_prescribed_probs_mod", num_clusters);
  auto l_cluster_presc_prob_modified =
      workspace_scope.view<Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
          Kokkos::WithoutInitializing, "l_cluster_presc_prob_modified", num_clusters);
  for (int icluster = 0; icluster < num_clusters; ++icluster)
    for (int icategory = 0; icategory < num_categories_in_cluster(icluster); ++icategory)
      cluster_real_probs(icluster) = cluster_real_probs(icluster) +
                                     category_real_probs[cluster_categories(icluster, icategory)];
  bool any_modified = false;
  for (int icluster = 0; icluster < num_clusters; ++icluster) {
    if (cluster_real_probs(icluster) < silhs.importance_prob_thresh) {
      cluster_prescribed_probs_mod(icluster)  = cluster_real_probs(icluster);
      l_cluster_presc_prob_modified(icluster) = true;
    } else {
      cluster_prescribed_probs_mod(icluster)  = cluster_prescribed_probs(icluster);
      l_cluster_presc_prob_modified(icluster) = false;
    }
    any_modified = any_modified || l_cluster_presc_prob_modified(icluster);
  }
  if (any_modified) {
    real_t nonzero_real_clust_sum = zero;
    for (int icluster = 0; icluster < num_clusters; ++icluster)
      if (!l_cluster_presc_prob_modified(icluster))
        nonzero_real_clust_sum = nonzero_real_clust_sum + cluster_real_probs(icluster);
    // Transfer each thresholded cluster in Native order. Summing all differences
    // first, as JAX does, would change the sequence of rounding operations.
    for (int icluster = 0; icluster < num_clusters; ++icluster) {
      if (l_cluster_presc_prob_modified(icluster)) {
        const real_t presc_prob_difference =
            cluster_prescribed_probs(icluster) - cluster_prescribed_probs_mod(icluster);
        for (int jcluster = 0; jcluster < num_clusters; ++jcluster)
          if (!l_cluster_presc_prob_modified(jcluster))
            cluster_prescribed_probs_mod(jcluster) =
                cluster_prescribed_probs_mod(jcluster) +
                (presc_prob_difference * cluster_real_probs(jcluster) / nonzero_real_clust_sum);
      }
    }
  }
  SilhsProbabilities category_prescribed_probs;
  for (int icluster = 0; icluster < num_clusters; ++icluster) {
    for (int icategory = 0; icategory < num_categories_in_cluster(icluster); ++icategory) {
      const int cat_idx = cluster_categories(icluster, icategory);
      if (l_cluster_presc_prob_modified(icluster))
        category_prescribed_probs[cat_idx] = category_real_probs[cat_idx];
      else
        category_prescribed_probs[cat_idx] =
            (category_real_probs[cat_idx] / cluster_real_probs(icluster)) *
            cluster_prescribed_probs_mod(icluster);
    }
  }
  return category_prescribed_probs;
}
} // namespace

namespace {

//-----------------------------------------------------------------------
SilhsProbabilities cloud_importance_sampling(const SilhsCategories& importance_categories,  // In
                                             const SilhsProbabilities& category_real_probs, // In
                                             real_t cloud_frac_1, real_t cloud_frac_2,      // In
                                             real_t mixt_frac                               // In
) {
  // Description:
  //   Applies cloud weighted sampling such that approximately half of all
  //   sample points land in cloud and half land out of cloud!
  // References:
  //   None :(
  // importance_categories: A list of importance categories
  // category_real_probs: The actual PDF probability for each category
  constexpr real_t cloud_frac_min_samp = real_t(.001), cloud_frac_max_samp = real_t(.5);
  const real_t cloud_frac = compute_mean_binormal(cloud_frac_1, cloud_frac_2, mixt_frac);
  SilhsProbabilities category_prescribed_probs;
  if (cloud_frac >= cloud_frac_min_samp && cloud_frac < cloud_frac_max_samp) {
    for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
      if (importance_categories[icategory].l_in_cloud)
        category_prescribed_probs[icategory] = category_real_probs[icategory] / (two * cloud_frac);
      else
        category_prescribed_probs[icategory] =
            category_real_probs[icategory] / (two * (one - cloud_frac));
    }
  } else
    category_prescribed_probs = category_real_probs;
  return category_prescribed_probs;
}
} // namespace

namespace {

//-----------------------------------------------------------------------
void importance_sampling_assertions(
    int num_samples, const SilhsCategories& importance_categories,                      // In
    const SilhsProbabilities& category_real_probs,                                      // In
    const SilhsProbabilities& category_prescribed_probs,                                // In
    const SilhsProbabilities& category_sample_weights, SilhsHostReal1D X_u_chi_one_lev, // In
    SilhsHostReal1D X_u_dp1_one_lev, SilhsHostReal1D X_u_dp2_one_lev,                   // In
    SilhsHostReal1D lh_sample_point_weights, SilhsHostInt1D int_sample_category,        // In
    real_t cloud_frac_1, real_t cloud_frac_2,                                           // In
    real_t mixt_frac,                                                                   // In
    real_t precip_frac_1, real_t precip_frac_2,                                         // In
    bool l_lh_normalize_weights,                                                        // In
    bool& l_error                                                                       // Out
) {
  // Description:
  //   Various assertion checks for importance sampling are performed here.
  // References:
  //   None
  // num_samples: Number of SILHS sample points
  // importance_categories: The defined importance categories
  // category_real_probs: The real PDF probabilities for each category
  // category_prescribed_probs: Prescribed probability for each category
  // category_sample_weights: Sample weight for each category
  // X_u_chi_one_lev: Samples of chi in uniform space
  // X_u_dp1_one_lev: Samples of the dp1 variate
  // X_u_dp2_one_lev: Samples of the dp2 variate
  // lh_sample_point_weights: Weights of samples
  // int_sample_category: An integer for each sample corresponding to the category picked for the
  //   sample
  // l_lh_normalize_weights: Normalize weights to sum to num_samples
  // l_error: True if the assertion check fails.
  l_error = false;
  const real_t tolerance =
      real_t(num_importance_categories) * std::numeric_limits<real_t>::epsilon();
  real_t category_sum = zero;
  for (int j = 0; j < 8; ++j)
    category_sum = category_sum + category_real_probs[j];
  if (std::abs(category_sum - one) > tolerance) {
    std::cerr << "The real category probabilities do not sum to one.\n"
              << "sum( category_real_probs ) = " << category_sum << '\n';
    l_error = true;
  }
  category_sum = zero;
  for (int j = 0; j < 8; ++j)
    category_sum = category_sum + category_prescribed_probs[j];
  if (std::abs(category_sum - one) > tolerance) {
    std::cerr << "The prescribed category probabilities do not sum to one.\n"
              << "sum( category_prescribed_probs ) = " << category_sum
              << "\ncategory_real_probs = ";
    for (int j = 0; j < 8; ++j)
      std::cerr << category_real_probs[j] << ' ';
    std::cerr << "\ncategory_prescribed_probs = ";
    for (int j = 0; j < 8; ++j)
      std::cerr << category_prescribed_probs[j] << ' ';
    std::cerr << '\n';
    l_error = true;
  }
  category_sum = zero;
  for (int j = 0; j < 8; ++j)
    category_sum = category_sum + category_sample_weights[j] * category_prescribed_probs[j];
  if (std::abs(category_sum - one) > tolerance) {
    std::cerr << "The weighted prescribed category probabilities do not sum to one.\n"
              << "sum( category_sample_weights * category_prescribed_probs ) = " << category_sum
              << '\n';
    l_error = true;
  }
  if (l_lh_normalize_weights) {
    const real_t num_samples_real = real_t(num_samples);
    real_t weights_sum            = zero;
    for (int j = 0; j < num_samples; ++j)
      weights_sum = weights_sum + lh_sample_point_weights(j);
    const real_t weights_avg = weights_sum / num_samples_real;
    if (std::abs(weights_avg - one) > num_samples_real * std::numeric_limits<real_t>::epsilon()) {
      std::cerr << "The sample point weights do not average to one.\n"
                << "num_samples = " << num_samples
                << "\nsum( lh_sample_point_weights(:) ) = " << weights_sum
                << "\navg( lh_sample_point_weights(:) ) = " << weights_avg << '\n';
      l_error = true;
    }
  }
  for (int isample = 0; isample < num_samples; ++isample) {
    const auto category = importance_categories[int_sample_category(isample)];
    real_t cloud_frac_i, precip_frac_i;
    if (category.l_in_component_1) {
      if (X_u_dp1_one_lev(isample) > mixt_frac) {
        std::cerr << "The component of a sample is incorrect.\n";
        l_error = true;
      }
      cloud_frac_i  = cloud_frac_1;
      precip_frac_i = precip_frac_1;
    } else {
      if (X_u_dp1_one_lev(isample) < mixt_frac) {
        std::cerr << "The component of a sample is incorrect.\n";
        l_error = true;
      }
      cloud_frac_i  = cloud_frac_2;
      precip_frac_i = precip_frac_2;
    }
    if (category.l_in_cloud) {
      if (X_u_chi_one_lev(isample) < one - cloud_frac_i) {
        std::cerr << "The chi element of a sample is incorrect.\n";
        l_error = true;
      }
    } else {
      if (X_u_chi_one_lev(isample) > one - cloud_frac_i) {
        std::cerr << "The chi element of a sample is incorrect.\n";
        l_error = true;
      }
    }
    if (category.l_in_precip) {
      if (X_u_dp2_one_lev(isample) > precip_frac_i) {
        std::cerr << "The in-precipitation status of a sample is incorrect.\n";
        l_error = true;
      }
    } else {
      if (X_u_dp2_one_lev(isample) < precip_frac_i) {
        std::cerr << "The in-precipitation status of a sample is incorrect.\n";
        l_error = true;
      }
    }
  }
}
} // namespace

//-----------------------------------------------------------------------
void cloud_weighted_sampling_driver(ScratchWorkspace& workspace, int num_samples,
                                    SilhsHostInt1D p_matrix_chi,
                                    SilhsHostInt1D p_matrix_dp1, // In
                                    real_t cloud_frac_1, real_t cloud_frac_2,
                                    real_t mixt_frac,                                 // In
                                    SilhsHostReal1D X_u_chi, SilhsHostReal1D X_u_dp1, // InOut
                                    SilhsHostReal1D lh_sample_point_weights,          // Out
                                    mt95& rng                                         // Case state
) {
  // Description:
  //   Performs importance sampling such that half of sample points are in cloud
  // References:
  //   None
  // num_samples: Number of SILHS sample points
  // p_matrix_chi: Permutation of the integers 0..num_samples-1 from p_matrix for chi
  // p_matrix_dp1: Elements from p_matrix for dp1 element
  // cloud_frac_1: Cloud fraction in PDF component 1
  // cloud_frac_2: Cloud fraction in PDF component 2
  // mixt_frac: Weight of first gaussian component
  // X_u_chi: Samples of chi in uniform space
  // X_u_dp1: Samples of the dp1 variate for determining mixture component
  // lh_sample_point_weights: Weight of SILHS sample points (these must be applied when averaging
  //   results from, e.g., calling microphysics
  auto workspace_scope                 = workspace.scope();
  constexpr real_t cloud_frac_min_samp = real_t(.001), cloud_frac_max_samp = real_t(.5);
  const real_t cloud_frac = compute_mean_binormal(cloud_frac_1, cloud_frac_2, mixt_frac);
  if (cloud_frac >= cloud_frac_min_samp && cloud_frac < cloud_frac_max_samp) {
    int n_cloudy_samples = 0, n_clear_samples = 0;
    SilhsHostInt2D mixt_permuted_cloud = workspace_scope.view<SilhsHostInt2D>(
        Kokkos::WithoutInitializing, "mixt_permuted_cloud", (num_samples / 2), 1);
    SilhsHostInt2D mixt_permuted_clear = workspace_scope.view<SilhsHostInt2D>(
        Kokkos::WithoutInitializing, "mixt_permuted_clear", (num_samples / 2), 1);
    SilhsHostReal2D mixt_rand_cloud = workspace_scope.view<SilhsHostReal2D>(
        Kokkos::WithoutInitializing, "mixt_rand_cloud", (num_samples / 2), 1);
    SilhsHostReal2D mixt_rand_clear = workspace_scope.view<SilhsHostReal2D>(
        Kokkos::WithoutInitializing, "mixt_rand_clear", (num_samples / 2), 1);
    for (int sample = 0; sample < num_samples; ++sample) {
      if (p_matrix_dp1(sample) < (num_samples / 2)) {
        ++n_cloudy_samples;
        mixt_permuted_cloud((n_cloudy_samples - 1), 0) = p_matrix_dp1(sample);
      } else if (p_matrix_dp1(sample) >= (num_samples / 2)) {
        ++n_clear_samples;
        mixt_permuted_clear((n_clear_samples - 1), 0) = p_matrix_dp1(sample) - (num_samples / 2);
      }
    }
    for (int sample = 0; sample < (num_samples / 2); ++sample) {
      mixt_rand_cloud(sample, 0) =
          choose_permuted_random((num_samples / 2), mixt_permuted_cloud(sample, 0), // In
                                 rng                                                // Case state
          );
      mixt_rand_clear(sample, 0) =
          choose_permuted_random((num_samples / 2), mixt_permuted_clear(sample, 0), // In
                                 rng                                                // Case state
          );
    }
    n_cloudy_samples                    = 0;
    n_clear_samples                     = 0;
    const real_t lh_sample_cloud_weight = two * cloud_frac;
    const real_t lh_sample_clear_weight = two - lh_sample_cloud_weight;
    for (int sample = 0; sample < num_samples; ++sample) {
      bool l_cloudy_sample;
      real_t mixt_rand_element;
      if (p_matrix_chi(sample) < (num_samples / 2)) {
        l_cloudy_sample                 = false;
        lh_sample_point_weights(sample) = lh_sample_clear_weight;
        ++n_clear_samples;
        mixt_rand_element = mixt_rand_clear((n_clear_samples - 1), 0);
      } else {
        l_cloudy_sample                 = true;
        lh_sample_point_weights(sample) = lh_sample_cloud_weight;
        ++n_cloudy_samples;
        mixt_rand_element = mixt_rand_cloud((n_cloudy_samples - 1), 0);
      }
      choose_X_u_scaled(l_cloudy_sample,                   // In
                        p_matrix_chi(sample), num_samples, // In
                        cloud_frac_1, cloud_frac_2,        // In
                        mixt_frac, mixt_rand_element,      // In
                        X_u_dp1(sample), X_u_chi(sample),  // Out
                        rng                                // Case state
      );
    }
  } else
    Kokkos::deep_copy(lh_sample_point_weights, one);
}

//-----------------------------------------------------------------------
SilhsHostReal1D generate_strat_uniform_variate(ScratchWorkspace& workspace, // InOut
                                               int num_samples,             // In
                                               mt95& rng                    // Case state
) {
  // Description:
  //   Generates a stratified uniform sample for a single variable
  //   The returned profile remains reserved in the caller's active workspace scope.
  // References:
  //   None
  // num_samples: Number of SILHS sample points
  SilhsHostInt1D pvect =
      workspace.view<SilhsHostInt1D>(Kokkos::WithoutInitializing, "pvect", num_samples);
  rand_permute(num_samples, pvect, // In / Out
               rng                 // Case state
  );
  SilhsHostReal1D stratified_variate = workspace.view<SilhsHostReal1D>(
      Kokkos::WithoutInitializing, "stratified_variate", num_samples);
  for (int sample = 0; sample < num_samples; ++sample)
    stratified_variate(sample) = real_t(choose_permuted_random(num_samples, pvect(sample), // In
                                                               rng // Case state
                                                               ));
  return stratified_variate;
}

namespace {

//-----------------------------------------------------------------------
void choose_X_u_scaled(bool l_cloudy_sample,                             // In
                       int p_matrix_element, int num_samples,            // In
                       real_t cloud_frac_1, real_t cloud_frac_2,         // In
                       real_t mixt_frac, real_t mixt_rand_element,       // In
                       real_t& X_u_dp1_element, real_t& X_u_chi_element, // Out
                       mt95& rng                                         // Case state
) {
  // Description:
  //   Find a clear or cloudy point for sampling.
  // References:
  //   None
  // l_cloudy_sample: Whether his is a cloudy or clear air sample point
  // p_matrix_element: Integer from 0..num_samples-1 for this sample
  // num_samples: Total number of calls to the microphysics
  // cloud_frac_1: Cloud fraction associated with mixture component 1     [-]
  // cloud_frac_2: Cloud fraction associated with mixture component 2     [-]
  // mixt_frac: Mixture fraction                                       [-]
  // mixt_rand_element: Random number (0,1) for determining mixture component
  // X_u_dp1_element: Elements from X_u (uniform dist.)
  // X_u_chi_element: Elements from X_u (uniform dist.)
  const real_t cld_comp1_frac = mixt_frac * cloud_frac_1;
  const real_t cld_comp2_frac = (one - mixt_frac) * cloud_frac_2;
  real_t conditional_mixt_frac;
  if (l_cloudy_sample)
    conditional_mixt_frac = cld_comp1_frac / (cld_comp1_frac + cld_comp2_frac);
  else {
    const real_t nocld_comp1_frac = mixt_frac - cld_comp1_frac;
    const real_t nocld_comp2_frac = (one - mixt_frac) - cld_comp2_frac;
    conditional_mixt_frac         = nocld_comp1_frac / (nocld_comp1_frac + nocld_comp2_frac);
  }
  int X_mixt_comp_one_lev;
  if (mixt_rand_element < conditional_mixt_frac)
    X_mixt_comp_one_lev = 1;
  else
    X_mixt_comp_one_lev = 2;
  real_t mixt_rand_element_scaled;
  if (X_mixt_comp_one_lev == 1) {
    mixt_rand_element_scaled = mixt_rand_element / conditional_mixt_frac;
    X_u_dp1_element          = mixt_rand_element_scaled * mixt_frac;
  } else if (X_mixt_comp_one_lev == 2) {
    mixt_rand_element_scaled =
        (mixt_rand_element - conditional_mixt_frac) / (one - conditional_mixt_frac);
    X_u_dp1_element = mixt_rand_element_scaled * (one - mixt_frac) + mixt_frac;
  } else
    throw std::runtime_error("Should not be here");
  int p_matrix_element_ranged;
  if (p_matrix_element >= (num_samples / 2))
    p_matrix_element_ranged = (-(num_samples / 2) + p_matrix_element);
  else
    p_matrix_element_ranged = p_matrix_element;
  const real_t chi_rand_element =
      choose_permuted_random((num_samples / 2), p_matrix_element_ranged, // In
                             rng                                         // Case state
      );
  real_t cloud_frac_i;
  if (X_mixt_comp_one_lev == 1)
    cloud_frac_i = cloud_frac_1;
  else
    cloud_frac_i = cloud_frac_2;
  if (l_cloudy_sample)
    X_u_chi_element = cloud_frac_i * chi_rand_element + (one - cloud_frac_i);
  else
    X_u_chi_element = (one - cloud_frac_i) * chi_rand_element;
}
} // namespace

//-----------------------------------------------------------------------
SilhsHostInt1D determine_sample_categories(ScratchWorkspace& workspace, int num_samples,
                                           int pdf_dim,
                                           const hm_metadata_type& hm_metadata, // In
                                           SilhsHostReal2D X_nl_one_lev,        // In
                                           SilhsHostInt1D X_mixt_comp_one_lev,
                                           const SilhsCategories& importance_categories // In
) {
  // Description:
  //   Determines the importance category of each sample.
  //   Returns slots 0..7 in the caller's active workspace scope.
  // References:
  //   None
  // num_samples: Number of SILHS sample points
  // pdf_dim: Number of variates in X_nl
  // X_nl_one_lev: SILHS sample vector at one height level
  // X_mixt_comp_one_lev: Category of each sample point
  SilhsHostInt1D int_sample_category = workspace.view<SilhsHostInt1D>(
      Kokkos::WithoutInitializing, "int_sample_category", num_samples);
  for (int isample = 0; isample < num_samples; ++isample) {
    importance_category_type sample_category;
    if (X_nl_one_lev(isample, hm_metadata.iiPDF_chi) < zero)
      sample_category.l_in_cloud = false;
    else
      sample_category.l_in_cloud = true;
    if (hm_metadata.iiPDF_rr == -2)
      throw std::runtime_error(
          "iiPDF_rr must be greater than zero for the category sampler to work.");
    if (X_nl_one_lev(isample, hm_metadata.iiPDF_rr) > zero)
      sample_category.l_in_precip = true;
    else
      sample_category.l_in_precip = false;
    if (X_mixt_comp_one_lev(isample) == 1)
      sample_category.l_in_component_1 = true;
    else
      sample_category.l_in_component_1 = false;
    int found_category_index = -1;
    for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
      if (importance_categories[icategory].l_in_cloud == sample_category.l_in_cloud &&
          importance_categories[icategory].l_in_precip == sample_category.l_in_precip &&
          importance_categories[icategory].l_in_component_1 == sample_category.l_in_component_1) {
        found_category_index = icategory;
        break;
      }
    }
    if (found_category_index == -1)
      throw std::runtime_error("Fatal error determining category in determine_sample_categories");
    int_sample_category(isample) = found_category_index;
  }
  return int_sample_category;
}
} // namespace clubb_kokkos
