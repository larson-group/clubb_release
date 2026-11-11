// Estimate microphysics spread within SILHS importance-sampling categories.
//
// silhs_category_variance_driver writes the category-conditioned rain-tendency
// root mean square into C++ stats. silhs_sample_category_variance fills the
// same fields without updating stats. Samples are (column,sample,level,variate),
// with column-first LayoutLeft storage and zero-based metadata slots.
// The eight categories come from the SILHS importance-sampling owner. Sample
// accumulation remains ordered; empty-probability categories receive -999.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "SILHS/silhs_importance_sample_module.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <stdexcept>
namespace clubb_kokkos {
using Int3D = Kokkos::View<int***, Kokkos::LayoutLeft>;
inline void silhs_sample_category_variance(ScratchWorkspace& workspace, int ngrdcol, int nzt,
                                           int num_samples, int pdf_dim, Real4D X_nl_all_levs,
                                           Int3D X_mixt_comp_all_levs, Real3D samples_all,
                                           Real3D lh_sample_point_weights,
                                           const pdf_parameter& pdf_params,
                                           const precipitation_fractions& precip_fracs,
                                           const hm_metadata_type& hm_metadata,
                                           Real3D root_weight_mean_sq_cat);

inline void
silhs_category_variance_driver(ScratchWorkspace& workspace, int ngrdcol, int nzt, int num_samples,
                               int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
                               Real4D X_nl_all_levs, Int3D X_mixt_comp_all_levs,
                               Real4D lh_hydromet_mc_all, Real3D lh_sample_point_weights,
                               const pdf_parameter& pdf_params,
                               const precipitation_fractions& precip_fracs, stats_type& stats) {
  // Compute the variance of a microphysics variable in each importance category
  // for all columns. Input hydrometeor tendencies are [(hydromet units)/s];
  // output statistics retain the source root-mean-square units and names.
  auto workspace_scope = workspace.scope();
  Real3D samples_all   = workspace_scope.view<Real3D>("silhs_category_variance_driver::samples_all",
                                                      ngrdcol, num_samples, nzt);
  Real3D root_weight_mean_sq_cat =
      workspace_scope.view<Real3D>("silhs_category_variance_driver::root_weight_mean_sq_cat",
                                   ngrdcol, nzt, num_importance_categories);
  // Sample rrm_mc from the hydrometeor tendencies generated for each SILHS point.
  Kokkos::deep_copy(samples_all, Kokkos::subview(lh_hydromet_mc_all, Kokkos::ALL, Kokkos::ALL,
                                                 Kokkos::ALL, hm_metadata.iirr));
  silhs_sample_category_variance(workspace, ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs,
                                 X_mixt_comp_all_levs, samples_all, lh_sample_point_weights,
                                 pdf_params, precip_fracs, hm_metadata, root_weight_mean_sq_cat);
  if (stats.l_sample) {
    for (int icat = 0; icat < num_importance_categories; ++icat) {
      const std::string cat_num  = std::to_string((((icat + 1))));
      const std::string cat_name = "silhs_var_cat_" + cat_num;
      stats_update(cat_name,
                   Kokkos::subview(root_weight_mean_sq_cat, Kokkos::ALL, Kokkos::ALL, icat), stats);
    }
  }
}

inline void silhs_sample_category_variance(ScratchWorkspace& workspace, int ngrdcol, int nzt,
                                           int num_samples, int pdf_dim, Real4D X_nl_all_levs,
                                           Int3D X_mixt_comp_all_levs, Real3D samples_all,
                                           Real3D lh_sample_point_weights,
                                           const pdf_parameter& pdf_params,
                                           const precipitation_fractions& precip_fracs,
                                           const hm_metadata_type& hm_metadata,
                                           Real3D root_weight_mean_sq_cat) {
  // Compute category-conditioned root-mean-square values in every column.
  // Sample weights and the category's actual PDF probability define the mean;
  // each column/level owns an independent ordered sample accumulation.
  auto workspace_scope = workspace.scope();
  if (hm_metadata.iiPDF_rr == -2)
    throw std::runtime_error(
        "iiPDF_rr must be greater than zero for the category sampler to work.");
  Kokkos::deep_copy(root_weight_mean_sq_cat, zero);
  const auto importance_categories = define_importance_categories();
  Int3D int_sample_category        = workspace_scope.view<Int3D>(
      "silhs_sample_category_variance::int_sample_category", ngrdcol, num_samples, nzt);
  const int iiPDF_chi = hm_metadata.iiPDF_chi, iiPDF_rr = hm_metadata.iiPDF_rr;
  const auto mixt_frac    = pdf_params.mixt_frac;
  const auto cloud_frac_1 = pdf_params.cloud_frac_1, cloud_frac_2 = pdf_params.cloud_frac_2;
  const auto precip_frac_1 = precip_fracs.precip_frac_1, precip_frac_2 = precip_fracs.precip_frac_2;
  // Keep sample classifications consistent with define_importance_categories.
  // Columns/levels are independent; each category's source sample order remains serial.
  Kokkos::parallel_for(
      "silhs_sample_category_variance::columns", ((ngrdcol * nzt)), KOKKOS_LAMBDA(int cell) {
        const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
        for (int isample = 0; isample < num_samples; ++isample) {
          const bool l_in_cloud       = X_nl_all_levs(i, isample, k, iiPDF_chi) >= zero;
          const bool l_in_precip      = X_nl_all_levs(i, isample, k, iiPDF_rr) > zero;
          const bool l_in_component_1 = X_mixt_comp_all_levs(i, isample, k) == 1;
          int found_category_index    = -1;
          for (int icat = 0; icat < num_importance_categories; ++icat) {
            const auto cat = importance_categories[icat];
            if (found_category_index == -1 && cat.l_in_cloud == l_in_cloud &&
                cat.l_in_precip == l_in_precip && cat.l_in_component_1 == l_in_component_1)
              found_category_index = icat;
          }
          // Native categories enumerate all eight triples, so every sample matches.
          if (found_category_index == -1)
            Kokkos::abort("Fatal error determining category in determine_sample_categories");
          int_sample_category(i, isample, k) = found_category_index;
        }
        real_t category_real_probs[num_importance_categories];
        for (int icat = 0; icat < num_importance_categories; ++icat) {
          const auto cat = importance_categories[icat];
          if (cat.l_in_component_1) {
            if (cat.l_in_cloud) {
              if (cat.l_in_precip)
                category_real_probs[icat] =
                    mixt_frac(i, k) * cloud_frac_1(i, k) * precip_frac_1(i, k);
              else
                category_real_probs[icat] =
                    mixt_frac(i, k) * cloud_frac_1(i, k) * (one - precip_frac_1(i, k));
            } else {
              if (cat.l_in_precip)
                category_real_probs[icat] =
                    mixt_frac(i, k) * (one - cloud_frac_1(i, k)) * precip_frac_1(i, k);
              else
                category_real_probs[icat] =
                    mixt_frac(i, k) * (one - cloud_frac_1(i, k)) * (one - precip_frac_1(i, k));
            }
          } else {
            if (cat.l_in_cloud) {
              if (cat.l_in_precip)
                category_real_probs[icat] =
                    (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * precip_frac_2(i, k);
              else
                category_real_probs[icat] =
                    (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * (one - precip_frac_2(i, k));
            } else {
              if (cat.l_in_precip)
                category_real_probs[icat] =
                    (one - mixt_frac(i, k)) * (one - cloud_frac_2(i, k)) * precip_frac_2(i, k);
              else
                category_real_probs[icat] = (one - mixt_frac(i, k)) * (one - cloud_frac_2(i, k)) *
                                            (one - precip_frac_2(i, k));
            }
          }
        }
        for (int isample = 0; isample < num_samples; ++isample) {
          const int icat = int_sample_category(i, isample, k);
          root_weight_mean_sq_cat(i, k, icat) =
              root_weight_mean_sq_cat(i, k, icat) +
              lh_sample_point_weights(i, isample, k) *
                  (samples_all(i, isample, k) * samples_all(i, isample, k));
        }
        for (int icat = 0; icat < num_importance_categories; ++icat) {
          root_weight_mean_sq_cat(i, k, icat) =
              root_weight_mean_sq_cat(i, k, icat) / static_cast<real_t>(num_samples);
          if (category_real_probs[icat] > zero)
            root_weight_mean_sq_cat(i, k, icat) =
                Kokkos::sqrt(root_weight_mean_sq_cat(i, k, icat) / category_real_probs[icat]);
          else
            root_weight_mean_sq_cat(i, k, icat) = real_t(-999.);
        }
      });
}
// Stage host or strided arrays only at the Native API boundary. The PDF records
// retain their device Views; formal outputs are restored to the caller allocation.
template <class X, class C, class H, class W>
inline void
silhs_category_variance_driver(ScratchWorkspace& workspace, int ngrdcol, int nzt, int num_samples,
                               int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
                               X X_nl_all_levs, C X_mixt_comp_all_levs, H lh_hydromet_mc_all,
                               W lh_sample_point_weights, const pdf_parameter& pdf_params,
                               const precipitation_fractions& precip_fracs, stats_type& stats) {
  auto workspace_scope = workspace.scope();
  auto x               = to_execution_space(workspace, X_nl_all_levs);
  auto c               = to_execution_space(workspace, X_mixt_comp_all_levs);
  auto h               = to_execution_space(workspace, lh_hydromet_mc_all);
  auto w               = to_execution_space(workspace, lh_sample_point_weights);
  silhs_category_variance_driver(workspace, ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim,
                                 hm_metadata, x, c, h, w, pdf_params, precip_fracs, stats);
}
template <class X, class C, class S, class W, class O>
inline void
silhs_sample_category_variance(ScratchWorkspace& workspace, int ngrdcol, int nzt, int num_samples,
                               int pdf_dim, X X_nl_all_levs, C X_mixt_comp_all_levs, S samples_all,
                               W lh_sample_point_weights, const pdf_parameter& pdf_params,
                               const precipitation_fractions& precip_fracs,
                               const hm_metadata_type& hm_metadata, O root_weight_mean_sq_cat) {
  auto workspace_scope = workspace.scope();
  auto x               = to_execution_space(workspace, X_nl_all_levs);
  auto c               = to_execution_space(workspace, X_mixt_comp_all_levs);
  auto sample          = to_execution_space(workspace, samples_all);
  auto w               = to_execution_space(workspace, lh_sample_point_weights);
  auto o               = to_execution_space(workspace, root_weight_mean_sq_cat, false);
  silhs_sample_category_variance(workspace, ngrdcol, nzt, num_samples, pdf_dim, x, c, sample, w,
                                 pdf_params, precip_fracs, hm_metadata, o);
  from_execution_space(workspace, root_weight_mean_sq_cat, o);
}
} // namespace clubb_kokkos
