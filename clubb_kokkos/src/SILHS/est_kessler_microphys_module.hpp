// Kessler autoconversion estimates from sampled and analytic cloud water.
// est_kessler_microphys_api returns gridbox means and standard deviations
// [kg/kg/s], plus the sampled cloud-water mean [kg/kg], for every column.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"

namespace clubb_kokkos {
void est_kessler_microphys_api(ScratchWorkspace& workspace,                        // InOut
                               int nzt, int num_samples, int pdf_dim, int ngrdcol, // In
                               Real4D X_nl_all_levs, const pdf_parameter& pdf_params, Real2D rcm,
                               Real2D cloud_frac,                                             // In
                               Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs, // In
                               Real3D lh_sample_point_weights,                                // In
                               bool l_lh_importance_sampling,                                 // In
                               Real2D lh_AKm, Real2D AKm, Real2D AKstd, Real2D AKstd_cld,     // Out
                               Real2D AKm_rcm, Real2D AKm_rcc, Real2D lh_rcm_avg              // Out
);

template <class Field_X_nl_all_levs, class Field_rcm, class Field_cloud_frac,
          class Field_X_mixt_comp_all_levs, class Field_lh_sample_point_weights, class Field_lh_AKm,
          class Field_AKm, class Field_AKstd, class Field_AKstd_cld, class Field_AKm_rcm,
          class Field_AKm_rcc, class Field_lh_rcm_avg>
void est_kessler_microphys_api(ScratchWorkspace& workspace,                        // InOut
                               int nzt, int num_samples, int pdf_dim, int ngrdcol, // In
                               Field_X_nl_all_levs X_nl_all_levs, const pdf_parameter& pdf_params,
                               Field_rcm rcm,                                         // In
                               Field_cloud_frac cloud_frac,                           // In
                               Field_X_mixt_comp_all_levs X_mixt_comp_all_levs,       // In
                               Field_lh_sample_point_weights lh_sample_point_weights, // In
                               bool l_lh_importance_sampling,                         // In
                               Field_lh_AKm lh_AKm, Field_AKm AKm, Field_AKstd AKstd,
                               Field_AKstd_cld AKstd_cld, // Out
                               Field_AKm_rcm AKm_rcm, Field_AKm_rcc AKm_rcc,
                               Field_lh_rcm_avg lh_rcm_avg // Out
) {
  // Shared-owner conversion at the host/strided public field boundary.
  // nzt: Number of vertical levels
  // num_samples: Number of sample points
  // pdf_dim: Number of variates
  // ngrdcol: Number of model columns
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // pdf_params: PDF parameters       [units vary]
  // rcm: Liquid water mixing ratio                [kg/kg]
  // cloud_frac: Cloud fraction           [-]
  // X_mixt_comp_all_levs: Whether we're in mixture component 1 or 2
  // lh_sample_point_weights: Weight for cloud weighted sampling
  // l_lh_importance_sampling: Do importance sampling (SILHS) [-]
  // lh_AKm: Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
  // AKm: Exact Kessler autoconversion, AKm,             [kg/kg/s]
  // AKstd: Exact standard deviation of gba Kessler        [kg/kg/s]
  // AKstd_cld: Exact w/in cloud std of gba Kessler            [kg/kg/s]
  // AKm_rcm: Exact local gba Kessler auto based on rcm      [kg/kg/s]
  // AKm_rcc: Exact local gba Kessler based on w/in cloud rc [kg/kg/s]
  // lh_rcm_avg: LH estimate of grid box avg liquid water [kg/kg]
  auto workspace_scope      = workspace.scope();
  auto lh_AKm_execution     = to_execution_space(workspace, lh_AKm, false);
  auto AKm_execution        = to_execution_space(workspace, AKm, false);
  auto AKstd_execution      = to_execution_space(workspace, AKstd, false);
  auto AKstd_cld_execution  = to_execution_space(workspace, AKstd_cld, false);
  auto AKm_rcm_execution    = to_execution_space(workspace, AKm_rcm, false);
  auto AKm_rcc_execution    = to_execution_space(workspace, AKm_rcc, false);
  auto lh_rcm_avg_execution = to_execution_space(workspace, lh_rcm_avg, false);
  est_kessler_microphys_api(
      workspace, nzt, num_samples, pdf_dim, ngrdcol, // In
      to_execution_space(workspace, X_nl_all_levs), pdf_params, to_execution_space(workspace, rcm),
      to_execution_space(workspace, cloud_frac), // In
      to_execution_space(workspace, X_mixt_comp_all_levs),
      to_execution_space(workspace, lh_sample_point_weights),                // In
      l_lh_importance_sampling,                                              // In
      lh_AKm_execution, AKm_execution, AKstd_execution, AKstd_cld_execution, // Out
      AKm_rcm_execution, AKm_rcc_execution, lh_rcm_avg_execution             // Out
  );
  from_execution_space(workspace, lh_AKm, lh_AKm_execution);
  from_execution_space(workspace, AKm, AKm_execution);
  from_execution_space(workspace, AKstd, AKstd_execution);
  from_execution_space(workspace, AKstd_cld, AKstd_cld_execution);
  from_execution_space(workspace, AKm_rcm, AKm_rcm_execution);
  from_execution_space(workspace, AKm_rcc, AKm_rcc_execution);
  from_execution_space(workspace, lh_rcm_avg, lh_rcm_avg_execution);
}

} // namespace clubb_kokkos
