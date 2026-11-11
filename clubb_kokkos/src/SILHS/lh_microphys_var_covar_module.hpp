// SILHS microphysical variance/covariance source terms in every column.
// lh_microphys_var_covar_driver_api includes optional finite-dt corrections
// while the instantaneous option agrees with the analytic KK source limit.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"

namespace clubb_kokkos {
void lh_microphys_var_covar_driver_api(
    ScratchWorkspace& workspace,                                                      // InOut
    int nzt, int num_samples, int ngrdcol, real_t dt, Real3D lh_sample_point_weights, // In
    const pdf_parameter& pdf_params, Real3D lh_rt_all, Real3D lh_thl_all,             // In
    Real3D lh_w_all,                                                                  // In
    Real3D lh_rcm_mc_all, Real3D lh_rvm_mc_all, Real3D lh_thlm_mc_all,                // In
    bool l_lh_instant_var_covar_src,                                                  // In
    Real2D lh_rtp2_mc_zt, Real2D lh_thlp2_mc_zt, Real2D lh_wprtp_mc_zt,               // Out
    Real2D lh_wpthlp_mc_zt, Real2D lh_rtpthlp_mc_zt                                   // Out
);

template <class Field_lh_sample_point_weights, class Field_lh_rt_all, class Field_lh_thl_all,
          class Field_lh_w_all, class Field_lh_rcm_mc_all, class Field_lh_rvm_mc_all,
          class Field_lh_thlm_mc_all, class Field_lh_rtp2_mc_zt, class Field_lh_thlp2_mc_zt,
          class Field_lh_wprtp_mc_zt, class Field_lh_wpthlp_mc_zt, class Field_lh_rtpthlp_mc_zt>
void lh_microphys_var_covar_driver_api(
    ScratchWorkspace& workspace,                                                   // InOut
    int nzt, int num_samples, int ngrdcol, real_t dt,                              // In
    Field_lh_sample_point_weights lh_sample_point_weights,                         // In
    const pdf_parameter& pdf_params, Field_lh_rt_all lh_rt_all,                    // In
    Field_lh_thl_all lh_thl_all, Field_lh_w_all lh_w_all,                          // In
    Field_lh_rcm_mc_all lh_rcm_mc_all, Field_lh_rvm_mc_all lh_rvm_mc_all,          // In
    Field_lh_thlm_mc_all lh_thlm_mc_all,                                           // In
    bool l_lh_instant_var_covar_src,                                               // In
    Field_lh_rtp2_mc_zt lh_rtp2_mc_zt, Field_lh_thlp2_mc_zt lh_thlp2_mc_zt,        // Out
    Field_lh_wprtp_mc_zt lh_wprtp_mc_zt,                                           // Out
    Field_lh_wpthlp_mc_zt lh_wpthlp_mc_zt, Field_lh_rtpthlp_mc_zt lh_rtpthlp_mc_zt // Out
) {
  // Shared-owner conversion at the host/strided public field boundary.
  // nzt: Number of vertical levels
  // num_samples: Number of SILHS sample points
  // ngrdcol: Number of model columns
  // dt: Model time step                             [s]
  // lh_sample_point_weights: Weight of SILHS sample points
  // pdf_params: The PDF parameters_silhs
  // lh_rt_all: SILHS samples of total water                [kg/kg]
  // lh_thl_all: SILHS samples of potential temperature      [K]
  // lh_w_all: SILHS samples of vertical velocity          [m/s]
  // lh_rcm_mc_all: SILHS microphys. tendency of rcm            [kg/kg/s]
  // lh_rvm_mc_all: SILHS microphys. tendency of rvm            [kg/kg/s]
  // lh_thlm_mc_all: SILHS microphys. tendency of thlm           [K/s]
  // l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies  [-]
  // lh_rtp2_mc_zt: SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
  // lh_thlp2_mc_zt: SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
  // lh_wprtp_mc_zt: SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
  // lh_wpthlp_mc_zt: SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
  // lh_rtpthlp_mc_zt: SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]
  auto workspace_scope = workspace.scope();
  auto lh_rtp2_mc_zt_execution =
      to_execution_space(workspace, Kokkos::WithoutInitializing, lh_rtp2_mc_zt);
  auto lh_thlp2_mc_zt_execution =
      to_execution_space(workspace, Kokkos::WithoutInitializing, lh_thlp2_mc_zt);
  auto lh_wprtp_mc_zt_execution =
      to_execution_space(workspace, Kokkos::WithoutInitializing, lh_wprtp_mc_zt);
  auto lh_wpthlp_mc_zt_execution =
      to_execution_space(workspace, Kokkos::WithoutInitializing, lh_wpthlp_mc_zt);
  auto lh_rtpthlp_mc_zt_execution =
      to_execution_space(workspace, Kokkos::WithoutInitializing, lh_rtpthlp_mc_zt);
  lh_microphys_var_covar_driver_api(
      workspace, nzt, num_samples, ngrdcol, dt,
      to_execution_space(workspace, lh_sample_point_weights), // In
      pdf_params, to_execution_space(workspace, lh_rt_all),
      to_execution_space(workspace, lh_thl_all), to_execution_space(workspace, lh_w_all), // In
      to_execution_space(workspace, lh_rcm_mc_all), to_execution_space(workspace, lh_rvm_mc_all),
      to_execution_space(workspace, lh_thlm_mc_all),                               // In
      l_lh_instant_var_covar_src,                                                  // In
      lh_rtp2_mc_zt_execution, lh_thlp2_mc_zt_execution, lh_wprtp_mc_zt_execution, // Out
      lh_wpthlp_mc_zt_execution, lh_rtpthlp_mc_zt_execution                        // Out
  );
  from_execution_space(workspace, lh_rtp2_mc_zt, lh_rtp2_mc_zt_execution);
  from_execution_space(workspace, lh_thlp2_mc_zt, lh_thlp2_mc_zt_execution);
  from_execution_space(workspace, lh_wprtp_mc_zt, lh_wprtp_mc_zt_execution);
  from_execution_space(workspace, lh_wpthlp_mc_zt, lh_wpthlp_mc_zt_execution);
  from_execution_space(workspace, lh_rtpthlp_mc_zt, lh_rtpthlp_mc_zt_execution);
}

} // namespace clubb_kokkos
