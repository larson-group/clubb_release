// Keep one case's microphysics settings and sedimentation flags.
//
// parameters_microphys mirrors the Native module's fields/defaults in declaration
// order. A ClubbDriver owns one value, restores it before each initialization and
// appends it to coordinator calls that read implicit Native module state. The
// l_hydromet_sed View is species indexed; atmospheric fields remain caller owned.
// No mutable model configuration is stored in C++ namespace globals.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include <string>
#include <cstdint>
namespace clubb_kokkos {
inline constexpr int lh_microphys_interactive     = 1; // Samples feed back into the model.
inline constexpr int lh_microphys_non_interactive = 2; // Diagnostic-only sampled tendencies.
inline constexpr int lh_microphys_disabled        = 3; // Disable sampled microphysics.
inline constexpr int morrison_no_aerosol          = 0;
inline constexpr int morrison_power_law           = 1;
inline constexpr int morrison_lognormal           = 2;
struct parameters_microphys {
  // Constant parameters are shared; this object's fields replace Native module state.
  bool l_cloud_sed       = false; // Cloud water sedimentation (KK/no microphysics).
  bool l_ice_microphys   = false; // Include ice species (Morrison).
  bool l_upwind_diff_sed = false; // Upwind differencing for sedimentation.
  bool l_graupel         = false; // Include graupel.
  bool l_hail           = false; // Morrison graupel/hail properties; outside supported JAX options.
  bool l_seifert_beheng = false; // Seifert-Beheng warm drizzle; outside JAX options.
  bool l_predict_Nc     = false; // Predict droplet concentration.
  bool l_subgrid_w      = true;  // Add subgrid vertical velocity in Morrison.
  bool l_arctic_nucl    = false; // Arctic nucleation; outside supported JAX options.
  bool l_fix_pgam       = false; // Fixed droplet gamma shape; outside supported JAX options.
  bool l_in_cloud_Nc_diff = true;  // Diffuse in-cloud Nc rather than grid-mean Nc.
  bool l_var_covar_src    = false; // Upscaled microphysics variance/covariance sources.
  // Mean-only clipping allows SILHS convergence to analytic KK (source ticket558).
  bool l_silhs_KK_convergence_adj_mean = false;
  bool l_cloud_edge_activation         = false; // Lateral activation; outside JAX options.
  bool l_local_kk                      = false; // Use local KK drizzle formulas.
  std::string specify_aerosol          = "morrison_lognormal";
  int lh_num_samples                   = 2;    // Number of SILHS microphysics calls per timestep.
  int lh_sequence_length               = 1;    // Steps before SILHS sequence repeats.
  std::int32_t lh_seed                 = 5489; // Native selected_int_kind(9) random seed.
  int lh_microphys_type                = 3;    // Interactive/non-interactive/disabled selector.
  std::string microphys_scheme         = "none";
  Bool1D l_hydromet_sed; // Whether each mean hydrometeor sediments.
  // Native supplies this uninitialized module flag during init_microphys; C++
  // initializes it to the value Native explicitly restores before namelist I/O.
  bool l_gfdl_activation      = false;
  double microphys_start_time = 0.;            // Delayed microphysics start [s].
  real_t Nc0_in_cloud         = real_t(1.0e8); // Initial in-cloud concentration [num/m^3].
  real_t sigma_g              = real_t(1.5);   // Geometric droplet spread for Stokes settling [-].
};
} // namespace clubb_kokkos
