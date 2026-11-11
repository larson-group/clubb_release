// Advance non-BUGS radiation and sample its diagnostics in Native call order.
// advance_clubb_radiation retains soil updates, radiation cadence, fatal returns
// and per-sample stats. radiation_driver dispatches simplified/BOMEX/LBA physics;
// update_radiation_variables records its case-owned fields. SILHS radiation stays
// an explicit unsupported API, matching the existing Native/JAX feature gate.
// Native formals come first; appended state replaces Native module variables.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "parameters_radiation.hpp"
#include "radiation_variables_module.hpp"
#include "simple_rad_module.hpp"
#include "soil_vegetation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

namespace clubb_kokkos {
void advance_clubb_radiation(
    ScratchWorkspace& workspace,                                                    // InOut
    const grid& gr, int ngrdcol, int hydromet_dim, int pdf_dim, int lh_num_samples, // In
    bool l_rad_itime, real_t dt_main, int day, int month, int year,                 // In
    real_t lat_vals, real_t lon_vals, double time_current, double time_initial,     // In
    Real2D rho, Real2D rho_zm, Real2D p_in_Pa, Real2D exner,                        // In
    Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D p_sfc,                              // In
    Real2D cloud_frac, Real2D ice_supersat_frac,                                    // In
    Real2D thlm, Real2D rtm, Real2D rcm,                                            // In
    Real4D X_nl_all_levs, Real3D lh_rt_clipped, Real3D lh_thl_clipped, Real3D lh_rc_clipped, // In
    Real3D lh_sample_point_weights, Real3D hydromet, const hm_metadata_type& hm_metadata,    // In
    stats_type& stats, err_info_type& err_info,                         // InOut
    Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K, Real1D veg_T_in_K, // InOut
    Real2D radht,                                                       // InOut
    const parameters_radiation_state& parameters_radiation,             // Native module In
    radiation_variables_state& radiation_variables,                     // Native module InOut
    const simple_rad_state& simple_rad,                                 // Native module In
    const soil_vegetation_state& soil_vegetation                        // Native module In
);

void radiation_driver(ScratchWorkspace& workspace, const grid& gr, double time_current,
                      double time_initial,
                      int hydromet_dim,                                                // In
                      int ngrdcol,                                                     // In
                      int day, int month, int year,                                    // In
                      real_t lat_vals, real_t lon_vals,                                // In
                      Real2D rho, Real2D rho_zm, Real2D p_in_Pa,                       // In
                      Real2D exner, Real2D cloud_frac, Real2D ice_supersat_frac,       // In
                      Real2D thlm, Real2D rtm, Real2D rcm, Real3D hydromet,            // In
                      const hm_metadata_type& hm_metadata,                             // In
                      stats_type& stats,                                               // InOut
                      err_info_type& err_info,                                         // InOut
                      Real2D radht, Real2D Frad, Real2D Frad_SW_up, Real2D Frad_LW_up, // Out
                      Real2D Frad_SW_down, Real2D Frad_LW_down,                        // Out
                      const parameters_radiation_state& parameters_radiation, // Native module In
                      radiation_variables_state& radiation_variables,         // Native module InOut
                      const simple_rad_state& simple_rad                      // Native module In
);

void update_radiation_variables(
    int ngrdcol, int nzm, int nzt, Real2D radht, Real2D Frad, Real2D Frad_SW_up, // In
    Real2D Frad_LW_up,                                                           // In
    Real2D Frad_SW_down, Real2D Frad_LW_down,                                    // In
    stats_type& stats,                                                           // InOut
    const parameters_radiation_state& parameters_radiation,                      // Native module In
    const radiation_variables_state& radiation_variables                         // Native module In
);

void silhs_radiation_driver(const grid& gr, int ngrdcol, int nzm, int nzt, int lh_num_samples,
                            int pdf_dim,                                           // In
                            int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                            int day, int month, int year,                          // In
                            real_t lat_vals, real_t lon_vals,                      // In
                            double time_current, double time_initial, Real2D rho, Real2D rho_zm,
                            Real2D p_in_Pa,                                                    // In
                            Real2D exner,                                                      // In
                            Real2D cloud_frac, Real2D ice_supersat_frac, Real4D X_nl_all_levs, // In
                            Real3D lh_rt_clipped, Real3D lh_thl_clipped, Real3D lh_rc_clipped, // In
                            Real3D lh_sample_point_weights, Real3D hydromet,                   // In
                            err_info_type& err_info, // InOut
                            Real2D radht, Real2D Frad, Real2D Frad_SW_up, Real2D Frad_LW_up,
                            Real2D Frad_SW_down, // Out
                            Real2D Frad_LW_down  // Out
);

template <class Field_rho, class Field_rho_zm, class Field_p_in_Pa, class Field_exner,
          class Field_wpthlp_sfc, class Field_wprtp_sfc, class Field_p_sfc, class Field_cloud_frac,
          class Field_ice_supersat_frac, class Field_thlm, class Field_rtm, class Field_rcm,
          class Field_X_nl_all_levs, class Field_lh_rt_clipped, class Field_lh_thl_clipped,
          class Field_lh_rc_clipped, class Field_lh_sample_point_weights, class Field_hydromet,
          class Field_deep_soil_T_in_K, class Field_sfc_soil_T_in_K, class Field_veg_T_in_K,
          class Field_radht>
void advance_clubb_radiation(
    ScratchWorkspace& workspace,                                                    // InOut
    const grid& gr, int ngrdcol, int hydromet_dim, int pdf_dim, int lh_num_samples, // In
    bool l_rad_itime, real_t dt_main, int day, int month, int year,                 // In
    real_t lat_vals, real_t lon_vals, double time_current, double time_initial,     // In
    Field_rho rho, Field_rho_zm rho_zm, Field_p_in_Pa p_in_Pa, Field_exner exner,   // In
    Field_wpthlp_sfc wpthlp_sfc, Field_wprtp_sfc wprtp_sfc, Field_p_sfc p_sfc,      // In
    Field_cloud_frac cloud_frac, Field_ice_supersat_frac ice_supersat_frac,         // In
    Field_thlm thlm, Field_rtm rtm, Field_rcm rcm,                                  // In
    Field_X_nl_all_levs X_nl_all_levs, Field_lh_rt_clipped lh_rt_clipped,           // In
    Field_lh_thl_clipped lh_thl_clipped, Field_lh_rc_clipped lh_rc_clipped,         // In
    Field_lh_sample_point_weights lh_sample_point_weights, Field_hydromet hydromet, // In
    const hm_metadata_type& hm_metadata,                                            // In
    stats_type& stats, err_info_type& err_info,                                     // InOut
    Field_deep_soil_T_in_K deep_soil_T_in_K, Field_sfc_soil_T_in_K sfc_soil_T_in_K, // InOut
    Field_veg_T_in_K veg_T_in_K,                                                    // InOut
    Field_radht radht,                                                              // InOut
    const parameters_radiation_state& parameters_radiation, // Native module In
    radiation_variables_state& radiation_variables,         // Native module InOut
    const simple_rad_state& simple_rad,                     // Native module In
    const soil_vegetation_state& soil_vegetation            // Native module In
) {
  // Use the shared owner only at host/strided API boundaries.
  auto workspace_scope            = workspace.scope();
  auto deep_soil_T_in_K_execution = to_execution_space(workspace, deep_soil_T_in_K);
  auto sfc_soil_T_in_K_execution  = to_execution_space(workspace, sfc_soil_T_in_K);
  auto veg_T_in_K_execution       = to_execution_space(workspace, veg_T_in_K);
  auto radht_execution            = to_execution_space(workspace, radht);
  advance_clubb_radiation(
      workspace, gr, ngrdcol, hydromet_dim, pdf_dim, lh_num_samples, l_rad_itime, dt_main, day,
      month, year, lat_vals, lon_vals, time_current, time_initial,
      to_execution_space(workspace, rho), to_execution_space(workspace, rho_zm),
      to_execution_space(workspace, p_in_Pa), to_execution_space(workspace, exner),
      to_execution_space(workspace, wpthlp_sfc), to_execution_space(workspace, wprtp_sfc),
      to_execution_space(workspace, p_sfc), to_execution_space(workspace, cloud_frac),
      to_execution_space(workspace, ice_supersat_frac), to_execution_space(workspace, thlm),
      to_execution_space(workspace, rtm), to_execution_space(workspace, rcm),
      to_execution_space(workspace, X_nl_all_levs), to_execution_space(workspace, lh_rt_clipped),
      to_execution_space(workspace, lh_thl_clipped), to_execution_space(workspace, lh_rc_clipped),
      to_execution_space(workspace, lh_sample_point_weights),
      to_execution_space(workspace, hydromet), hm_metadata, stats, err_info,
      deep_soil_T_in_K_execution, sfc_soil_T_in_K_execution, veg_T_in_K_execution, radht_execution,
      parameters_radiation, radiation_variables, simple_rad, soil_vegetation);
  from_execution_space(workspace, deep_soil_T_in_K, deep_soil_T_in_K_execution);
  from_execution_space(workspace, sfc_soil_T_in_K, sfc_soil_T_in_K_execution);
  from_execution_space(workspace, veg_T_in_K, veg_T_in_K_execution);
  from_execution_space(workspace, radht, radht_execution);
}

template <class Field_rho, class Field_rho_zm, class Field_p_in_Pa, class Field_exner,
          class Field_cloud_frac, class Field_ice_supersat_frac, class Field_thlm, class Field_rtm,
          class Field_rcm, class Field_hydromet, class Field_radht, class Field_Frad,
          class Field_Frad_SW_up, class Field_Frad_LW_up, class Field_Frad_SW_down,
          class Field_Frad_LW_down>
void radiation_driver(ScratchWorkspace& workspace, const grid& gr, double time_current,
                      double time_initial,
                      int hydromet_dim,                                          // In
                      int ngrdcol,                                               // In
                      int day, int month, int year,                              // In
                      real_t lat_vals, real_t lon_vals,                          // In
                      Field_rho rho, Field_rho_zm rho_zm, Field_p_in_Pa p_in_Pa, // In
                      Field_exner exner, Field_cloud_frac cloud_frac,
                      Field_ice_supersat_frac ice_supersat_frac,                              // In
                      Field_thlm thlm, Field_rtm rtm, Field_rcm rcm, Field_hydromet hydromet, // In
                      const hm_metadata_type& hm_metadata,                                    // In
                      stats_type& stats,                                                // InOut
                      err_info_type& err_info,                                          // InOut
                      Field_radht radht, Field_Frad Frad, Field_Frad_SW_up Frad_SW_up,  // Out
                      Field_Frad_LW_up Frad_LW_up,                                      // Out
                      Field_Frad_SW_down Frad_SW_down, Field_Frad_LW_down Frad_LW_down, // Out
                      const parameters_radiation_state& parameters_radiation, // Native module In
                      radiation_variables_state& radiation_variables,         // Native module InOut
                      const simple_rad_state& simple_rad                      // Native module In
) {
  // Use the shared owner only at host/strided API boundaries.
  auto workspace_scope        = workspace.scope();
  auto radht_execution        = to_execution_space(workspace, radht, false);
  auto Frad_execution         = to_execution_space(workspace, Frad, false);
  auto Frad_SW_up_execution   = to_execution_space(workspace, Frad_SW_up, false);
  auto Frad_LW_up_execution   = to_execution_space(workspace, Frad_LW_up, false);
  auto Frad_SW_down_execution = to_execution_space(workspace, Frad_SW_down, false);
  auto Frad_LW_down_execution = to_execution_space(workspace, Frad_LW_down, false);
  radiation_driver(workspace, gr, time_current, time_initial, hydromet_dim, ngrdcol, day, month,
                   year, lat_vals, lon_vals, to_execution_space(workspace, rho),
                   to_execution_space(workspace, rho_zm), to_execution_space(workspace, p_in_Pa),
                   to_execution_space(workspace, exner), to_execution_space(workspace, cloud_frac),
                   to_execution_space(workspace, ice_supersat_frac),
                   to_execution_space(workspace, thlm), to_execution_space(workspace, rtm),
                   to_execution_space(workspace, rcm), to_execution_space(workspace, hydromet),
                   hm_metadata, stats, err_info, radht_execution, Frad_execution,
                   Frad_SW_up_execution, Frad_LW_up_execution, Frad_SW_down_execution,
                   Frad_LW_down_execution, parameters_radiation, radiation_variables, simple_rad);
  from_execution_space(workspace, radht, radht_execution);
  from_execution_space(workspace, Frad, Frad_execution);
  from_execution_space(workspace, Frad_SW_up, Frad_SW_up_execution);
  from_execution_space(workspace, Frad_LW_up, Frad_LW_up_execution);
  from_execution_space(workspace, Frad_SW_down, Frad_SW_down_execution);
  from_execution_space(workspace, Frad_LW_down, Frad_LW_down_execution);
}

template <class Field_radht, class Field_Frad, class Field_Frad_SW_up, class Field_Frad_LW_up,
          class Field_Frad_SW_down, class Field_Frad_LW_down>
void update_radiation_variables(
    int ngrdcol, int nzm, int nzt, Field_radht radht, Field_Frad Frad, // In
    Field_Frad_SW_up Frad_SW_up, Field_Frad_LW_up Frad_LW_up,          // In
    Field_Frad_SW_down Frad_SW_down, Field_Frad_LW_down Frad_LW_down,  // In
    stats_type& stats,                                                 // InOut
    const parameters_radiation_state& parameters_radiation,            // Native module In
    const radiation_variables_state& radiation_variables               // Native module In
) {
  // Use the shared owner only at host/strided API boundaries.
  update_radiation_variables(ngrdcol, nzm, nzt, to_execution_space(radht), to_execution_space(Frad),
                             to_execution_space(Frad_SW_up), to_execution_space(Frad_LW_up),
                             to_execution_space(Frad_SW_down), to_execution_space(Frad_LW_down),
                             stats, parameters_radiation, radiation_variables);
}

} // namespace clubb_kokkos
