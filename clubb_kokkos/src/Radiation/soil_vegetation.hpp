// Advance soil/vegetation temperatures using surface radiation and heat fluxes.
// initialize_soil_veg sets the GABLS3 reference temperatures and disables the
// scheme; advance_soil_veg applies the surface energy budget for each column.
// Temperatures are [K], radiation and soil heat fluxes are [W/m^2], dt is [s].
// Mutable configuration belongs to the case's soil_vegetation_state.
// Reference: Duynkerke (1991), Monthly Weather Review119,324–341.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

namespace clubb_kokkos {

struct soil_vegetation_state {
  bool l_soil_veg{}; // Enable the soil/vegetation surface energy scheme.
};

void advance_soil_veg(ScratchWorkspace& workspace,                       // InOut
                      int ngrdcol, real_t dt, Real1D rho_sfc,            // In
                      Real1D Frad_SW_up_sfc, Real1D Frad_SW_down_sfc,    // In
                      Real1D Frad_LW_down_sfc,                           // In
                      Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D p_sfc, // In
                      stats_type& stats,                                 // InOut
                      Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K,   // InOut
                      Real1D veg_T_in_K                                  // InOut
);

void initialize_soil_veg(int ngrdcol,                                                        // In
                         Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K, Real1D veg_T_in_K, // Out
                         soil_vegetation_state& soil_vegetation // Native module Out
);

template <class Rho, class SWUp, class SWDown, class LWDown, class ThlFlux, class RtFlux,
          class Pressure, class Deep, class Surface, class Vegetation>
void advance_soil_veg(ScratchWorkspace& workspace,                          // InOut
                      int ngrdcol, real_t dt, Rho rho_sfc,                  // In
                      SWUp Frad_SW_up_sfc, SWDown Frad_SW_down_sfc,         // In
                      LWDown Frad_LW_down_sfc,                              // In
                      ThlFlux wpthlp_sfc, RtFlux wprtp_sfc, Pressure p_sfc, // In
                      stats_type& stats,                                    // InOut
                      Deep deep_soil_T_in_K, Surface sfc_soil_T_in_K,       // InOut
                      Vegetation veg_T_in_K                                 // InOut
) {
  // Shared array-owner conversion is needed only at host/strided boundaries.
  auto workspace_scope = workspace.scope();
  auto deep            = to_execution_space(workspace, deep_soil_T_in_K);
  auto surface         = to_execution_space(workspace, sfc_soil_T_in_K);
  auto vegetation      = to_execution_space(workspace, veg_T_in_K);
  advance_soil_veg(workspace, ngrdcol, dt, to_execution_space(workspace, rho_sfc), // In
                   to_execution_space(workspace, Frad_SW_up_sfc),
                   to_execution_space(workspace, Frad_SW_down_sfc), // In
                   to_execution_space(workspace, Frad_LW_down_sfc), // In
                   to_execution_space(workspace, wpthlp_sfc),
                   to_execution_space(workspace, wprtp_sfc), // In
                   to_execution_space(workspace, p_sfc),     // In
                   stats, deep, surface, vegetation);        // InOut
  from_execution_space(workspace, deep_soil_T_in_K, deep);
  from_execution_space(workspace, sfc_soil_T_in_K, surface);
  from_execution_space(workspace, veg_T_in_K, vegetation);
}

template <class Deep, class Surface, class Vegetation>
void initialize_soil_veg(int ngrdcol, // In
                         Deep deep_soil_T_in_K, Surface sfc_soil_T_in_K,
                         Vegetation veg_T_in_K,                 // Out
                         soil_vegetation_state& soil_vegetation // Native module Out
) {
  auto deep       = to_execution_space(deep_soil_T_in_K, false);
  auto surface    = to_execution_space(sfc_soil_T_in_K, false);
  auto vegetation = to_execution_space(veg_T_in_K, false);
  initialize_soil_veg(ngrdcol, deep, surface, vegetation, soil_vegetation);
  from_execution_space(deep_soil_T_in_K, deep);
  from_execution_space(sfc_soil_T_in_K, surface);
  from_execution_space(veg_T_in_K, vegetation);
}

} // namespace clubb_kokkos
