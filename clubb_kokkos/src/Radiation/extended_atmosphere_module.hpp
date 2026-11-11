// Construct and release the atmosphere extension shared by sounding/radiation.
//
// extended_atmosphere_state owns the Native module arrays for one CLUBB case.
// determine_extended_atmos_bounds selects source-profile bounds and builds the
// complete thermodynamic/momentum grids; finalize_extended_atm releases arrays.
// Heights increase with index. Returned source bounds are zero-based.
// Reference: McClatchey et al. (1972), Environmental Research Papers411,p.94.
#pragma once

#include "CLUBB_core/array_types.hpp"
#include "parameters_radiation.hpp"

namespace clubb_kokkos {

struct extended_atmosphere_state {
  int extended_atmos_dim{}; // Number of levels in the extension profile [-].
  int total_atmos_dim{};    // Computational + buffer + extended levels [-].

  Real1D complete_alt;      // Complete thermodynamic-grid altitudes [m].
  Real1D complete_momentum; // Complete momentum-grid altitudes [m].

  Real1D extended_alt;      // Extension altitudes, increasing with index [m].
  Real1D extended_T_in_K;   // Extension temperature [K].
  Real1D extended_sp_hmdty; // Extension specific humidity [kg/kg].
  Real1D extended_p_in_mb;  // Extension pressure [mb].
  Real1D extended_o3l;      // Extension ozone mass fraction [kg/kg].
};

void determine_extended_atmos_bounds(
    ScratchWorkspace& workspace,                               // InOut
    int grid_size, Real1D zt_grid,                             // In
    Real1D zm_grid, Real1D zm_grid_spacing, Real1D p_in_Pa_zm, // In
    int& extended_atmos_bottom_level,                          // Out
    int& extended_atmos_top_level,                             // Out
    int& extended_atmos_range_size,                            // Out
    int& lin_int_buffer_size,                                  // Out
    const parameters_radiation_state& parameters_radiation,    // Native module input
    extended_atmosphere_state& extended_atmosphere             // Native module InOut
);

void finalize_extended_atm(extended_atmosphere_state& extended_atmosphere);

template <class Zt, class Zm, class Spacing, class Pressure>
void determine_extended_atmos_bounds(
    ScratchWorkspace& workspace,                              // InOut
    int grid_size, Zt zt_grid,                                // In
    Zm zm_grid, Spacing zm_grid_spacing, Pressure p_in_Pa_zm, // In
    int& extended_atmos_bottom_level,                         // Out
    int& extended_atmos_top_level,                            // Out
    int& extended_atmos_range_size,                           // Out
    int& lin_int_buffer_size,                                 // Out
    const parameters_radiation_state& parameters_radiation,   // Native module input
    extended_atmosphere_state& extended_atmosphere            // Native module InOut
) {
  // Pack host/strided input profiles at this public storage boundary; callers
  // already using the shared device array owner retain their allocations.
  auto workspace_scope = workspace.scope();
  determine_extended_atmos_bounds(
      workspace, grid_size, to_execution_space(workspace, zt_grid),                           // In
      to_execution_space(workspace, zm_grid), to_execution_space(workspace, zm_grid_spacing), // In
      to_execution_space(workspace, p_in_Pa_zm),                                              // In
      extended_atmos_bottom_level, extended_atmos_top_level,                                  // Out
      extended_atmos_range_size, lin_int_buffer_size,                                         // Out
      parameters_radiation, extended_atmosphere);
}

} // namespace clubb_kokkos
