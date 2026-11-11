// Simplified longwave, BOMEX and LBA radiation on CLUBB's bottom-up grid.
// simple_rad computes liquid-water-path flux/heating and optional inversion
// diagnostics; simple_rad_bomex applies the GCSS prescribed tendency.
// simple_rad_lba interpolates case-owned tables loaded by simple_rad_lba_init.
// Fluxes are on momentum levels [W/m^2], heating on thermodynamic levels [K/s].
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "parameters_radiation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

namespace clubb_kokkos {

struct simple_rad_state {
  Real1D lba_zrad; // LBA table altitudes [m],33 levels.
  Real2D lba_krad; // LBA heating [K/s],(height,time),33x36 entries.
};

void simple_rad(ScratchWorkspace& workspace,                            // InOut
                const grid& gr, int ngrdcol, Real2D rho, Real2D rho_zm, // In
                Real2D rtm, Real2D rcm, Real2D exner,                   // In
                stats_type& stats,                                      // InOut
                Real2D Frad_LW, Real2D radht_LW,                        // Out
                const parameters_radiation_state& parameters_radiation  // Native module In
);
void simple_rad_bomex(ScratchWorkspace& workspace,                                         // InOut
                      const grid& gr, int ngrdcol,                                         // In
                      Real2D radht);                                                       // Out
void simple_rad_lba(ScratchWorkspace& workspace,                                           // InOut
                    const grid& gr, int ngrdcol, double time_current, double time_initial, // In
                    Real2D radht,                                                          // Out
                    const simple_rad_state& simple_rad // Native module In
);
void simple_rad_lba_init(int iunit, const std::string& file_path, // In
                         simple_rad_state& simple_rad);           // Native module Out

template <class Rho, class RhoZm, class Rt, class Rc, class Exner, class Flux, class Heating>
void simple_rad(ScratchWorkspace& workspace,                           // InOut
                const grid& gr, int ngrdcol, Rho rho, RhoZm rho_zm,    // In
                Rt rtm, Rc rcm, Exner exner,                           // In
                stats_type& stats,                                     // InOut
                Flux Frad_LW, Heating radht_LW,                        // Out
                const parameters_radiation_state& parameters_radiation // Native module In
) {
  // Shared owner conversions preserve host/strided callers at the API boundary.
  auto workspace_scope = workspace.scope();
  auto flux            = to_execution_space(workspace, Kokkos::WithoutInitializing, Frad_LW);
  auto heating         = to_execution_space(workspace, Kokkos::WithoutInitializing, radht_LW);
  simple_rad(workspace, gr, ngrdcol, to_execution_space(workspace, rho),
             to_execution_space(workspace, rho_zm),                                  // In
             to_execution_space(workspace, rtm), to_execution_space(workspace, rcm), // In
             to_execution_space(workspace, exner),                                   // In
             stats,                                                                  // InOut
             flux, heating, parameters_radiation); // Out/module In
  from_execution_space(workspace, Frad_LW, flux);
  from_execution_space(workspace, radht_LW, heating);
}

template <class Heating>
void simple_rad_bomex(ScratchWorkspace& workspace, // InOut
                      const grid& gr, int ngrdcol, // In
                      Heating radht) {             // Out
  auto staging_scope = workspace.scope();
  auto heating       = to_execution_space(workspace, Kokkos::WithoutInitializing, radht);
  simple_rad_bomex(workspace, gr, ngrdcol, heating);
  from_execution_space(workspace, radht, heating);
}

template <class Heating>
void simple_rad_lba(ScratchWorkspace& workspace,                                           // InOut
                    const grid& gr, int ngrdcol, double time_current, double time_initial, // In
                    Heating radht,                                                         // Out
                    const simple_rad_state& simple_rad // Native module In
) {
  auto workspace_scope = workspace.scope();
  auto heating         = to_execution_space(workspace, Kokkos::WithoutInitializing, radht);
  simple_rad_lba(workspace, gr, ngrdcol, time_current, time_initial, heating, simple_rad);
  from_execution_space(workspace, radht, heating);
}

} // namespace clubb_kokkos
