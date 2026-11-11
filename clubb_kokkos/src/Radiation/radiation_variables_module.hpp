// Allocate, reset and release radiation arrays owned by one CLUBB case.
// setup/reset/cleanup_radiation_variables manage the regular flux [W/m^2] and
// heating [K/s] fields. The corresponding bugsrad routines manage optional
// extended-grid diagnostics; allocation is not an implementation of BUGSrad.
// Array shapes retain Native (column,level) order with explicit LayoutLeft.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include "extended_atmosphere_module.hpp"

namespace clubb_kokkos {
using RadiationDouble2D = Kokkos::View<double**, Kokkos::LayoutLeft>;

struct radiation_variables_state {
  Real2D Frad;                             // Momentum-level radiative flux [W/m^2].
  Real2D Frad_SW_up;                       // Momentum-level radiative flux [W/m^2].
  Real2D Frad_LW_up;                       // Momentum-level radiative flux [W/m^2].
  Real2D Frad_SW_down;                     // Momentum-level radiative flux [W/m^2].
  Real2D Frad_LW_down;                     // Momentum-level radiative flux [W/m^2].
  Real2D radht_LW;                         // Thermodynamic heating [K/s].
  Real2D radht_SW;                         // Thermodynamic heating [K/s].
  Real2D Frad_SW;                          // Momentum-level radiative flux [W/m^2].
  Real2D Frad_LW;                          // Momentum-level radiative flux [W/m^2].
  RadiationDouble2D fdswcl;                // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D fuswcl;                // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D fdlwcl;                // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D fulwcl;                // Extended-grid Native double-precision diagnostic.
  int extended_atmos_bottom_level = -1;    // Lowest zero-based source level; -1 before setup.
  int extended_atmos_top_level    = -1;    // Highest zero-based source level; -1 before setup.
  int extended_atmos_range_size{};         // Number of extension/buffer levels.
  int lin_int_buffer{};                    // Number of extension/buffer levels.
  RadiationDouble2D T_in_K;                // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D rcil;                  // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D o3l;                   // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D rsm_rad;               // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D rcm_in_cloud_rad;      // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D cloud_frac_rad;        // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D ice_supersat_frac_rad; // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D radht_SW_rad;          // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D radht_LW_rad;          // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D p_in_mb;               // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D sp_humidity;           // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D Frad_uLW;              // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D Frad_dLW;              // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D Frad_uSW;              // Extended-grid Native double-precision diagnostic.
  RadiationDouble2D Frad_dSW;              // Extended-grid Native double-precision diagnostic.
};

void setup_radiation_variables(int ngrdcol, int nzm, int nzt,                   // In
                               radiation_variables_state& radiation_variables); // Native module Out
void reset_radiation_variables(radiation_variables_state& radiation_variables);
void setup_bugsrad_variables(
    ScratchWorkspace& workspace,                            // InOut
    int ngrdcol, int nzmax, Real1D zt_grid, Real1D zm_grid, // In
    Real1D zm_grid_spacing, Real1D p_in_Pa_zm,              // In
    radiation_variables_state& radiation_variables,         // Native module Out
    const parameters_radiation_state& parameters_radiation, // Native module In
    extended_atmosphere_state& extended_atmosphere          // Native module InOut
);
void reset_bugsrad_variables(radiation_variables_state& radiation_variables);
void cleanup_radiation_variables(radiation_variables_state& radiation_variables);
void cleanup_bugsrad_variables(radiation_variables_state& radiation_variables);

template <class Zt, class Zm, class Spacing, class Pressure>
void setup_bugsrad_variables(
    ScratchWorkspace& workspace,                            // InOut
    int ngrdcol, int nzmax, Zt zt_grid, Zm zm_grid,         // In
    Spacing zm_grid_spacing, Pressure p_in_Pa_zm,           // In
    radiation_variables_state& radiation_variables,         // Native module Out
    const parameters_radiation_state& parameters_radiation, // Native module In
    extended_atmosphere_state& extended_atmosphere          // Native module InOut
) {
  // The Native profile arguments may be strided slices of a column array.
  auto workspace_scope = workspace.scope();
  setup_bugsrad_variables(workspace, ngrdcol, nzmax, to_execution_space(workspace, zt_grid),
                          to_execution_space(workspace, zm_grid), // In
                          to_execution_space(workspace, zm_grid_spacing),
                          to_execution_space(workspace, p_in_Pa_zm), // In
                          radiation_variables, parameters_radiation, extended_atmosphere);
}

} // namespace clubb_kokkos
