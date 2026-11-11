// Integrate hydrostatic pressure/density and inverse sounding altitude.
//
// hydrostatic fills pressure [Pa], Exner[-] and density [kg/m^3] on the two
// staggered grids from virtual potential temperature [K] and surface pressure.
// inverse_hydrostatic finds sounding altitudes [m] from surface pressure [Pa],
// surface altitude [m], virtual potential temperature [K] and Exner[-].
//
// Both retain Native argument order and use column-first Kokkos Views. The
// forward calculation reuses the permanent pressure and grid operators. The
// inverse preserves the ordered sounding recurrence and surface bracketing,
// with a fatal diagnostic when the entire sounding lies below the surface.
// Scratch uses the caller workspace; no persistent case state or Native science is retained.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"

namespace clubb_kokkos {

void hydrostatic(ScratchWorkspace& workspace,                            // InOut
                 int ngrdcol, const grid& gr, Real2D thvm, Real1D p_sfc, // In
                 Real2D p_in_Pa, Real2D p_in_Pa_zm,                      // Out
                 Real2D exner, Real2D exner_zm,                          // Out
                 Real2D rho, Real2D rho_zm);                             // Out

void inverse_hydrostatic(ScratchWorkspace& workspace,               // InOut
                         real_t p_sfc, real_t zm_init, int nlevels, // In
                         Real1D thvm, Real1D exner,                 // In
                         Real1D z);                                 // Out

// Host/strided callers cross the existing array owner once at this API boundary.
template <class Thvm, class Surface, class Pressure, class MomentumPressure, class Exner,
          class MomentumExner, class Density, class MomentumDensity>
void hydrostatic(ScratchWorkspace& workspace,                           // InOut
                 int ngrdcol, const grid& gr, Thvm thvm, Surface p_sfc, // In
                 Pressure p_in_Pa, MomentumPressure p_in_Pa_zm,         // Out
                 Exner exner, MomentumExner exner_zm,                   // Out
                 Density rho, MomentumDensity rho_zm)                   // Out
{
  auto workspace_scope = workspace.scope();
  auto thvm_exec       = to_execution_space(workspace, thvm);
  auto p_sfc_exec      = to_execution_space(workspace, p_sfc);
  auto p_exec          = to_execution_space(workspace, p_in_Pa, false);
  auto p_zm_exec       = to_execution_space(workspace, p_in_Pa_zm, false);
  auto exner_exec      = to_execution_space(workspace, exner, false);
  auto exner_zm_exec   = to_execution_space(workspace, exner_zm, false);
  auto rho_exec        = to_execution_space(workspace, rho, false);
  auto rho_zm_exec     = to_execution_space(workspace, rho_zm, false);
  hydrostatic(workspace, ngrdcol, gr, thvm_exec, p_sfc_exec, // In
              p_exec, p_zm_exec,                             // Out
              exner_exec, exner_zm_exec,                     // Out
              rho_exec, rho_zm_exec);                        // Out
  from_execution_space(workspace, p_in_Pa, p_exec);
  from_execution_space(workspace, p_in_Pa_zm, p_zm_exec);
  from_execution_space(workspace, exner, exner_exec);
  from_execution_space(workspace, exner_zm, exner_zm_exec);
  from_execution_space(workspace, rho, rho_exec);
  from_execution_space(workspace, rho_zm, rho_zm_exec);
}

template <class Thvm, class Exner, class Height>
void inverse_hydrostatic(ScratchWorkspace& workspace,               // InOut
                         real_t p_sfc, real_t zm_init, int nlevels, // In
                         Thvm thvm, Exner exner,                    // In
                         Height z)                                  // Out
{
  auto workspace_scope = workspace.scope();
  auto thvm_exec       = to_execution_space(workspace, thvm);
  auto exner_exec      = to_execution_space(workspace, exner);
  auto z_exec          = to_execution_space(workspace, z, false);
  inverse_hydrostatic(workspace, p_sfc, zm_init, nlevels, thvm_exec, exner_exec, // In
                      z_exec);                                                   // Out
  from_execution_space(workspace, z, z_exec);
}

} // namespace clubb_kokkos
