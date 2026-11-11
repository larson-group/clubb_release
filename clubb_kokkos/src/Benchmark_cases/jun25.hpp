// Interpolate supplied heat fluxes for the June 25 altocumulus case.
//
// jun25_altocu_read_t_dependent returns sensible and latent heat fluxes [W/m^2]
// at time [s] from the case-owned surface series. The bottom model level is
// above the surface, so this file does not provide a surface-layer routine.
// The final context replaces the Native module's implicit time-series state.
#pragma once

#include "time_dependent_input.hpp"
#include "CLUBB_core/interpolation.hpp"

namespace clubb_kokkos {

//-----------------------------------------------------------------------
inline void jun25_altocu_read_t_dependent(ScratchWorkspace& workspace, double time, real_t& sens_ht,
                                          real_t& latent_ht,
                                          const time_dependent_input_state& input) {
  // Description:
  // Interpolate the case's sensible/latent heat input at the requested time
  // and return both scalar fluxes in W/m^2. Bracketing errors follow time_select.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;

  // Interpolate the values read from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto sensible = to_host_space(workspace, input.sens_ht_given);
  const auto latent   = to_host_space(workspace, input.latent_ht_given);
  sens_ht   = linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
  latent_ht = linear_interp_factor(time_frac, latent(after_time), latent(before_time));
}

} // namespace clubb_kokkos
