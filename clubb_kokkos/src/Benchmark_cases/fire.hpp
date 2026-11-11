// Compute aerodynamic surface fluxes for the GCSS FIRE case.
//
// fire_sfclyr interpolates sea-surface temperature from the supplied case data,
// computes saturation and bulk heat/water exchange, and prescribes ustar=0.3 m/s.
// Inputs include time [s], wind speed [m/s], pressure [Pa], near-surface theta_l
// [K], total water [kg/kg], Exner [-] and the saturation formula selector. Outputs
// are heat [K m/s], water [(kg/kg) m/s], ustar [m/s] and surface temperature [K].
// The appended context replaces implicit Native time-series module storage.
// Reference: Moeng et al. (1996), Bull. Amer. Meteor. Soc. 77, 261–278.
#pragma once

#include "scratch_workspace.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"
#include "CLUBB_core/saturation.hpp"

namespace clubb_kokkos {

//======================================================================
inline void fire_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D ubar,
                        Real1D p_sfc, Real1D thlm_sfc, Real1D rtm_sfc, Real1D exner_sfc,
                        int saturation_formula, Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar,
                        Real1D T_sfc, const time_dependent_input_state& input) {
  // Description:
  // Compute aerodynamic surface heat/moisture fluxes with the prescribed FIRE
  // exchange coefficient and interpolated surface temperature.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;

  // Interpolate variables from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "fire_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        T_sfc(i) = T_sfc_interp;
        const real_t Cz   = real_t(0.0013);
        ustar(i) = real_t(0.3);
        const real_t rsat = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);
        // Consume the column-local coefficients without intermediate fields.
        wpthlp_sfc(i) = -Cz * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
        wprtp_sfc(i)  = -Cz * ubar(i) * (rtm_sfc(i) - rsat);
      });
}

} // namespace clubb_kokkos
