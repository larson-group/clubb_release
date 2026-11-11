// Compute bulk surface heat and moisture exchange for the January 2006 TWP-ICE case.
//
// twp_ice_sfclyr interpolates surface temperature [K] at time [s], scales
// RICO 20-m exchange coefficients to the supplied height [m], and supplies heat
// [K m/s], water[(kg/kg) m/s], ustar [m/s] and surface temperature [K]. Inputs
// also include theta_l [K], total water [kg/kg], wind [m/s] and pressure [Pa].
// The appended context owns the surface series replacing Native module state.
// Supply surface Exner [-] explicitly. Reference: GCSS TWP-ICE specification;
// http://users.monash.edu.au/~ladavies/gcss.html; Stevens et al. (2000), eq. 3.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/saturation.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void twp_ice_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D z,
                           Real1D exner_sfc, Real1D thlm_sfc, Real1D ubar, Real1D rtm, Real1D p_sfc,
                           int saturation_formula, Real1D wpthlp_sfc, Real1D wprtp_sfc,
                           Real1D ustar, Real1D T_sfc, const time_dependent_input_state& input) {
  // Description:
  // Set surface exchange using time-dependent surface temperature, RICO bulk
  // coefficients and ATEX roughness scaling to the lowest model level.
  auto workspace_scope               = workspace.scope();
  constexpr real_t C_h_20            = real_t(0.001094); // RICO heat coefficient.
  constexpr real_t C_q_20            = real_t(0.001133); // RICO moisture coefficient.
  constexpr real_t z0                = real_t(0.00015);  // ATEX roughness length [m].
  constexpr real_t standard_flux_alt = real_t(20.0);     // Reference exchange height [m].
  real_t time_frac;
  int before_time, after_time;

  // Interpolate surface temperature from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "twp_ice_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        T_sfc(i) = T_sfc_interp;
        // Declare the prescribed value of ustar.
        ustar(i) = real_t(0.3);

        // Stevens et al. (2000), eq. 3: modify exchange for the lowest model height.
        const real_t Ch = C_h_20 * (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z(i) / z0)) *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z(i) / z0));
        const real_t Cq = C_q_20 * (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z(i) / z0)) *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z(i) / z0));
        const real_t rsat = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);
        // Consume the column-local coefficients without intermediate fields.
        wpthlp_sfc(i) = -Ch * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
        wprtp_sfc(i)  = -Cq * ubar(i) * (rtm(i) - rsat);
      });
}

} // namespace clubb_kokkos
