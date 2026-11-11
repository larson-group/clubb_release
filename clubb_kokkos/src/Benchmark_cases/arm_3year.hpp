// Compute horizontal-momentum, heat and moisture surface fluxes for ARM three-year.
//
// arm_3year_sfclyr interpolates the supplied heat fluxes, converts them using
// each column's surface density, and diagnoses friction velocity. Supply time
// [s], height [m], density [kg/m^3], theta_l [K] and wind speed [m/s]. Returned
// heat/water fluxes have units [K m/s]/[(kg/kg) m/s]; ustar is [m/s]. The appended
// case context replaces Native's implicit surface time-series module state.
// Reference: None specified in arm_3year.F90.
#pragma once

#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void arm_3year_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D z,
                             Real1D rho_sfc, Real1D thlm_sfc, Real1D ubar, Real1D wpthlp_sfc,
                             Real1D wprtp_sfc, Real1D ustar,
                             const time_dependent_input_state& input) {
  // Description:
  // Compute surface momentum, heat and moisture fluxes according to the ARM
  // specification, using the complete supplied surface time series.
  auto workspace_scope = workspace.scope();
  constexpr real_t z0 = real_t(0.035); // ARM cumulus momentum roughness height [m].
  real_t heat_flx, moisture_flx;       // Interpolated surface heat fluxes [W/m^2].

  // Compute heat and moisture fluxes from ARM data in W/m^2.
  compute_ht_mostr_flux(workspace, time, static_cast<int>(input.time_sfc_given.extent(0)),
                        input.latent_ht_given, input.sens_ht_given, input.time_sfc_given, heat_flx,
                        moisture_flx);

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "arm_3year_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Convert W/m^2 into the kinematic heat and water flux units.
        wpthlp_sfc(i) = convert_sens_ht_to_km_s(heat_flx, rho_sfc(i));
        wprtp_sfc(i)  = convert_latent_ht_to_m_s(moisture_flx, rho_sfc(i));

        // Compute momentum fluxes using ARM cumulus formulae.
        const real_t bflx = grav / thlm_sfc(i) * wpthlp_sfc(i);
        ustar(i)          = diag_ustar(z(i), bflx, ubar(i), z0);
      });
}

} // namespace clubb_kokkos
