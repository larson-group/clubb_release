// Compute surface fluxes for the GCSS ARM shallow-cumulus case.
//
// arm_sfclyr uses seven supplied heat-flux samples and a fixed density of
// 1.1 kg/m^3, then diagnoses friction velocity from surface similarity. Supply
// time [s], height [m], theta_l [K] and wind speed [m/s]; outputs are heat flux
// [K m/s], water flux [(kg/kg) m/s] and ustar [m/s]. The final input context owns
// the surface time series that Native reads from implicit module storage.
// Reference: Brown et al. (2002), Q. J. R. Meteorol. Soc. 128, 1075–1093.
#pragma once

#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void arm_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D z,
                       Real1D thlm_sfc, Real1D ubar, Real1D wpthlp_sfc, Real1D wprtp_sfc,
                       Real1D ustar, const time_dependent_input_state& input) {
  // Description:
  // Compute surface fluxes of horizontal momentum, heat and moisture according
  // to GCSS ARM specifications. Heat and moisture are supplied in W/m^2.
  auto workspace_scope = workspace.scope();
  constexpr real_t z0  = real_t(0.035); // ARM cumulus momentum roughness height [m].
  constexpr int ntimes = 7;
  real_t heat_flx, moisture_flx;

  // Compute heat and moisture fluxes from ARM data in W/m^2.
  compute_ht_mostr_flux(workspace, time, ntimes, input.latent_ht_given, input.sens_ht_given,
                        input.time_sfc_given, heat_flx, moisture_flx);
  const real_t rho_sfc = real_t(1.1); // Native prescribed surface density [kg/m^3].

  // Convert heat_flx and moisture_flx to natural units.
  const real_t heat_flx2     = convert_sens_ht_to_km_s(heat_flx, rho_sfc);
  const real_t moisture_flx2 = convert_latent_ht_to_m_s(moisture_flx, rho_sfc);

  // Compute momentum fluxes.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "arm_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Buoyancy flux [m^2/s^3] is needed by diag_ustar.
        const real_t bflx = grav / thlm_sfc(i) * heat_flx2;
        ustar(i)          = diag_ustar(z(i), bflx, ubar(i), z0);

        // Assign fluxes.
        wpthlp_sfc(i) = heat_flx2;
        wprtp_sfc(i)  = moisture_flx2;
      });
}

} // namespace clubb_kokkos
