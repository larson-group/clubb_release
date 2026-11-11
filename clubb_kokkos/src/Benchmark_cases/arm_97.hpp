// Compute surface fluxes for the ARM IOP A case of July 26–30, 1997.
//
// arm_97_sfclyr interpolates the case heat fluxes when l_t_dependent is true;
// otherwise its output Views are left unchanged, matching Native. Supply time
// [s], height [m], density [kg/m^3], theta_l [K] and wind speed [m/s]. Outputs
// are heat flux [K m/s], water flux [(kg/kg) m/s] and ustar [m/s]. The appended
// input context replaces Native module access, without persistent local state.
// Reference: Native GCSS working group 4 case 3 specification.
#pragma once

#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void arm_97_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D z,
                          Real1D rho_sfc, Real1D thlm_sfc, Real1D ubar, Real1D wpthlp_sfc,
                          Real1D wprtp_sfc, Real1D ustar, const time_dependent_input_state& input) {
  // Description:
  // Compute surface momentum, heat and moisture fluxes according to the ARM
  // specification, conditional on time-dependent input being enabled.
  auto workspace_scope = workspace.scope();
  constexpr real_t z0 = real_t(0.035); // ARM cumulus momentum roughness height [m].
  if (input.l_t_dependent) {
    real_t time_frac = -real_t(1.0);
    int before_time, after_time;
    time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
                input.time_sfc_given, before_time, after_time, time_frac);
    const auto sens_ht   = to_host_space(workspace, input.sens_ht_given);
    const auto latent_ht = to_host_space(workspace, input.latent_ht_given);
    const real_t heat_flx =
        linear_interp_factor(time_frac, sens_ht(after_time), sens_ht(before_time));
    const real_t moisture_flx =
        linear_interp_factor(time_frac, latent_ht(after_time), latent_ht(before_time));

    using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
    Kokkos::parallel_for(
        "arm_97_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          // Convert W/m^2 into kinematic heat/water flux units.
          wpthlp_sfc(i) = convert_sens_ht_to_km_s(heat_flx, rho_sfc(i));
          wprtp_sfc(i)  = convert_latent_ht_to_m_s(moisture_flx, rho_sfc(i));

          // Compute momentum fluxes using ARM cumulus formulae.
          const real_t bflx = grav / thlm_sfc(i) * wpthlp_sfc(i);
          ustar(i)          = diag_ustar(z(i), bflx, ubar(i), z0);
        });
  }
}

} // namespace clubb_kokkos
