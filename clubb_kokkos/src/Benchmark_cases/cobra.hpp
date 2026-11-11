// Compute heat, water and carbon-dioxide surface exchange for the COBRA case.
//
// cobra_sfclyr interpolates 49 surface samples at time [s] and converts supplied
// sensible/latent heat [W/m^2] and CO2[umol/(m^2 s)] using local air density
// [kg/m^3]. It returns heat [K m/s], water[(kg/kg) m/s], CO2[ppmv m/s] scalar
// fluxes, temperature [K] and friction velocity [m/s]. Height [m], theta_l [K] and
// wind [m/s] enter similarity. The appended case context owns all time series.
// Selected passive/eddy scalar copies follow the Native explicit-diffusion
// guards; unspecified scalar fields retain their supplied values.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//-----------------------------------------------------------------------
inline void cobra_sfclyr(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                         const sclr_idx_type& sclr_idx, bool l_do_expldiff_rtm_thlm, double time,
                         Real1D z, Real1D rho_sfc, Real1D thlm_sfc, Real1D ubar, Real1D wpthlp_sfc,
                         Real1D wprtp_sfc, Real1D ustar, Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc,
                         Real1D T_sfc, const time_dependent_input_state& input) {
  // Description:
  // Compute ARM-format surface momentum/heat/moisture fluxes with data altered
  // for COBRA CO2, preserving scalar copies and the no-double-forcing guards.
  auto workspace_scope  = workspace.scope();
  constexpr real_t M_da = real_t(0.02897); // Dry-air molecular weight [kg/mol].
  constexpr int ntimes  = 49;
  constexpr real_t z0   = real_t(1.75); // COBRA momentum roughness height [m].
  real_t time_frac;
  int before_time, after_time;

  // Compute heat/moisture fluxes from ARM-format data [W/m^2]. Identify the time
  // interval first, then interpolate heat, water, CO2 and temperature in order.
  time_select(time, ntimes, input.time_sfc_given, before_time, after_time, time_frac);
  const auto sensible = to_host_space(workspace, input.sens_ht_given);
  const real_t heat_flx =
      linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
  const auto latent = to_host_space(workspace, input.latent_ht_given);
  const real_t moisture_flx =
      linear_interp_factor(time_frac, latent(after_time), latent(before_time));
  const auto co2         = to_host_space(workspace, input.CO2_sfc_given);
  const real_t CO2_flx   = linear_interp_factor(time_frac, co2(after_time), co2(before_time));
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_calc =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "cobra_sfclyr::fluxes", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        T_sfc(i) = T_sfc_calc;

        // Convert heat_flx and moisture_flx to natural units.
        const real_t heat_flx2     = convert_sens_ht_to_km_s(heat_flx, rho_sfc(i));
        const real_t moisture_flx2 = convert_latent_ht_to_m_s(moisture_flx, rho_sfc(i));

        // CO2 is supplied in micromoles/(m^2 s); concentration is ppmv, or 10^6
        // times mol(CO2)/mol(dry air). w'CO2'=flux*M(dry air)/rho(dry air).
        const real_t CO2_flx2 = CO2_flx * (M_da / rho_sfc(i));

        // Buoyancy flux [m^2/s^3] is needed by diag_ustar.
        const real_t bflx = grav / thlm_sfc(i) * heat_flx2;
        ustar(i)          = diag_ustar(z(i), bflx, ubar(i), z0);

        // Assign fluxes.
        wpthlp_sfc(i) = heat_flx2;
        wprtp_sfc(i)  = moisture_flx2;
        if (sclr_dim > 0) {
          if (sclr_idx.iisclr_CO2 > -1)
            wpsclrp_sfc(i, sclr_idx.iisclr_CO2) = CO2_flx2;
          if (sclr_idx.iisclr_thl > -1)
            wpsclrp_sfc(i, sclr_idx.iisclr_thl) = wpthlp_sfc(i);
          if (sclr_idx.iisclr_rt > -1)
            wpsclrp_sfc(i, sclr_idx.iisclr_rt) = wprtp_sfc(i);
        }
        if (edsclr_dim > 0) {
          if (sclr_idx.iiedsclr_CO2 > -1)
            wpedsclrp_sfc(i, sclr_idx.iiedsclr_CO2) = CO2_flx2;

          // Explicit diffusion puts copies of theta_l/rt in these eddy scalars.
          // Their original fields already get the surface flux in advance_xm_wpxp;
          // do not apply those fluxes again to the copies.
          if (l_do_expldiff_rtm_thlm && edsclr_dim > 1) {
            if (sclr_idx.iiedsclr_thl > -1)
              wpedsclrp_sfc(i, sclr_idx.iiedsclr_thl) = zero;
            if (sclr_idx.iiedsclr_rt > -1)
              wpedsclrp_sfc(i, sclr_idx.iiedsclr_rt) = zero;
          } else {
            if (sclr_idx.iiedsclr_thl > -1)
              wpedsclrp_sfc(i, sclr_idx.iiedsclr_thl) = wpthlp_sfc(i);
            if (sclr_idx.iiedsclr_rt > -1)
              wpedsclrp_sfc(i, sclr_idx.iiedsclr_rt) = wprtp_sfc(i);
          }
        }
      });
}

} // namespace clubb_kokkos
