// Set subsidence and surface exchange for the ASTEX Khairoutdinov–Kogan case.
//
// astex_a209_tndcy supplies one profile of subsidence [m/s], zero theta_l/water
// forcing and selected scalar copies. astex_a209_sfclyr interpolates temperature
// from 41 supplied times, evaluates saturation and scales bulk coefficients to
// the lowest model height. Outputs are heat flux [K m/s], water flux
// [(kg/kg) m/s], ustar [m/s] and temperature [K], one per column. The final input
// context replaces Native's implicit surface time-series module storage.
// Reference: ASTEX Lagrangian specification,
// http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/saturation.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void astex_a209_tndcy(int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                             const grid& gr, Real1D wm_zt, Real1D wm_zm,
                             Real1D thlm_forcing, Real1D rtm_forcing,
                             Real2D sclrm_forcing, Real2D edsclrm_forcing) {
  // Description:
  // Set theta/water tendencies and linear large-scale subsidence for ASTEX KK.
  // These are single profiles on the first grid column, as in Native.
  (void)sclr_dim;
  (void)edsclr_dim;
  const int nzt = gr.nzt;
  const int nzm = gr.nzm;
  const auto zt = gr.zt;
  using Policy  = Kokkos::RangePolicy<typename Real1D::execution_space>;

  // Compute large-scale subsidence.
  Kokkos::parallel_for(
      "astex_a209_tndcy::subsidence", Policy(0, nzt),
      KOKKOS_LAMBDA(int i) { wm_zt(i) = -real_t(5.e-6) * zt(0, i); });

  // Interpolate to momentum levels.
  Kokkos::deep_copy(wm_zm, zt2zm_api(gr, wm_zt));

  // Boundary conditions on zm: at the surface and model top.
  Kokkos::parallel_for(
      "astex_a209_tndcy::boundaries", Policy(0, 1), KOKKOS_LAMBDA(int) {
        wm_zm(0)       = real_t(0.0);
        wm_zm((nzm - 1)) = real_t(0.0);
      });

  // Radiative theta-l tendency, then large-scale advective moisture tendency.
  Kokkos::deep_copy(thlm_forcing, real_t(0.0));
  Kokkos::deep_copy(rtm_forcing, real_t(0.0));

  // Test scalars with theta_l and rt if desired.
  if (sclr_idx.iisclr_thl > -1)
    Kokkos::parallel_for(
        "astex_a209_tndcy::thl", Policy(0, nzt),
        KOKKOS_LAMBDA(int i) { sclrm_forcing(i, sclr_idx.iisclr_thl) = thlm_forcing(i); });
  if (sclr_idx.iisclr_rt > -1)
    Kokkos::parallel_for(
        "astex_a209_tndcy::rt", Policy(0, nzt),
        KOKKOS_LAMBDA(int i) { sclrm_forcing(i, sclr_idx.iisclr_rt) = rtm_forcing(i); });
  if (sclr_idx.iiedsclr_thl > -1)
    Kokkos::parallel_for(
        "astex_a209_tndcy::eddy_thl", Policy(0, nzt),
        KOKKOS_LAMBDA(int i) { edsclrm_forcing(i, sclr_idx.iiedsclr_thl) = thlm_forcing(i); });
  if (sclr_idx.iiedsclr_rt > -1)
    Kokkos::parallel_for(
        "astex_a209_tndcy::eddy_rt", Policy(0, nzt),
        KOKKOS_LAMBDA(int i) { edsclrm_forcing(i, sclr_idx.iiedsclr_rt) = rtm_forcing(i); });
}

//----------------------------------------------------------------------
inline void astex_a209_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D ubar,
                              Real1D rtm, Real1D thlm, Real1D lowestlevel, Real1D exner_sfc,
                              Real1D p_sfc, int saturation_formula, Real1D wpthlp_sfc,
                              Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc,
                              const time_dependent_input_state& input) {
  // Description:
  // Compute ASTEX surface fluxes with the Khairoutdinov–Kogan alteration.
  // Scale RICO 20-m drag coefficients to the lowest height with ATEX roughness.
  auto workspace_scope               = workspace.scope();
  constexpr int ntimes               = 41;
  constexpr real_t C_h_20            = real_t(0.001094); // RICO heat exchange coefficient.
  constexpr real_t C_q_20            = real_t(0.001133); // RICO water exchange coefficient.
  constexpr real_t z0                = real_t(0.00015);  // ATEX roughness length [m].
  constexpr real_t standard_flux_alt = real_t(20.0);     // Reference exchange height [m].
  real_t time_frac;
  int before_time, after_time;

  // Find time indices and the fraction for linear interpolation.
  time_select(time, ntimes, input.time_sfc_given, before_time, after_time, time_frac);

  // Interpolate the value for surface temperature.
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "astex_a209_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        T_sfc(i) = T_sfc_interp;

        // Set ustar as in RICO.
        ustar(i) = real_t(0.155);
        const real_t rsat = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);

        // Stevens et al. (2000), eq. 3: adjust the reference exchange height using
        // ATEX roughness. Retain both ratios and their multiplication order.
        const real_t Ch = C_h_20 *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowestlevel(i) / z0)) *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowestlevel(i) / z0));
        const real_t Cq = C_q_20 *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowestlevel(i) / z0)) *
                          (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowestlevel(i) / z0));
        // Consume the column-local coefficients without intermediate fields.
        wpthlp_sfc(i) = -Ch * ubar(i) * (thlm(i) - T_sfc(i) / exner_sfc(i));
        wprtp_sfc(i)  = -Cq * ubar(i) * (rtm(i) - rsat);
      });

  // Momentum fluxes are computed elsewhere.
}

} // namespace clubb_kokkos
