// Supply subsidence, tendencies and surface fluxes for MPACE-B.
//
// mpace_b_tndcy computes pressure-dependent vertical motion [m/s] and prescribed
// temperature/water advection [K/s,(kg/kg)/s] with selected scalar copies.
// mpace_b_sfclyr interpolates supplied heat fluxes [W/m^2] at time [s] and converts
// them with each surface density [kg/m^3], returning heat [K m/s], water[(kg/kg) m/s]
// and ustar [m/s]. Model fields use column-first Views. The appended surface
// context replaces Native time-series module state; all scratch is call-owned.
// References: Michael Falk, July–November 2006; Liou, Wallace and Hobbs, Shettle and Weinman;
// http://science.arm.gov/wg/cpm/scm/scmic5/index.html
#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void mpace_b_tndcy(int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                          const grid& gr, Real2D p_in_Pa, Real2D thvm,
                          Real2D wm_zt, Real2D wm_zm, Real2D thlm_forcing, Real2D rtm_forcing,
                          Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Compute MPACE-B large-scale subsidence and temperature/water tendencies
  // from pressure [Pa] and virtual potential temperature [K], then copy scalars.
  constexpr real_t D     = real_t(5.8e-6);   // Divergence[1/s].
  constexpr real_t p_sfc = real_t(101000.0); // Prescribed surface pressure [Pa].
  constexpr real_t pinv  = real_t(85000.0);  // Inversion pressure [Pa].
  const int nzt          = gr.nzt;
  const int nzm          = gr.nzm;
  using Policy           = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Compute vertical motion.
  Kokkos::parallel_for(
      "mpace_b_tndcy::motion", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i                 = ((index % ngrdcol));
        const int k                 = ((index / ngrdcol));
        const real_t velocity_omega = Kokkos::min(D * (p_sfc - p_in_Pa(i, k)), D * (p_sfc - pinv));
        wm_zt(i, k)                 = -velocity_omega * Rd * thvm(i, k) / p_in_Pa(i, k) / grav;
      });

  // Interpolate to momentum levels, then enforce surface/top boundaries.
  zt2zm_api(nzm, nzt, ngrdcol, gr, wm_zt, wm_zm);
  Kokkos::parallel_for(
      "mpace_b_tndcy::boundaries", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wm_zm(i, 0)       = real_t(0.0);
        wm_zm(i, ((nzm - 1))) = real_t(0.0);
      });

  // Compute large-scale tendencies. Temperature begins in K/day; water begins
  // in g/(kg day), with each conversion retained in its Native order.
  Kokkos::parallel_for(
      "mpace_b_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i             = ((index % ngrdcol));
        const int k             = ((index / ngrdcol));
        const real_t t_tendency = Kokkos::min(
            -real_t(4.0),
            -real_t(15.0) * (real_t(1.0) - ((p_sfc - p_in_Pa(i, k)) / real_t(21818.0))));
        thlm_forcing(i, k) =
            (t_tendency * Kokkos::pow(p_sfc / p_in_Pa(i, k), Rd / Cp)) / sec_per_day;
        rtm_forcing(i, k) =
            Kokkos::min(real_t(0.164), -real_t(3.0) * (real_t(1.0) - ((p_sfc - p_in_Pa(i, k)) /
                                                                      real_t(15171.0)))) /
            g_per_kg / sec_per_day;
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "mpace_b_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          // Test scalars with theta_l and rt if desired.
          if (sclr_idx.iisclr_thl > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iisclr_rt > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "mpace_b_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//----------------------------------------------------------------------
inline void mpace_b_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D rho_sfc,
                           Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar,
                           const time_dependent_input_state& input) {
  // Description:
  // Interpolate the MPACE-B sensible/latent heat input [W/m^2], convert it to
  // kinematic fluxes with local density, and prescribe ustar = 0.25 m/s.
  auto workspace_scope = workspace.scope();
  real_t time_frac;
  int before_time, after_time;
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);

  // Get sensible and latent heat from the input, in Native call order.
  const auto sensible = to_host_space(workspace, input.sens_ht_given);
  const real_t sensible_heat_flx =
      linear_interp_factor(time_frac, sensible(after_time), sensible(before_time));
  const auto latent = to_host_space(workspace, input.latent_ht_given);
  const real_t latent_heat_flx =
      linear_interp_factor(time_frac, latent(after_time), latent(before_time));
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "mpace_b_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Declare ustar, then compute heat and moisture fluxes.
        ustar(i)      = real_t(0.25);
        wpthlp_sfc(i) = convert_sens_ht_to_km_s(sensible_heat_flx, rho_sfc(i));
        wprtp_sfc(i)  = convert_latent_ht_to_m_s(latent_heat_flx, rho_sfc(i));
      });
}

} // namespace clubb_kokkos
