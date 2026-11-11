// Apply large-scale and surface forcing for the GABLS2 case.
//
// gabls2_tndcy supplies subsidence after 26 h with zero theta_l/water advection
// and selected scalar copies. gabls2_sfclyr computes the piecewise diurnal sea-
// surface temperature, bulk fluxes and similarity friction velocity. Time is [s],
// height [m], pressure [Pa], theta_l/temperature [K], water [kg/kg], wind [m/s] and
// Exner [-]. Outputs are heat [K m/s], water[(kg/kg) m/s] and ustar [m/s]. Model
// fields use (column, level[,scalar]) Views; no persistent state is stored.
// Reference: GABLS2 specification, http://people.su.se/~gsven/gabls/
// Original case implementation: Michael Falk, 29 December 2006.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/saturation.hpp"
#include "diag_ustar_module.hpp"
#include "sfc_flux.hpp"

namespace clubb_kokkos {

//-------------------------------------------------------------------------------
inline void gabls2_tndcy(int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                         const grid& gr, double time, double time_initial,
                         Real2D wm_zt, Real2D wm_zm, Real2D thlm_forcing, Real2D rtm_forcing,
                         Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Apply GABLS2 large-scale motion after 26 h, zero heat/water advection and
  // copy selected scalar tendencies. Michael Falk, 29 December 2006.
  const int nzt = gr.nzt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;

  const int nzm = gr.nzm;
  const auto zt = gr.zt;

  // Compute vertical motion after 93600 s (26 h) of simulation time.
  if (time > time_initial + double(93600.0)) {
    Kokkos::parallel_for(
        "gabls2_tndcy::subsidence", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (zt(i, k) <= real_t(1000.0))
            wm_zt(i, k) = -real_t(0.005) * (zt(i, k) / real_t(1000.0));
          else
            wm_zt(i, k) = -real_t(0.005);
        });
  } else {
    Kokkos::deep_copy(wm_zt, real_t(0.0));
  }
  zt2zm_api(nzm, nzt, ngrdcol, gr, wm_zt, wm_zm);

  // Boundary conditions on vertical motion: surface and model top.
  Kokkos::parallel_for(
      "gabls2_tndcy::boundaries", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wm_zm(i, 0)       = real_t(0.0);
        wm_zm(i, (nzm - 1)) = real_t(0.0);
      });

  Kokkos::parallel_for(
      "gabls2_tndcy::profiles", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i = (index % ngrdcol);
        const int k = (index / ngrdcol);
        // Large-scale temperature and advective moisture tendencies.
        thlm_forcing(i, k) = real_t(0.0);
        rtm_forcing(i, k)  = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "gabls2_tndcy::passive", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          // Test scalars with theta_l and rt if desired.
          if (sclr_idx.iisclr_thl > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iisclr_rt > -1)
            sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        });
  }
  if (edsclr_dim > 0) {
    Kokkos::parallel_for(
        "gabls2_tndcy::eddy", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//-------------------------------------------------------------------------------
inline void gabls2_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time,
                          double time_initial, Real1D lowest_level, Real1D p_sfc, Real1D ubar,
                          Real1D thlm, Real1D rtm, Real1D exner_sfc, int saturation_formula,
                          Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc) {
  // Description:
  // Surface forcing for GABLS2: use the specified piecewise temperature in
  // local hours, ATEX height scaling, and 2.5% of potential latent heat flux.
  auto workspace_scope               = workspace.scope();
  constexpr real_t standard_flux_alt = real_t(10.0);   // Exchange reference height [m].
  constexpr real_t C_10              = real_t(0.0013); // ATEX exchange coefficient.
  constexpr real_t z0                = real_t(0.03);   // GABLS2 roughness length [m].
  Real1D Cz                          = workspace_scope.view<Real1D>("gabls2_sfclyr::Cz", ngrdcol);
  Real1D rsat                        = workspace_scope.view<Real1D>("gabls2_sfclyr::rsat", ngrdcol);

  // The experiment starts at 14 local (19 UTC). Cast only after the double time
  // expression, matching Native time_precision and core_rknd conversion order.
  const double time_in_hours_init = 14.0;
  const real_t time_in_hours      = static_cast<real_t>(
      (time - time_initial) / static_cast<double>(sec_per_hr) + time_in_hours_init);
  real_t T_sfc_calc;
  // Compute sea-surface temperature in Celsius with the specified constants.
  if (time_in_hours <= real_t(17.4)) {
    T_sfc_calc =
        -real_t(10.0) - real_t(25.0) * Kokkos::cos(time_in_hours * real_t(0.22) + real_t(0.2));
  } else if (time_in_hours <= real_t(30.0)) {
    T_sfc_calc = (-real_t(0.54) * time_in_hours) + real_t(15.2);
  } else if (time_in_hours <= real_t(41.9)) {
    T_sfc_calc =
        -real_t(7.0) - real_t(25.0) * Kokkos::cos(time_in_hours * real_t(0.21) + real_t(1.8));
  } else if (time_in_hours <= real_t(53.3)) {
    T_sfc_calc = (-real_t(0.37) * time_in_hours) + real_t(18.0);
  } else if (time_in_hours <= real_t(65.6)) {
    T_sfc_calc =
        -real_t(4.0) - real_t(25.0) * Kokkos::cos(time_in_hours * real_t(0.22) + real_t(2.5));
  } else {
    T_sfc_calc = real_t(4.4);
  }

  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "gabls2_sfclyr::temperature", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Scale to the lowest level when it differs from 10 m (Stevens et al., 2000, eq. 3).
        Cz(i) = C_10 * (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowest_level(i) / z0)) *
                (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(lowest_level(i) / z0));
        T_sfc(i) = T_sfc_calc + real_t(273.15);
        rsat(i)  = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);
      });

  // Compute heat and moisture fluxes.
  compute_wpthlp_sfc(ngrdcol, Cz, ubar, thlm, T_sfc, exner_sfc, wpthlp_sfc);
  compute_wprtp_sfc(ngrdcol, Cz, ubar, rtm, rsat, wprtp_sfc);
  Kokkos::parallel_for(
      "gabls2_sfclyr::momentum", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        // Latent heat flux is 2.5% of its potential value.
        wprtp_sfc(i) = wprtp_sfc(i) * real_t(0.025);

        // Compute momentum fluxes.
        const real_t sstheta = T_sfc(i) * Kokkos::pow(p0 / p_sfc(i), Rd / Cp);
        const real_t bflx    = wpthlp_sfc(i) * grav / sstheta;
        ustar(i)             = diag_ustar(lowest_level(i), bflx, ubar(i), z0);
      });
}

} // namespace clubb_kokkos
