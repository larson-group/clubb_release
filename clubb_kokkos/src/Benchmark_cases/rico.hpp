// Prescribe large-scale tendencies and bulk surface fluxes for the RICO case.
//
// rico_tndcy supplies piecewise temperature/water advection, converts specific-
// humidity forcing to mixing ratio, and copies selected scalars. rico_sfclyr
// interpolates surface temperature and applies RICO exchange coefficients with
// ATEX height scaling. Fields use (column, level[,scalar]) LayoutLeft Views.
// Height [m], pressure [Pa], theta_l/temperature [K], water [kg/kg], wind [m/s] and
// Exner [-] determine tendencies [K/s,(kg/kg)/s], surface fluxes [K m/s,(kg/kg)m/s,
// (m/s)^2] and ustar [m/s]. The appended context supplies surface input series.
// References: Michael Falk, December 2006; RICO 3D LES and ATEX specifications,
// http://www.knmi.nl/samenw/rico/setup3d.html
// http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/saturation.hpp"
#include "sfc_flux.hpp"
#include "spec_hum_to_mixing_ratio.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
inline void rico_tndcy(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                       const sclr_idx_type& sclr_idx, const grid& gr, Real2D rtm, Real2D exner,
                       Real2D thlm_forcing, Real2D rtm_forcing, Real3D sclrm_forcing,
                       Real3D edsclrm_forcing) {
  // Description:
  // Apply the RICO short-composite large-scale temperature/radiation and
  // specific-humidity tendencies, converting them to CLUBB field units.
  auto workspace_scope = workspace.scope();
  const int nzt = gr.nzt;
  const auto zt = gr.zt;
  Real2D qtm_forcing   = workspace_scope.view<Real2D>("rico_tndcy::qtm_forcing", ngrdcol, nzt);
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Horizontal temperature advection and radiation(M. Falk, 15 December 2006).
  Kokkos::parallel_for(
      "rico_tndcy::temperature", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i = (index % ngrdcol);
        const int k = (index / ngrdcol);
        real_t t_tendency;
        if (zt(i, k) < real_t(4000.0)) {
          t_tendency = -real_t(2.51) / real_t(86400.0) + (-real_t(2.18) + real_t(2.51)) /
                                                             (real_t(86400.0) * real_t(4000.0)) *
                                                             zt(i, k);
        } else if (zt(i, k) < real_t(5000.0)) {
          t_tendency = -real_t(2.18) / real_t(86400.0) +
                       real_t(2.18) / (real_t(86400.0) * (real_t(5000.0) - real_t(4000.0))) *
                           (zt(i, k) - real_t(4000.0));
        } else {
          t_tendency = real_t(0.0);
        }
        // Convert absolute-temperature forcing to potential-temperature forcing.
        thlm_forcing(i, k) = t_tendency / exner(i, k);
      });

  // Large-scale horizontal moisture advection starts in g/(kg s).
  Kokkos::parallel_for(
      "rico_tndcy::water", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i = (index % ngrdcol);
        const int k = (index / ngrdcol);
        if (zt(i, k) < real_t(3000.0)) {
          qtm_forcing(i, k) =
              -real_t(1.0) / real_t(86400.0) +
              (real_t(0.345) + real_t(1.0)) / (real_t(86400.0) * real_t(3000.0)) * zt(i, k);
        } else if (zt(i, k) < real_t(4000.0)) {
          qtm_forcing(i, k) = real_t(0.345) / real_t(86400.0);
        } else if (zt(i, k) < real_t(5000.0)) {
          qtm_forcing(i, k) = real_t(0.345) / real_t(86400.0) +
                              (-real_t(0.345)) /
                                  (real_t(86400.0) * (real_t(5000.0) - real_t(4000.0))) *
                                  (zt(i, k) - real_t(4000.0));
        } else {
          qtm_forcing(i, k) = real_t(0.0);
        }
        qtm_forcing(i, k) = qtm_forcing(i, k) / g_per_kg; // g/(kg s) -> kg/(kg s).
      });

  // Convert total-water specific-humidity forcing to mixing-ratio forcing.
  force_spec_hum_to_mixing_ratio(ngrdcol, nzt, rtm, qtm_forcing, rtm_forcing);
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "rico_tndcy::passive", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
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
        "rico_tndcy::eddy", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//----------------------------------------------------------------------
inline void rico_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D um_sfc,
                        Real1D vm_sfc, Real1D thlm, Real1D rtm, Real1D z_bot, Real1D p_sfc,
                        Real1D exner_sfc, int saturation_formula, Real1D upwp_sfc, Real1D vpwp_sfc,
                        Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc,
                        const time_dependent_input_state& input) {
  // Description:
  // Compute RICO surface forcing using the RICO 3D LES surface formulations.
  // Retain the source's inactive old-ATEX branch and ordered per-column formulas.
  auto workspace_scope               = workspace.scope();
  constexpr real_t C_10              = real_t(0.0013);   // ATEX drag coefficient.
  constexpr real_t C_m_20            = real_t(0.001229); // RICO momentum coefficient.
  constexpr real_t C_h_20            = real_t(0.001094); // RICO heat coefficient.
  constexpr real_t C_q_20            = real_t(0.001133); // RICO water coefficient.
  constexpr real_t z0                = real_t(0.00015);  // ATEX roughness length [m].
  constexpr real_t standard_flux_alt = real_t(20.0);     // Reference height [m].
  Real1D rsat = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "rico_sfclyr::rsat",
                                             ngrdcol),
         ubar = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "rico_sfclyr::ubar",
                                             ngrdcol);
  real_t time_frac;
  int before_time, after_time;

  // Interpolate variables from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  // Choose the scheme, then define variable values.
  const bool l_use_old_atex = false;
  compute_ubar(ngrdcol, um_sfc, vm_sfc, ubar);
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "rico_sfclyr::temperature", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        T_sfc(i) = T_sfc_interp;
        ustar(i) = real_t(0.3);
        rsat(i)  = sat_mixrat_liq_api(p_sfc(i), T_sfc(i), saturation_formula);
      });

  // Compute heat and moisture fluxes.
  if (l_use_old_atex) {
    Kokkos::parallel_for(
        "rico_sfclyr::atex", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          // Stevens et al. (2000), eq. 3: height adjustment for the ATEX specification.
          const real_t Cz = C_10 * (Kokkos::log(real_t(10.0) / z0) / Kokkos::log(z_bot(i) / z0)) *
                            (Kokkos::log(real_t(10.0) / z0) / Kokkos::log(z_bot(i) / z0));
          // Consume the column-local coefficients without intermediate fields.
          wpthlp_sfc(i) = -Cz * ubar(i) * (thlm(i) - T_sfc(i) / exner_sfc(i));
          wprtp_sfc(i)  = -Cz * ubar(i) * (rtm(i) - rsat(i));
        });

    compute_momentum_flux(ngrdcol, um_sfc, vm_sfc, ubar, ustar, upwp_sfc, vpwp_sfc);
  } else {
    Kokkos::parallel_for(
        "rico_sfclyr::rico", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          // ATEX height adjustment applied to each RICO reference coefficient.
          const real_t Cm = C_m_20 *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0)) *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0));
          const real_t Ch = C_h_20 *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0)) *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0));
          const real_t Cq = C_q_20 *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0)) *
                            (Kokkos::log(standard_flux_alt / z0) / Kokkos::log(z_bot(i) / z0));
          // Consume the column-local coefficients without intermediate fields.
          wpthlp_sfc(i) = -Ch * ubar(i) * (thlm(i) - T_sfc(i) / exner_sfc(i));
          wprtp_sfc(i)  = -Cq * ubar(i) * (rtm(i) - rsat(i));

          upwp_sfc(i) = -um_sfc(i) * Cm * ubar(i);
          vpwp_sfc(i) = -vm_sfc(i) * Cm * ubar(i);
        });
  }
}

} // namespace clubb_kokkos
