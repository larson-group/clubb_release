// Prescribe tendencies, subsidence and surface fluxes for the long ATEX case.
//
// atex_long_tndcy supplies piecewise heat/water tendencies [K/s,(kg/kg)/s] and
// vertical motion [m/s], with 12 h of spin up and selected scalar copies. Its private
// calc_forcings retains the specified profiles. atex_long_sfclyr interpolates
// temperature [K] at time [s] and supplies bulk heat [K m/s], water[(kg/kg) m/s] and
// ustar [m/s]. The source's inactive prescribed-flux alternative is also retained.
// Model fields are column-first; the final surface context owns the input series.
// References: Hing Ong, The nontraditional Coriolis terms and trade-wind cumuli;
// Stevens et al. (2000), Simulations of trade-wind cumuli under a strong inversion,
// J. Atmos. Sci. 58,1870–1891.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"

namespace clubb_kokkos {

//----------------------------------------------------------------------
static inline void calc_forcings(int ngrdcol, const grid& gr, Real2D thlm_forcing, Real2D rtm_forcing) {
  // Description:
  // Fill the specified long-ATEX theta_l and water forcing profiles from each
  // column's height [m], retaining the source intervals and constants.
  const int nzt = gr.nzt;
  const auto zt = gr.zt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;
  Kokkos::parallel_for(
      "atex_long::calc_forcings", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        // Theta_l tendency(Hing's specified constants).
        if (zt(i, k) >= real_t(0.0) && zt(i, k) < real_t(1400.0)) {
          thlm_forcing(i, k) = -real_t(3.5805e-5);
        } else if (zt(i, k) >= real_t(1400.0) && zt(i, k) < real_t(1650.0)) {
          thlm_forcing(i, k) =
              -real_t(3.5805e-5) + real_t(1.1935e-5) * (zt(i, k) - real_t(1400.0)) * real_t(0.004);
        } else if (zt(i, k) >= real_t(1650.0) && zt(i, k) < real_t(2990.0)) {
          thlm_forcing(i, k) =
              -real_t(2.3870e-5) - real_t(0.1155e-5) * (zt(i, k) - real_t(1650.0)) / real_t(1350.0);
        } else {
          thlm_forcing(i, k) = real_t(0.0);
        }
        // Moisture tendency(Brian's specified constants).
        if (zt(i, k) >= real_t(0.0) && zt(i, k) < real_t(1050.0))
          rtm_forcing(i, k) = -real_t(1.58e-8) * (real_t(1.0) - zt(i, k) / real_t(1050.0));
        else
          rtm_forcing(i, k) = real_t(0.0);
      });
}

//======================================================================
inline void atex_long_tndcy(int ngrdcol, int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx,
                            const grid& gr, double time,
                            Real2D wm_zt, Real2D wm_zm,
                            Real2D thlm_forcing, Real2D rtm_forcing,
                            Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set theta_l/water tendencies and subsidence for long ATEX, scaling the
  // first 43200 s linearly through spin-up before interpolating vertical motion.
  const int nzt = gr.nzt;
  const int nzm = gr.nzm;
  const auto zt = gr.zt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // Large-scale subsidence(Hing's specified constants).
  Kokkos::parallel_for(
      "atex_long_tndcy::motion", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i = ((index % ngrdcol));
        const int k = ((index / ngrdcol));
        if (zt(i, k) >= real_t(0.0) && zt(i, k) < real_t(1050.0)) {
          wm_zt(i, k) = -real_t(0.00636) * zt(i, k) / real_t(1050.0);
        } else if (zt(i, k) >= real_t(1050.0) && zt(i, k) < real_t(1650.0)) {
          wm_zt(i, k) =
              -real_t(0.00636) - real_t(0.00079) * (zt(i, k) - real_t(1050.0)) / real_t(600.0);
        } else {
          wm_zt(i, k) = -real_t(0.00715);
        }
      });
  calc_forcings(ngrdcol, gr, thlm_forcing, rtm_forcing);

  // Spin-up period; the mixed time/core precision follows Native arithmetic.
  const real_t spinup = real_t(43200.0);
  if (time < spinup) {
    Kokkos::parallel_for(
        "atex_long_tndcy::spinup", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i        = ((index % ngrdcol));
          const int k        = ((index / ngrdcol));
          thlm_forcing(i, k) = thlm_forcing(i, k) * time / spinup;
          rtm_forcing(i, k)  = rtm_forcing(i, k) * time / spinup;
          wm_zt(i, k)        = wm_zt(i, k) * time / spinup;
        });
  }
  zt2zm_api(nzm, nzt, ngrdcol, gr, wm_zt, wm_zm);

  // Boundary conditions on vertical motion.
  Kokkos::parallel_for(
      "atex_long_tndcy::boundaries", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wm_zm(i, 0)       = real_t(0.0);
        wm_zm(i, ((nzm - 1))) = real_t(0.0);
      });
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "atex_long_tndcy::passive", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
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
        "atex_long_tndcy::eddy", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
          const int i = ((index % ngrdcol));
          const int k = ((index / ngrdcol));
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//======================================================================
inline void atex_long_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D ubar,
                             Real1D thlm_sfc, Real1D rtm_sfc, Real1D exner_sfc, Real1D rho_sfc,
                             Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc,
                             const time_dependent_input_state& input) {
  // Description:
  // Compute GCSS long-ATEX surface heat/water exchange with the source's
  // prescribed adjustment, retaining its inactive fixed-flux spin-up branch.
  auto workspace_scope      = workspace.scope();
  const bool l_compute_flux = true;
  int before_time, after_time;
  real_t time_frac;

  // Interpolate surface temperature from time_dependent_input.
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));

  real_t sens_ht{}, latent_ht{}; // Assigned in the inactive prescribed-flux branch.
  if (!l_compute_flux) {
    sens_ht             = real_t(10.70);
    latent_ht           = real_t(154.0);
    const real_t spinup = real_t(86400.0);
    if (time < spinup) {
      sens_ht   = sens_ht - real_t(2.20) * (spinup - time) / spinup;
      latent_ht = latent_ht + real_t(40.0) * (spinup - time) / spinup;
    }
  }
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "atex_long_sfclyr::coefficients", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        const real_t C_10 = real_t(0.0013);
        // Source value approximates 0.981*sat_mixrat_liq(101600 Pa, 298 K, Flatau).
        const real_t adjustment = real_t(0.0194664);
        ustar(i)      = real_t(0.3);
        T_sfc(i)      = T_sfc_interp;
        if (!l_compute_flux) {
          wpthlp_sfc(i) = convert_sens_ht_to_km_s(sens_ht, rho_sfc(i));
          wprtp_sfc(i)  = convert_latent_ht_to_m_s(latent_ht, rho_sfc(i));
        }
        if (l_compute_flux) {
          wpthlp_sfc(i) = -C_10 * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
          wprtp_sfc(i)  = -C_10 * ubar(i) * (rtm_sfc(i) - adjustment);
        }
      });
}

} // namespace clubb_kokkos
