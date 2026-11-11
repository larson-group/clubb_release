// Prescribe inversion-based forcing and surface fluxes for the GCSS ATEX case.
//
// atex_tndcy applies forcing after 90 min, locates the 6.5 g/kg water threshold,
// prescribes subsidence and computes heat/water tendencies. An existing dycore
// grid can supply the forcing profiles through the permanent interpolation and
// conservative-remapping owners. atex_sfclyr interpolates temperature [K] at
// time [s] and supplies heat [K m/s], water[(kg/kg) m/s] and ustar [m/s]. Model fields
// are (column, level[,scalar]); scratch reuses the explicit case workspace. The final surface
// context replaces Native time-series module state. Missing inversions set the
// relevant column's fatal status at debug level 2, preserving source branch order.
// Reference: Stevens et al. (2000), Simulations of trade-wind cumuli under a strong
// inversion, J. Atmos. Sci. 58,1870–1891.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/remapping_module.hpp"
#include "sfc_flux.hpp"
#include "time_dependent_input.hpp"
#include <iostream>

namespace clubb_kokkos {

//----------------------------------------------------------------------
static inline void calc_forcings(int ngrdcol, const grid& gr, Real1D z_inversion,
                                 Real2D thlm_forcing, Real2D rtm_forcing) {
  // Description:
  // Fill ATEX temperature and moisture forcing from the diagnosed inversion
  // height [m], retaining the strict source bounds at zero and the inversion.
  const int nzt = gr.nzt;
  const auto zt = gr.zt;
  using Policy  = Kokkos::RangePolicy<typename Real2D::execution_space>;
  Kokkos::parallel_for(
      "atex::calc_forcings", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i = (index % ngrdcol);
        const int k = (index / ngrdcol);
        // Theta_l tendency from the specified ATEX constants.
        if (zt(i, k) > real_t(0.0) && zt(i, k) < z_inversion(i)) {
          thlm_forcing(i, k) = -real_t(1.1575e-5) * (real_t(3.0) - zt(i, k) / z_inversion(i));
        } else if (zt(i, k) > z_inversion(i) && zt(i, k) <= z_inversion(i) + real_t(300.0)) {
          thlm_forcing(i, k) =
              -real_t(2.315e-5) * (real_t(1.0) - (zt(i, k) - z_inversion(i)) / real_t(300.0));
        } else {
          thlm_forcing(i, k) = real_t(0.0);
        }
        // Moisture tendency(Brian's specified constants).
        if (zt(i, k) > real_t(0.0) && zt(i, k) < z_inversion(i))
          rtm_forcing(i, k) = -real_t(1.58e-8) * (real_t(1.0) - zt(i, k) / z_inversion(i));
        else
          rtm_forcing(i, k) = real_t(0.0);
      });
}

//======================================================================
inline void atex_tndcy(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                       const sclr_idx_type& sclr_idx, const grid& gr, double time,
                       double time_initial, Real2D rtm, bool l_add_dycore_grid,
                       int grid_remap_method, const grid& gr_dycore, int total_idx_rho_lin_spline,
                       Real2D rho_lin_spline_vals, Real2D rho_lin_spline_levels, Real1D p_sfc,
                       err_info_type& err_info, Real2D wm_zt, Real2D wm_zm, Real2D thlm_forcing,
                       Real2D rtm_forcing, Real3D sclrm_forcing, Real3D edsclrm_forcing) {
  // Description:
  // Set ATEX theta_l/water forcing and vertical motion after 5400 s from the
  // current inversion height, optionally evaluating the forcing on a dycore
  // grid before conservative remapping. Return early on missing-inversion error.
  auto workspace_scope   = workspace.scope();
  constexpr int iv_other = 1;
  const int nzt          = gr.nzt;
  const int nzm          = gr.nzm;
  const int dycore_nzt   = gr_dycore.nzt;
  const auto zt          = gr.zt;
  const auto dycore_zt   = gr_dycore.zt;
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  if (time >= time_initial + double(5400.0)) {
    auto z_lev = workspace.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
        Kokkos::WithoutInitializing, "atex_tndcy::z_lev", ngrdcol);
    auto z_inversion =
        workspace.view<Real1D>(Kokkos::WithoutInitializing, "atex_tndcy::z_inversion", ngrdcol);
    // Identify the height of the 6.5 g/kg moisture level in each column.
    Kokkos::parallel_for(
        "atex_tndcy::level", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          z_lev(i) = 0;
          while (z_lev(i) < nzt && rtm(i, z_lev(i)) > real_t(6.5e-3))
            z_lev(i) = z_lev(i) + 1;
        });
    if (clubb_at_least_debug_level_api(2)) {
      const auto levels = to_host_space(workspace, z_lev);
      for (int i = 0; i < ngrdcol; ++i) {
        if (levels(i) == nzt || levels(i) == 0) {
          const auto water   = to_host_space(workspace, rtm);
          const auto headers = to_host_space(workspace, err_info.err_header);
          auto codes         = to_host_space(workspace, err_info.err_code);
          std::cerr.write(headers.data() + (headers.extent(0) * i), headers.extent(0));
          std::cerr << "\nIdentification of 6.5 g/kg level failed\n";
          std::cerr << "Subroutine: atex_tndcy. File: atex.F\n";
          std::cerr << "k = " << levels(i) + 1 << " i = " << ((i + 1)) << '\n';
          // Native prints rtm(nzt+1) after an unsuccessful top search. That
          // diagnostic read is undefined; print the value only for a valid
          // index, while preserving its fatal status and early return.
          if (levels(i) < nzt)
            std::cerr << "rtm(k) = " << water(i, levels(i)) << '\n';
          else
            std::cerr << "rtm(k) = <outside sounding>\n";
          codes(i) = clubb_fatal_error;
          Kokkos::deep_copy(err_info.err_code, codes);
          return;
        }
      }
    }
    Kokkos::parallel_for(
        "atex_tndcy::inversion", Policy(0, ngrdcol),
        KOKKOS_LAMBDA(int i) { z_inversion(i) = zt(i, z_lev(i) - 1); });

    // Large-scale subsidence from the diagnosed inversion height.
    Kokkos::parallel_for(
        "atex_tndcy::motion", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (zt(i, k) > real_t(0.0) && zt(i, k) <= z_inversion(i)) {
            wm_zt(i, k) = -real_t(0.0065) * zt(i, k) / z_inversion(i);
          } else if (zt(i, k) > z_inversion(i) && zt(i, k) <= z_inversion(i) + real_t(300.0)) {
            wm_zt(i, k) =
                -real_t(0.0065) * (real_t(1.0) - (zt(i, k) - z_inversion(i)) / real_t(300.0));
          } else {
            wm_zt(i, k) = real_t(0.0);
          }
        });
    zt2zm_api(nzm, nzt, ngrdcol, gr, wm_zt, wm_zm);
    // Boundary conditions on momentum-grid vertical motion.
    Kokkos::parallel_for(
        "atex_tndcy::boundaries", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          wm_zm(i, 0)       = real_t(0.0);
          wm_zm(i, (nzm - 1)) = real_t(0.0);
        });

    if (l_add_dycore_grid) {
      auto z_lev_dycore = workspace.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "atex_tndcy::z_lev_dycore", ngrdcol);
      auto z_inversion_dycore  = workspace.view<Real1D>(Kokkos::WithoutInitializing,
                                                        "atex_tndcy::z_inversion_dycore", ngrdcol);
      auto thlm_forcing_dycore = workspace.view<Real2D>(
          Kokkos::WithoutInitializing, "atex_tndcy::thlm_forcing_dycore", ngrdcol, dycore_nzt);
      auto rtm_forcing_dycore = workspace.view<Real2D>(
          Kokkos::WithoutInitializing, "atex_tndcy::rtm_forcing_dycore", ngrdcol, dycore_nzt);
      auto rtm_dycore = workspace.view<Real2D>(Kokkos::WithoutInitializing,
                                               "atex_tndcy::rtm_dycore", ngrdcol, dycore_nzt);
      for (int i = 0; i < ngrdcol; ++i) {
        lin_interp_between_grids(
            workspace, dycore_nzt, nzt, Kokkos::subview(dycore_zt, i, Kokkos::ALL()),
            Kokkos::subview(zt, i, Kokkos::ALL()), Kokkos::subview(rtm, i, Kokkos::ALL()),
            Kokkos::subview(rtm_dycore, i, Kokkos::ALL()));
      }
      // Locate the inversion again using the interpolated dycore water field.
      Kokkos::parallel_for(
          "atex_tndcy::dycore_level", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
            z_lev_dycore(i) = 0;
            while (z_lev_dycore(i) < dycore_nzt && rtm_dycore(i, z_lev_dycore(i)) > real_t(6.5e-3))
              z_lev_dycore(i) = z_lev_dycore(i) + 1;
          });
      if (clubb_at_least_debug_level_api(2)) {
        const auto levels = to_host_space(workspace, z_lev_dycore);
        for (int i = 0; i < ngrdcol; ++i) {
          if (levels(i) == dycore_nzt || levels(i) == 0) {
            const auto water   = to_host_space(workspace, rtm_dycore);
            const auto headers = to_host_space(workspace, err_info.err_header);
            auto codes         = to_host_space(workspace, err_info.err_code);
            std::cerr.write(headers.data() + (headers.extent(0) * i), headers.extent(0));
            std::cerr << "\nIdentification of 6.5 g/kg level failed on dycore grid\n";
            std::cerr << "Subroutine: atex_tndcy. File: atex.F\n";
            std::cerr << "k = " << levels(i) + 1 << " i = " << ((i + 1)) << '\n';
            // Same bounded diagnostic adaptation as the physics-grid error.
            if (levels(i) < dycore_nzt)
              std::cerr << "rtm_dycore(k) = " << water(i, levels(i)) << '\n';
            else
              std::cerr << "rtm_dycore(k) = <outside sounding>\n";
            codes(i) = clubb_fatal_error;
            Kokkos::deep_copy(err_info.err_code, codes);
            return;
          }
        }
      }
      Kokkos::parallel_for(
          "atex_tndcy::dycore_inversion", Policy(0, ngrdcol),
          KOKKOS_LAMBDA(int i) { z_inversion_dycore(i) = dycore_zt(i, z_lev_dycore(i) - 1); });
      calc_forcings(ngrdcol, gr_dycore, z_inversion_dycore,
                     thlm_forcing_dycore, rtm_forcing_dycore);
      const bool l_zt_variable = true;
      remap_vals_to_target(workspace, ngrdcol, gr_dycore, gr, dycore_nzt, thlm_forcing_dycore, nzt,
                           total_idx_rho_lin_spline, rho_lin_spline_vals, rho_lin_spline_levels,
                           iv_other, p_sfc, grid_remap_method, l_zt_variable, thlm_forcing);
      remap_vals_to_target(workspace, ngrdcol, gr_dycore, gr, dycore_nzt, rtm_forcing_dycore, nzt,
                           total_idx_rho_lin_spline, rho_lin_spline_vals, rho_lin_spline_levels,
                           iv_other, p_sfc, grid_remap_method, l_zt_variable, rtm_forcing);
    } else {
      calc_forcings(ngrdcol, gr, z_inversion, thlm_forcing, rtm_forcing);
    }
  } else {
    // Forcing is applied only after 5400 s from initialization.
    Kokkos::parallel_for(
        "atex_tndcy::before_start", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i        = (index % ngrdcol);
          const int k        = (index / ngrdcol);
          wm_zt(i, k)        = real_t(0.0);
          thlm_forcing(i, k) = real_t(0.0);
          rtm_forcing(i, k)  = real_t(0.0);
        });
    Kokkos::deep_copy(wm_zm, real_t(0.0));
  }
  if (sclr_dim > 0) {
    Kokkos::parallel_for(
        "atex_tndcy::passive", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
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
        "atex_tndcy::eddy", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
          const int i = (index % ngrdcol);
          const int k = (index / ngrdcol);
          if (sclr_idx.iiedsclr_thl > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
          if (sclr_idx.iiedsclr_rt > -1)
            edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
        });
  }
}

//======================================================================
inline void atex_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D ubar,
                        Real1D thlm_sfc, Real1D rtm_sfc, Real1D exner_sfc, Real1D wpthlp_sfc,
                        Real1D wprtp_sfc, Real1D ustar, Real1D T_sfc,
                        const time_dependent_input_state& input) {
  // Description:
  // Compute GCSS ATEX heat/water fluxes from the interpolated temperature and
  // specified bulk coefficient/moisture adjustment, with ustar = 0.3 m/s.
  auto workspace_scope = workspace.scope();
  int before_time, after_time;
  real_t time_frac;
  time_select(time, static_cast<int>(input.time_sfc_given.extent(0)),
              input.time_sfc_given, before_time, after_time, time_frac);
  const auto temperature = to_host_space(workspace, input.T_sfc_given);
  const real_t T_sfc_interp =
      linear_interp_factor(time_frac, temperature(after_time), temperature(before_time));
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  // The bulk coefficients are constant; compute both fluxes in the same column kernel.
  Kokkos::parallel_for(
      "atex_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        ustar(i)      = real_t(0.3);
        T_sfc(i)      = T_sfc_interp;
        wpthlp_sfc(i) = -real_t(0.0013) * ubar(i) * (thlm_sfc(i) - T_sfc(i) / exner_sfc(i));
        wprtp_sfc(i)  = -real_t(0.0013) * ubar(i) * (rtm_sfc(i) - real_t(0.0198293));
      });
}

} // namespace clubb_kokkos
