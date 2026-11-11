// Accumulate sampled model diagnostics; see stats_clubb_utilities.hpp.
#include "scratch_workspace.hpp"
#include "CLUBB_core/stats_clubb_utilities.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/numerical_check.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/Skx_module.hpp"
#include "CLUBB_core/T_in_K_module.hpp"

#include <cmath>
#include <limits>
#include <string>

namespace clubb_kokkos {

//===============================================================================
void stats_accumulate(ScratchWorkspace& workspace,                                 // InOut
                      int nzm, int nzt, int ngrdcol, int sclr_dim, int edsclr_dim, // In
                      const grid& gr, real_t dt,                                   // In
                      bool l_implemented, bool l_host_applies_sfc_fluxes,          // In
                      bool l_stability_correct_tau_zm,                             // In
                      StatsHostReal2D clubb_params,                                // In
                      StatsHostReal2D um, StatsHostReal2D vm, StatsHostReal2D upwp,
                      StatsHostReal2D vpwp, StatsHostReal2D up2, StatsHostReal2D vp2, // In
                      StatsHostReal2D thlm, StatsHostReal2D rtm, StatsHostReal2D thlm_before,
                      StatsHostReal2D rtm_before, StatsHostReal2D thlm_forcing,
                      StatsHostReal2D rtm_forcing, // In
                      StatsHostReal1D wpthlp_sfc, StatsHostReal1D wprtp_sfc, StatsHostReal2D wprtp,
                      StatsHostReal2D wpthlp, // In
                      StatsHostReal2D wp2, StatsHostReal2D wp3, StatsHostReal2D rtp2,
                      StatsHostReal2D rtp3, StatsHostReal2D thlp2, StatsHostReal2D thlp3, // In
                      StatsHostReal2D rtpthlp,                                            // In
                      StatsHostReal2D p_in_Pa, StatsHostReal2D exner, StatsHostReal2D rho,
                      StatsHostReal2D rho_zm, // In
                      StatsHostReal2D rho_ds_zm, StatsHostReal2D rho_ds_zt,
                      StatsHostReal2D thv_ds_zm,
                      StatsHostReal2D thv_ds_zt,                                           // In
                      StatsHostReal2D wm_zt, StatsHostReal2D wm_zm, StatsHostReal2D rcm,   // In
                      StatsHostReal2D cloud_frac,                                          // In
                      StatsHostReal2D thvm, StatsHostReal2D ug, StatsHostReal2D vg,        // In
                      StatsHostReal2D ddzt_umvm_sqd, StatsHostReal2D stability_correction, // In
                      StatsHostReal2D Kh_zt,                                               // In
                      StatsHostReal2D rsat,                                                // In
                      StatsHostReal2D Kh_zm,                                               // In
                      StatsHostReal2D em,                                                  // In
                      StatsHostReal3D sclrm, StatsHostReal3D sclrp2,                       // In
                      StatsHostReal3D sclrprtp, StatsHostReal3D sclrpthlp,
                      StatsHostReal3D sclrm_forcing, // In
                      StatsHostReal3D wpsclrp, StatsHostReal3D wpedsclrp,
                      StatsHostReal3D edsclrm,         // In
                      StatsHostReal3D edsclrm_forcing, // In
                      int saturation_formula,          // In
                      stats_type& stats                // InOut
) {
  // Description:
  //   Sample the model's thermodynamic/momentum profiles and column diagnostics.
  //   Field units follow the corresponding Native stats registry entries.
  //   stats retains the sampled values; all scientific inputs remain unchanged.
  // References: None.

  // Local Variables
  auto workspace_scope = workspace.scope();
  int sclr, edsclr, k, i;
  StatsHostReal2D T_in_K, rsati,
      rcm_in_cloud;         // Temperature [K], ice saturation and cloudy water [kg/kg].
  StatsHostReal2D shear;    // Wind shear production term [m^2/s^3].
  StatsHostReal2D dzm, dzt; // Signed momentum/thermodynamic-grid layer thickness [m].
  StatsHostReal2D wp3_on_wp2, wp3_on_wp2_zt; // Smoothed ratio of third to second w moment [m/s].
  Real2D wp3_zm, Skw_zm, gamma_Skw_fnc; // Momentum third moment [m^3/s^3], skewness and gamma [-].
  StatsHostReal1D xtmp;
  StatsHostReal1D rtm_integral_before, rtm_integral_after, rtm_integral_forcing;
  StatsHostReal1D rtm_flux_top, rtm_flux_sfc, rtm_spur_src;
  StatsHostReal1D thlm_integral_before, thlm_integral_after, thlm_integral_forcing;
  StatsHostReal1D thlm_flux_top, thlm_flux_sfc, thlm_spur_src;
  std::string sclr_idx; // Native ADJUSTL/TRIM of an integer written into character(len=50).

  // ---- Begin Code ----

  // Sample fields
  if (stats.l_sample) {
    dzm = workspace_scope.view<StatsHostReal2D>(Kokkos::WithoutInitializing,
                                                "stats_accumulate::dzm", ngrdcol, nzm);
    dzt = workspace_scope.view<StatsHostReal2D>(Kokkos::WithoutInitializing,
                                                "stats_accumulate::dzt", ngrdcol, nzt);
    const auto gr_dzm = to_host_space(workspace, gr.dzm);
    const auto gr_dzt = to_host_space(workspace, gr.dzt);
    for (k = 0; k < nzm; ++k)
      for (i = 0; i < ngrdcol; ++i)
        dzm(i, k) = gr.grid_dir * gr_dzm(i, k);
    for (k = 0; k < nzt; ++k)
      for (i = 0; i < ngrdcol; ++i)
        dzt(i, k) = gr.grid_dir * gr_dzt(i, k);

    // The ratio is used only by these requested diagnostics; no model state depends on it.
    const bool sample_wp3_ratio = var_on_stats_list(stats, "wp3_on_wp2") ||
                                  var_on_stats_list(stats, "wp3_on_wp2_zt") ||
                                  var_on_stats_list(stats, "wp3_on_wp2_cfl_num");
    if (sample_wp3_ratio) {
      // The existing smoothing owner uses execution Views. Cross that API once;
      // the sampled diagnostics below use its returned host snapshots.
      Real2D wp3_on_wp2_exec = workspace_scope.view<Real2D>(
          Kokkos::WithoutInitializing, "stats_accumulate::wp3_on_wp2", ngrdcol, nzm);
      Real2D wp3_on_wp2_zt_exec = workspace_scope.view<Real2D>(
          Kokkos::WithoutInitializing, "stats_accumulate::wp3_on_wp2_zt", ngrdcol, nzt);
      calc_wp3_on_wp2(workspace, nzm, nzt, ngrdcol, gr, to_execution_space(workspace, wp2),
                      to_execution_space(workspace, wp3), wp3_on_wp2_exec, wp3_on_wp2_zt_exec);
      wp3_on_wp2    = to_host_space(workspace, wp3_on_wp2_exec);
      wp3_on_wp2_zt = to_host_space(workspace, wp3_on_wp2_zt_exec);
    }

    // Native updates its sampled arrays on the host at this boundary. Generic
    // entry overloads already own those snapshots, including zero-sized scalar
    // fields, which transfer no values when their Native dimensions are zero.
    if (var_on_stats_list(stats, "T_in_K") || var_on_stats_list(stats, "rsati")) {
      // These inputs already belong to the host statistics boundary. Retain the
      // Native 2D temperature expression without a host/device round trip.
      T_in_K = workspace_scope.view<StatsHostReal2D>(Kokkos::WithoutInitializing,
                                                     "stats_accumulate::T_in_K", ngrdcol, nzt);
      for (k = 0; k < nzt; ++k)
        for (i = 0; i < ngrdcol; ++i)
          T_in_K(i, k) = thlm(i, k) * exner(i, k) + Lv * rcm(i, k) / Cp;
      stats_update("T_in_K", T_in_K, stats);
    }
    if (var_on_stats_list(stats, "rsati")) {
      rsati = workspace_scope.view<StatsHostReal2D>(Kokkos::WithoutInitializing,
                                                    "stats_accumulate::rsati", ngrdcol, nzt);
      for (k = 0; k < nzt; ++k)
        for (i = 0; i < ngrdcol; ++i)
          rsati(i, k) = sat_mixrat_ice_k(p_in_Pa(i, k), T_in_K(i, k), saturation_formula);
      stats_update("rsati", rsati, stats);
    }
    if (var_on_stats_list(stats, "rcm_in_cloud")) {
      rcm_in_cloud = workspace_scope.view<StatsHostReal2D>(
          Kokkos::WithoutInitializing, "stats_accumulate::rcm_in_cloud", ngrdcol, nzt);
      for (k = 0; k < nzt; ++k)
        for (i = 0; i < ngrdcol; ++i)
          if (cloud_frac(i, k) > cloud_frac_min)
            rcm_in_cloud(i, k) = rcm(i, k) / cloud_frac(i, k);
          else
            rcm_in_cloud(i, k) = rcm(i, k);
      stats_update("rcm_in_cloud", rcm_in_cloud, stats);
    }
    if (var_on_stats_list(stats, "shear")) {
      shear = workspace_scope.view<StatsHostReal2D>(Kokkos::WithoutInitializing,
                                                    "stats_accumulate::shear", ngrdcol, nzm);
      const auto gr_invrs_dzm = to_host_space(workspace, gr.invrs_dzm);
      for (i = 0; i < ngrdcol; ++i)
        shear(i, 0) = real_t(0.0);
      for (k = 1; k < (nzm - 1); ++k)
        for (i = 0; i < ngrdcol; ++i)
          shear(i, k) = -upwp(i, k) * (um(i, k) - um(i, (k - 1))) * gr_invrs_dzm(i, k) -
                        vpwp(i, k) * (vm(i, k) - vm(i, (k - 1))) * gr_invrs_dzm(i, k);
      for (i = 0; i < ngrdcol; ++i)
        shear(i, (nzm - 1)) = real_t(0.0);
      stats_update("shear", shear, stats);
    }

    // stats_zt variables
    stats_update("thlm", thlm, stats);
    stats_update("thvm", thvm, stats);
    stats_update("rtm", rtm, stats);
    stats_update("rcm", rcm, stats);
    stats_update("um", um, stats);
    stats_update("vm", vm, stats);
    stats_update("wm_zt", wm_zt, stats);
    stats_update("ug", ug, stats);
    stats_update("vg", vg, stats);
    stats_update("p_in_Pa", p_in_Pa, stats);
    stats_update("exner", exner, stats);
    stats_update("rho_ds_zt", rho_ds_zt, stats);
    stats_update("thv_ds_zt", thv_ds_zt, stats);
    stats_update("wp3", wp3, stats);
    stats_update("Kh_zt", Kh_zt, stats);
    stats_update("rho", rho, stats);
    stats_update("rsat", rsat, stats);
    stats_update("thlp3", thlp3, stats);
    stats_update("rtp3", rtp3, stats);
    stats_update("wp3_on_wp2_zt", wp3_on_wp2_zt, stats);
    if (sclr_dim > 0) {
      for (sclr = 0; sclr < sclr_dim; ++sclr) {
        sclr_idx = std::to_string(sclr + 1);
        stats_update("sclr" + sclr_idx + "m",
                     Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("sclr" + sclr_idx + "m_f",
                     Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr), stats);
      }
    }
    if (edsclr_dim > 0) {
      for (edsclr = 0; edsclr < edsclr_dim; ++edsclr) {
        sclr_idx = std::to_string(edsclr + 1);
        stats_update("edsclr" + sclr_idx + "m",
                     Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr), stats);
        stats_update("edsclr" + sclr_idx + "m_f",
                     Kokkos::subview(edsclrm_forcing, Kokkos::ALL, Kokkos::ALL, edsclr), stats);
      }
    }

    // stats_zm variables
    stats_update("wm_zm", wm_zm, stats);
    stats_update("ddzt_umvm_sqd", ddzt_umvm_sqd, stats);
    stats_update("wp2", wp2, stats);
    stats_update("rtp2", rtp2, stats);
    stats_update("thlp2", thlp2, stats);
    stats_update("rtpthlp", rtpthlp, stats);
    stats_update("wprtp", wprtp, stats);
    stats_update("wpthlp", wpthlp, stats);
    if (l_stability_correct_tau_zm)
      stats_update("stability_correction", stability_correction, stats);
    stats_update("Kh_zm", Kh_zm, stats);
    stats_update("upwp", upwp, stats);
    stats_update("vpwp", vpwp, stats);
    stats_update("vp2", vp2, stats);
    stats_update("up2", up2, stats);
    stats_update("rho_zm", rho_zm, stats);
    stats_update("rho_ds_zm", rho_ds_zm, stats);
    stats_update("thv_ds_zm", thv_ds_zm, stats);
    stats_update("em", em, stats);
    stats_update("wp3_on_wp2", wp3_on_wp2, stats);
    if (var_on_stats_list(stats, "wp3_on_wp2_cfl_num")) {
      StatsHostReal2D wp3_on_wp2_cfl_num = workspace_scope.view<StatsHostReal2D>(
          Kokkos::WithoutInitializing, "stats_accumulate::wp3_on_wp2_cfl_num", ngrdcol, nzm);
      for (k = 0; k < nzm; ++k)
        for (i = 0; i < ngrdcol; ++i)
          wp3_on_wp2_cfl_num(i, k) = wp3_on_wp2(i, k) * dt / dzm(i, k);
      stats_update("wp3_on_wp2_cfl_num", wp3_on_wp2_cfl_num, stats);
    }

    if (var_on_stats_list(stats, "gamma_Skw_fnc")) {
      wp3_zm        = to_execution_space(workspace, zt2zm_api(nzm, nzt, ngrdcol, gr, wp3));
      Skw_zm        = workspace_scope.view<Real2D>("stats_accumulate::Skw_zm", ngrdcol, nzm);
      gamma_Skw_fnc = workspace_scope.view<Real2D>("stats_accumulate::gamma_Skw_fnc", ngrdcol, nzm);
      const auto clubb_params_exec = to_execution_space(workspace, clubb_params);
      Skx_func(nzm, ngrdcol, to_execution_space(workspace, wp2), wp3_zm, w_tol, clubb_params_exec,
               Skw_zm);
      compute_gamma_Skw(nzm, ngrdcol, Skw_zm, clubb_params_exec, // In
                        gamma_Skw_fnc);                          // Out
      stats_update("gamma_Skw_fnc", gamma_Skw_fnc, stats);
    }
    if (sclr_dim > 0) {
      for (sclr = 0; sclr < sclr_dim; ++sclr) {
        sclr_idx = std::to_string(sclr + 1);
        stats_update("sclr" + sclr_idx + "p2",
                     Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("sclr" + sclr_idx + "prtp",
                     Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("sclr" + sclr_idx + "pthlp",
                     Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("wpsclr" + sclr_idx + "p",
                     Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
      }
    }
    if (edsclr_dim > 0) {
      for (edsclr = 0; edsclr < edsclr_dim; ++edsclr) {
        sclr_idx = std::to_string(edsclr + 1);
        stats_update("wpedsclr" + sclr_idx + "p",
                     Kokkos::subview(wpedsclrp, Kokkos::ALL, Kokkos::ALL, edsclr), stats);
      }
    }

    // stats_sfc variables
    if (var_on_stats_list(stats, "cc")) {
      for (i = 0; i < ngrdcol; ++i) {
        real_t cc = -std::numeric_limits<real_t>::max();
        for (k = 0; k < nzt; ++k)
          if (cloud_frac(i, k) > cc)
            cc = cloud_frac(i, k);
        stats_update("cc", cc, stats, i);
      }
    }

    xtmp = workspace_scope.view<StatsHostReal1D>("stats_accumulate::xtmp", ngrdcol);
    if (var_on_stats_list(stats, "z_cloud_base")) {
      const auto gr_zt = to_host_space(workspace, gr.zt);
      for (i = 0; i < ngrdcol; ++i) {
        k = 0;
        while (rcm(i, k) < rc_tol && k < (nzt - 1))
          ++k;
        if (k == 0)
          xtmp(i) = gr_zt(i, 0);
        else if (k > 0 && k < (nzt - 1))
          xtmp(i) = lin_interpolate_two_points(rc_tol, rcm(i, k), rcm(i, (k - 1)), gr_zt(i, k),
                                               gr_zt(i, (k - 1)));
        else
          xtmp(i) = real_t(-10.0);
      }
      stats_update("z_cloud_base", xtmp, stats);
    }
    if (var_on_stats_list(stats, "lwp")) {
      for (i = 0; i < ngrdcol; ++i)
        xtmp(i) = vertical_integral(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                    Kokkos::subview(rcm, i, Kokkos::ALL),
                                    Kokkos::subview(dzt, i, Kokkos::ALL));
      stats_update("lwp", xtmp, stats);
    }
    if (var_on_stats_list(stats, "vwp")) {
      StatsHostReal2D vapor = workspace_scope.view<StatsHostReal2D>(
          Kokkos::WithoutInitializing, "stats_accumulate::vapor", ngrdcol, nzt);
      for (k = 0; k < nzt; ++k)
        for (i = 0; i < ngrdcol; ++i)
          vapor(i, k) = rtm(i, k) - rcm(i, k);
      for (i = 0; i < ngrdcol; ++i)
        xtmp(i) = vertical_integral(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                    Kokkos::subview(vapor, i, Kokkos::ALL),
                                    Kokkos::subview(dzt, i, Kokkos::ALL));
      stats_update("vwp", xtmp, stats);
    }
    for (i = 0; i < ngrdcol; ++i) {
      if (var_on_stats_list(stats, "thlm_vert_avg")) {
        stats_update("thlm_vert_avg",
                     vertical_avg(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(thlm, i, Kokkos::ALL),
                                  Kokkos::subview(dzt, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "rtm_vert_avg")) {
        stats_update("rtm_vert_avg",
                     vertical_avg(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(rtm, i, Kokkos::ALL),
                                  Kokkos::subview(dzt, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "um_vert_avg")) {
        stats_update("um_vert_avg",
                     vertical_avg(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(um, i, Kokkos::ALL),
                                  Kokkos::subview(dzt, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "vm_vert_avg")) {
        stats_update("vm_vert_avg",
                     vertical_avg(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(vm, i, Kokkos::ALL),
                                  Kokkos::subview(dzt, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "wp2_vert_avg")) {
        stats_update("wp2_vert_avg",
                     vertical_avg(workspace, nzm, Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),
                                  Kokkos::subview(wp2, i, Kokkos::ALL),
                                  Kokkos::subview(dzm, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "up2_vert_avg")) {
        stats_update("up2_vert_avg",
                     vertical_avg(workspace, nzm, Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),
                                  Kokkos::subview(up2, i, Kokkos::ALL),
                                  Kokkos::subview(dzm, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "vp2_vert_avg")) {
        stats_update("vp2_vert_avg",
                     vertical_avg(workspace, nzm, Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),
                                  Kokkos::subview(vp2, i, Kokkos::ALL),
                                  Kokkos::subview(dzm, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "rtp2_vert_avg")) {
        stats_update("rtp2_vert_avg",
                     vertical_avg(workspace, nzm, Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),
                                  Kokkos::subview(rtp2, i, Kokkos::ALL),
                                  Kokkos::subview(dzm, i, Kokkos::ALL)),
                     stats, i);
      }
      if (var_on_stats_list(stats, "thlp2_vert_avg")) {
        stats_update("thlp2_vert_avg",
                     vertical_avg(workspace, nzm, Kokkos::subview(rho_ds_zm, i, Kokkos::ALL),
                                  Kokkos::subview(thlp2, i, Kokkos::ALL),
                                  Kokkos::subview(dzm, i, Kokkos::ALL)),
                     stats, i);
      }
    }

    if (var_on_stats_list(stats, "tot_vartn_normlzd_rtm")) {
      for (i = 0; i < ngrdcol; ++i) {
        if (std::abs(rtm(i, (nzt - 1)) - rtm(i, 0)) < eps)
          xtmp(i) = real_t(-999);
        else {
          xtmp(i) = real_t(0);
          for (k = 1; k < nzt; ++k)
            xtmp(i) = xtmp(i) + std::abs(rtm(i, k) - rtm(i, (k - 1))) /
                                    std::abs(rtm(i, (nzt - 1)) - rtm(i, 0));
        }
      }
      stats_update("tot_vartn_normlzd_rtm", xtmp, stats);
    }

    if (var_on_stats_list(stats, "tot_vartn_normlzd_thlm")) {
      for (i = 0; i < ngrdcol; ++i) {
        if (std::abs(thlm(i, (nzt - 1)) - thlm(i, 0)) < eps)
          xtmp(i) = real_t(-999);
        else {
          xtmp(i) = real_t(0);
          for (k = 1; k < nzt; ++k)
            xtmp(i) = xtmp(i) + std::abs(thlm(i, k) - thlm(i, (k - 1))) /
                                    std::abs(thlm(i, (nzt - 1)) - thlm(i, 0));
        }
      }
      stats_update("tot_vartn_normlzd_thlm", xtmp, stats);
    }

    if (var_on_stats_list(stats, "tot_vartn_normlzd_wprtp")) {
      for (i = 0; i < ngrdcol; ++i) {
        if (std::abs(wprtp(i, (nzm - 1)) - wprtp(i, 0)) < eps)
          xtmp(i) = real_t(-999);
        else {
          xtmp(i) = real_t(0);
          for (k = 1; k < nzm; ++k)
            xtmp(i) = xtmp(i) + std::abs(wprtp(i, k) - wprtp(i, (k - 1))) /
                                    std::abs(wprtp(i, (nzm - 1)) - wprtp(i, 0));
        }
      }
      stats_update("tot_vartn_normlzd_wprtp", xtmp, stats);
    }

    if (var_on_stats_list(stats, "rtm_spur_src") || var_on_stats_list(stats, "thlm_spur_src")) {
      // Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
      // Therefore, wm must be zero or l_implemented must be true.
      rtm_integral_before =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_integral_before", ngrdcol);
      rtm_integral_after =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_integral_after", ngrdcol);
      rtm_integral_forcing =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_integral_forcing", ngrdcol);
      rtm_flux_top =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_flux_top", ngrdcol);
      rtm_flux_sfc =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_flux_sfc", ngrdcol);
      rtm_spur_src =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::rtm_spur_src", ngrdcol);
      thlm_integral_before =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_integral_before", ngrdcol);
      thlm_integral_after =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_integral_after", ngrdcol);
      thlm_integral_forcing =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_integral_forcing", ngrdcol);
      thlm_flux_top =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_flux_top", ngrdcol);
      thlm_flux_sfc =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_flux_sfc", ngrdcol);
      thlm_spur_src =
          workspace_scope.view<StatsHostReal1D>("stats_accumulate::thlm_spur_src", ngrdcol);
      for (i = 0; i < ngrdcol; ++i) {
        bool all_wm_zt = true, all_wm_zm = true;
        for (k = 0; k < nzt; ++k)
          all_wm_zt = all_wm_zt && std::abs(wm_zt(i, k)) < eps;
        for (k = 0; k < nzm; ++k)
          all_wm_zm = all_wm_zm && std::abs(wm_zm(i, k)) < eps;
        if (l_implemented || (all_wm_zt && all_wm_zm)) {
          rtm_flux_top(i) = rho_ds_zm(i, gr.k_ub_zm) * wprtp(i, gr.k_ub_zm);
          if (!l_host_applies_sfc_fluxes)
            rtm_flux_sfc(i) = rho_ds_zm(i, gr.k_lb_zm) * wprtp_sfc(i);
          else
            rtm_flux_sfc(i) = real_t(0.0);

          // Get the vertical integral of rtm before this function begins
          // so that spurious source can be calculated.
          rtm_integral_before(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(rtm_before, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          rtm_integral_after(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(rtm, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          rtm_integral_forcing(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(rtm_forcing, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          // Calculate the spurious source for rtm.
          rtm_spur_src(i) = calculate_spurious_source(rtm_integral_after(i), rtm_integral_before(i),
                                                      rtm_flux_top(i), rtm_flux_sfc(i),
                                                      rtm_integral_forcing(i), dt);
          thlm_flux_top(i) = rho_ds_zm(i, gr.k_ub_zm) * wpthlp(i, gr.k_ub_zm);
          if (!l_host_applies_sfc_fluxes)
            thlm_flux_sfc(i) = rho_ds_zm(i, gr.k_lb_zm) * wpthlp_sfc(i);
          else
            thlm_flux_sfc(i) = real_t(0.0);

          // Get the vertical integral of thlm before this function begins
          // so that spurious source can be calculated.
          thlm_integral_before(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(thlm_before, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          thlm_integral_after(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(thlm, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          thlm_integral_forcing(i) = vertical_integral(
              workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
              Kokkos::subview(thlm_forcing, i, Kokkos::ALL), Kokkos::subview(dzt, i, Kokkos::ALL));
          // Calculate the spurious source for thlm.
          thlm_spur_src(i) = calculate_spurious_source(
              thlm_integral_after(i), thlm_integral_before(i), thlm_flux_top(i), thlm_flux_sfc(i),
              thlm_integral_forcing(i), dt);
        } else { // If l_implemented is false, we don't want spurious source output.
          rtm_spur_src(i)  = real_t(-9999.0);
          thlm_spur_src(i) = real_t(-9999.0);
        }
      }
      stats_update("rtm_spur_src", rtm_spur_src, stats);
      stats_update("thlm_spur_src", thlm_spur_src, stats);
    }
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void stats_accumulate_hydromet_api(ScratchWorkspace& workspace,                   // InOut
                                   const grid& gr, int ngrdcol, int hydromet_dim, // In
                                   const hm_metadata_type& hm_metadata,
                                   StatsHostReal3D hydromet,  // In
                                   StatsHostReal2D rho_ds_zt, // In
                                   stats_type& stats          // InOut
) {
  // Description:
  //   Compute hydrometeor statistics for all model columns.
  //   hydromet contains all species except cloud water in Native metadata order;
  //   rho_ds_zt is dry static density [kg/m^3]. stats retains profiles/water paths.
  auto workspace_scope = workspace.scope();
  (void)hydromet_dim; // Species count is carried by the Native array contract.
  StatsHostReal1D xtmp;
  int i;

  // ---- Begin Code ----
  if (!stats.l_sample)
    return;
  if (hm_metadata.iirr > -1) {
    stats_update("rrm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirr),
                 stats);
  }
  if (hm_metadata.iirs > -1) {
    stats_update("rsm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirs),
                 stats);
  }
  if (hm_metadata.iiri > -1) {
    stats_update("rim", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiri),
                 stats);
  }
  if (hm_metadata.iirg > -1) {
    stats_update("rgm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirg),
                 stats);
  }
  if (hm_metadata.iiNi > -1) {
    stats_update("Nim", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNi),
                 stats);
  }
  if (hm_metadata.iiNr > -1) {
    stats_update("Nrm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNr),
                 stats);
  }
  if (hm_metadata.iiNs > -1) {
    stats_update("Nsm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNs),
                 stats);
  }
  if (hm_metadata.iiNg > -1) {
    stats_update("Ngm", Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNg),
                 stats);
  }
  xtmp = workspace_scope.view<StatsHostReal1D>("stats_accumulate_hydromet_api::xtmp", ngrdcol);
  const auto gr_dzt = to_host_space(workspace, gr.dzt);
  // Snow Water Path
  if (var_on_stats_list(stats, "swp") && hm_metadata.iirs > -1) {
    for (i = 0; i < ngrdcol; ++i)
      xtmp(i) = vertical_integral(workspace, gr.nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(hydromet, i, Kokkos::ALL, hm_metadata.iirs),
                                  Kokkos::subview(gr_dzt, i, Kokkos::ALL));
    stats_update("swp", xtmp, stats);
  }
  // Ice Water Path
  if (var_on_stats_list(stats, "iwp") && hm_metadata.iiri > -1) {
    for (i = 0; i < ngrdcol; ++i)
      xtmp(i) = vertical_integral(workspace, gr.nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(hydromet, i, Kokkos::ALL, hm_metadata.iiri),
                                  Kokkos::subview(gr_dzt, i, Kokkos::ALL));
    stats_update("iwp", xtmp, stats);
  }
  // Rain Water Path
  if (var_on_stats_list(stats, "rwp") && hm_metadata.iirr > -1) {
    for (i = 0; i < ngrdcol; ++i)
      xtmp(i) = vertical_integral(workspace, gr.nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                  Kokkos::subview(hydromet, i, Kokkos::ALL, hm_metadata.iirr),
                                  Kokkos::subview(gr_dzt, i, Kokkos::ALL));
    stats_update("rwp", xtmp, stats);
  }
}

//------------------------------------------------------------------------------
void stats_accumulate_lh_tend(ScratchWorkspace& workspace, // InOut
                              const grid& gr, int ngrdcol, int hydromet_dim,
                              const hm_metadata_type& hm_metadata,                       // In
                              StatsHostReal3D lh_hydromet_mc, StatsHostReal2D lh_Ncm_mc, // In
                              StatsHostReal2D lh_thlm_mc, StatsHostReal2D lh_rvm_mc,
                              StatsHostReal2D lh_rcm_mc, // In
                              StatsHostReal2D lh_AKm, StatsHostReal2D AKm, StatsHostReal2D AKstd,
                              StatsHostReal2D AKstd_cld, // In
                              StatsHostReal2D lh_rcm_avg, StatsHostReal2D AKm_rcm,
                              StatsHostReal2D AKm_rcc, // In
                              stats_type& stats        // InOut
) {
  // Description:
  //   Compute full-field stats for the tendency of Latin hypercube sample points.
  //   Fields retain Native (column, thermodynamic-level[, hydrometeor]) shapes;
  //   stats receives tendencies and analytic/Monte Carlo Kessler diagnostics.
  // References: None.
  (void)gr;
  (void)ngrdcol;
  (void)hydromet_dim;
  if (stats.l_sample) {
    // Native synchronizes its host microphysics results with device copies here.
    // The generic C++ entry snapshots those current fields before any update.
    stats_update("lh_thlm_mc", lh_thlm_mc, stats);
    stats_update("lh_rcm_mc", lh_rcm_mc, stats);
    stats_update("lh_rvm_mc", lh_rvm_mc, stats);
    stats_update("lh_Ncm_mc", lh_Ncm_mc, stats);
    if (hm_metadata.iirr > -1) {
      stats_update("lh_rrm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirr),
                   stats);
    }
    if (hm_metadata.iirs > -1) {
      stats_update("lh_rsm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirs),
                   stats);
    }
    if (hm_metadata.iiri > -1) {
      stats_update("lh_rim_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiri),
                   stats);
    }
    if (hm_metadata.iirg > -1) {
      stats_update("lh_rgm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirg),
                   stats);
    }
    if (hm_metadata.iiNi > -1) {
      stats_update("lh_Nim_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNi),
                   stats);
    }
    if (hm_metadata.iiNr > -1) {
      stats_update("lh_Nrm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNr),
                   stats);
    }
    if (hm_metadata.iiNs > -1) {
      stats_update("lh_Nsm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNs),
                   stats);
    }
    if (hm_metadata.iiNg > -1) {
      stats_update("lh_Ngm_mc",
                   Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNg),
                   stats);
    }
    stats_update("AKm", AKm, stats);
    stats_update("lh_AKm", lh_AKm, stats);
    stats_update("lh_rcm_avg", lh_rcm_avg, stats);
    stats_update("AKstd", AKstd, stats);
    stats_update("AKstd_cld", AKstd_cld, stats);
    stats_update("AKm_rcm", AKm_rcm, stats);
    stats_update("AKm_rcc", AKm_rcc, stats);
  }
}
} // namespace clubb_kokkos
