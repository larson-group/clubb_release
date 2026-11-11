// Coordinate sampled microphysics and select interactive SILHS feedback.
//
// lh_microphys_driver computes SILHS tendencies and sampled variance/covariance
// sources, writes the sampled fields and, in interactive mode, copies means into
// the ordinary model tendencies. Non-interactive mode retains diagnostic output.
// Input samples are (column,sample,level[,variate]); fields/tendencies use
// (column,level[,species]), all column-first LayoutLeft. Appended case records
// replace implicit Native settings; the callback retains its Native arguments.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/error_code.hpp"
#include "Microphys/estimate_scm_microphys_module.hpp"
#include "SILHS/math_utilities.hpp"
#include "SILHS/latin_hypercube_driver_module.hpp"
#include "SILHS/lh_microphys_var_covar_module.hpp"
#include "SILHS/est_kessler_microphys_module.hpp"
#include "Microphys/microphys_host_types.hpp"
#include "Microphys/KK_microphys/KK_microphys.hpp"
#include "Microphys/advance_microphys_module.hpp"
#include "Microphys/Morrison_microphys/engine_contracts.hpp"
namespace clubb_kokkos::microphys_host {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Microphysics>
inline void lh_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 lh_sample_point_weights, const pdf_parameter& pdf_params,
    const precipitation_fractions& precip_fracs, Array2 p_in_Pa, Array3 exner, Array4 rho,
    Array5 rcm, Array6 delta_zt, Array7 cloud_frac, Array8 hydromet, Array9 X_mixt_comp_all_levs,
    Array10 lh_rt_clipped, Array11 lh_thl_clipped, Array12 lh_rc_clipped, Array13 lh_rv_clipped,
    Array14 lh_Nc_clipped, bool l_lh_importance_sampling, bool l_lh_instant_var_covar_src,
    int saturation_formula, stats_type& stats, Array15 lh_hydromet_mc, Array16 lh_hydromet_vel,
    Array17 lh_Ncm_mc, Array18 lh_rcm_mc, Array19 lh_rvm_mc, Array20 lh_thlm_mc, Array21 lh_rtp2_mc,
    Array22 lh_thlp2_mc, Array23 lh_wprtp_mc, Array24 lh_wpthlp_mc, Array25 lh_rtpthlp_mc,
    Array26 lh_AKm, Array27 AKm, Array28 AKstd, Array29 AKstd_cld, Array30 lh_rcm_avg,
    Array31 AKm_rcm, Array32 AKm_rcc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Microphysics>
inline void lh_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 lh_sample_point_weights, const pdf_parameter& pdf_params,
    const precipitation_fractions& precip_fracs, Array2 p_in_Pa, Array3 exner, Array4 rho,
    Array5 rcm, Array6 delta_zt, Array7 cloud_frac, Array8 hydromet, Array9 X_mixt_comp_all_levs,
    Array10 lh_rt_clipped, Array11 lh_thl_clipped, Array12 lh_rc_clipped, Array13 lh_rv_clipped,
    Array14 lh_Nc_clipped, bool l_lh_importance_sampling, bool l_lh_instant_var_covar_src,
    int saturation_formula, stats_type& stats, Array15 lh_hydromet_mc, Array16 lh_hydromet_vel,
    Array17 lh_Ncm_mc, Array18 lh_rcm_mc, Array19 lh_rvm_mc, Array20 lh_thlm_mc, Array21 lh_rtp2_mc,
    Array22 lh_thlp2_mc, Array23 lh_wprtp_mc, Array24 lh_wpthlp_mc, Array25 lh_rtpthlp_mc,
    Array26 lh_AKm, Array27 AKm, Array28 AKstd, Array29 AKstd_cld, Array30 lh_rcm_avg,
    Array31 AKm_rcm, Array32 AKm_rcc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Description:
  //   Computes an estimate of the change due to microphysics given a set of
  //   subcolumns of thlm, rtm, et cetera from the subcolumn generator.
  //
  // References:
  //   None
  //---------------------------------------------------------------------------
  // Interface block
  // Input Variables
  // Input Variables
  // Output Variables
  // Output Variables

  // ---- Begin Code ----

  // Perform LH and analytic microphysical calculations
  // As a test of SILHS, compute an estimate of Kessler microphysics
  auto workspace_scope = workspace.scope();
  if (clubb_at_least_debug_level_api(2)) {
    est_kessler_microphys_api(workspace, nzt, num_samples, pdf_dim, ngrdcol, X_nl_all_levs,
                              pdf_params, rcm, cloud_frac, X_mixt_comp_all_levs,
                              lh_sample_point_weights, l_lh_importance_sampling, lh_AKm, AKm, AKstd,
                              AKstd_cld, AKm_rcm, AKm_rcc, lh_rcm_avg);
  } else {
    for (int index_1 = 0; index_1 < lh_AKm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_AKm.extent_int(0); ++index_0) {
        lh_AKm(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < AKm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < AKm.extent_int(0); ++index_0) {
        AKm(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < AKstd.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < AKstd.extent_int(0); ++index_0) {
        AKstd(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < AKstd_cld.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < AKstd_cld.extent_int(0); ++index_0) {
        AKstd_cld(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < lh_rcm_avg.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_rcm_avg.extent_int(0); ++index_0) {
        lh_rcm_avg(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < AKm_rcm.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < AKm_rcm.extent_int(0); ++index_0) {
        AKm_rcm(index_0, index_1) = real_t(0.);
      }
    }
    for (int index_1 = 0; index_1 < AKm_rcc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < AKm_rcc.extent_int(0); ++index_0) {
        AKm_rcc(index_0, index_1) = real_t(0.);
      }
    }
  }

  // Call the latin hypercube microphysics driver for microphys_sub
  microphys_host::est_silhs_tndcy(
      workspace, gr, ngrdcol, dt, nzt, nzm, num_samples, pdf_dim, hydromet_dim, hm_metadata,
      X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights, pdf_params, precip_fracs,
      p_in_Pa, exner, rho, delta_zt, hydromet, rcm, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,
      lh_rv_clipped, lh_Nc_clipped, l_lh_instant_var_covar_src, saturation_formula, stats,
      lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, lh_rtp2_mc,
      lh_thlp2_mc, lh_wprtp_mc, lh_wpthlp_mc, lh_rtpthlp_mc, microphys_sub, microphys, kk);

  return;
}
} // namespace clubb_kokkos::microphys_host

namespace clubb_kokkos {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Microphysics>
inline void lh_microphys_driver(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 lh_sample_point_weights, const pdf_parameter& pdf_params,
    const precipitation_fractions& precip_fracs, Array2 p_in_Pa, Array3 exner, Array4 rho,
    Array5 rcm, Array6 delta_zt, Array7 cloud_frac, Array8 hydromet, Array9 X_mixt_comp_all_levs,
    Array10 lh_rt_clipped, Array11 lh_thl_clipped, Array12 lh_rc_clipped, Array13 lh_rv_clipped,
    Array14 lh_Nc_clipped, bool l_lh_importance_sampling, bool l_lh_instant_var_covar_src,
    int saturation_formula, stats_type& stats, Array15 lh_hydromet_mc, Array16 lh_hydromet_vel,
    Array17 lh_Ncm_mc, Array18 lh_rcm_mc, Array19 lh_rvm_mc, Array20 lh_thlm_mc, Array21 lh_rtp2_mc,
    Array22 lh_thlp2_mc, Array23 lh_wprtp_mc, Array24 lh_wpthlp_mc, Array25 lh_rtpthlp_mc,
    Array26 lh_AKm, Array27 AKm, Array28 AKstd, Array29 AKstd_cld, Array30 lh_rcm_avg,
    Array31 AKm_rcm, Array32 AKm_rcc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  // In: gr
  // In: ngrdcol — Number of model columns
  // In: dt — Model timestep       [s]
  // In: nzt — Number of thermodynamic vertical model levels
  // In: nzm — Number of momentum vertical model levels
  // In: num_samples — Number of calls to microphysics per timestep (normally=2)
  // In: pdf_dim — Number of variables to sample
  // In: hydromet_dim
  // In: hm_metadata
  // In: X_nl_all_levs (ngrdcol,num_samples,nzt,pdf_dim) — Sample that is transformed ultimately to normal-lognormal
  // In: lh_sample_point_weights (ngrdcol,num_samples,nzt) — Weight given the individual sample points
  // In: pdf_params — PDF parameters       [units vary]
  // In: precip_fracs — Precipitation fractions      [-]
  // In: p_in_Pa (ngrdcol,nzt) — Pressure                     [Pa]
  // In: exner (ngrdcol,nzt) — Exner function               [-]
  // In: rho (ngrdcol,nzt) — Density on thermo. grid      [kg/m^3]
  // In: rcm (ngrdcol,nzt) — Liquid water mixing ratio    [kg/kg]
  // In: delta_zt (ngrdcol,nzt) — Change in meters with height [m]
  // In: cloud_frac (ngrdcol,nzt) — Cloud fraction               [-]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor species    [units vary]
  // In: X_mixt_comp_all_levs (ngrdcol,num_samples,nzt) — Which mixture component we're in
  // In: lh_rt_clipped (ngrdcol,num_samples,nzt) — rt generated from silhs sample points
  // In: lh_thl_clipped (ngrdcol,num_samples,nzt) — thl generated from silhs sample points
  // In: lh_rc_clipped (ngrdcol,num_samples,nzt) — rc generated from silhs sample points
  // In: lh_rv_clipped (ngrdcol,num_samples,nzt) — rv generated from silhs sample points
  // In: lh_Nc_clipped (ngrdcol,num_samples,nzt) — Nc generated from silhs sample points
  // In: l_lh_importance_sampling — Do importance sampling (SILHS) [-]
  // In: l_lh_instant_var_covar_src — Produce instantaneous var/covar tendencies [-]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // Inout: stats
  // Out: lh_hydromet_mc (ngrdcol,nzt,hydromet_dim) — LH estimate of hydrometeor time tendency          [(units vary)/s]
  // Out: lh_hydromet_vel (ngrdcol,nzt,hydromet_dim) — LH estimate of hydrometeor sedimentation velocity [m/s]
  // Out: lh_Ncm_mc (ngrdcol,nzt) — LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
  // Out: lh_rcm_mc (ngrdcol,nzt) — LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
  // Out: lh_rvm_mc (ngrdcol,nzt) — LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
  // Out: lh_thlm_mc (ngrdcol,nzt) — LH estimate of time tndcy. of liquid potential temp.  [K/s]
  // Out: lh_rtp2_mc (ngrdcol,nzm) — LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
  // Out: lh_thlp2_mc (ngrdcol,nzm) — LH microphysics tendency for <thl'^2>                 [K^2/s]
  // Out: lh_wprtp_mc (ngrdcol,nzm) — LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
  // Out: lh_wpthlp_mc (ngrdcol,nzm) — LH microphysics tendency for <w'thl'>                 [m*K/s^2]
  // Out: lh_rtpthlp_mc (ngrdcol,nzm) — LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]
  // Out: lh_AKm (ngrdcol,nzt) — Kessler ac estimate                 [kg/kg/s]
  // Out: AKm (ngrdcol,nzt) — Exact Kessler ac                    [kg/kg/s]
  // Out: AKstd (ngrdcol,nzt) — St dev of exact Kessler ac          [kg/kg/s]
  // Out: AKstd_cld (ngrdcol,nzt) — Stdev of exact w/in cloud ac        [kg/kg/s]
  // Out: lh_rcm_avg (ngrdcol,nzt) — Monte Carlo rcm estimate            [kg/kg]
  // Out: AKm_rcm (ngrdcol,nzt) — Kessler ac based on rcm             [kg/kg/s]
  // Out: AKm_rcc (ngrdcol,nzt) — Kessler ac based on rcm/cloud_frac  [kg/kg/s]
  // In: microphys_sub
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: kk replaces the Native parameters_KK module.

  auto workspace_scope              = workspace.scope();
  auto X_nl_all_levs_host           = to_host_space(workspace, X_nl_all_levs);
  auto lh_sample_point_weights_host = to_host_space(workspace, lh_sample_point_weights);
  auto p_in_Pa_host                 = to_host_space(workspace, p_in_Pa);
  auto exner_host                   = to_host_space(workspace, exner);
  auto rho_host                     = to_host_space(workspace, rho);
  auto rcm_host                     = to_host_space(workspace, rcm);
  auto delta_zt_host                = to_host_space(workspace, delta_zt);
  auto cloud_frac_host              = to_host_space(workspace, cloud_frac);
  auto hydromet_host                = to_host_space(workspace, hydromet);
  auto X_mixt_comp_all_levs_host    = to_host_space(workspace, X_mixt_comp_all_levs);
  auto lh_rt_clipped_host           = to_host_space(workspace, lh_rt_clipped);
  auto lh_thl_clipped_host          = to_host_space(workspace, lh_thl_clipped);
  auto lh_rc_clipped_host           = to_host_space(workspace, lh_rc_clipped);
  auto lh_rv_clipped_host           = to_host_space(workspace, lh_rv_clipped);
  auto lh_Nc_clipped_host           = to_host_space(workspace, lh_Nc_clipped);
  auto lh_hydromet_mc_host          = to_host_space(workspace, lh_hydromet_mc);
  auto lh_hydromet_vel_host         = to_host_space(workspace, lh_hydromet_vel);
  auto lh_Ncm_mc_host               = to_host_space(workspace, lh_Ncm_mc);
  auto lh_rcm_mc_host               = to_host_space(workspace, lh_rcm_mc);
  auto lh_rvm_mc_host               = to_host_space(workspace, lh_rvm_mc);
  auto lh_thlm_mc_host              = to_host_space(workspace, lh_thlm_mc);
  auto lh_rtp2_mc_host              = to_host_space(workspace, lh_rtp2_mc);
  auto lh_thlp2_mc_host             = to_host_space(workspace, lh_thlp2_mc);
  auto lh_wprtp_mc_host             = to_host_space(workspace, lh_wprtp_mc);
  auto lh_wpthlp_mc_host            = to_host_space(workspace, lh_wpthlp_mc);
  auto lh_rtpthlp_mc_host           = to_host_space(workspace, lh_rtpthlp_mc);
  auto lh_AKm_host                  = to_host_space(workspace, lh_AKm);
  auto AKm_host                     = to_host_space(workspace, AKm);
  auto AKstd_host                   = to_host_space(workspace, AKstd);
  auto AKstd_cld_host               = to_host_space(workspace, AKstd_cld);
  auto lh_rcm_avg_host              = to_host_space(workspace, lh_rcm_avg);
  auto AKm_rcm_host                 = to_host_space(workspace, AKm_rcm);
  auto AKm_rcc_host                 = to_host_space(workspace, AKm_rcc);
  microphys_host::lh_microphys_driver(
      workspace, gr, ngrdcol, dt, nzt, nzm, num_samples, pdf_dim, hydromet_dim, hm_metadata,
      X_nl_all_levs_host, lh_sample_point_weights_host, pdf_params, precip_fracs, p_in_Pa_host,
      exner_host, rho_host, rcm_host, delta_zt_host, cloud_frac_host, hydromet_host,
      X_mixt_comp_all_levs_host, lh_rt_clipped_host, lh_thl_clipped_host, lh_rc_clipped_host,
      lh_rv_clipped_host, lh_Nc_clipped_host, l_lh_importance_sampling, l_lh_instant_var_covar_src,
      saturation_formula, stats, lh_hydromet_mc_host, lh_hydromet_vel_host, lh_Ncm_mc_host,
      lh_rcm_mc_host, lh_rvm_mc_host, lh_thlm_mc_host, lh_rtp2_mc_host, lh_thlp2_mc_host,
      lh_wprtp_mc_host, lh_wpthlp_mc_host, lh_rtpthlp_mc_host, lh_AKm_host, AKm_host, AKstd_host,
      AKstd_cld_host, lh_rcm_avg_host, AKm_rcm_host, AKm_rcc_host, microphys_sub, microphys, kk);
  from_execution_space(workspace, lh_hydromet_mc, lh_hydromet_mc_host);
  from_execution_space(workspace, lh_hydromet_vel, lh_hydromet_vel_host);
  from_execution_space(workspace, lh_Ncm_mc, lh_Ncm_mc_host);
  from_execution_space(workspace, lh_rcm_mc, lh_rcm_mc_host);
  from_execution_space(workspace, lh_rvm_mc, lh_rvm_mc_host);
  from_execution_space(workspace, lh_thlm_mc, lh_thlm_mc_host);
  from_execution_space(workspace, lh_rtp2_mc, lh_rtp2_mc_host);
  from_execution_space(workspace, lh_thlp2_mc, lh_thlp2_mc_host);
  from_execution_space(workspace, lh_wprtp_mc, lh_wprtp_mc_host);
  from_execution_space(workspace, lh_wpthlp_mc, lh_wpthlp_mc_host);
  from_execution_space(workspace, lh_rtpthlp_mc, lh_rtpthlp_mc_host);
  from_execution_space(workspace, lh_AKm, lh_AKm_host);
  from_execution_space(workspace, AKm, AKm_host);
  from_execution_space(workspace, AKstd, AKstd_host);
  from_execution_space(workspace, AKstd_cld, AKstd_cld_host);
  from_execution_space(workspace, lh_rcm_avg, lh_rcm_avg_host);
  from_execution_space(workspace, AKm_rcm, AKm_rcm_host);
  from_execution_space(workspace, AKm_rcc, AKm_rcc_host);
}
} // namespace clubb_kokkos
