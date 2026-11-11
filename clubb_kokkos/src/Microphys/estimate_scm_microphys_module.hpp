// Estimate grid-mean microphysics tendencies from SILHS samples.
//
// est_silhs_tndcy calls the selected microphysics for every sample, computes
// weighted means and variance/covariance sources, and writes sampled diagnostics.
// adjust_KK_src_means applies the source mean-only clipping option for convergence
// comparisons with analytic KK. Samples are (column,sample,level[,variate]);
// atmospheric fields are (column,level[,species]), all with LayoutLeft storage.
// A generic callable keeps the full Native microphysics callback argument order;
// appended case records replace the implicit settings read by the source helpers.
#pragma once
#include "scratch_workspace.hpp"
#include "Microphys/KK_microphys_module.hpp"
#include "Microphys/silhs_category_variance_module.hpp"
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
          class Array24, class Microphysics>
inline void est_silhs_tndcy(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 X_mixt_comp_all_levs, Array2 lh_sample_point_weights,
    const pdf_parameter& pdf_params, const precipitation_fractions& precip_fracs, Array3 p_in_Pa,
    Array4 exner, Array5 rho, Array6 dzq, Array7 hydromet, Array8 rcm, Array9 lh_rt_clipped,
    Array10 lh_thl_clipped, Array11 lh_rc_clipped, Array12 lh_rv_clipped, Array13 lh_Nc_clipped,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats,
    Array14 lh_hydromet_mc, Array15 lh_hydromet_vel, Array16 lh_Ncm_mc, Array17 lh_rvm_mc,
    Array18 lh_rcm_mc, Array19 lh_thlm_mc, Array20 lh_rtp2_mc, Array21 lh_thlp2_mc,
    Array22 lh_wprtp_mc, Array23 lh_wpthlp_mc, Array24 lh_rtpthlp_mc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20>
inline void adjust_KK_src_means(ScratchWorkspace& workspace, real_t dt, int nzt, int ngrdcol,
                                Array0 exner, Array1 rcm, Array2 rrm, Array3 Nrm, Array4 hydromet,
                                int hydromet_dim, int iiri, Array5 rrm_auto, Array6 rrm_accr,
                                Array7 rrm_evap, Array8 Nrm_auto, Array9 Nrm_evap, Array10 lh_Vrr,
                                Array11 lh_VNr, Array12 rrm_mc, Array13 Nrm_mc, Array14 rvm_mc,
                                Array15 rcm_mc, Array16 thlm_mc, Array17 rrm_src_adj,
                                Array18 Nrm_src_adj, Array19 rrm_evap_adj, Array20 Nrm_evap_adj,
                                const parameters_KK& kk);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Microphysics>
inline void est_silhs_tndcy(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 X_mixt_comp_all_levs, Array2 lh_sample_point_weights,
    const pdf_parameter& pdf_params, const precipitation_fractions& precip_fracs, Array3 p_in_Pa,
    Array4 exner, Array5 rho, Array6 dzq, Array7 hydromet, Array8 rcm, Array9 lh_rt_clipped,
    Array10 lh_thl_clipped, Array11 lh_rc_clipped, Array12 lh_rv_clipped, Array13 lh_Nc_clipped,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats,
    Array14 lh_hydromet_mc, Array15 lh_hydromet_vel, Array16 lh_Ncm_mc, Array17 lh_rvm_mc,
    Array18 lh_rcm_mc, Array19 lh_thlm_mc, Array20 lh_rtp2_mc, Array21 lh_thlp2_mc,
    Array22 lh_wprtp_mc, Array23 lh_wpthlp_mc, Array24 lh_rtpthlp_mc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Description:
  //   Estimate the tendency of a microphysics scheme via latin hypercube sampling
  //
  // References:
  //   None
  //-------------------------------------------------------------------------------
  // External
  // Constant parameters
  // Input Variables
  // Output Variables
  // Local Variables
  // These parameters are not used by the microphysics scheme when SILHS is
  // turned on.
  // Just to avoid typing hm_metadata%iiPDF_x everywhere
  auto workspace_scope = workspace.scope();
  // Reuse one execution-space snapshot for the unchanged sample weights.
  const auto sample_weights_execution = to_execution_space(workspace, lh_sample_point_weights);
  const auto hm_metadata_iirr         = hm_metadata.iirr;
  const auto hm_metadata_iiri         = hm_metadata.iiri;
  const auto hm_metadata_iiNr         = hm_metadata.iiNr;
  const auto hm_metadata_iiPDF_chi = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_w      = hm_metadata.iiPDF_w;
  constexpr bool l_latin_hypercube = true;
  HostReal4D lh_hydromet_mc_all       = workspace_scope.view<HostReal4D>(
      "est_silhs_tndcy::lh_hydromet_mc_all", ngrdcol, num_samples, nzt,
      hydromet_dim); // LH est of hydrometeor time tendency          [(units vary)/s]
  HostReal4D lh_hydromet_vel_all = workspace_scope.view<HostReal4D>(
      "est_silhs_tndcy::lh_hydromet_vel_all", ngrdcol, num_samples, nzt,
      hydromet_dim); // LH est of hydrometeor sedimentation velocity [m/s]
  HostReal3D lh_Ncm_mc_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::lh_Ncm_mc_all", ngrdcol, num_samples,
      nzt); // LH est of time tendency of cloud droplet concentration  [#/kg/s]
  HostReal3D lh_rcm_mc_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::lh_rcm_mc_all", ngrdcol, num_samples,
      nzt); // LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
  HostReal3D lh_rvm_mc_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::lh_rvm_mc_all", ngrdcol, num_samples,
      nzt); // LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
  HostReal3D lh_thlm_mc_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::lh_thlm_mc_all", ngrdcol, num_samples,
      nzt); // LH est of time tendency of liquid potential temperature     [K/s]
  HostReal4D hydromet_all_points = workspace_scope.view<HostReal4D>(
      "est_silhs_tndcy::hydromet_all_points", ngrdcol, num_samples, nzt,
      hydromet_dim); // Hydrometeor species                    [units vary]
  HostReal3D Ncn_all_points = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::Ncn_all_points", ngrdcol, num_samples,
      nzt); // Cloud Nuclei conc. (simplified); Nc=Ncn*H(chi) [#/kg]
  HostReal3D chi_all_points = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::chi_all_points", ngrdcol, num_samples,
      nzt); // 's' (Mellor 1977)                              [kg/kg]
  HostReal3D w_all_points =
      workspace_scope.view<HostReal3D>("est_silhs_tndcy::w_all_points", ngrdcol, num_samples,
                                       nzt); // Vertical velocity                              [m/s]
  HostReal2D cloud_frac_unused =
      workspace_scope.view<HostReal2D>("est_silhs_tndcy::cloud_frac_unused", ngrdcol, nzt);
  HostReal2D w_std_dev_unused =
      workspace_scope.view<HostReal2D>("est_silhs_tndcy::w_std_dev_unused", ngrdcol, nzt);
  HostReal3D rrm_auto_diag_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::rrm_auto_diag_all", ngrdcol, num_samples, nzt);
  HostReal3D rrm_accr_diag_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::rrm_accr_diag_all", ngrdcol, num_samples, nzt);
  HostReal3D rrm_evap_diag_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::rrm_evap_diag_all", ngrdcol, num_samples, nzt);
  HostReal3D Nrm_auto_diag_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::Nrm_auto_diag_all", ngrdcol, num_samples, nzt);
  HostReal3D Nrm_evap_diag_all = workspace_scope.view<HostReal3D>(
      "est_silhs_tndcy::Nrm_evap_diag_all", ngrdcol, num_samples, nzt);
  int ivar;
  int sample;
  int iiNr;
  int iirr;
  int iiPDF_chi;
  int iiPDF_w;

  // ---- Begin Code ----

  iiNr      = hm_metadata_iiNr;
  iirr      = hm_metadata_iirr;
  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_w   = hm_metadata_iiPDF_w;

  for (int index_2 = 0; index_2 < w_all_points.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < w_all_points.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < w_all_points.extent_int(0); ++index_0) {
        w_all_points(index_0, index_1, index_2) =
            static_cast<real_t>(X_nl_all_levs(index_0, index_1, index_2, iiPDF_w));
      }
    }
  }
  for (int index_2 = 0; index_2 < chi_all_points.extent_int(2); ++index_2) {
    for (int index_1 = 0; index_1 < chi_all_points.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < chi_all_points.extent_int(0); ++index_0) {
        chi_all_points(index_0, index_1, index_2) =
            static_cast<real_t>(X_nl_all_levs(index_0, index_1, index_2, iiPDF_chi));
      }
    }
  }

  copy_X_nl_into_hydromet_all_pts(nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs, hydromet_dim,
                                  hm_metadata, hydromet, hydromet_all_points, Ncn_all_points);

  for (int index_1 = 0; index_1 < cloud_frac_unused.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < cloud_frac_unused.extent_int(0); ++index_0) {
      cloud_frac_unused(index_0, index_1) = unused_var;
    }
  }
  for (int index_1 = 0; index_1 < w_std_dev_unused.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < w_std_dev_unused.extent_int(0); ++index_0) {
      w_std_dev_unused(index_0, index_1) = unused_var;
    }
  }

  for (sample = 0; sample < num_samples; sample += 1) {
    // Call the microphysics scheme to obtain a sample point
    microphys_sub(
        gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube,
        Kokkos::subview(lh_thl_clipped, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(w_all_points, Kokkos::ALL, sample, Kokkos::ALL), p_in_Pa, exner, rho,
        cloud_frac_unused, w_std_dev_unused, dzq,
        Kokkos::subview(lh_rc_clipped, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(lh_Nc_clipped, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(chi_all_points, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(lh_rv_clipped, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(hydromet_all_points, Kokkos::ALL, sample, Kokkos::ALL, Kokkos::ALL),
        saturation_formula,
        Kokkos::subview(lh_sample_point_weights, Kokkos::ALL, sample, Kokkos::ALL), stats,
        Kokkos::subview(lh_hydromet_mc_all, Kokkos::ALL, sample, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(lh_hydromet_vel_all, Kokkos::ALL, sample, Kokkos::ALL, Kokkos::ALL),
        Kokkos::subview(lh_Ncm_mc_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(lh_rcm_mc_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(lh_rvm_mc_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(lh_thlm_mc_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(rrm_auto_diag_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(rrm_accr_diag_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(rrm_evap_diag_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(Nrm_auto_diag_all, Kokkos::ALL, sample, Kokkos::ALL),
        Kokkos::subview(Nrm_evap_diag_all, Kokkos::ALL, sample, Kokkos::ALL));
  }

  if (microphys.l_var_covar_src) {
    HostReal2D lh_rtp2_mc_zt = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rtp2_mc_zt", ngrdcol,
        nzt); // LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
    HostReal2D lh_thlp2_mc_zt = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_thlp2_mc_zt", ngrdcol,
        nzt); // LH microphysics tendency for <thl'^2>                 [K^2/s]
    HostReal2D lh_wprtp_mc_zt = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_wprtp_mc_zt", ngrdcol,
        nzt); // LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
    HostReal2D lh_wpthlp_mc_zt = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_wpthlp_mc_zt", ngrdcol,
        nzt); // LH microphysics tendency for <w'thl'>                 [m*K/s^2]
    HostReal2D lh_rtpthlp_mc_zt = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rtpthlp_mc_zt", ngrdcol,
        nzt); // LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    lh_microphys_var_covar_driver_api(
        workspace, nzt, num_samples, ngrdcol, dt, sample_weights_execution, pdf_params,
        lh_rt_clipped, lh_thl_clipped, w_all_points, lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all,
        l_lh_instant_var_covar_src, lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, lh_wpthlp_mc_zt,
        lh_rtpthlp_mc_zt);

    // Convert from the zt grid to the zm grid.
    zt2zm_api(nzm, nzt, ngrdcol, gr, lh_rtp2_mc_zt, lh_rtp2_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, lh_thlp2_mc_zt, lh_thlp2_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, lh_wprtp_mc_zt, lh_wprtp_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, lh_wpthlp_mc_zt, lh_wpthlp_mc);
    zt2zm_api(nzm, nzt, ngrdcol, gr, lh_rtpthlp_mc_zt, lh_rtpthlp_mc);

    // This section has not been fully GPUized yet, but these variables have GPU copies.
    // Copy the GPU interpolation results back for the remaining CPU code.

    // Stats sampling for LH variance/covariance tendencies.
    if (stats.l_sample) {
      stats_update(std::string("lh_rtp2_mc"), lh_rtp2_mc, stats);
      stats_update(std::string("lh_thlp2_mc"), lh_thlp2_mc, stats);
      stats_update(std::string("lh_wprtp_mc"), lh_wprtp_mc, stats);
      stats_update(std::string("lh_wpthlp_mc"), lh_wpthlp_mc, stats);
      stats_update(std::string("lh_rtpthlp_mc"), lh_rtpthlp_mc, stats);
    }
  } else {
    for (int index_1 = 0; index_1 < lh_rtp2_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_rtp2_mc.extent_int(0); ++index_0) {
        lh_rtp2_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < lh_thlp2_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_thlp2_mc.extent_int(0); ++index_0) {
        lh_thlp2_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < lh_wprtp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_wprtp_mc.extent_int(0); ++index_0) {
        lh_wprtp_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < lh_wpthlp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_wpthlp_mc.extent_int(0); ++index_0) {
        lh_wpthlp_mc(index_0, index_1) = zero;
      }
    }
    for (int index_1 = 0; index_1 < lh_rtpthlp_mc.extent_int(1); ++index_1) {
      for (int index_0 = 0; index_0 < lh_rtpthlp_mc.extent_int(0); ++index_0) {
        lh_rtpthlp_mc(index_0, index_1) = zero;
      }
    }
  }

  // Grid box average.
  for (ivar = 0; ivar < hydromet_dim; ivar += 1) {
    compute_sample_mean(
        workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
        Kokkos::subview(lh_hydromet_vel_all, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, ivar),
        Kokkos::subview(lh_hydromet_vel, Kokkos::ALL, Kokkos::ALL, ivar));
    compute_sample_mean(
        workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
        Kokkos::subview(lh_hydromet_mc_all, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, ivar),
        Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, ivar));
  }

  compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution, lh_Ncm_mc_all,
                      lh_Ncm_mc);
  compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution, lh_rcm_mc_all,
                      lh_rcm_mc);
  compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution, lh_rvm_mc_all,
                      lh_rvm_mc);
  compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                      lh_thlm_mc_all, lh_thlm_mc);

  // Adjust the mean if l_silhs_KK_convergence_adj_mean is true
  if (microphys.l_silhs_KK_convergence_adj_mean) {
    // These means and adjustment diagnostics have no consumers outside this option.
    HostReal2D rrm_auto_diag_avg = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::rrm_auto_diag_avg", ngrdcol, nzt);
    HostReal2D rrm_accr_diag_avg = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::rrm_accr_diag_avg", ngrdcol, nzt);
    HostReal2D rrm_evap_diag_avg = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::rrm_evap_diag_avg", ngrdcol, nzt);
    HostReal2D Nrm_auto_diag_avg = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::Nrm_auto_diag_avg", ngrdcol, nzt);
    HostReal2D Nrm_evap_diag_avg = workspace_scope.view<HostReal2D>(
        Kokkos::WithoutInitializing, "est_silhs_tndcy::Nrm_evap_diag_avg", ngrdcol, nzt);
    HostReal2D lh_rrm_src_adj =
        workspace_scope.view<HostReal2D>("est_silhs_tndcy::lh_rrm_src_adj", ngrdcol, nzt);
    HostReal2D lh_Nrm_src_adj =
        workspace_scope.view<HostReal2D>("est_silhs_tndcy::lh_Nrm_src_adj", ngrdcol, nzt);
    HostReal2D lh_rrm_evap_adj =
        workspace_scope.view<HostReal2D>("est_silhs_tndcy::lh_rrm_evap_adj", ngrdcol, nzt);
    HostReal2D lh_Nrm_evap_adj =
        workspace_scope.view<HostReal2D>("est_silhs_tndcy::lh_Nrm_evap_adj", ngrdcol, nzt);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                        rrm_auto_diag_all, rrm_auto_diag_avg);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                        rrm_accr_diag_all, rrm_accr_diag_avg);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                        rrm_evap_diag_all, rrm_evap_diag_avg);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                        Nrm_auto_diag_all, Nrm_auto_diag_avg);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, sample_weights_execution,
                        Nrm_evap_diag_all, Nrm_evap_diag_avg);

    microphys_host::adjust_KK_src_means(
        workspace, dt, nzt, ngrdcol, exner, rcm,
        Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, iirr),
        Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, iiNr), hydromet, hydromet_dim,
        hm_metadata_iiri, rrm_auto_diag_avg, rrm_accr_diag_avg, rrm_evap_diag_avg,
        Nrm_auto_diag_avg, Nrm_evap_diag_avg,
        Kokkos::subview(lh_hydromet_vel, Kokkos::ALL, Kokkos::ALL, iirr),
        Kokkos::subview(lh_hydromet_vel, Kokkos::ALL, Kokkos::ALL, iiNr),
        Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, iirr),
        Kokkos::subview(lh_hydromet_mc, Kokkos::ALL, Kokkos::ALL, iiNr), lh_rvm_mc, lh_rcm_mc,
        lh_thlm_mc, lh_rrm_src_adj, lh_Nrm_src_adj, lh_rrm_evap_adj, lh_Nrm_evap_adj, kk);

    // Statistical sampling (moved here with the adjustment diagnostics).
    if (stats.l_sample) {
      stats_update(std::string("lh_rrm_src_adj"), lh_rrm_src_adj, stats);
      stats_update(std::string("lh_Nrm_src_adj"), lh_Nrm_src_adj, stats);
      stats_update(std::string("lh_rrm_evap_adj"), lh_rrm_evap_adj, stats);
      stats_update(std::string("lh_Nrm_evap_adj"), lh_Nrm_evap_adj, stats);
    }
  }

  if (stats.l_sample) {
    // Invoke the SILHS category variance sampler (if desired by user)!!
    if (var_on_stats_list(stats, std::string("silhs_var_cat_1"))) {
      silhs_category_variance_driver(workspace, ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim,
                                     hm_metadata, X_nl_all_levs, X_mixt_comp_all_levs,
                                     lh_hydromet_mc_all, sample_weights_execution, pdf_params,
                                     precip_fracs, stats);
    }
  }
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20>
inline void adjust_KK_src_means(ScratchWorkspace& workspace, real_t dt, int nzt, int ngrdcol,
                                Array0 exner, Array1 rcm, Array2 rrm, Array3 Nrm, Array4 hydromet,
                                int hydromet_dim, int iiri, Array5 rrm_auto, Array6 rrm_accr,
                                Array7 rrm_evap, Array8 Nrm_auto, Array9 Nrm_evap, Array10 lh_Vrr,
                                Array11 lh_VNr, Array12 rrm_mc, Array13 Nrm_mc, Array14 rvm_mc,
                                Array15 rcm_mc, Array16 thlm_mc, Array17 rrm_src_adj,
                                Array18 Nrm_src_adj, Array19 rrm_evap_adj, Array20 Nrm_evap_adj,
                                const parameters_KK& kk) {
  // Description:
  //   Adjusts the means of microphysics terms for KK microphysics by calling the
  //   KK microphysics adjustment subroutine for every model column.
  //
  // References:
  //   clubb:ticket:558
  //-----------------------------------------------------------------------------

  auto workspace_scope              = workspace.scope();
  constexpr bool l_src_adj_enabled  = true;
  constexpr bool l_evap_adj_enabled = true;
  // In: dt — Model timestep [s]
  // In: nzt — Number of thermodynamic vertical levels
  // In: ngrdcol — Number of model columns
  // In: exner (ngrdcol,nzt) — Exner function                            [-]
  // In: rcm (ngrdcol,nzt) — Mean liquid water mixing ratio            [kg/kg]
  // In: rrm (ngrdcol,nzt) — Rain water mixing ratio                   [kg/kg]
  // In: Nrm (ngrdcol,nzt) — Rain drop concentration                   [num/kg]
  // In: hydromet (ngrdcol,nzt,hydromet_dim) — Hydrometeor species    [units vary]
  // In: hydromet_dim
  // In: iiri
  // In: rrm_auto (ngrdcol,nzt) — Mean change in rain due to autoconversion [(kg/kg)/s]
  // In: rrm_accr (ngrdcol,nzt) — Mean change in rain due to accretion      [(kg/kg)/s]
  // In: rrm_evap (ngrdcol,nzt) — Mean change in rain due to evap           [(kg/kg)/s]
  // In: Nrm_auto (ngrdcol,nzt) — Mean change in Nrm due to autoconversion  [(num/kg)/s]
  // In: Nrm_evap (ngrdcol,nzt) — Mean change in Nrm due to evaporation     [(num/kg)/s]
  // Inout: lh_Vrr (ngrdcol,nzt) — Mean sedimentation velocity of < r_r > [m/s]
  // Inout: lh_VNr (ngrdcol,nzt) — Mean sedimentation velocity of < N_r > [m/s]
  // Out: rrm_mc (ngrdcol,nzt) — Mean change in rain due to microphysics [(kg/kg)/s]
  // Out: Nrm_mc (ngrdcol,nzt) — Mean change in Nrm due to microphysics  [(kg/kg)/s]
  // Out: rvm_mc (ngrdcol,nzt) — Time tendency of rvm                    [(kg/kg)/s]
  // Out: rcm_mc (ngrdcol,nzt) — Time tendency of rcm                    [(kg/kg)/s]
  // Out: thlm_mc (ngrdcol,nzt) — Time tendency of thlm                   [K/s]
  // Out: rrm_src_adj (ngrdcol,nzt)
  // Out: Nrm_src_adj (ngrdcol,nzt)
  // Out: rrm_evap_adj (ngrdcol,nzt)
  // Out: Nrm_evap_adj (ngrdcol,nzt)
  // Case state: kk replaces the Native parameters_KK module.

  HostAdj2D adj_terms("adjust_KK_src_means::adj_terms", ngrdcol,
                      nzt); // Adjustment terms returned from the adjustment routine
  HostInt1D cloud_top_level =
      workspace_scope.view<HostInt1D>("adjust_KK_src_means::cloud_top_level", ngrdcol);
  int i;
  int k;

  //-------------------------- Begin Code --------------------------

  // Initialize output
  for (int index_1 = 0; index_1 < rrm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_mc.extent_int(0); ++index_0) {
      rrm_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < Nrm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_mc.extent_int(0); ++index_0) {
      Nrm_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < rvm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rvm_mc.extent_int(0); ++index_0) {
      rvm_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < rcm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rcm_mc.extent_int(0); ++index_0) {
      rcm_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < thlm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < thlm_mc.extent_int(0); ++index_0) {
      thlm_mc(index_0, index_1) = zero;
    }
  }

  for (k = 0; k < nzt; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      // We call KK_microphys_adjust to adjust the means of the mc terms
      microphys_host::KK_microphys_adjust(
          dt, exner(i, k), rcm(i, k), rrm(i, k), Nrm(i, k), rrm_evap(i, k), rrm_auto(i, k),
          rrm_accr(i, k), Nrm_evap(i, k), Nrm_auto(i, k), l_src_adj_enabled, l_evap_adj_enabled,
          rrm_mc(i, k), Nrm_mc(i, k), rvm_mc(i, k), rcm_mc(i, k), thlm_mc(i, k), adj_terms(i, k),
          kk);
    }
  }

  // Clip positive values of Vrr and VNr
  for (k = 0; k < nzt - 1; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      if ((lh_Vrr(i, k) > zero)) {
        lh_Vrr(i, k) = zero;
      }
      if ((lh_VNr(i, k) > zero)) {
        lh_VNr(i, k) = zero;
      }
    }
  }

  for (int index_0 = 0; index_0 < cloud_top_level.extent_int(0); ++index_0) {
    cloud_top_level(index_0) = 1;
  }
  if ((iiri > -1)) {
    for (k = (nzt - 1); k >= 0; k += (-1)) {
      for (i = 0; i < ngrdcol; i += 1) {
        if ((cloud_top_level(i) == 1)) {
          if (((rcm(i, k) > rc_tol) || (hydromet(i, k, iiri) > ri_tol))) {
            cloud_top_level(i) = ((k + 1));
          }
        }
      }
    }
  } else {
    for (k = (nzt - 1); k >= 0; k += (-1)) {
      for (i = 0; i < ngrdcol; i += 1) {
        if (((cloud_top_level(i) == 1) && (rcm(i, k) > rc_tol))) {
          cloud_top_level(i) = ((k + 1));
        }
      }
    }
  }

  //!! Mean sedimentation above cloud top should have a value of 0.
  for (k = 1; k < nzt - 1; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      if ((((k + 1)) > cloud_top_level(i))) {
        lh_Vrr(i, k) = zero;
        lh_VNr(i, k) = zero;
      }
    }
  }

  // Set boundary conditions
  for (int index_0 = 0; index_0 < rrm_mc.extent_int(0); ++index_0) {
    rrm_mc(index_0, ((nzt - 1))) = zero;
  }
  for (int index_0 = 0; index_0 < Nrm_mc.extent_int(0); ++index_0) {
    Nrm_mc(index_0, ((nzt - 1))) = zero;
  }
  for (int index_0 = 0; index_0 < rvm_mc.extent_int(0); ++index_0) {
    rvm_mc(index_0, ((nzt - 1))) = zero;
  }
  for (int index_0 = 0; index_0 < rcm_mc.extent_int(0); ++index_0) {
    rcm_mc(index_0, ((nzt - 1))) = zero;
  }
  for (int index_0 = 0; index_0 < thlm_mc.extent_int(0); ++index_0) {
    thlm_mc(index_0, ((nzt - 1))) = zero;
  }

  for (int index_1 = 0; index_1 < rrm_src_adj.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_src_adj.extent_int(0); ++index_0) {
      rrm_src_adj(index_0, index_1) = adj_terms(index_0, index_1).rrm_src_adj;
    }
  }
  for (int index_1 = 0; index_1 < Nrm_src_adj.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_src_adj.extent_int(0); ++index_0) {
      Nrm_src_adj(index_0, index_1) = adj_terms(index_0, index_1).Nrm_src_adj;
    }
  }
  for (int index_1 = 0; index_1 < rrm_evap_adj.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rrm_evap_adj.extent_int(0); ++index_0) {
      rrm_evap_adj(index_0, index_1) = adj_terms(index_0, index_1).rrm_evap_adj;
    }
  }
  for (int index_1 = 0; index_1 < Nrm_evap_adj.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Nrm_evap_adj.extent_int(0); ++index_0) {
      Nrm_evap_adj(index_0, index_1) = adj_terms(index_0, index_1).Nrm_evap_adj;
    }
  }
}
} // namespace clubb_kokkos::microphys_host

namespace clubb_kokkos {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Microphysics>
inline void est_silhs_tndcy(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int nzt, int nzm,
    int num_samples, int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata,
    Array0 X_nl_all_levs, Array1 X_mixt_comp_all_levs, Array2 lh_sample_point_weights,
    const pdf_parameter& pdf_params, const precipitation_fractions& precip_fracs, Array3 p_in_Pa,
    Array4 exner, Array5 rho, Array6 dzq, Array7 hydromet, Array8 rcm, Array9 lh_rt_clipped,
    Array10 lh_thl_clipped, Array11 lh_rc_clipped, Array12 lh_rv_clipped, Array13 lh_Nc_clipped,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats,
    Array14 lh_hydromet_mc, Array15 lh_hydromet_vel, Array16 lh_Ncm_mc, Array17 lh_rvm_mc,
    Array18 lh_rcm_mc, Array19 lh_thlm_mc, Array20 lh_rtp2_mc, Array21 lh_thlp2_mc,
    Array22 lh_wprtp_mc, Array23 lh_wpthlp_mc, Array24 lh_rtpthlp_mc, Microphysics microphys_sub,
    const parameters_microphys& microphys, const parameters_KK& kk) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  auto workspace_scope              = workspace.scope();
  auto X_nl_all_levs_host           = to_host_space(workspace, X_nl_all_levs);
  auto X_mixt_comp_all_levs_host    = to_host_space(workspace, X_mixt_comp_all_levs);
  auto lh_sample_point_weights_host = to_host_space(workspace, lh_sample_point_weights);
  auto p_in_Pa_host                 = to_host_space(workspace, p_in_Pa);
  auto exner_host                   = to_host_space(workspace, exner);
  auto rho_host                     = to_host_space(workspace, rho);
  auto dzq_host                     = to_host_space(workspace, dzq);
  auto hydromet_host                = to_host_space(workspace, hydromet);
  auto rcm_host                     = to_host_space(workspace, rcm);
  auto lh_rt_clipped_host           = to_host_space(workspace, lh_rt_clipped);
  auto lh_thl_clipped_host          = to_host_space(workspace, lh_thl_clipped);
  auto lh_rc_clipped_host           = to_host_space(workspace, lh_rc_clipped);
  auto lh_rv_clipped_host           = to_host_space(workspace, lh_rv_clipped);
  auto lh_Nc_clipped_host           = to_host_space(workspace, lh_Nc_clipped);
  auto lh_hydromet_mc_host =
      workspace_scope.view<decltype(to_host_space(workspace, lh_hydromet_mc))>(
          Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_hydromet_mc_host",
          lh_hydromet_mc.extent(0), lh_hydromet_mc.extent(1), lh_hydromet_mc.extent(2));
  auto lh_hydromet_vel_host =
      workspace_scope.view<decltype(to_host_space(workspace, lh_hydromet_vel))>(
          Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_hydromet_vel_host",
          lh_hydromet_vel.extent(0), lh_hydromet_vel.extent(1), lh_hydromet_vel.extent(2));
  auto lh_Ncm_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_Ncm_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_Ncm_mc_host", lh_Ncm_mc.extent(0),
      lh_Ncm_mc.extent(1));
  auto lh_rvm_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_rvm_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rvm_mc_host", lh_rvm_mc.extent(0),
      lh_rvm_mc.extent(1));
  auto lh_rcm_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_rcm_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rcm_mc_host", lh_rcm_mc.extent(0),
      lh_rcm_mc.extent(1));
  auto lh_thlm_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_thlm_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_thlm_mc_host", lh_thlm_mc.extent(0),
      lh_thlm_mc.extent(1));
  auto lh_rtp2_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_rtp2_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rtp2_mc_host", lh_rtp2_mc.extent(0),
      lh_rtp2_mc.extent(1));
  auto lh_thlp2_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_thlp2_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_thlp2_mc_host", lh_thlp2_mc.extent(0),
      lh_thlp2_mc.extent(1));
  auto lh_wprtp_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_wprtp_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_wprtp_mc_host", lh_wprtp_mc.extent(0),
      lh_wprtp_mc.extent(1));
  auto lh_wpthlp_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_wpthlp_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_wpthlp_mc_host", lh_wpthlp_mc.extent(0),
      lh_wpthlp_mc.extent(1));
  auto lh_rtpthlp_mc_host = workspace_scope.view<decltype(to_host_space(workspace, lh_rtpthlp_mc))>(
      Kokkos::WithoutInitializing, "est_silhs_tndcy::lh_rtpthlp_mc_host", lh_rtpthlp_mc.extent(0),
      lh_rtpthlp_mc.extent(1));
  microphys_host::est_silhs_tndcy(
      workspace, gr, ngrdcol, dt, nzt, nzm, num_samples, pdf_dim, hydromet_dim, hm_metadata,
      X_nl_all_levs_host, X_mixt_comp_all_levs_host, lh_sample_point_weights_host, pdf_params,
      precip_fracs, p_in_Pa_host, exner_host, rho_host, dzq_host, hydromet_host, rcm_host,
      lh_rt_clipped_host, lh_thl_clipped_host, lh_rc_clipped_host, lh_rv_clipped_host,
      lh_Nc_clipped_host, l_lh_instant_var_covar_src, saturation_formula, stats,
      lh_hydromet_mc_host, lh_hydromet_vel_host, lh_Ncm_mc_host, lh_rvm_mc_host, lh_rcm_mc_host,
      lh_thlm_mc_host, lh_rtp2_mc_host, lh_thlp2_mc_host, lh_wprtp_mc_host, lh_wpthlp_mc_host,
      lh_rtpthlp_mc_host, microphys_sub, microphys, kk);
  from_execution_space(workspace, lh_hydromet_mc, lh_hydromet_mc_host);
  from_execution_space(workspace, lh_hydromet_vel, lh_hydromet_vel_host);
  from_execution_space(workspace, lh_Ncm_mc, lh_Ncm_mc_host);
  from_execution_space(workspace, lh_rvm_mc, lh_rvm_mc_host);
  from_execution_space(workspace, lh_rcm_mc, lh_rcm_mc_host);
  from_execution_space(workspace, lh_thlm_mc, lh_thlm_mc_host);
  from_execution_space(workspace, lh_rtp2_mc, lh_rtp2_mc_host);
  from_execution_space(workspace, lh_thlp2_mc, lh_thlp2_mc_host);
  from_execution_space(workspace, lh_wprtp_mc, lh_wprtp_mc_host);
  from_execution_space(workspace, lh_wpthlp_mc, lh_wpthlp_mc_host);
  from_execution_space(workspace, lh_rtpthlp_mc, lh_rtpthlp_mc_host);
}
} // namespace clubb_kokkos
