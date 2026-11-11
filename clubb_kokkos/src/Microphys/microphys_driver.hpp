#include "scratch_workspace.hpp"
#include "CLUBB_core/stats_clubb_utilities.hpp"
// Dispatch supported microphysics and combine its CLUBB source tendencies.
//
// calc_microphys_scheme_tendcies initializes outputs, selects none/KK/Morrison,
// coordinates SILHS and optional cloud sedimentation, then updates variance and
// covariance tendencies. Its Native argument prefix is unchanged; appended
// records own microphysics settings, KK parameters and Morrison engine state.
// Arrays are column-first LayoutLeft (column,level[,species]); SILHS samples add
// a sample dimension after column. Public host staging copies bits only and
// preserves output allocations. The same explicitly gated options as JAX are
// excluded locally at their Native branch; numerical engines are C++.
#pragma once
#include "Microphys/KK_microphys_module.hpp"
#include "Microphys/morrison_microphys_module.hpp"
#include "Microphys/lh_microphys_driver_module.hpp"
#include "Microphys/cloud_sed_module.hpp"
#include "CLUBB_core/advance_xp2_xpyp_module.hpp"
#include "Microphys/microphys_host_types.hpp"
#include "Microphys/KK_microphys/KK_microphys.hpp"
#include "Microphys/advance_microphys_module.hpp"
#include "Microphys/Morrison_microphys/module_mp_graupel.hpp"
namespace clubb_kokkos::microphys_host {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
inline void calc_microphys_scheme_tendcies(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, double time_current,
    int pdf_dim, int hydromet_dim, const std::string& runtype, Array0 thlm, Array1 p_in_Pa,
    Array2 exner, Array3 rho, Array4 rho_zm, Array5 rtm, Array6 rcm, Array7 cloud_frac,
    Array8 wm_zt, Array9 wm_zm, Array10 wp2, Array11 wp3, Array12 clubb_params, Array13 hydromet,
    Array14 Nc_in_cloud, const hm_metadata_type& hm_metadata, const pdf_parameter& pdf_params,
    const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array15 X_nl_all_levs,
    Array16 X_mixt_comp_all_levs, Array17 lh_sample_point_weights, Array18 mu_x_1_n,
    Array19 mu_x_2_n, Array20 sigma_x_1_n, Array21 sigma_x_2_n, Array22 corr_array_1_n,
    Array23 corr_array_2_n, Array24 lh_rt_clipped, Array25 lh_thl_clipped, Array26 lh_rc_clipped,
    Array27 lh_rv_clipped, Array28 lh_Nc_clipped, bool l_lh_importance_sampling,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats, Array29 Nccnm,
    Array30 hydromet_mc, Array31 Ncm_mc, Array32 rcm_mc, Array33 rvm_mc, Array34 thlm_mc,
    Array35 hydromet_vel_zt, Array36 hydromet_vel_covar_zt_impc, Array37 hydromet_vel_covar_zt_expc,
    Array38 wprtp_mc, Array39 wpthlp_mc, Array40 rtp2_mc, Array41 thlp2_mc, Array42 rtpthlp_mc,
    Array43 Skw_zm_smooth, const parameters_microphys& microphys, const parameters_KK& kk,
    morrison_microphysics_state& morrison);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
inline void calc_microphys_scheme_tendcies(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, double time_current,
    int pdf_dim, int hydromet_dim, const std::string& runtype, Array0 thlm, Array1 p_in_Pa,
    Array2 exner, Array3 rho, Array4 rho_zm, Array5 rtm, Array6 rcm, Array7 cloud_frac,
    Array8 wm_zt, Array9 wm_zm, Array10 wp2, Array11 wp3, Array12 clubb_params, Array13 hydromet,
    Array14 Nc_in_cloud, const hm_metadata_type& hm_metadata, const pdf_parameter& pdf_params,
    const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array15 X_nl_all_levs,
    Array16 X_mixt_comp_all_levs, Array17 lh_sample_point_weights, Array18 mu_x_1_n,
    Array19 mu_x_2_n, Array20 sigma_x_1_n, Array21 sigma_x_2_n, Array22 corr_array_1_n,
    Array23 corr_array_2_n, Array24 lh_rt_clipped, Array25 lh_thl_clipped, Array26 lh_rc_clipped,
    Array27 lh_rv_clipped, Array28 lh_Nc_clipped, bool l_lh_importance_sampling,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats, Array29 Nccnm,
    Array30 hydromet_mc, Array31 Ncm_mc, Array32 rcm_mc, Array33 rvm_mc, Array34 thlm_mc,
    Array35 hydromet_vel_zt, Array36 hydromet_vel_covar_zt_impc, Array37 hydromet_vel_covar_zt_expc,
    Array38 wprtp_mc, Array39 wpthlp_mc, Array40 rtp2_mc, Array41 thlp2_mc, Array42 rtpthlp_mc,
    Array43 Skw_zm_smooth, const parameters_microphys& microphys, const parameters_KK& kk,
    morrison_microphysics_state& morrison) {
  // Description:
  // Call a microphysics scheme and output microphysics tendencies for the
  // predictive variables.
  // References:
  // H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
  // moment microphysics scheme for application in cloud and
  // climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.
  //
  // Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
  // large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
  // 229-243, 2000.
  //-----------------------------------------------------------------------
  // Input Variables
  //------------------------ Input/Output Variables ------------------------
  // Note:
  // For COAMPS Nccnm is initialized and Nim & Ncm are computed within
  // subroutine adjtg.
  // Output Variables
  // Local Variables
  //Input aerosol mass concentration: the unit is 10^12 ug/m3.
  //For example, aeromass=2.25e-12 means that the aerosol mass
  //concentration is 2.25 ug/m3.
  //This value of aeromass was recommended by Huan Guo
  // Just to avoid typing hm_metadata%iixx everywhere
  // In: gr
  // In: ngrdcol — Number of grid columns [-]
  // In: dt — Timestep         [s]
  // In: time_current — Current time     [s]
  // In: pdf_dim — Number of variables in the correlation arrays
  // In: hydromet_dim — Number of hydrometeor species
  // In: runtype — Name of the run, for case specific effects.
  // In: thlm (ngrdcol,gr%nzt) — Liquid potential temp.                    [K]
  // In: p_in_Pa (ngrdcol,gr%nzt) — Pressure                                  [Pa]
  // In: exner (ngrdcol,gr%nzt) — Exner function                            [-]
  // In: rho (ngrdcol,gr%nzt) — Density on thermodynamic levels           [kg/m^3]
  // In: rho_zm (ngrdcol,gr%nzm) — Density on momentum levels                [kg/m^3]
  // In: rtm (ngrdcol,gr%nzt) — Total water mixing ratio                  [kg/kg]
  // In: rcm (ngrdcol,gr%nzt) — Liquid water mixing ratio                 [kg/kg]
  // In: cloud_frac (ngrdcol,gr%nzt) — Cloud fraction                            [-]
  // In: wm_zt (ngrdcol,gr%nzt) — w wind component on thermodynamic levels  [m/s]
  // In: wm_zm (ngrdcol,gr%nzm) — w wind component on momentum levels       [m/s]
  // In: wp2 (ngrdcol,gr%nzm) — w'^2 on the momentum grid                 [m^2/s^2]
  // In: wp3 (ngrdcol,gr%nzt) — w'^3 on the thermo. grid                  [m^3/s^3]
  // In: clubb_params (ngrdcol,nparams) — CLUBB's tunable parameters                 [units vary]
  // In: hydromet (ngrdcol,gr%nzt,hydromet_dim) — Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
  // In: Nc_in_cloud (ngrdcol,gr%nzt) — Mean (in-cloud) cloud droplet concentration      [num/kg]
  // In: hm_metadata
  // In: pdf_params — PDF parameters
  // In: hydromet_pdf_params (ngrdcol,gr%nzt) — PDF parameters
  // In: precip_fracs — Precipitation fractions      [-]
  // In: X_nl_all_levs (ngrdcol,lh_num_samples,gr%nzt,pdf_dim) — Normally and lognormally distributed hydrometeors and other variables
  // In: X_mixt_comp_all_levs (ngrdcol,lh_num_samples,gr%nzt) — Which mixture component the sample is in
  // In: lh_sample_point_weights (ngrdcol,lh_num_samples,gr%nzt) — Weights for cloud weighted sampling
  // In: mu_x_1_n (ngrdcol,gr%nzt,pdf_dim) — Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // In: mu_x_2_n (ngrdcol,gr%nzt,pdf_dim) — Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // In: sigma_x_1_n (ngrdcol,gr%nzt,pdf_dim) — Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // In: sigma_x_2_n (ngrdcol,gr%nzt,pdf_dim) — Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // In: corr_array_1_n (ngrdcol,gr%nzt,pdf_dim,pdf_dim) — Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // In: corr_array_2_n (ngrdcol,gr%nzt,pdf_dim,pdf_dim) — Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // In: lh_rt_clipped (ngrdcol,lh_num_samples,gr%nzt) — rt generated from silhs sample points
  // In: lh_thl_clipped (ngrdcol,lh_num_samples,gr%nzt) — thl generated from silhs sample points
  // In: lh_rc_clipped (ngrdcol,lh_num_samples,gr%nzt) — rc generated from silhs sample points
  // In: lh_rv_clipped (ngrdcol,lh_num_samples,gr%nzt) — rv generated from silhs sample points
  // In: lh_Nc_clipped (ngrdcol,lh_num_samples,gr%nzt) — Nc generated from silhs sample points
  // In: l_lh_importance_sampling — Do importance sampling (SILHS) [-]
  // In: l_lh_instant_var_covar_src — Produce instantaneous var/covar tendencies [-]
  // In: saturation_formula — Integer that stores the saturation formula to be used
  // Inout: stats
  // Inout: Nccnm (ngrdcol,gr%nzt) — Cloud condensation nuclei concentration (COAMPS)  [num/kg]
  // Out: hydromet_mc (ngrdcol,gr%nzt,hydromet_dim) — Change in hydrometeors due to microphysics  [units/s]
  // Out: Ncm_mc (ngrdcol,gr%nzt) — Change in Ncm due to microphysics  [num/kg/s]
  // Out: rcm_mc (ngrdcol,gr%nzt) — Microphysics contributions to liquid water          [kg/kg/s]
  // Out: rvm_mc (ngrdcol,gr%nzt) — Microphysics contributions to vapor water           [kg/kg/s]
  // Out: thlm_mc (ngrdcol,gr%nzt) — Microphysics contributions to theta-l               [K/s]
  // Out: hydromet_vel_zt (ngrdcol,gr%nzt,hydromet_dim) — Mean hydrometeor sed. velocity on thermo. levs. [m/s]
  // Out: hydromet_vel_covar_zt_impc (ngrdcol,gr%nzt,hydromet_dim) — Imp. comp. <V_xx'x_x'> t-levs [m/s]
  // Out: hydromet_vel_covar_zt_expc (ngrdcol,gr%nzt,hydromet_dim) — Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]
  // Out: wprtp_mc (ngrdcol,gr%nzm) — Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
  // Out: wpthlp_mc (ngrdcol,gr%nzm) — Microphysics tendency for <w'thl'>  [m*K/s^2]
  // Out: rtp2_mc (ngrdcol,gr%nzm) — Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
  // Out: thlp2_mc (ngrdcol,gr%nzm) — Microphysics tendency for <thl'^2>  [K^2/s]
  // Out: rtpthlp_mc (ngrdcol,gr%nzm) — Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
  // Out: Skw_zm_smooth (ngrdcol,gr%nzm) — Smoothed vertical velocity skewness [-]
  // Case state: microphys replaces the Native parameters_microphys module.
  // Case state: kk replaces the Native parameters_KK module.
  // Case state: morrison replaces the Native morrison_microphysics_state module.

  auto workspace_scope            = workspace.scope();
  const auto gr_nzm               = gr.nzm;
  const auto gr_nzt               = gr.nzt;
  const auto gr_dzm               = to_host_space(workspace, gr.dzm);
  const auto hm_metadata_iirr     = hm_metadata.iirr;
  const auto hm_metadata_iirs     = hm_metadata.iirs;
  const auto hm_metadata_iiri     = hm_metadata.iiri;
  const auto hm_metadata_iirg     = hm_metadata.iirg;
  const auto hm_metadata_iiNr     = hm_metadata.iiNr;
  const auto hm_metadata_iiNi     = hm_metadata.iiNi;
  const auto pdf_params_rt_1      = to_host_space(workspace, pdf_params.rt_1);
  const auto pdf_params_rt_2      = to_host_space(workspace, pdf_params.rt_2);
  const auto pdf_params_thl_1     = to_host_space(workspace, pdf_params.thl_1);
  const auto pdf_params_thl_2     = to_host_space(workspace, pdf_params.thl_2);
  const auto pdf_params_chi_1     = to_host_space(workspace, pdf_params.chi_1);
  const auto pdf_params_chi_2     = to_host_space(workspace, pdf_params.chi_2);
  const auto pdf_params_rc_1      = to_host_space(workspace, pdf_params.rc_1);
  const auto pdf_params_rc_2      = to_host_space(workspace, pdf_params.rc_2);
  const auto pdf_params_mixt_frac = to_host_space(workspace, pdf_params.mixt_frac);
  HostReal2D delta_zt             = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::delta_zt", ngrdcol,
      gr_nzt); // Difference in thermo. height levels     [m]
  HostReal2D wp2_zt = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::wp2_zt", ngrdcol,
      gr_nzt); // w'^2 on the thermo. grid [m^2/s^2]
  HostReal2D wp3_zm = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::wp3_zm", ngrdcol,
      gr_nzm); // w'^3 (momentum levels)                         [m^3/s^3]
  HostReal2D Skw_zm = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::Skw_zm", ngrdcol,
      gr_nzm); // Skewness of w on momentum levels  [-]
  HostReal2D T_in_K = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::T_in_K", ngrdcol,
      gr_nzt); // Temperature                                [K]
  HostReal2D rvm = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::rvm", ngrdcol,
      gr_nzt); // Vapor water mixing ratio                   [kg/kg]
  HostReal2D thlm_morr = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::thlm_morr", ngrdcol,
      gr_nzt); // Thlm fed into morrison microphysics        [K]
  HostReal2D Ncm_microphys = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::Ncm_microphys", ngrdcol,
      gr_nzt); // Mean cloud droplet concentration, <N_c>    [num/kg]
  HostReal2D wtmp = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::wtmp", ngrdcol,
      gr_nzt); // Standard dev. of w                   [m/s]
  HostReal2D chi = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::chi", ngrdcol,
      gr_nzt); // The variable 's' in Mellor (1977)    [kg/kg]
  HostReal2D rrm_auto_diag = workspace_scope.view<HostReal2D>(
      "calc_microphys_scheme_tendcies::rrm_auto_diag", ngrdcol, gr_nzt);
  HostReal2D rrm_accr_diag = workspace_scope.view<HostReal2D>(
      "calc_microphys_scheme_tendcies::rrm_accr_diag", ngrdcol, gr_nzt);
  HostReal2D rrm_evap_diag =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::rrm_evap_diag", ngrdcol,
                                       gr_nzt); // Evaporation rate for rrm         [kg/kg/s]
  HostReal2D Nrm_auto_diag = workspace_scope.view<HostReal2D>(
      "calc_microphys_scheme_tendcies::Nrm_auto_diag", ngrdcol, gr_nzt);
  HostReal2D Nrm_evap_diag = workspace_scope.view<HostReal2D>(
      "calc_microphys_scheme_tendcies::Nrm_evap_diag", ngrdcol, gr_nzt);
  int i; // Zero-based column index
  int k; // Zero-based level index
  HostReal2D lh_AKm =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::lh_AKm", ngrdcol,
                                       gr_nzt); // Kessler ac estimate                 [kg/kg/s]
  HostReal2D AKm =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::AKm", ngrdcol,
                                       gr_nzt); // Exact Kessler ac                    [kg/kg/s]
  HostReal2D AKstd =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::AKstd", ngrdcol,
                                       gr_nzt); // St dev of exact Kessler ac          [kg/kg/s]
  HostReal2D AKstd_cld =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::AKstd_cld", ngrdcol,
                                       gr_nzt); // Stdev of exact w/in cloud ac        [kg/kg/s]
  HostReal2D lh_rcm_avg =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::lh_rcm_avg", ngrdcol,
                                       gr_nzt); // Monte Carlo rcm estimate            [kg/kg]
  HostReal2D AKm_rcm =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::AKm_rcm", ngrdcol,
                                       gr_nzt); // Kessler ac based on rcm             [kg/kg/s]
  HostReal2D AKm_rcc =
      workspace_scope.view<HostReal2D>("calc_microphys_scheme_tendcies::AKm_rcc", ngrdcol,
                                       gr_nzt); // Kessler ac based on rcm/cloud_frac  [kg/kg/s]
  bool l_latin_hypercube_input;
  HostReal2D unit_sample_weight_2d = workspace_scope.view<HostReal2D>(
      Kokkos::WithoutInitializing, "calc_microphys_scheme_tendcies::unit_sample_weight_2d", ngrdcol,
      gr_nzt);
  int iirr;
  int iirs;
  int iiri;
  int iirg;
  int iiNr;
  int iiNi;

  //------------------------ Begin Code ------------------------

  iirr = hm_metadata_iirr;
  iirs = hm_metadata_iirs;
  iiri = hm_metadata_iiri;
  iirg = hm_metadata_iirg;
  iiNr = hm_metadata_iiNr;
  iiNi = hm_metadata_iiNi;
  for (int index_1 = 0; index_1 < unit_sample_weight_2d.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < unit_sample_weight_2d.extent_int(0); ++index_0) {
      unit_sample_weight_2d(index_0, index_1) = one;
    }
  }

  // We must initialize intent(out) variables. If we do not, and they do not
  // otherwise get assigned (e.g. microphys_scheme=="none"), then we are not
  // within the realms of Fortran standards compliance.
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
  for (int index_1 = 0; index_1 < wprtp_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < wprtp_mc.extent_int(0); ++index_0) {
      wprtp_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < wpthlp_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < wpthlp_mc.extent_int(0); ++index_0) {
      wpthlp_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < rtp2_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rtp2_mc.extent_int(0); ++index_0) {
      rtp2_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < thlp2_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < thlp2_mc.extent_int(0); ++index_0) {
      thlp2_mc(index_0, index_1) = zero;
    }
  }
  for (int index_1 = 0; index_1 < rtpthlp_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < rtpthlp_mc.extent_int(0); ++index_0) {
      rtpthlp_mc(index_0, index_1) = zero;
    }
  }

  // Calculate Skw_zm for use in advance_microphys.
  zt2zm_api(workspace, gr_nzm, gr_nzt, ngrdcol, gr, wp3, wp3_zm);

  {
    auto staging_scope = workspace.scope();
    auto field_2       = to_execution_space(workspace, wp2, true);
    auto field_3       = to_execution_space(workspace, wp3_zm, true);
    auto field_5       = to_execution_space(workspace, clubb_params, true);
    auto field_6       = to_execution_space(workspace, Kokkos::WithoutInitializing, Skw_zm);
    Skx_func(gr_nzm, ngrdcol, field_2, field_3, w_tol, field_5, field_6);
    from_execution_space(workspace, Skw_zm, field_6);
  }

  // This field is smoothed by interpolating to thermodynamic levels and then
  // interpolating back to momentum levels.
  zm2zt2zm(workspace, gr_nzm, gr_nzt, ngrdcol, gr, Skw_zm, Skw_zm_smooth);
  zm2zt_api(workspace, gr_nzm, gr_nzt, ngrdcol, gr, wp2, wp2_zt, w_tol_sqd);

  // Initialize predicitive hydrometeor tendencies and sedimentation
  // velocities.
  if ((hydromet_dim > 0)) {
    for (int index_2 = 0; index_2 < hydromet_mc.extent_int(2); ++index_2) {
      for (int index_1 = 0; index_1 < hydromet_mc.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_mc.extent_int(0); ++index_0) {
          hydromet_mc(index_0, index_1, index_2) = zero;
        }
      }
    }
    for (int index_2 = 0; index_2 < hydromet_vel_zt.extent_int(2); ++index_2) {
      for (int index_1 = 0; index_1 < hydromet_vel_zt.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_vel_zt.extent_int(0); ++index_0) {
          hydromet_vel_zt(index_0, index_1, index_2) = zero;
        }
      }
    }
  }
  // Initialize predicitive Ncm tendency.
  for (int index_1 = 0; index_1 < Ncm_mc.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_mc.extent_int(0); ++index_0) {
      Ncm_mc(index_0, index_1) = zero;
    }
  }

  // Initialize the values of the implicit and explicit components used to
  // calculate the covariances of hydrometeor sedimentation velocities and
  // their associated hydrometeors (for example, <V_rr'r_r'> and <V_Nr'N_r'>)
  // to 0.
  if ((hydromet_dim > 0)) {
    for (int index_2 = 0; index_2 < hydromet_vel_covar_zt_impc.extent_int(2); ++index_2) {
      for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_impc.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_impc.extent_int(0); ++index_0) {
          hydromet_vel_covar_zt_impc(index_0, index_1, index_2) = zero;
        }
      }
    }
    for (int index_2 = 0; index_2 < hydromet_vel_covar_zt_expc.extent_int(2); ++index_2) {
      for (int index_1 = 0; index_1 < hydromet_vel_covar_zt_expc.extent_int(1); ++index_1) {
        for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_expc.extent_int(0); ++index_0) {
          hydromet_vel_covar_zt_expc(index_0, index_1, index_2) = zero;
        }
      }
    }
  }

  // Return if there is delay between the model start time and start of the
  // microphysics
  if ((time_current < microphys.microphys_start_time)) {
    return;
  }

  // Make some compiler warnings go away for external users
  if ((runtype == std::string(""))) {
    throw std::runtime_error(std::string("Runtype is null, which should not happen."));
  }

  // Calculate the updated value of mean cloud droplet concentration, based on
  // Nc_in_cloud and the updated value of cloud_frac.
  for (int index_1 = 0; index_1 < Ncm_microphys.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < Ncm_microphys.extent_int(0); ++index_0) {
      Ncm_microphys(index_0, index_1) =
          (Nc_in_cloud(index_0, index_1) * cloud_frac(index_0, index_1));
    }
  }

  // Determine 's' from Mellor (1977)
  for (int index_1 = 0; index_1 < chi.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < chi.extent_int(0); ++index_0) {
      chi(index_0, index_1) =
          ((pdf_params_mixt_frac(index_0, index_1) * pdf_params_chi_1(index_0, index_1)) +
           ((one - pdf_params_mixt_frac(index_0, index_1)) * pdf_params_chi_2(index_0, index_1)));
    }
  }

  // Compute standard deviation of vertical velocity in the grid column
  for (int index_1 = 0; index_1 < wtmp.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < wtmp.extent_int(0); ++index_0) {
      wtmp(index_0, index_1) = Kokkos::sqrt(wp2_zt(index_0, index_1));
    }
  }

  // Compute difference in thermodynamic height levels.
  // The variable delta_zt is on the Morrison Microphysics grid, which
  // identifies delta_zt at level k as being the difference between zt(k+1) and
  // zt(k). On CLUBB's grid, the difference between zt(k+1) and zt(k) is
  // identified as gr%dzm at level k+1 (since momentum level k+1 is between
  // thermodynamic levels k and k+1).
  for (k = 1; k < gr_nzm; k += 1) {
    for (i = 0; i < ngrdcol; i += 1) {
      delta_zt(i, (k - 1)) = gr_dzm(i, k);
    }
  }

  // Calculate T_in_K
  for (int index_1 = 0; index_1 < T_in_K.extent_int(1); ++index_1) {
    for (int index_0 = 0; index_0 < T_in_K.extent_int(0); ++index_0) {
      T_in_K(index_0, index_1) =
          thlm2T_in_K_api(thlm(index_0, index_1), exner(index_0, index_1), rcm(index_0, index_1));
    }
  }

  // Begin by calling Brian Griffin's implementation of the
  // Khairoutdinov and Kogan microphysics (analytic or local formulas),
  // the COAMPS implementation of Rutlege and Hobbes, or the Morrison
  // microphysics.
  // Note: COAMPS appears to have some K&K elements to it as well.

  if (microphys_detail::trim(microphys.microphys_scheme) == std::string("coamps")) {
    // TODO: This scheme is excluded by the actual JAX initialization gate.
    throw std::runtime_error(
        std::string("Requested scheme is outside the supported JAX/Kokkos options"));
  } else if (microphys_detail::trim(microphys.microphys_scheme) == std::string("morrison")) {

    if ((microphys.lh_microphys_type != lh_microphys_disabled)) {
      // Procedure
      microphys_host::lh_microphys_driver(
          workspace, gr, ngrdcol, dt, gr_nzt, gr_nzm, microphys.lh_num_samples, pdf_dim,
          hydromet_dim, hm_metadata, X_nl_all_levs, lh_sample_point_weights, pdf_params,
          precip_fracs, p_in_Pa, exner, rho, rcm, delta_zt, cloud_frac, hydromet,
          X_mixt_comp_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped,
          lh_Nc_clipped, l_lh_importance_sampling, l_lh_instant_var_covar_src, saturation_formula,
          stats, hydromet_mc, hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc, rtp2_mc, thlp2_mc,
          wprtp_mc, wpthlp_mc, rtpthlp_mc, lh_AKm, AKm, AKstd, AKstd_cld, lh_rcm_avg, AKm_rcm,
          AKm_rcc,
          [&](auto&&... values) {
            microphys_host::morrison_microphys_driver(
                workspace, std::forward<decltype(values)>(values)..., microphys, morrison);
          },
          microphys, kk);

      stats_accumulate_lh_tend(workspace, gr, ngrdcol, hydromet_dim, hm_metadata, hydromet_mc,
                               Ncm_mc, thlm_mc, rvm_mc, rcm_mc, lh_AKm, AKm, AKstd, AKstd_cld,
                               lh_rcm_avg, AKm_rcm, AKm_rcc, stats);
    }

    // Call the microphysics if we don't want to have feedback effects from the
    // latin hypercube result (above)
    if ((microphys.lh_microphys_type != lh_microphys_interactive)) {
      l_latin_hypercube_input = false;

      for (k = 0; k < gr_nzt; k += 1) {
        for (i = 0; i < ngrdcol; i += 1) {
          if (l_morr_xp2_mc) {
            //Use the moister rt_1/rt_2 rather than rtm in morrison microphys
            if ((pdf_params_rt_1(i, k) > pdf_params_rt_2(i, k))) {
              rvm(i, k) = (pdf_params_rt_1(i, k) - pdf_params_rc_1(i, k));
            } else {
              rvm(i, k) = (pdf_params_rt_2(i, k) - pdf_params_rc_2(i, k));
            }

            //Also use the colder of thl_1/thl_2
            if ((pdf_params_thl_1(i, k) < pdf_params_thl_2(i, k))) {
              thlm_morr(i, k) = pdf_params_thl_1(i, k);
            } else {
              thlm_morr(i, k) = pdf_params_thl_2(i, k);
            }
          } else {
            rvm(i, k)       = (rtm(i, k) - rcm(i, k));
            thlm_morr(i, k) = thlm(i, k);
          }
        }
      }

      if ((!l_morr_xp2_mc)) {
        // Set microphysics tendencies for model variances to 0.
        for (k = 0; k < gr_nzm; k += 1) {
          for (i = 0; i < ngrdcol; i += 1) {
            rtp2_mc(i, k)    = zero;
            thlp2_mc(i, k)   = zero;
            wprtp_mc(i, k)   = zero;
            wpthlp_mc(i, k)  = zero;
            rtpthlp_mc(i, k) = zero;
          }
        }
      }

      microphys_host::morrison_microphys_driver(
          workspace, gr, ngrdcol, dt, gr_nzt, hydromet_dim, hm_metadata, l_latin_hypercube_input,
          thlm_morr, wm_zt, p_in_Pa, exner, rho, cloud_frac, wtmp, delta_zt, rcm, Ncm_microphys,
          chi, rvm, hydromet, saturation_formula, unit_sample_weight_2d, stats, hydromet_mc,
          hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc, rrm_auto_diag, rrm_accr_diag,
          rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag, microphys, morrison);

      if (l_morr_xp2_mc) {
        {
          auto staging_scope = workspace.scope();
          auto field_5       = to_execution_space(workspace, cloud_frac, true);
          auto field_6       = to_execution_space(workspace, rcm, true);
          auto field_7       = to_execution_space(workspace, rvm, true);
          auto field_8       = to_execution_space(workspace, thlm_morr, true);
          auto field_9       = to_execution_space(workspace, wm_zt, true);
          auto field_10      = to_execution_space(workspace, exner, true);
          auto field_11      = to_execution_space(workspace, rrm_evap_diag, true);
          auto field_13      = to_execution_space(workspace, rtp2_mc, true);
          auto field_14      = to_execution_space(workspace, thlp2_mc, true);
          auto field_15      = to_execution_space(workspace, wprtp_mc, true);
          auto field_16      = to_execution_space(workspace, wpthlp_mc, true);
          auto field_17      = to_execution_space(workspace, rtpthlp_mc, true);
          update_xp2_mc(workspace, gr, gr_nzm, gr_nzt, ngrdcol, dt, field_5, field_6, field_7,
                        field_8, field_9, field_10, field_11, pdf_params, field_13, field_14,
                        field_15, field_16, field_17);
          from_execution_space(workspace, rtp2_mc, field_13);
          from_execution_space(workspace, thlp2_mc, field_14);
          from_execution_space(workspace, wprtp_mc, field_15);
          from_execution_space(workspace, wpthlp_mc, field_16);
          from_execution_space(workspace, rtpthlp_mc, field_17);
        }
      }
    }

    if (stats.l_sample) {
      // Output rain sedimentation velocity

      const auto array_actual_1 =
          zt2zm_api(workspace, gr_nzm, gr_nzt, ngrdcol, gr,
                    Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, iirr));
      stats_update(std::string("Vrr"), array_actual_1, stats);
    }

  } else if (microphys_detail::trim(microphys.microphys_scheme) ==
             std::string("khairoutdinov_kogan")) {

    if ((microphys.lh_microphys_type != lh_microphys_disabled)) {

      // Procedure
      microphys_host::lh_microphys_driver(
          workspace, gr, ngrdcol, dt, gr_nzt, gr_nzm, microphys.lh_num_samples, pdf_dim,
          hydromet_dim, hm_metadata, X_nl_all_levs, lh_sample_point_weights, pdf_params,
          precip_fracs, p_in_Pa, exner, rho, rcm, delta_zt, cloud_frac, hydromet,
          X_mixt_comp_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, lh_rv_clipped,
          lh_Nc_clipped, l_lh_importance_sampling, l_lh_instant_var_covar_src, saturation_formula,
          stats, hydromet_mc, hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc, rtp2_mc, thlp2_mc,
          wprtp_mc, wpthlp_mc, rtpthlp_mc, lh_AKm, AKm, AKstd, AKstd_cld, lh_rcm_avg, AKm_rcm,
          AKm_rcc,
          [&](auto&&... values) {
            microphys_host::KK_local_microphys_driver(
                workspace, std::forward<decltype(values)>(values)..., microphys, kk);
          },
          microphys, kk);

      stats_accumulate_lh_tend(workspace, gr, ngrdcol, hydromet_dim, hm_metadata, hydromet_mc,
                               Ncm_mc, thlm_mc, rvm_mc, rcm_mc, lh_AKm, AKm, AKstd, AKstd_cld,
                               lh_rcm_avg, AKm_rcm, AKm_rcc, stats);

      if (stats.l_sample) {
        // Latin hypercube estimate for sedimentation velocities
        stats_update(std::string("lh_Vrr"),
                     Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, iirr), stats);
        stats_update(std::string("lh_VNr"),
                     Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, iiNr), stats);
      }
    }

    // Call the microphysics if we don't want to have feedback effects from the
    // latin hypercube result (above)
    if ((microphys.lh_microphys_type != lh_microphys_interactive)) {

      l_latin_hypercube_input = false;
      for (k = 0; k < gr_nzt; k += 1) {
        for (i = 0; i < ngrdcol; i += 1) {
          rvm(i, k) = (rtm(i, k) - rcm(i, k));
        }
      }

      if (microphys.l_local_kk) {
        microphys_host::KK_local_microphys_driver(
            workspace, gr, ngrdcol, dt, gr_nzt, hydromet_dim, hm_metadata, l_latin_hypercube_input,
            thlm, wm_zt, p_in_Pa, exner, rho, cloud_frac, wtmp, delta_zt, rcm, Ncm_microphys, chi,
            rvm, hydromet, saturation_formula, unit_sample_weight_2d, stats, hydromet_mc,
            hydromet_vel_zt, Ncm_mc, rcm_mc, rvm_mc, thlm_mc, rrm_auto_diag, rrm_accr_diag,
            rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag, microphys, kk);
      } else {
        microphys_host::KK_upscaled_microphys(
            workspace, gr, ngrdcol, dt, gr_nzt, gr_nzm, pdf_dim, hydromet_dim, hm_metadata, wm_zt,
            rtm, thlm, p_in_Pa, exner, rho, rcm, pdf_params, hydromet_pdf_params, precip_fracs,
            hydromet, mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n, corr_array_1_n, corr_array_2_n,
            saturation_formula, stats, hydromet_mc, hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc,
            hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, wprtp_mc, wpthlp_mc, rtp2_mc,
            thlp2_mc, rtpthlp_mc, microphys, kk);
      }

      if (((!microphys.l_local_kk) && microphys.l_silhs_KK_convergence_adj_mean)) {
        // SILHS cannot currently predict the turbulent sedimentation
        // terms, so zero them out for a convergence test.
        for (k = 0; k < gr_nzt; k += 1) {
          for (i = 0; i < ngrdcol; i += 1) {
            for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_impc.extent_int(2); ++index_0) {
              hydromet_vel_covar_zt_impc(i, k, index_0) = zero;
            }
            for (int index_0 = 0; index_0 < hydromet_vel_covar_zt_expc.extent_int(2); ++index_0) {
              hydromet_vel_covar_zt_expc(i, k, index_0) = zero;
            }
          }
        }
      }
    }

    if (stats.l_sample) {
      // Output rain sedimentation velocity

      const auto array_actual_2 =
          zt2zm_api(workspace, gr_nzm, gr_nzt, ngrdcol, gr,
                    Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, iirr));
      stats_update(std::string("Vrr"), array_actual_2, stats);
      // Output rain droplet number sedimentation velocity

      const auto array_actual_3 =
          zt2zm_api(workspace, gr_nzm, gr_nzt, ngrdcol, gr,
                    Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, iiNr));
      stats_update(std::string("VNr"), array_actual_3, stats);
    }

  } else if (microphys_detail::trim(microphys.microphys_scheme) == std::string("simplified_ice")) {
    // TODO: This scheme is excluded by the actual JAX initialization gate.
    throw std::runtime_error(
        std::string("Requested scheme is outside the supported JAX/Kokkos options"));
  } else {
    // Do nothing
    ;

    // microphys_scheme
  }

  // Sample microphysics variables if necessary
  if (stats.l_sample) {
    stats_update(std::string("Nccnm"), Nccnm, stats);
  }

  // Call GFDL activation code
  if (microphys.l_gfdl_activation) {
    // TODO: GFDL activation is excluded by the actual JAX initialization gate.
    throw std::runtime_error(
        std::string("GFDL activation is outside the supported JAX/Kokkos options"));
  }

  // Cloud water sedimentation.
  // Note:  it would be very easy to upscale the cloud water sedimentation
  //        flux, so we should look into adding an upscaled option.
  if (microphys.l_cloud_sed) {
    {
      auto staging_scope = workspace.scope();
      auto field_2       = to_execution_space(workspace, rcm, true);
      auto field_3       = to_execution_space(workspace, Ncm_microphys, true);
      auto field_4       = to_execution_space(workspace, rho_zm, true);
      auto field_5       = to_execution_space(workspace, rho, true);
      auto field_6       = to_execution_space(workspace, exner, true);
      auto field_9       = to_execution_space(workspace, rcm_mc, true);
      auto field_10      = to_execution_space(workspace, thlm_mc, true);
      cloud_drop_sed(workspace, gr, ngrdcol, field_2, field_3, field_4, field_5, field_6,
                     microphys.sigma_g, stats, field_9, field_10);
      from_execution_space(workspace, rcm_mc, field_9);
      from_execution_space(workspace, thlm_mc, field_10);
    }
  }

  return;
}
} // namespace clubb_kokkos::microphys_host

namespace clubb_kokkos {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
inline void calc_microphys_scheme_tendcies(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, double time_current,
    int pdf_dim, int hydromet_dim, const std::string& runtype, Array0 thlm, Array1 p_in_Pa,
    Array2 exner, Array3 rho, Array4 rho_zm, Array5 rtm, Array6 rcm, Array7 cloud_frac,
    Array8 wm_zt, Array9 wm_zm, Array10 wp2, Array11 wp3, Array12 clubb_params, Array13 hydromet,
    Array14 Nc_in_cloud, const hm_metadata_type& hm_metadata, const pdf_parameter& pdf_params,
    const hydromet_pdf_parameter_collection& hydromet_pdf_params,
    const precipitation_fractions& precip_fracs, Array15 X_nl_all_levs,
    Array16 X_mixt_comp_all_levs, Array17 lh_sample_point_weights, Array18 mu_x_1_n,
    Array19 mu_x_2_n, Array20 sigma_x_1_n, Array21 sigma_x_2_n, Array22 corr_array_1_n,
    Array23 corr_array_2_n, Array24 lh_rt_clipped, Array25 lh_thl_clipped, Array26 lh_rc_clipped,
    Array27 lh_rv_clipped, Array28 lh_Nc_clipped, bool l_lh_importance_sampling,
    bool l_lh_instant_var_covar_src, int saturation_formula, stats_type& stats, Array29 Nccnm,
    Array30 hydromet_mc, Array31 Ncm_mc, Array32 rcm_mc, Array33 rvm_mc, Array34 thlm_mc,
    Array35 hydromet_vel_zt, Array36 hydromet_vel_covar_zt_impc, Array37 hydromet_vel_covar_zt_expc,
    Array38 wprtp_mc, Array39 wpthlp_mc, Array40 rtp2_mc, Array41 thlp2_mc, Array42 rtpthlp_mc,
    Array43 Skw_zm_smooth, const parameters_microphys& microphys, const parameters_KK& kk,
    morrison_microphysics_state& morrison) {
  // Stage the Native engine call boundary into call-local host Views.
  // The C++ science implementation above retains Native float casts/order;
  // model arrays and every InOut/Out allocation remain caller owned.
  auto workspace_scope                 = workspace.scope();
  auto thlm_host                       = to_host_space(workspace, thlm);
  auto p_in_Pa_host                    = to_host_space(workspace, p_in_Pa);
  auto exner_host                      = to_host_space(workspace, exner);
  auto rho_host                        = to_host_space(workspace, rho);
  auto rho_zm_host                     = to_host_space(workspace, rho_zm);
  auto rtm_host                        = to_host_space(workspace, rtm);
  auto rcm_host                        = to_host_space(workspace, rcm);
  auto cloud_frac_host                 = to_host_space(workspace, cloud_frac);
  auto wm_zt_host                      = to_host_space(workspace, wm_zt);
  auto wm_zm_host                      = to_host_space(workspace, wm_zm);
  auto wp2_host                        = to_host_space(workspace, wp2);
  auto wp3_host                        = to_host_space(workspace, wp3);
  auto clubb_params_host               = to_host_space(workspace, clubb_params);
  auto hydromet_host                   = to_host_space(workspace, hydromet);
  auto Nc_in_cloud_host                = to_host_space(workspace, Nc_in_cloud);
  auto X_nl_all_levs_host              = to_host_space(workspace, X_nl_all_levs);
  auto X_mixt_comp_all_levs_host       = to_host_space(workspace, X_mixt_comp_all_levs);
  auto lh_sample_point_weights_host    = to_host_space(workspace, lh_sample_point_weights);
  auto mu_x_1_n_host                   = to_host_space(workspace, mu_x_1_n);
  auto mu_x_2_n_host                   = to_host_space(workspace, mu_x_2_n);
  auto sigma_x_1_n_host                = to_host_space(workspace, sigma_x_1_n);
  auto sigma_x_2_n_host                = to_host_space(workspace, sigma_x_2_n);
  auto corr_array_1_n_host             = to_host_space(workspace, corr_array_1_n);
  auto corr_array_2_n_host             = to_host_space(workspace, corr_array_2_n);
  auto lh_rt_clipped_host              = to_host_space(workspace, lh_rt_clipped);
  auto lh_thl_clipped_host             = to_host_space(workspace, lh_thl_clipped);
  auto lh_rc_clipped_host              = to_host_space(workspace, lh_rc_clipped);
  auto lh_rv_clipped_host              = to_host_space(workspace, lh_rv_clipped);
  auto lh_Nc_clipped_host              = to_host_space(workspace, lh_Nc_clipped);
  auto Nccnm_host                      = to_host_space(workspace, Nccnm);
  auto hydromet_mc_host                = to_host_space(workspace, hydromet_mc);
  auto Ncm_mc_host                     = to_host_space(workspace, Ncm_mc);
  auto rcm_mc_host                     = to_host_space(workspace, rcm_mc);
  auto rvm_mc_host                     = to_host_space(workspace, rvm_mc);
  auto thlm_mc_host                    = to_host_space(workspace, thlm_mc);
  auto hydromet_vel_zt_host            = to_host_space(workspace, hydromet_vel_zt);
  auto hydromet_vel_covar_zt_impc_host = to_host_space(workspace, hydromet_vel_covar_zt_impc);
  auto hydromet_vel_covar_zt_expc_host = to_host_space(workspace, hydromet_vel_covar_zt_expc);
  auto wprtp_mc_host                   = to_host_space(workspace, wprtp_mc);
  auto wpthlp_mc_host                  = to_host_space(workspace, wpthlp_mc);
  auto rtp2_mc_host                    = to_host_space(workspace, rtp2_mc);
  auto thlp2_mc_host                   = to_host_space(workspace, thlp2_mc);
  auto rtpthlp_mc_host                 = to_host_space(workspace, rtpthlp_mc);
  auto Skw_zm_smooth_host              = to_host_space(workspace, Skw_zm_smooth);
  microphys_host::calc_microphys_scheme_tendcies(
      workspace, gr, ngrdcol, dt, time_current, pdf_dim, hydromet_dim, runtype, thlm_host,
      p_in_Pa_host, exner_host, rho_host, rho_zm_host, rtm_host, rcm_host, cloud_frac_host,
      wm_zt_host, wm_zm_host, wp2_host, wp3_host, clubb_params_host, hydromet_host,
      Nc_in_cloud_host, hm_metadata, pdf_params, hydromet_pdf_params, precip_fracs,
      X_nl_all_levs_host, X_mixt_comp_all_levs_host, lh_sample_point_weights_host, mu_x_1_n_host,
      mu_x_2_n_host, sigma_x_1_n_host, sigma_x_2_n_host, corr_array_1_n_host, corr_array_2_n_host,
      lh_rt_clipped_host, lh_thl_clipped_host, lh_rc_clipped_host, lh_rv_clipped_host,
      lh_Nc_clipped_host, l_lh_importance_sampling, l_lh_instant_var_covar_src, saturation_formula,
      stats, Nccnm_host, hydromet_mc_host, Ncm_mc_host, rcm_mc_host, rvm_mc_host, thlm_mc_host,
      hydromet_vel_zt_host, hydromet_vel_covar_zt_impc_host, hydromet_vel_covar_zt_expc_host,
      wprtp_mc_host, wpthlp_mc_host, rtp2_mc_host, thlp2_mc_host, rtpthlp_mc_host,
      Skw_zm_smooth_host, microphys, kk, morrison);
  from_execution_space(workspace, Nccnm, Nccnm_host);
  from_execution_space(workspace, hydromet_mc, hydromet_mc_host);
  from_execution_space(workspace, Ncm_mc, Ncm_mc_host);
  from_execution_space(workspace, rcm_mc, rcm_mc_host);
  from_execution_space(workspace, rvm_mc, rvm_mc_host);
  from_execution_space(workspace, thlm_mc, thlm_mc_host);
  from_execution_space(workspace, hydromet_vel_zt, hydromet_vel_zt_host);
  from_execution_space(workspace, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_impc_host);
  from_execution_space(workspace, hydromet_vel_covar_zt_expc, hydromet_vel_covar_zt_expc_host);
  from_execution_space(workspace, wprtp_mc, wprtp_mc_host);
  from_execution_space(workspace, wpthlp_mc, wpthlp_mc_host);
  from_execution_space(workspace, rtp2_mc, rtp2_mc_host);
  from_execution_space(workspace, thlp2_mc, thlp2_mc_host);
  from_execution_space(workspace, rtpthlp_mc, rtpthlp_mc_host);
  from_execution_space(workspace, Skw_zm_smooth, Skw_zm_smooth_host);
}
} // namespace clubb_kokkos
