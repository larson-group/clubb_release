// Host/device entries to the ported CLUBB core.
//
// advance_clubb_core_api mirrors advance_clubb_core_api_multi_col in
// src/CLUBB_core/clubb_api_module.F90. It stages host fields into column-first
// execution-space Views, calls the C++ core and copies InOut/Out fields back.
// Compatible device callers share their storage. Statistics run in C++; optional
// BFB PDF matrices retain f2k. This entry never calls the native core routine.
// Public entries:
//   advance_clubb_core_api: stage fields and advance the multicolumn C++ core.
//   calculate_thlp2_rad_api: update radiative temperature-variance forcing.
//   gregorian2julian_day_api/leap_year_api: query Gregorian day and leap rules.
//   vertical_integral_api: integrate one density-weighted vertical profile.
//   setup_grid_heights_api: refresh heights/weights for one or many columns.
//   setup_grid_api: allocate and initialize one or many column grids.
//   lin_interpolate_two_points_api: interpolate a scalar between two heights.
//   calc_derived_params_api: update diffusivities, minimum length and mixture limit.
//   init_pdf_params_api: allocate and initialize caller-owned PDF parameters.
//   calculate_spurious_source_api: evaluate a column conservation residual.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/advance_clubb_core_module.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/calendar.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/numerical_check.hpp"
#include "CLUBB_core/error_code.hpp"
#include "stats_netcdf.hpp"
#include "grid_class.hpp"
#include "parameters_tunable.hpp"
#include <iostream>

namespace clubb_kokkos {

template <
    class Array6, class Array7, class Array8, class Array11, class Array14, class Array15,
    class Array16, class Array17, class Array18, class Array19, class Array20, class Array21,
    class Array22, class Array23, class Array24, class Array25, class Array26, class Array27,
    class Array28, class Array29, class Array30, class Array31, class Array32, class Array33,
    class Array34, class Array35, class Array36, class Array37, class Array38, class Array39,
    class Array40, class Array41, class Array42, class Array43, class Array44, class Array45,
    class Array46, class Array47, class Array48, class Array49, class Array50, class Array51,
    class Array52, class Array53, class Array54, class Array55, class Array56, class Array57,
    class Array58, class Array59, class Array60, class Array62, class Array63, class Array70,
    class Array71, class Array72, class Array73, class Array74, class Array75, class Array76,
    class Array77, class Array78, class Array79, class Array80, class Array81, class Array82,
    class Array83, class Array84, class Array85, class Array86, class Array87, class Array88,
    class Array89, class Array90, class Array91, class Array92, class Array93, class Array94,
    class Array95, class Array97, class Array98, class Array99, class Array100, class Array101,
    class Array102, class Array103, class Array104, class Array105, class Array106, class Array107,
    class Array108, class Array109, class Array110, class Array111, class Array112, class Array113,
    class Array114, class Array115, class Array116, class Array117, class Array118, class Array119,
    class Array123, class Array124, class Array125, class Array126, class Array127, class Array128,
    class Array129, class Array130, class Array131, class Array132, class Array133, class Array134>
void advance_clubb_core_api(
    ScratchWorkspace& workspace,                                                        // InOut
    const grid& gr, int nzm, int nzt, int ngrdcol,                                      // In
    bool l_implemented, real_t dt, const Array6& fcor, const Array7& fcor_y,            // In
    const Array8& sfc_elevation,                                                        // In
    int hydromet_dim,                                                                   // In
    int sclr_dim, const Array11& sclr_tol, int edsclr_dim,                              // In
    const sclr_idx_type& sclr_idx,                                                      // In
    const Array14& thlm_forcing, const Array15& rtm_forcing, const Array16& um_forcing, // In
    const Array17& vm_forcing,                                                          // In
    const Array18& sclrm_forcing, const Array19& edsclrm_forcing,                       // In
    const Array20& wprtp_forcing,                                                       // In
    const Array21& wpthlp_forcing, const Array22& rtp2_forcing,                         // In
    const Array23& thlp2_forcing,                                                       // In
    const Array24& rtpthlp_forcing, const Array25& wm_zm, const Array26& wm_zt,         // In
    const Array27& wpthlp_sfc, const Array28& wprtp_sfc, const Array29& upwp_sfc,       // In
    const Array30& vpwp_sfc, const Array31& p_sfc,                                      // In
    const Array32& wpsclrp_sfc, const Array33& wpedsclrp_sfc,                           // In
    const Array34& upwp_sfc_pert, const Array35& vpwp_sfc_pert,                         // In
    const Array36& rtm_ref, const Array37& thlm_ref, const Array38& um_ref,             // In
    const Array39& vm_ref, const Array40& ug, const Array41& vg,                        // In
    Array42 p_in_Pa,                                                                    // InOut
    const Array43& rho_zm, const Array44& rho,                                          // In
    Array45 exner,                                                                      // InOut
    const Array46& rho_ds_zm, const Array47& rho_ds_zt, const Array48& invrs_rho_ds_zm, // In
    const Array49& invrs_rho_ds_zt, const Array50& thv_ds_zm, const Array51& thv_ds_zt, // In
    const Array52& l_mix_rat_hm,                                                        // In
    const Array53& rfrzm,                                                               // In
    const Array54& wphydrometp, const Array55& wp2hmp, const Array56& rtphmp_zt,        // In
    const Array57& thlphmp_zt,                                                          // In
    const Array58& host_dx, const Array59& host_dy,                                     // In
    const Array60& clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,            // In
    const Array62& lmin,                                                                // In
    const Array63& mixt_frac_max_mag, real_t T0, real_t ts_nudge,                       // In
    real_t rtm_min, real_t rtm_nudge_max_altitude,                                      // In
    const clubb_config_flags_type& clubb_config_flags,                                  // In
    stats_type& stats,                                                                  // InOut
    Array70 um, Array71 vm, Array72 upwp, Array73 vpwp, Array74 up2, Array75 vp2,       // InOut
    Array76 up3, Array77 vp3,                                                           // InOut
    Array78 thlm, Array79 rtm, Array80 wprtp, Array81 wpthlp,                           // InOut
    Array82 wp2, Array83 wp3, Array84 rtp2, Array85 rtp3, Array86 thlp2, Array87 thlp3, // InOut
    Array88 rtpthlp,                                                                    // InOut
    Array89 sclrm,                                                                      // InOut
    Array90 sclrp2, Array91 sclrp3, Array92 sclrprtp, Array93 sclrpthlp,                // InOut
    Array94 wpsclrp, Array95 edsclrm, err_info_type& err_info_api,                      // InOut
    Array97 rcm, Array98 cloud_frac,                                                    // InOut
    Array99 wpthvp, Array100 wp2thvp, Array101 wp2up, Array102 rtpthvp,                 // InOut
    Array103 thlpthvp,                                                                  // InOut
    Array104 sclrpthvp,                                                                 // InOut
    Array105 wp2rtp, Array106 wp2thlp, Array107 uprcp, Array108 vprcp,                  // InOut
    Array109 rc_coef_zm, Array110 wp4,                                                  // InOut
    Array111 wpup2, Array112 wpvp2, Array113 wp2up2, Array114 wp2vp2,                   // InOut
    Array115 ice_supersat_frac,                                                         // InOut
    Array116 um_pert, Array117 vm_pert, Array118 upwp_pert, Array119 vpwp_pert,         // InOut
    pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm,                            // InOut
    implicit_coefs_terms& pdf_implicit_coefs_terms,                                     // InOut
    Array123 Kh_zm, Array124 Kh_zt,                                                     // Out
    Array125 thlprcp, Array126 wprcp, Array127 w_up_in_cloud, Array128 w_down_in_cloud, // Out
    Array129 cloudy_updraft_frac, Array130 cloudy_downdraft_frac,                       // Out
    Array131 rcm_in_layer, Array132 cloud_cover, Array133 invrs_tau_zm,                 // Out
    Array134 Lscale) {                                                                  // Out
  // Advance one core timestep, restoring InOut/Out fields and err_info to the
  // caller. Native API data-region ownership: host inputs are staged once here.
  // The driver keeps runtime fields resident; compatible device Views bypass
  // allocation/copy here. Other callers retain the host-facing API boundary.
  auto workspace_scope               = workspace.scope();
  auto fcor_execution                = to_execution_space(workspace, fcor, true);
  auto fcor_y_execution              = to_execution_space(workspace, fcor_y, true);
  auto sfc_elevation_execution       = to_execution_space(workspace, sfc_elevation, true);
  auto sclr_tol_execution            = to_execution_space(workspace, sclr_tol, true);
  auto thlm_forcing_execution        = to_execution_space(workspace, thlm_forcing, true);
  auto rtm_forcing_execution         = to_execution_space(workspace, rtm_forcing, true);
  auto um_forcing_execution          = to_execution_space(workspace, um_forcing, true);
  auto vm_forcing_execution          = to_execution_space(workspace, vm_forcing, true);
  auto sclrm_forcing_execution       = to_execution_space(workspace, sclrm_forcing, true);
  auto edsclrm_forcing_execution     = to_execution_space(workspace, edsclrm_forcing, true);
  auto wprtp_forcing_execution       = to_execution_space(workspace, wprtp_forcing, true);
  auto wpthlp_forcing_execution      = to_execution_space(workspace, wpthlp_forcing, true);
  auto rtp2_forcing_execution        = to_execution_space(workspace, rtp2_forcing, true);
  auto thlp2_forcing_execution       = to_execution_space(workspace, thlp2_forcing, true);
  auto rtpthlp_forcing_execution     = to_execution_space(workspace, rtpthlp_forcing, true);
  auto wm_zm_execution               = to_execution_space(workspace, wm_zm, true);
  auto wm_zt_execution               = to_execution_space(workspace, wm_zt, true);
  auto wpthlp_sfc_execution          = to_execution_space(workspace, wpthlp_sfc, true);
  auto wprtp_sfc_execution           = to_execution_space(workspace, wprtp_sfc, true);
  auto upwp_sfc_execution            = to_execution_space(workspace, upwp_sfc, true);
  auto vpwp_sfc_execution            = to_execution_space(workspace, vpwp_sfc, true);
  auto p_sfc_execution               = to_execution_space(workspace, p_sfc, true);
  auto wpsclrp_sfc_execution         = to_execution_space(workspace, wpsclrp_sfc, true);
  auto wpedsclrp_sfc_execution       = to_execution_space(workspace, wpedsclrp_sfc, true);
  auto upwp_sfc_pert_execution       = to_execution_space(workspace, upwp_sfc_pert, true);
  auto vpwp_sfc_pert_execution       = to_execution_space(workspace, vpwp_sfc_pert, true);
  auto rtm_ref_execution             = to_execution_space(workspace, rtm_ref, true);
  auto thlm_ref_execution            = to_execution_space(workspace, thlm_ref, true);
  auto um_ref_execution              = to_execution_space(workspace, um_ref, true);
  auto vm_ref_execution              = to_execution_space(workspace, vm_ref, true);
  auto ug_execution                  = to_execution_space(workspace, ug, true);
  auto vg_execution                  = to_execution_space(workspace, vg, true);
  auto p_in_Pa_execution             = to_execution_space(workspace, p_in_Pa, true);
  auto rho_zm_execution              = to_execution_space(workspace, rho_zm, true);
  auto rho_execution                 = to_execution_space(workspace, rho, true);
  auto exner_execution               = to_execution_space(workspace, exner, true);
  auto rho_ds_zm_execution           = to_execution_space(workspace, rho_ds_zm, true);
  auto rho_ds_zt_execution           = to_execution_space(workspace, rho_ds_zt, true);
  auto invrs_rho_ds_zm_execution     = to_execution_space(workspace, invrs_rho_ds_zm, true);
  auto invrs_rho_ds_zt_execution     = to_execution_space(workspace, invrs_rho_ds_zt, true);
  auto thv_ds_zm_execution           = to_execution_space(workspace, thv_ds_zm, true);
  auto thv_ds_zt_execution           = to_execution_space(workspace, thv_ds_zt, true);
  auto l_mix_rat_hm_execution        = to_execution_space(workspace, l_mix_rat_hm, true);
  auto rfrzm_execution               = to_execution_space(workspace, rfrzm, true);
  auto wphydrometp_execution         = to_execution_space(workspace, wphydrometp, true);
  auto wp2hmp_execution              = to_execution_space(workspace, wp2hmp, true);
  auto rtphmp_zt_execution           = to_execution_space(workspace, rtphmp_zt, true);
  auto thlphmp_zt_execution          = to_execution_space(workspace, thlphmp_zt, true);
  auto host_dx_execution             = to_execution_space(workspace, host_dx, true);
  auto host_dy_execution             = to_execution_space(workspace, host_dy, true);
  auto clubb_params_execution        = to_execution_space(workspace, clubb_params, true);
  auto lmin_execution                = to_execution_space(workspace, lmin, true);
  auto mixt_frac_max_mag_execution   = to_execution_space(workspace, mixt_frac_max_mag, true);
  auto um_execution                  = to_execution_space(workspace, um, true);
  auto vm_execution                  = to_execution_space(workspace, vm, true);
  auto upwp_execution                = to_execution_space(workspace, upwp, true);
  auto vpwp_execution                = to_execution_space(workspace, vpwp, true);
  auto up2_execution                 = to_execution_space(workspace, up2, true);
  auto vp2_execution                 = to_execution_space(workspace, vp2, true);
  auto up3_execution                 = to_execution_space(workspace, up3, true);
  auto vp3_execution                 = to_execution_space(workspace, vp3, true);
  auto thlm_execution                = to_execution_space(workspace, thlm, true);
  auto rtm_execution                 = to_execution_space(workspace, rtm, true);
  auto wprtp_execution               = to_execution_space(workspace, wprtp, true);
  auto wpthlp_execution              = to_execution_space(workspace, wpthlp, true);
  auto wp2_execution                 = to_execution_space(workspace, wp2, true);
  auto wp3_execution                 = to_execution_space(workspace, wp3, true);
  auto rtp2_execution                = to_execution_space(workspace, rtp2, true);
  auto rtp3_execution                = to_execution_space(workspace, rtp3, true);
  auto thlp2_execution               = to_execution_space(workspace, thlp2, true);
  auto thlp3_execution               = to_execution_space(workspace, thlp3, true);
  auto rtpthlp_execution             = to_execution_space(workspace, rtpthlp, true);
  auto sclrm_execution               = to_execution_space(workspace, sclrm, true);
  auto sclrp2_execution              = to_execution_space(workspace, sclrp2, true);
  auto sclrp3_execution              = to_execution_space(workspace, sclrp3, true);
  auto sclrprtp_execution            = to_execution_space(workspace, sclrprtp, true);
  auto sclrpthlp_execution           = to_execution_space(workspace, sclrpthlp, true);
  auto wpsclrp_execution             = to_execution_space(workspace, wpsclrp, true);
  auto edsclrm_execution             = to_execution_space(workspace, edsclrm, true);
  auto rcm_execution                 = to_execution_space(workspace, rcm, true);
  auto cloud_frac_execution          = to_execution_space(workspace, cloud_frac, true);
  auto wpthvp_execution              = to_execution_space(workspace, wpthvp, true);
  auto wp2thvp_execution             = to_execution_space(workspace, wp2thvp, true);
  auto wp2up_execution               = to_execution_space(workspace, wp2up, true);
  auto rtpthvp_execution             = to_execution_space(workspace, rtpthvp, true);
  auto thlpthvp_execution            = to_execution_space(workspace, thlpthvp, true);
  auto sclrpthvp_execution           = to_execution_space(workspace, sclrpthvp, true);
  auto wp2rtp_execution              = to_execution_space(workspace, wp2rtp, true);
  auto wp2thlp_execution             = to_execution_space(workspace, wp2thlp, true);
  auto uprcp_execution               = to_execution_space(workspace, uprcp, true);
  auto vprcp_execution               = to_execution_space(workspace, vprcp, true);
  auto rc_coef_zm_execution          = to_execution_space(workspace, rc_coef_zm, true);
  auto wp4_execution                 = to_execution_space(workspace, wp4, true);
  auto wpup2_execution               = to_execution_space(workspace, wpup2, true);
  auto wpvp2_execution               = to_execution_space(workspace, wpvp2, true);
  auto wp2up2_execution              = to_execution_space(workspace, wp2up2, true);
  auto wp2vp2_execution              = to_execution_space(workspace, wp2vp2, true);
  auto ice_supersat_frac_execution   = to_execution_space(workspace, ice_supersat_frac, true);
  auto um_pert_execution             = to_execution_space(workspace, um_pert, true);
  auto vm_pert_execution             = to_execution_space(workspace, vm_pert, true);
  auto upwp_pert_execution           = to_execution_space(workspace, upwp_pert, true);
  auto vpwp_pert_execution           = to_execution_space(workspace, vpwp_pert, true);
  auto Kh_zm_execution               = to_execution_space(workspace, Kh_zm, false);
  auto Kh_zt_execution               = to_execution_space(workspace, Kh_zt, false);
  auto thlprcp_execution             = to_execution_space(workspace, thlprcp, false);
  auto wprcp_execution               = to_execution_space(workspace, wprcp, false);
  auto w_up_in_cloud_execution       = to_execution_space(workspace, w_up_in_cloud, false);
  auto w_down_in_cloud_execution     = to_execution_space(workspace, w_down_in_cloud, false);
  auto cloudy_updraft_frac_execution = to_execution_space(workspace, cloudy_updraft_frac, false);
  auto cloudy_downdraft_frac_execution =
      to_execution_space(workspace, cloudy_downdraft_frac, false);
  auto rcm_in_layer_execution = to_execution_space(workspace, rcm_in_layer, false);
  auto cloud_cover_execution  = to_execution_space(workspace, cloud_cover, false);
  auto invrs_tau_zm_execution = to_execution_space(workspace, invrs_tau_zm, false);
  auto Lscale_execution       = to_execution_space(workspace, Lscale, false);

  advance_clubb_core(
      workspace, gr, nzm, nzt, ngrdcol,                                               // In
      l_implemented, dt, fcor_execution, fcor_y_execution, sfc_elevation_execution,   // In
      hydromet_dim,                                                                   // In
      sclr_dim, sclr_tol_execution, edsclr_dim, sclr_idx,                             // In
      thlm_forcing_execution, rtm_forcing_execution, um_forcing_execution,            // In
      vm_forcing_execution,                                                           // In
      sclrm_forcing_execution, edsclrm_forcing_execution, wprtp_forcing_execution,    // In
      wpthlp_forcing_execution, rtp2_forcing_execution, thlp2_forcing_execution,      // In
      rtpthlp_forcing_execution, wm_zm_execution, wm_zt_execution,                    // In
      wpthlp_sfc_execution, wprtp_sfc_execution, upwp_sfc_execution,                  // In
      vpwp_sfc_execution, p_sfc_execution,                                            // In
      wpsclrp_sfc_execution, wpedsclrp_sfc_execution,                                 // In
      upwp_sfc_pert_execution, vpwp_sfc_pert_execution,                               // In
      rtm_ref_execution, thlm_ref_execution, um_ref_execution, vm_ref_execution,      // In
      ug_execution, vg_execution,                                                     // In
      p_in_Pa_execution,                                                              // InOut
      rho_zm_execution, rho_execution,                                                // In
      exner_execution,                                                                // InOut
      rho_ds_zm_execution, rho_ds_zt_execution, invrs_rho_ds_zm_execution,            // In
      invrs_rho_ds_zt_execution, thv_ds_zm_execution, thv_ds_zt_execution,            // In
      l_mix_rat_hm_execution,                                                         // In
      rfrzm_execution,                                                                // In
      wphydrometp_execution, wp2hmp_execution, rtphmp_zt_execution,                   // In
      thlphmp_zt_execution,                                                           // In
      host_dx_execution, host_dy_execution,                                           // In
      clubb_params_execution, nu_vert_res_dep, lmin_execution,                        // In
      mixt_frac_max_mag_execution, T0, ts_nudge,                                      // In
      rtm_min, rtm_nudge_max_altitude,                                                // In
      clubb_config_flags,                                                             // In
      stats,                                                                          // InOut
      um_execution, vm_execution, upwp_execution, vpwp_execution, up2_execution,      // InOut
      vp2_execution, up3_execution, vp3_execution,                                    // InOut
      thlm_execution, rtm_execution, wprtp_execution, wpthlp_execution,               // InOut
      wp2_execution, wp3_execution, rtp2_execution, rtp3_execution, thlp2_execution,  // InOut
      thlp3_execution, rtpthlp_execution,                                             // InOut
      sclrm_execution,                                                                // InOut
      sclrp2_execution, sclrp3_execution, sclrprtp_execution, sclrpthlp_execution,    // InOut
      wpsclrp_execution, edsclrm_execution,                                           // InOut
      rcm_execution, cloud_frac_execution,                                            // InOut
      wpthvp_execution, wp2thvp_execution, wp2up_execution, rtpthvp_execution,        // InOut
      thlpthvp_execution,                                                             // InOut
      sclrpthvp_execution,                                                            // InOut
      wp2rtp_execution, wp2thlp_execution, uprcp_execution, vprcp_execution,          // InOut
      rc_coef_zm_execution, wp4_execution,                                            // InOut
      wpup2_execution, wpvp2_execution, wp2up2_execution, wp2vp2_execution,           // InOut
      ice_supersat_frac_execution,                                                    // InOut
      um_pert_execution, vm_pert_execution, upwp_pert_execution, vpwp_pert_execution, // InOut
      pdf_params, pdf_params_zm,                                                      // InOut
      pdf_implicit_coefs_terms,                                                       // InOut
      err_info_api,                                                                   // InOut
      Kh_zm_execution, Kh_zt_execution,                                               // Out
      thlprcp_execution, wprcp_execution, w_up_in_cloud_execution,                    // Out
      w_down_in_cloud_execution,                                                      // Out
      cloudy_updraft_frac_execution, cloudy_downdraft_frac_execution,                 // Out
      rcm_in_layer_execution, cloud_cover_execution, invrs_tau_zm_execution,          // Out
      Lscale_execution);                                                              // Out

  from_execution_space(workspace, p_in_Pa, p_in_Pa_execution);
  from_execution_space(workspace, exner, exner_execution);
  from_execution_space(workspace, um, um_execution);
  from_execution_space(workspace, vm, vm_execution);
  from_execution_space(workspace, upwp, upwp_execution);
  from_execution_space(workspace, vpwp, vpwp_execution);
  from_execution_space(workspace, up2, up2_execution);
  from_execution_space(workspace, vp2, vp2_execution);
  from_execution_space(workspace, up3, up3_execution);
  from_execution_space(workspace, vp3, vp3_execution);
  from_execution_space(workspace, thlm, thlm_execution);
  from_execution_space(workspace, rtm, rtm_execution);
  from_execution_space(workspace, wprtp, wprtp_execution);
  from_execution_space(workspace, wpthlp, wpthlp_execution);
  from_execution_space(workspace, wp2, wp2_execution);
  from_execution_space(workspace, wp3, wp3_execution);
  from_execution_space(workspace, rtp2, rtp2_execution);
  from_execution_space(workspace, rtp3, rtp3_execution);
  from_execution_space(workspace, thlp2, thlp2_execution);
  from_execution_space(workspace, thlp3, thlp3_execution);
  from_execution_space(workspace, rtpthlp, rtpthlp_execution);
  from_execution_space(workspace, sclrm, sclrm_execution);
  from_execution_space(workspace, sclrp2, sclrp2_execution);
  from_execution_space(workspace, sclrp3, sclrp3_execution);
  from_execution_space(workspace, sclrprtp, sclrprtp_execution);
  from_execution_space(workspace, sclrpthlp, sclrpthlp_execution);
  from_execution_space(workspace, wpsclrp, wpsclrp_execution);
  from_execution_space(workspace, edsclrm, edsclrm_execution);
  from_execution_space(workspace, rcm, rcm_execution);
  from_execution_space(workspace, cloud_frac, cloud_frac_execution);
  from_execution_space(workspace, wpthvp, wpthvp_execution);
  from_execution_space(workspace, wp2thvp, wp2thvp_execution);
  from_execution_space(workspace, wp2up, wp2up_execution);
  from_execution_space(workspace, rtpthvp, rtpthvp_execution);
  from_execution_space(workspace, thlpthvp, thlpthvp_execution);
  from_execution_space(workspace, sclrpthvp, sclrpthvp_execution);
  from_execution_space(workspace, wp2rtp, wp2rtp_execution);
  from_execution_space(workspace, wp2thlp, wp2thlp_execution);
  from_execution_space(workspace, uprcp, uprcp_execution);
  from_execution_space(workspace, vprcp, vprcp_execution);
  from_execution_space(workspace, rc_coef_zm, rc_coef_zm_execution);
  from_execution_space(workspace, wp4, wp4_execution);
  from_execution_space(workspace, wpup2, wpup2_execution);
  from_execution_space(workspace, wpvp2, wpvp2_execution);
  from_execution_space(workspace, wp2up2, wp2up2_execution);
  from_execution_space(workspace, wp2vp2, wp2vp2_execution);
  from_execution_space(workspace, ice_supersat_frac, ice_supersat_frac_execution);
  from_execution_space(workspace, um_pert, um_pert_execution);
  from_execution_space(workspace, vm_pert, vm_pert_execution);
  from_execution_space(workspace, upwp_pert, upwp_pert_execution);
  from_execution_space(workspace, vpwp_pert, vpwp_pert_execution);
  from_execution_space(workspace, Kh_zm, Kh_zm_execution);
  from_execution_space(workspace, Kh_zt, Kh_zt_execution);
  from_execution_space(workspace, thlprcp, thlprcp_execution);
  from_execution_space(workspace, wprcp, wprcp_execution);
  from_execution_space(workspace, w_up_in_cloud, w_up_in_cloud_execution);
  from_execution_space(workspace, w_down_in_cloud, w_down_in_cloud_execution);
  from_execution_space(workspace, cloudy_updraft_frac, cloudy_updraft_frac_execution);
  from_execution_space(workspace, cloudy_downdraft_frac, cloudy_downdraft_frac_execution);
  from_execution_space(workspace, rcm_in_layer, rcm_in_layer_execution);
  from_execution_space(workspace, cloud_cover, cloud_cover_execution);
  from_execution_space(workspace, invrs_tau_zm, invrs_tau_zm_execution);
  from_execution_space(workspace, Lscale, Lscale_execution);

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info_api.err_code, clubb_fatal_error)) {
      std::cerr << err_info_api.err_header_global << '\n';
      std::cerr << "Fatal error advance_clubb_core in API module\n";
    }
  }
}

//=============================================================================
template <class Rcm, class Thlprcp, class Radht, class Parameters, class Forcing>
void calculate_thlp2_rad_api(ScratchWorkspace& workspace,                   // InOut
                             int ngrdcol, int nzm, int nzt, const grid& gr, // In
                             Rcm rcm, Thlprcp thlprcp, Radht radht,
                             Parameters clubb_params, // In
                             Forcing thlp2_forcing) {
  auto workspace_scope = workspace.scope(); // InOut
  // Add the radiative contribution to temperature variance forcing [K^2/s].
  // rcm [kg/kg] and radht [K/s] are on thermodynamic levels; thlprcp [K kg/kg]
  // and the forcing are on momentum levels. Fields are column-first Views.
  // This boundary mirrors the Native API's input/inout data region.
  auto rcm_execution           = to_execution_space(workspace, rcm);
  auto thlprcp_execution       = to_execution_space(workspace, thlprcp);
  auto radht_execution         = to_execution_space(workspace, radht);
  auto clubb_params_execution  = to_execution_space(workspace, clubb_params);
  auto thlp2_forcing_execution = to_execution_space(workspace, thlp2_forcing);

  calculate_thlp2_rad(workspace, ngrdcol, nzm, nzt, gr, // In
                      rcm_execution, thlprcp_execution, radht_execution,
                      clubb_params_execution,   // In
                      thlp2_forcing_execution); // InOut

  from_execution_space(workspace, thlp2_forcing, thlp2_forcing_execution);
}

//=============================================================================
inline int gregorian2julian_day_api(int day, int month, int year) { // In
  // Return the Gregorian day of year (1–366) for day/month/year.
  return gregorian2julian_day(day, month, year);
}

//=============================================================================
inline bool leap_year_api(int year) { // In
  // Determine whether the Gregorian year has a leap day.
  return leap_year(year);
}

//=============================================================================
template <class Rho, class Field, class Dz>
inline real_t vertical_integral_api(ScratchWorkspace& workspace, // InOut
                                    int total_idx, Rho rho_ds,   // In
                                    Field field, Dz dz) {        // In
  // Integrate total_idx levels of field [units vary], weighted by density
  // rho_ds [kg/m^3] and layer thickness dz [m], in source level order.
  return vertical_integral(workspace, total_idx, rho_ds, // In
                           field, dz);                   // In
}

// Both APIs keep Native error checks after the core calls.
template <class Momentum, class Thermodynamic>
void setup_grid_heights_api(bool l_implemented, bool l_ascending_grid,    // In
                            int grid_type, real_t deltaz, real_t zm_init, // In
                            Momentum momentum_heights,                    // In
                            grid& gr,                                     // InOut
                            Thermodynamic thermodynamic_heights,          // In
                            err_info_type& err_info_api) {                // InOut
  // Description:
  //   Reset an allocated single-column grid's heights and interpolation weights.
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D deltaz_col("setup_grid_heights_api::deltaz_col", 1);
  Host1D zm_init_col("setup_grid_heights_api::zm_init_col", 1);
  Real2D momentum_heights_col("setup_grid_heights_api::momentum_heights_col", 1, gr.nzm);
  Real2D thermodynamic_heights_col("setup_grid_heights_api::thermodynamic_heights_col", 1, gr.nzt);
  deltaz_col(0)  = deltaz;
  zm_init_col(0) = zm_init;
  Kokkos::deep_copy(Kokkos::subview(momentum_heights_col, 0, Kokkos::ALL),
                    to_execution_space(momentum_heights));
  Kokkos::deep_copy(Kokkos::subview(thermodynamic_heights_col, 0, Kokkos::ALL),
                    to_execution_space(thermodynamic_heights));
  setup_grid_heights(gr.nzm, gr.nzt, 1,                               // In
                     l_implemented, l_ascending_grid,                 // In
                     grid_type, deltaz_col, zm_init_col,              // In
                     momentum_heights_col, thermodynamic_heights_col, // In
                     gr, err_info_api);                               // InOut
  if (any_equal(err_info_api.err_code, clubb_fatal_error)) {
    std::cerr << err_info_api.err_header_global << '\n';
    std::cerr << "Error in CLUBB calling setup_grid_heights\n";
  }
}

template <class Spacing, class Start, class Momentum, class Thermodynamic>
void setup_grid_heights_api(int nzm, int nzt, int ngrdcol,                // In
                            bool l_implemented, bool l_ascending_grid,    // In
                            int grid_type, Spacing deltaz, Start zm_init, // In
                            Momentum momentum_heights,                    // In
                            grid& gr,                                     // InOut
                            Thermodynamic thermodynamic_heights,          // In
                            err_info_type& err_info_api) {                // InOut
  // Description:
  //   Reset allocated column grids' heights and interpolation weights.
  setup_grid_heights(nzm, nzt, ngrdcol,                       // In
                     l_implemented, l_ascending_grid,         // In
                     grid_type, deltaz, zm_init,              // In
                     momentum_heights, thermodynamic_heights, // In
                     gr, err_info_api);                       // InOut
  if (any_equal(err_info_api.err_code, clubb_fatal_error)) {
    std::cerr << err_info_api.err_header_global << '\n';
    std::cerr << "Error in CLUBB calling setup_grid_heights\n";
  }
}

template <class Momentum, class Thermodynamic>
void setup_grid_api(int nzmax, real_t sfc_elevation, bool l_implemented,            // In
                    bool l_ascending_grid, int grid_type,                           // In
                    real_t deltaz, real_t zm_init, real_t zm_top,                   // In
                    Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                    grid& gr, err_info_type& err_info_api) {                        // InOut
  // Description:
  //   Allocate and initialize the grid for a single-column caller.
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D sfc_elevation_col("setup_grid_api::sfc_elevation_col", 1);
  Host1D deltaz_col("setup_grid_api::deltaz_col", 1);
  Host1D zm_init_col("setup_grid_api::zm_init_col", 1);
  Host1D zm_top_col("setup_grid_api::zm_top_col", 1);
  Real2D momentum_heights_col("setup_grid_api::momentum_heights_col", 1, nzmax);
  Real2D thermodynamic_heights_col("setup_grid_api::thermodynamic_heights_col", 1, nzmax - 1);
  sfc_elevation_col(0) = sfc_elevation;
  deltaz_col(0)        = deltaz;
  zm_init_col(0)       = zm_init;
  zm_top_col(0)        = zm_top;
  Kokkos::deep_copy(Kokkos::subview(momentum_heights_col, 0, Kokkos::ALL),
                    to_execution_space(momentum_heights));
  Kokkos::deep_copy(Kokkos::subview(thermodynamic_heights_col, 0, Kokkos::ALL),
                    to_execution_space(thermodynamic_heights));
  setup_grid(nzmax, 1, sfc_elevation_col, l_implemented,                       // In
             l_ascending_grid, grid_type, deltaz_col, zm_init_col, zm_top_col, // In
             momentum_heights_col, thermodynamic_heights_col,                  // In
             gr, err_info_api);                                                // InOut
  if (any_equal(err_info_api.err_code, clubb_fatal_error)) {
    std::cerr << err_info_api.err_header_global << '\n';
    std::cerr << "Error in CLUBB calling setup_grid\n";
  }
}

template <class Surface, class Spacing, class Start, class Top, class Momentum, class Thermodynamic>
void setup_grid_api(int nzmax, int ngrdcol, Surface sfc_elevation, bool l_implemented, // In
                    bool l_ascending_grid, int grid_type, Spacing deltaz, Start zm_init,
                    Top zm_top,                                                     // In
                    Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                    grid& gr, err_info_type& err_info_api) {                        // InOut
  // Description:
  //   Allocate and initialize the grid for a column chunk.
  setup_grid(nzmax, ngrdcol, sfc_elevation, l_implemented,         // In
             l_ascending_grid, grid_type, deltaz, zm_init, zm_top, // In
             momentum_heights, thermodynamic_heights,              // In
             gr, err_info_api);                                    // InOut
  if (any_equal(err_info_api.err_code, clubb_fatal_error)) {
    std::cerr << err_info_api.err_header_global << '\n';
    std::cerr << "Error in CLUBB calling setup_grid\n";
  }
}

//=============================================================================
KOKKOS_INLINE_FUNCTION real_t lin_interpolate_two_points_api(real_t height_int, real_t height_high,
                                                             real_t height_low, // In
                                                             real_t var_high,
                                                             real_t var_low) { // In
  // Interpolate var [units vary] to height_int between height_low/high [m].
  return lin_interpolate_two_points(height_int, height_high, height_low, // In
                                    var_high, var_low);                  // In
}

template <class Vector>
void calc_derived_params_api(const grid& gr, int grid_type, real_t deltaz,       // In
                             Vector clubb_params,                                // In
                             bool l_prescribed_avg_deltaz,                       // In
                             nu_vertical_res_dep& nu_vert_res_dep, real_t& lmin, // InOut
                             real_t& mixt_frac_max_mag                           // InOut
) {
  // Description:
  //   Adjust background eddy diffusivity for vertical spacing in one column.
  // gr: Grid definition [-]. grid_type: Uniform(1), thermodynamic(2), momentum(3) [-].
  // deltaz: Change per height level [m]. clubb_params: Tunable parameter vector [-].
  // l_prescribed_avg_deltaz: Use deltaz as the average spacing when true [-].
  // nu_vert_res_dep: Resolution-dependent diffusivities [m^2/s].
  // lmin: Minimum mixing length [m]. mixt_frac_max_mag: Mixture-fraction limit [-].
  using Host1D = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host1D deltaz_col("calc_derived_params_api::deltaz_col", 1);
  Host1D lmin_col(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_derived_params_api::lmin_col"), 1);
  Host1D mixt_frac_max_mag_col(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                  "calc_derived_params_api::mixt_frac_max_mag_col"),
                               1);
  Real2D clubb_params_col("calc_derived_params_api::clubb_params_col", 1, nparams);
  deltaz_col(0) = deltaz;
  Kokkos::deep_copy(Kokkos::subview(clubb_params_col, 0, Kokkos::ALL),
                    to_execution_space(clubb_params));
  // Native passes its original rank-1 vector by sequence association. C++ uses
  // the corresponding one-column array, with identical values and ordering.
  calc_derived_params(gr, 1, grid_type, deltaz_col, // In
                      clubb_params_col,             // In
                      l_prescribed_avg_deltaz,      // In
                      nu_vert_res_dep, lmin_col,    // InOut
                      mixt_frac_max_mag_col);       // InOut
  lmin              = lmin_col(0);
  mixt_frac_max_mag = mixt_frac_max_mag_col(0);
}

template <class Spacing, class Matrix, class MinimumLength, class MixtureLimit>
void calc_derived_params_api(const grid& gr, int ngrdcol, int grid_type, Spacing deltaz, // In
                             Matrix clubb_params,                                        // In
                             bool l_prescribed_avg_deltaz,                               // In
                             nu_vertical_res_dep& nu_vert_res_dep, MinimumLength lmin,   // InOut
                             MixtureLimit mixt_frac_max_mag                              // InOut
) {
  // Description:
  //   Adjust background eddy diffusivity for vertical spacing in a column chunk.
  // gr: Grid definition [-]. ngrdcol: Number of columns [#].
  // grid_type: Uniform(1), thermodynamic(2), momentum(3) input grid [-].
  // deltaz: Change per height level [m]. clubb_params: Tunable parameter matrix [-].
  // l_prescribed_avg_deltaz: Use deltaz as the average spacing when true [-].
  // nu_vert_res_dep: Resolution-dependent diffusivities [m^2/s].
  // lmin: Minimum mixing length [m]. mixt_frac_max_mag: Mixture-fraction limit [-].
  calc_derived_params(gr, ngrdcol, grid_type, deltaz, // In
                      clubb_params,                   // In
                      l_prescribed_avg_deltaz,        // In
                      nu_vert_res_dep, lmin,          // InOut
                      mixt_frac_max_mag);             // InOut
}

//=============================================================================
inline void init_pdf_params_api(int nz, int ngrdcol,       // In
                                pdf_parameter& pdf_params) // Out
{
  // Description: Allocate and initialize the caller-owned PDF parameters.
  // nz: Number of vertical grid levels [-].
  // ngrdcol: Number of grid columns [-].
  // pdf_params: Caller-owned PDF parameters [units vary].
  init_pdf_params(nz, ngrdcol, // In
                  pdf_params); // Out
}

//=============================================================================
KOKKOS_INLINE_FUNCTION real_t calculate_spurious_source_api(real_t integral_after,
                                                            real_t integral_before,           // In
                                                            real_t flux_top, real_t flux_sfc, // In
                                                            real_t integral_forcing,
                                                            real_t dt) { // In
  // Return the change in a vertically integrated quantity [units vary] after
  // accounting for the boundary fluxes and integrated forcing over dt [s].
  return calculate_spurious_source(integral_after, integral_before, // In
                                   flux_top, flux_sfc,              // In
                                   integral_forcing, dt);           // In
}
} // namespace clubb_kokkos
