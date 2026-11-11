#include "scratch_workspace.hpp"
#include "stats_clubb_utilities.hpp"
// Advance CLUBB's equations one timestep with column-first Kokkos Views.
//
// advance_clubb_core retains the native solver/PDF/statistics/error sequence.
// Its own array operations run on the default execution space. Unported lower
// routines use f2k conversions; they never retain these local scratch Views.
// The private surface helper retains its native position below the main routine.
//
// References: Golaz et al. (2002), JAS 59, 3540-3551;
// https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
//
//                         Copyright Notice:
//
//   This code and the source code it references are (C) 2006-2020.
//
//   The distribution of this code and derived works thereof
//                   should include this notice.
//
//   Portions of this code derived from other sources (Hugh Morrison,
//   ACM TOMS, et cetera) are the intellectual
//   property of their respective authors as noted and are also subject
//   to copyright.
//
//
//
// Cloud Layers Unified By Binormals (CLUBB) user license
// agreement.
//
// Thank you for your interest in CLUBB. We work hard to create a
// code that implements the best software engineering practices,
// is supported to the extent allowed by our limited resources,
// and is available without cost to non-commercial users. You may
// use CLUBB if, in return, you abide by these conditions:
//
// 1. Please cite CLUBB in presentations and publications that
//  contain results obtained using CLUBB.
//
// 2. You may not use any part of CLUBB to create or modify
//  another single-column (1D) model that is not called CLUBB.
//  However, you may modify or augment CLUBB or parts of CLUBB if
//  you include "CLUBB" in the name of the resulting single-column
//  model. For example, a user at MIT might modify CLUBB and call
//  the modified version "CLUBB-MIT." Or, for example, a user of
//  the CLM land-surface model might interface CLM to CLUBB and
//  call it "CLM-CLUBB." This naming convention recognizes the
//  contributions of both sets of developers.
//
// 3. You may implement CLUBB as a parameterization in a large-
//  scale host model that has 2 or 3 spatial dimensions without
//  including "CLUBB" in the combined model name, but please
//  acknowledge in presentations and publications that CLUBB has
//  been included as a parameterization.
//
// 4. You may not provide all or part of CLUBB to anyone without
//  prior permission from Vincent Larson (vlarson@uwm.edu). If
//  you wish to share CLUBB with your collaborators without
//  seeking permission, please ask your collaborators to register
//  as CLUBB users at https://carson.math.uwm.edu/larson-group/clubb_site/ and to
//  download CLUBB from there.
//
// 5. You may not use CLUBB for commercial purposes unless you
//  receive permission from Vincent Larson.
//
// 6. You may not re-license all or any part of CLUBB.
//
// 7. CLUBB is provided "as is" and without warranty.
//
// We hope that CLUBB will develop into a community resource. We
// encourage users to contribute their CLUBB modifications or
// extensions to the CLUBB development group. We will then
// consider them for inclusion in CLUBB. Such contributions will
// benefit all CLUBB users. We would be pleased to acknowledge
// contributors and list their CLUBB-related papers on our "About
// CLUBB" webpage (https://carson.math.uwm.edu/larson-group/clubb_site/about.html) for
// those contributors who so desire.
//
// Thanks so much and best wishes for your research!
//
// The CLUBB Development Group
// (Present and past contributors to the source code include
// Vincent Larson, Chris Golaz, David Schanen, Brian Griffin,
// Joshua Fasching, Adam Smith, and Michael Falk).
//-----------------------------------------------------------------------
//
#include "clip_explicit.hpp"
#include "CLUBB_core/calc_pressure.hpp"
#include "CLUBB_core/advance_clubb_core_module.hpp"
#include "saturation.hpp"
#include "T_in_K_module.hpp"
#include "sigma_sqd_w_module.hpp"
#include "numerical_check.hpp"
#include "pdf_closure_module.hpp"
#include "advance_helper_module.hpp"
#include "mixing_length.hpp"
#include "sfc_varnce_module.hpp"
#include "advance_xp3_module.hpp"
#include "advance_xm_wpxp_module.hpp"
#include "advance_xp2_xpyp_module.hpp"
#include "advance_windm_edsclrm_module.hpp"
#include "advance_wp2_wp3_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/error_code.hpp"
#include "stats_netcdf.hpp"
#include <iostream>

namespace clubb_kokkos {
namespace {
void set_sfc_value_of_flux_profiles(int nzm, int ngrdcol, int sclr_dim, int edsclr_dim,
                                    const grid& gr, // In
                                    bool l_host_applies_sfc_fluxes,
                                    bool l_linearize_pbl_winds, // In
                                    const Real1D& wpthlp_sfc, const Real1D& wprtp_sfc,
                                    const Real1D& upwp_sfc,                                   // In
                                    const Real1D& vpwp_sfc,                                   // In
                                    const Real1D& upwp_sfc_pert, const Real1D& vpwp_sfc_pert, // In
                                    const Real2D& wpsclrp_sfc, const Real2D& wpedsclrp_sfc,   // In
                                    Real2D wpthlp, Real2D wprtp, Real2D upwp, Real2D vpwp, // InOut
                                    Real2D upwp_pert, Real2D vpwp_pert,                    // InOut
                                    Real3D wpsclrp, Real3D wpedsclrp);                     // InOut
} // namespace

void advance_clubb_core(
    ScratchWorkspace& workspace,                                                        // InOut
    const grid& gr, int nzm, int nzt, int ngrdcol,                                      // In
    bool l_implemented, real_t dt, const Real1D& fcor, const Real1D& fcor_y,            // In
    const Real1D& sfc_elevation,                                                        // In
    int hydromet_dim,                                                                   // In
    int sclr_dim, const Real1D& sclr_tol, int edsclr_dim,                               // In
    const sclr_idx_type& sclr_idx,                                                      // In
    const Real2D& thlm_forcing, const Real2D& rtm_forcing, const Real2D& um_forcing,    // In
    const Real2D& vm_forcing,                                                           // In
    const Real3D& sclrm_forcing, const Real3D& edsclrm_forcing,                         // In
    const Real2D& wprtp_forcing,                                                        // In
    const Real2D& wpthlp_forcing, const Real2D& rtp2_forcing,                           // In
    const Real2D& thlp2_forcing,                                                        // In
    const Real2D& rtpthlp_forcing, const Real2D& wm_zm, const Real2D& wm_zt,            // In
    const Real1D& wpthlp_sfc, const Real1D& wprtp_sfc, const Real1D& upwp_sfc,          // In
    const Real1D& vpwp_sfc, const Real1D& p_sfc,                                        // In
    const Real2D& wpsclrp_sfc, const Real2D& wpedsclrp_sfc,                             // In
    const Real1D& upwp_sfc_pert, const Real1D& vpwp_sfc_pert,                           // In
    const Real2D& rtm_ref, const Real2D& thlm_ref, const Real2D& um_ref,                // In
    const Real2D& vm_ref, const Real2D& ug, const Real2D& vg,                           // In
    Real2D p_in_Pa,                                                                     // InOut
    const Real2D& rho_zm, const Real2D& rho,                                            // In
    Real2D exner,                                                                       // InOut
    const Real2D& rho_ds_zm, const Real2D& rho_ds_zt, const Real2D& invrs_rho_ds_zm,    // In
    const Real2D& invrs_rho_ds_zt, const Real2D& thv_ds_zm, const Real2D& thv_ds_zt,    // In
    const Bool1D& l_mix_rat_hm,                                                         // In
    const Real2D& rfrzm,                                                                // In
    const Real3D& wphydrometp, const Real3D& wp2hmp, const Real3D& rtphmp_zt,           // In
    const Real3D& thlphmp_zt,                                                           // In
    const Real1D& host_dx, const Real1D& host_dy,                                       // In
    const Real2D& clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,             // In
    const Real1D& lmin,                                                                 // In
    const Real1D& mixt_frac_max_mag, real_t T0, real_t ts_nudge,                        // In
    real_t rtm_min, real_t rtm_nudge_max_altitude,                                      // In
    const clubb_config_flags_type& clubb_config_flags,                                  // In
    stats_type& stats,                                                                  // InOut
    Real2D um, Real2D vm, Real2D upwp, Real2D vpwp, Real2D up2, Real2D vp2, Real2D up3, // InOut
    Real2D vp3,                                                                         // InOut
    Real2D thlm, Real2D rtm, Real2D wprtp, Real2D wpthlp,                               // InOut
    Real2D wp2, Real2D wp3, Real2D rtp2, Real2D rtp3, Real2D thlp2, Real2D thlp3,       // InOut
    Real2D rtpthlp,                                                                     // InOut
    Real3D sclrm,                                                                       // InOut
    Real3D sclrp2, Real3D sclrp3, Real3D sclrprtp, Real3D sclrpthlp,                    // InOut
    Real3D wpsclrp, Real3D edsclrm,                                                     // InOut
    Real2D rcm, Real2D cloud_frac,                                                      // InOut
    Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D rtpthvp, Real2D thlpthvp,       // InOut
    Real3D sclrpthvp,                                                                   // InOut
    Real2D wp2rtp, Real2D wp2thlp, Real2D uprcp, Real2D vprcp, Real2D rc_coef_zm,       // InOut
    Real2D wp4,                                                                         // InOut
    Real2D wpup2, Real2D wpvp2, Real2D wp2up2, Real2D wp2vp2, Real2D ice_supersat_frac, // InOut
    Real2D um_pert, Real2D vm_pert, Real2D upwp_pert, Real2D vpwp_pert,                 // InOut
    pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm,                            // InOut
    implicit_coefs_terms& pdf_implicit_coefs_terms,                                     // InOut
    err_info_type& err_info,                                                            // InOut
    Real2D Kh_zm, Real2D Kh_zt,                                                         // Out
    Real2D thlprcp, Real2D wprcp, Real2D w_up_in_cloud, Real2D w_down_in_cloud,         // Out
    Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac,                           // Out
    Real2D rcm_in_layer, Real2D cloud_cover, Real2D invrs_tau_zm,                       // Out
    Real2D Lscale) {                                                                    // Out
  // Advance CLUBB one timestep. Prognostic/PDF/error fields are updated in
  // place; diffusivity, mixing-length and cloud diagnostics fill the outputs.
  auto workspace_scope       = workspace.scope();
  constexpr real_t tau_const = real_t(1000.0); // Constant dissipation time [s]
  constexpr real_t ufmin     = real_t(0.01);   // Minimum friction velocity [m/s]
  // Used only with tau-based mixing length, as in the native source.
  constexpr bool l_use_invrs_tau_N2_iso = false;

  // Capture scalar values: device min/max operands cannot refer to the host
  // storage of namespace-scope constexpr objects when passed by reference.
  const real_t em_min         = clubb_kokkos::em_min;
  const real_t zero_threshold = clubb_kokkos::zero_threshold;

  // --------------------------- Local Variables ---------------------------
  real_t dt_advance{}; // Advance-step duration; doubled for LMM stepping [s]
  // Virtual potential temperature [K]
  Real2D thvm = workspace_scope.view<Real2D>("advance_clubb_core::thvm", ngrdcol, nzt);
  // Saturation mixing ratio [kg/kg]
  Real2D rsat = workspace_scope.view<Real2D>("advance_clubb_core::rsat", ngrdcol, nzt);
  // PDF width parameter on thermodynamic levels [-]
  Real2D sigma_sqd_w_zt =
      workspace_scope.view<Real2D>("advance_clubb_core::sigma_sqd_w_zt", ngrdcol, nzt);
  // w'thl'^2 [m K^2/s]
  Real2D wpthlp2 = workspace_scope.view<Real2D>("advance_clubb_core::wpthlp2", ngrdcol, nzt);
  // w'rt'^2 [m (kg/kg)^2/s]
  Real2D wprtp2 = workspace_scope.view<Real2D>("advance_clubb_core::wprtp2", ngrdcol, nzt);
  // w'rt'thl' [m K kg/(kg s)]
  Real2D wprtpthlp = workspace_scope.view<Real2D>("advance_clubb_core::wprtpthlp", ngrdcol, nzt);
  // Upward mixing length [m]
  Real2D Lscale_up = workspace_scope.view<Real2D>("advance_clubb_core::Lscale_up", ngrdcol, nzt);
  // Downward mixing length [m]
  Real2D Lscale_down =
      workspace_scope.view<Real2D>("advance_clubb_core::Lscale_down", ngrdcol, nzt);
  // Momentum-level mixing length [m]
  Real2D Lscale_zm = workspace_scope.view<Real2D>("advance_clubb_core::Lscale_zm", ngrdcol, nzm);
  // Turbulent kinetic energy [m^2/s^2]
  Real2D em = workspace_scope.view<Real2D>("advance_clubb_core::em", ngrdcol, nzm);
  // Momentum-level dissipation time [s]
  Real2D tau_zm = workspace_scope.view<Real2D>("advance_clubb_core::tau_zm", ngrdcol, nzm);
  // Eddy-diffused scalar flux [m/s * scalar units]
  Real3D wpedsclrp =
      workspace_scope.view<Real3D>("advance_clubb_core::wpedsclrp", ngrdcol, nzm, edsclr_dim);
  // w'^2 sclr' [(m/s)^2 * scalar units]
  Real3D wp2sclrp =
      workspace_scope.view<Real3D>("advance_clubb_core::wp2sclrp", ngrdcol, nzt, sclr_dim);
  // w'sclr'^2 [m/s * scalar units^2]
  Real3D wpsclrp2 =
      workspace_scope.view<Real3D>("advance_clubb_core::wpsclrp2", ngrdcol, nzt, sclr_dim);
  // w'sclr'rt' [m/s * scalar units * kg/kg]
  Real3D wpsclrprtp =
      workspace_scope.view<Real3D>("advance_clubb_core::wpsclrprtp", ngrdcol, nzt, sclr_dim);
  // w'sclr'thl' [m/s * scalar units * K]
  Real3D wpsclrpthlp =
      workspace_scope.view<Real3D>("advance_clubb_core::wpsclrpthlp", ngrdcol, nzt, sclr_dim);
  // Momentum-level PDF width parameter [-]
  Real2D sigma_sqd_w =
      workspace_scope.view<Real2D>("advance_clubb_core::sigma_sqd_w", ngrdcol, nzm);
  // Thermodynamic-level turbulent kinetic energy [m^2/s^2]
  Real2D em_zt = workspace_scope.view<Real2D>("advance_clubb_core::em_zt", ngrdcol, nzt);
  // Mean vertical velocity of PDF component 1 [m/s]
  Real2D w_1_zm = workspace_scope.view<Real2D>("advance_clubb_core::w_1_zm", ngrdcol, nzm);
  // Mean vertical velocity of PDF component 2 [m/s]
  Real2D w_2_zm = workspace_scope.view<Real2D>("advance_clubb_core::w_2_zm", ngrdcol, nzm);
  // Vertical-velocity variance of PDF component 1 [m^2/s^2]
  Real2D varnce_w_1_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::varnce_w_1_zm", ngrdcol, nzm);
  // Vertical-velocity variance of PDF component 2 [m^2/s^2]
  Real2D varnce_w_2_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::varnce_w_2_zm", ngrdcol, nzm);
  // PDF component-1 mixture fraction [-]
  Real2D mixt_frac_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::mixt_frac_zm", ngrdcol, nzm);
  // Cloud-water variance [(kg/kg)^2]
  Real2D rcp2_zt = workspace_scope.view<Real2D>("advance_clubb_core::rcp2_zt", ngrdcol, nzt);
  // Total water before the core advance [kg/kg]
  Real2D rtm_before = workspace_scope.view<Real2D>("advance_clubb_core::rtm_before", ngrdcol, nzt);
  // Liquid potential temperature before advance [K]
  Real2D thlm_before =
      workspace_scope.view<Real2D>("advance_clubb_core::thlm_before", ngrdcol, nzt);
  // Relative humidity over liquid [-]
  Real2D rel_humidity =
      workspace_scope.view<Real2D>("advance_clubb_core::rel_humidity", ngrdcol, nzt);
  // Stability damping multiplier [-]
  Real2D stability_correction =
      workspace_scope.view<Real2D>("advance_clubb_core::stability_correction", ngrdcol, nzm);
  // Inverse dissipation time for C6 pressure terms [1/s]
  Real2D invrs_tau_C6_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_C6_zm", ngrdcol, nzm);
  // Inverse dissipation time for C1 damping terms [1/s]
  Real2D invrs_tau_C1_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_C1_zm", ngrdcol, nzm);
  // Inverse dissipation time for scalar variances [1/s]
  Real2D invrs_tau_xp2_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_xp2_zm", ngrdcol, nzm);
  // Inverse dissipation time for C4 terms [1/s]
  Real2D invrs_tau_C4_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_C4_zm", ngrdcol, nzm);
  // Inverse dissipation time for C14 terms [1/s]
  Real2D invrs_tau_C14_zm =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_C14_zm", ngrdcol, nzm);
  // Richardson-number closure factor [-]
  Real2D Cx_fnc_Richardson =
      workspace_scope.view<Real2D>("advance_clubb_core::Cx_fnc_Richardson", ngrdcol, nzm);
  // Squared norm of vertical horizontal-wind gradient [1/s^2]
  Real2D ddzt_umvm_sqd =
      workspace_scope.view<Real2D>("advance_clubb_core::ddzt_umvm_sqd", ngrdcol, nzm);
  // Buoyancy frequency squared [1/s^2]
  Real2D brunt_vaisala_freq_sqd =
      workspace_scope.view<Real2D>("advance_clubb_core::brunt_vaisala_freq_sqd", ngrdcol, nzm);
  // Dry/moist mixture of buoyancy frequency squared [1/s^2]
  Real2D brunt_vaisala_freq_sqd_mixed = workspace_scope.view<Real2D>(
      "advance_clubb_core::brunt_vaisala_freq_sqd_mixed", ngrdcol, nzm);
  // Vertically smoothed dry/moist buoyancy frequency squared [1/s^2]
  Real2D brunt_vaisala_freq_sqd_smth =
      workspace_scope.view<Real2D>("advance_clubb_core::brunt_vaisala_freq_sqd_smth", ngrdcol, nzm);
  // Thermodynamic-level inverse dissipation time [1/s]
  Real2D invrs_tau_zt =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_zt", ngrdcol, nzt);
  // Inverse dissipation time for third moment of w [1/s]
  Real2D invrs_tau_wp3_zt =
      workspace_scope.view<Real2D>("advance_clubb_core::invrs_tau_wp3_zt", ngrdcol, nzt);
  // Maximum momentum-level dissipation time [s]
  Real2D tau_max_zm = workspace_scope.view<Real2D>("advance_clubb_core::tau_max_zm", ngrdcol, nzm);
  // Maximum thermodynamic-level dissipation time [s]
  Real2D tau_max_zt = workspace_scope.view<Real2D>("advance_clubb_core::tau_max_zt", ngrdcol, nzt);
  bool l_samp_stats_in_pdf_call{}; // Whether this PDF call samples stats [-]
  // LHS coefficient of wp2 splatting [1/s]
  Real2D lhs_splat_wp2 =
      workspace_scope.view<Real2D>("advance_clubb_core::lhs_splat_wp2", ngrdcol, nzm);
  // LHS coefficient of wp3 splatting [1/s]
  Real2D lhs_splat_wp3 =
      workspace_scope.view<Real2D>("advance_clubb_core::lhs_splat_wp3", ngrdcol, nzt);

  // Variables associated with upgradient momentum contributions due to cumuli
  // (inactive in the native source).
  // Real1D Km_Skw_factor; // Factor < 1 reducing momentum diffusivity in skewed layers
  // constexpr real_t Km_Skw_thresh = zero_threshold; // Sk_w threshold for correction
  // constexpr real_t Km_Skw_factor_efold = real_t(0.5); // Exponential e-folding rate
  // constexpr real_t Km_Skw_factor_min = real_t(0.2); // Minimum correction factor

  int advance_order_loop_iter{}; // Advance-stage loop index [-]
  int wprtp_cl_num{};            // Covariance clipping count [-]
  int wpthlp_cl_num{};           // Covariance clipping count [-]
  int upwp_cl_num{};             // Covariance clipping count [-]
  int vpwp_cl_num{};             // Covariance clipping count [-]

  // ----------------------------- Begin Code ------------------------------

  // Select the same full or leapfrog/LMM advance duration as Fortran.
  if (clubb_config_flags.l_lmm_stepping) {
    dt_advance = two * dt;
  } else {
    dt_advance = dt;
  }

  // ----------------------------------------------------------------
  // Test input variables.
  // ----------------------------------------------------------------
  if (clubb_at_least_debug_level_api(2)) {
    parameterization_check(nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,               // In
                           thlm_forcing, rtm_forcing, um_forcing,                 // In
                           vm_forcing, wm_zm, wm_zt, p_in_Pa,                     // In
                           rho_zm, rho, exner, rho_ds_zm,                         // In
                           rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,           // In
                           thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, // In
                           vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,         // In
                           rtm, wprtp, thlm, wpthlp, wp2, wp3,                    // In
                           rtp2, thlp2, rtpthlp,                                  // In
                           "beginning of ",                                       // In
                           wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,            // In
                           sclrp2,                                                // In
                           sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,           // In
                           edsclrm_forcing,                                       // In
                           err_info);                                             // InOut
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << "Fatal error detected in parameterization_check when testing input\n";
      return;
    }
  }

  // C++ Views have no Fortran array-expression temporary. Materialize the
  // native field/dt expression on the execution space at the stats boundary.
  const auto per_dt = [&workspace, dt](const Real2D& field) {
    Real2D values     = workspace.view<Real2D>(Kokkos::WithoutInitializing, field.label() + "/dt",
                                               field.extent(0), field.extent(1));
    const int columns = field.extent_int(0);
    Kokkos::parallel_for(
        "advance_clubb_core::budget_field_per_dt", Kokkos::RangePolicy<>(0, field.size()),
        KOKKOS_LAMBDA(const int cell) {
          const int i  = cell % columns;
          const int k  = cell / columns;
          values(i, k) = field(i, k) / dt;
        });
    return values;
  };

  if (stats.l_sample) {
    auto budget_scope = workspace.scope();
    Kokkos::parallel_for(
        "advance_clubb_core::before_state", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
        KOKKOS_LAMBDA(const int cell) {
          const int i       = cell % ngrdcol;
          const int k       = cell / ngrdcol;
          rtm_before(i, k)  = rtm(i, k);
          thlm_before(i, k) = thlm(i, k);
        });

    stats_update("rfrzm", rfrzm, stats);

    // Set up budget stats variables, preserving source order and arithmetic.
    stats_begin_budget("wp2_bt", per_dt(wp2), stats);
    stats_begin_budget("vp2_bt", per_dt(vp2), stats);
    stats_begin_budget("up2_bt", per_dt(up2), stats);
    stats_begin_budget("wprtp_bt", per_dt(wprtp), stats);
    stats_begin_budget("wpthlp_bt", per_dt(wpthlp), stats);

    if (clubb_config_flags.l_predict_upwp_vpwp) {
      stats_begin_budget("upwp_bt", per_dt(upwp), stats);
      stats_begin_budget("vpwp_bt", per_dt(vpwp), stats);
    }

    stats_begin_budget("rtp2_bt", per_dt(rtp2), stats);
    stats_begin_budget("thlp2_bt", per_dt(thlp2), stats);
    stats_begin_budget("rtpthlp_bt", per_dt(rtpthlp), stats);

    stats_begin_budget("rtm_bt", per_dt(rtm), stats);
    stats_begin_budget("thlm_bt", per_dt(thlm), stats);
    stats_begin_budget("um_bt", per_dt(um), stats);
    stats_begin_budget("vm_bt", per_dt(vm), stats);

    stats_begin_budget("wp3_bt", per_dt(wp3), stats);

    // Set clipping counts to zero.
    wprtp_cl_num  = 0;
    wpthlp_cl_num = 0;
    upwp_cl_num   = 0;
    vpwp_cl_num   = 0;
  }

  set_sfc_value_of_flux_profiles(nzm, ngrdcol, sclr_dim, edsclr_dim, gr,       // In
                                 clubb_config_flags.l_host_applies_sfc_fluxes, // In
                                 clubb_config_flags.l_linearize_pbl_winds,     // In
                                 wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,    // In
                                 upwp_sfc_pert, vpwp_sfc_pert,                 // In
                                 wpsclrp_sfc, wpedsclrp_sfc,                   // In
                                 wpthlp, wprtp, upwp, vpwp,                    // InOut
                                 upwp_pert, vpwp_pert,                         // InOut
                                 wpsclrp, wpedsclrp);                          // InOut

  compute_sigma_sqd_w(workspace, nzm, nzt, ngrdcol, gr,       // In
                      wp3, wp2, thlp2, rtp2,                  // In
                      up2, vp2, wpthlp, wprtp, upwp, vpwp,    // In
                      clubb_params,                           // In
                      clubb_config_flags.l_predict_upwp_vpwp, // In
                      sigma_sqd_w);                           // Out

  if (clubb_config_flags.ipdf_call_placement == ipdf_pre_advance_fields ||
      clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields) {
    // Sample stats in this PDF call for both pre and pre/post placement.
    if (clubb_config_flags.ipdf_call_placement == ipdf_pre_advance_fields) {
      l_samp_stats_in_pdf_call = true;
    } else if (clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields) {
      l_samp_stats_in_pdf_call = true;
    }

    // ########################################################################
    // Call CLUBB's PDF and output PDF parameters and integrated quantities.
    // ########################################################################
    pdf_closure_driver(workspace, gr, nzm, nzt, ngrdcol,            // In
                       dt, hydromet_dim, sclr_dim, sclr_tol,        // In
                       wprtp, thlm, wpthlp, rtp2, rtp3,             // In
                       thlp2, thlp3, rtpthlp, wp2,                  // In
                       wp3, wm_zm, wm_zt,                           // In
                       um, up2, upwp, up3,                          // In
                       vm, vp2, vpwp, vp3,                          // In
                       p_in_Pa, exner,                              // In
                       thv_ds_zm, thv_ds_zt, rtm_ref,               // In
                       wphydrometp,                                 // In
                       wp2hmp, rtphmp_zt, thlphmp_zt,               // In
                       sclrm, wpsclrp, sclrp2,                      // In
                       sclrprtp, sclrpthlp, sclrp3,                 // In
                       p_sfc, l_samp_stats_in_pdf_call,             // In
                       mixt_frac_max_mag, ts_nudge,                 // In
                       rtm_min, rtm_nudge_max_altitude,             // In
                       clubb_params,                                // In
                       clubb_config_flags.iiPDF_type,               // In
                       clubb_config_flags.saturation_formula,       // In
                       clubb_config_flags.l_rtm_nudge,              // In
                       clubb_config_flags.l_trapezoidal_rule_zt,    // In
                       clubb_config_flags.l_trapezoidal_rule_zm,    // In
                       clubb_config_flags.l_call_pdf_closure_twice, // In
                       clubb_config_flags.l_use_cloud_cover,        // In
                       clubb_config_flags.l_rcm_supersat_adj,       // In
                       l_mix_rat_hm,                                // In
                       stats,                                       // InOut
                       rtm, sigma_sqd_w,                            // InOut
                       pdf_implicit_coefs_terms,                    // InOut
                       pdf_params, pdf_params_zm, err_info,         // InOut
                       rcm, cloud_frac,                             // Out
                       ice_supersat_frac, wprcp,                    // Out
                       wpthvp, wp2thvp, wp2up,                      // Out
                       rtpthvp, thlpthvp,                           // Out
                       rcm_in_layer, cloud_cover,                   // Out
                       rcp2_zt, thlprcp,                            // Out
                       rc_coef_zm, sclrpthvp,                       // Out
                       wpup2, wpvp2,                                // Out
                       wp2up2, wp2vp2, wp4,                         // Out
                       wp2rtp, wprtp2, wp2thlp,                     // Out
                       wpthlp2, wprtpthlp,                          // Out
                       uprcp, vprcp,                                // Out
                       w_up_in_cloud, w_down_in_cloud,              // Out
                       cloudy_updraft_frac,                         // Out
                       cloudy_downdraft_frac,                       // Out
                       wp2sclrp, wpsclrp2,                          // Out
                       wpsclrprtp, wpsclrpthlp);                    // Out

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << '\n';
        std::cerr << "Error calling pdf_closure_driver in advance_clubb_core\n";
        return;
      }
    }
  }

  // This feels like an awkward place to have this, but this keeps things BFB
  // with where it was before.
  if (stats.l_sample) {
    if (clubb_config_flags.iiPDF_type == iiPDF_ADG1) {
      zm2zt_api(nzm, nzt, ngrdcol, gr, sigma_sqd_w, sigma_sqd_w_zt, zero_threshold);
    } else {
      Kokkos::deep_copy(sigma_sqd_w_zt, zero);
    }
    stats_update("sigma_sqd_w_zt", sigma_sqd_w_zt, stats);
  }

  calculate_thvm(nzt, ngrdcol,                     // In
                 thlm, rtm, rcm, exner, thv_ds_zt, // In
                 thvm);                            // Out

  // Compute TKE (turbulent kinetic energy).
  Kokkos::parallel_for(
      "advance_clubb_core::tke", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        em(i, k)    = real_t(0.5) * (wp2(i, k) + vp2(i, k) + up2(i, k));
      });

  if (clubb_config_flags.l_call_pdf_closure_twice) {
    // Capture only the five component Views, not the entire owning PDF type.
    const auto pdf_w_1        = pdf_params_zm.w_1;
    const auto pdf_w_2        = pdf_params_zm.w_2;
    const auto pdf_varnce_w_1 = pdf_params_zm.varnce_w_1;
    const auto pdf_varnce_w_2 = pdf_params_zm.varnce_w_2;
    const auto pdf_mixt_frac  = pdf_params_zm.mixt_frac;
    Kokkos::parallel_for(
        "advance_clubb_core::momentum_pdf", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i         = cell % ngrdcol;
          const int k         = cell / ngrdcol;
          w_1_zm(i, k)        = pdf_w_1(i, k);
          w_2_zm(i, k)        = pdf_w_2(i, k);
          varnce_w_1_zm(i, k) = pdf_varnce_w_1(i, k);
          varnce_w_2_zm(i, k) = pdf_varnce_w_2(i, k);
          mixt_frac_zm(i, k)  = pdf_mixt_frac(i, k);
        });
  } else {
    zt2zm_api(nzm, nzt, ngrdcol, gr, pdf_params.w_1, w_1_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, pdf_params.w_2, w_2_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, pdf_params.varnce_w_1, varnce_w_1_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, pdf_params.varnce_w_2, varnce_w_2_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, pdf_params.mixt_frac, mixt_frac_zm);
  }

  if (stats.l_sample) {
    // Form vapor mixing ratio for the sampled statistics boundary.
    Real2D rvm = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                              "advance_clubb_core::rvm", ngrdcol, nzt);
    Kokkos::parallel_for(
        "advance_clubb_core::rvm", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          rvm(i, k)   = rtm(i, k) - rcm(i, k);
        });
    stats_update("rvm", rvm, stats);

    // Output relative humidity (q/q* where q* is saturation over liquid).
    // If neither output is registered, rsat is not computed. Test both names
    // so writing rel_humidity never uses an uninitialized saturation field.
    if (var_on_stats_list(stats, "rel_humidity") || var_on_stats_list(stats, "rsat")) {
      auto T_in_K = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                 "advance_clubb_core::T_in_K", ngrdcol, nzt);
      thlm2T_in_K_api(nzt, ngrdcol, thlm, exner, rcm, T_in_K);
      sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, T_in_K, clubb_config_flags.saturation_formula,
                         rsat);
      // Recompute rsat and relative humidity; they may have changed.
      Kokkos::parallel_for(
          "advance_clubb_core::relative_humidity", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
          KOKKOS_LAMBDA(const int cell) {
            const int i        = cell % ngrdcol;
            const int k        = cell / ngrdcol;
            rel_humidity(i, k) = (rtm(i, k) - rcm(i, k)) / rsat(i, k);
          });
      stats_update("rel_humidity", rel_humidity, stats);
    }
  }

  calc_brunt_vaisala_freq_sqd(workspace, nzm, nzt, ngrdcol, gr, thlm,                    // In
                              exner, rtm, rcm, p_in_Pa, thvm,                            // In
                              ice_supersat_frac,                                         // In
                              clubb_config_flags.saturation_formula,                     // In
                              clubb_config_flags.l_brunt_vaisala_freq_moist,             // In
                              clubb_config_flags.l_use_thvm_in_bv_freq,                  // In
                              clubb_config_flags.l_modify_limiters_for_cnvg_test,        // In
                              Kokkos::subview(clubb_params, Kokkos::ALL, ibv_efold), T0, // In
                              brunt_vaisala_freq_sqd,                                    // Out
                              brunt_vaisala_freq_sqd_mixed,                              // Out
                              brunt_vaisala_freq_sqd_smth,                               // Out
                              &stats); // Optional InOut
  // ----------------------------------------------------------------
  // Compute mixing length and dissipation time.
  // ----------------------------------------------------------------
  calc_ddzt_umvm_sqd(workspace, nzm, nzt, ngrdcol, gr, um, vm,                   // In
                     ddzt_umvm_sqd);                                             // Out
  calc_Lscale(workspace, nzm, nzt, ngrdcol, gr, l_implemented, host_dx, host_dy, // In
              p_in_Pa, exner, rtm, thlm, thvm,                                   // In
              thlp2, rtp2, rtpthlp,                                              // In
              pdf_params, em, thv_ds_zt, lmin,                                   // In
              upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac,              // In
              ufmin, tau_const, sfc_elevation, clubb_params,                     // In
              clubb_config_flags.saturation_formula,                             // In
              clubb_config_flags.l_Lscale_plume_centered,                        // In
              clubb_config_flags.Lscale_method,                                  // In
              clubb_config_flags.l_e3sm_config,                                  // In
              clubb_config_flags.l_smooth_Heaviside_tau_wpxp,                    // In
              clubb_config_flags.l_modify_limiters_for_cnvg_test,                // In
              l_use_invrs_tau_N2_iso,                                            // In
              brunt_vaisala_freq_sqd_smth,                                       // In
              stats, err_info,                                                   // InOut
              invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm,                      // Out
              invrs_tau_wp3_zt,                                                  // Out
              invrs_tau_C1_zm, invrs_tau_C4_zm,                                  // Out
              invrs_tau_C6_zm, invrs_tau_C14_zm,                                 // Out
              tau_max_zm, tau_max_zt, tau_zm,                                    // Out
              Lscale, Lscale_zm, Lscale_up, Lscale_down);                        // Out

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Error calling calc_Lscale in advance_clubb_core\n";
      return;
    }
  }

  // Determine whether to use tau_zm or tau_N2_zm, stability-corrected tau.
  // -dschanen 7 Nov 2014
  if (clubb_config_flags.l_stability_correct_tau_zm) {
    calc_stability_correction(
        nzm, ngrdcol, brunt_vaisala_freq_sqd,                                // In
        Lscale_zm, em,                                                       // In
        Kokkos::subview(clubb_params, Kokkos::ALL, ilambda0_stability_coef), // In
        stability_correction);                                               // Out
    // Stronger damping where the Brunt-Vaisala frequency (N^2) is large.
    Kokkos::parallel_for(
        "advance_clubb_core::stability_damping", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i           = cell % ngrdcol;
          const int k           = cell / ngrdcol;
          invrs_tau_C6_zm(i, k) = invrs_tau_zm(i, k) * stability_correction(i, k);
          invrs_tau_C1_zm(i, k) = invrs_tau_C6_zm(i, k);
        });
  }

  // ----------------------------------------------------------------
  // Eddy diffusivity coefficient.
  // ----------------------------------------------------------------
  // c_K is usually 0.548 (Duynkerke and Driedonks 1987). CLUBB uses a
  // smaller value to fit empirical data. Diffusivity is length times velocity.
  zm2zt_api(nzm, nzt, ngrdcol, gr, em, em_zt, em_min);
  Kokkos::parallel_for(
      "advance_clubb_core::Kh_zt", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        Kh_zt(i, k) = clubb_params(i, ic_K) * Lscale(i, k) * Kokkos::sqrt(em_zt(i, k));
      });

  Kokkos::parallel_for(
      "advance_clubb_core::Kh_zm", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        Kh_zm(i, k) = clubb_params(i, ic_K) * Kokkos::max(Lscale_zm(i, k), zero_threshold) *
                      Kokkos::sqrt(Kokkos::max(em(i, k), em_min));
      });

  wp23_term_splat_lhs(workspace, nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(clubb_params, Kokkos::ALL, iC_wp2_splat), // In
                      brunt_vaisala_freq_sqd_mixed, Lscale_zm, rho_ds_zm,       // In
                      stats,                                                    // InOut
                      lhs_splat_wp2, lhs_splat_wp3);                            // Out

  // ----------------------------------------------------------------
  // Set surface variances before advance_xp2_xpyp or advance_wp2_wp3.
  // Surface effects do not apply when the lowest level is above the ground.
  // Brian Griffin, December 22, 2005.
  // ----------------------------------------------------------------
  calc_sfc_varnce(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_idx,          // In
                  gr, dt, sfc_elevation,                                     // In
                  upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc,                     // In
                  um, vm, Lscale_up, wpsclrp_sfc,                            // In
                  lhs_splat_wp2, tau_zm,                                     // In
                  clubb_config_flags.l_vary_convect_depth,                   // In
                  T0,                                                        // In
                  Kokkos::subview(clubb_params, Kokkos::ALL, iup2_sfc_coef), // In
                  Kokkos::subview(clubb_params, Kokkos::ALL, ia_const),      // In
                  stats,                                                     // InOut
                  wp2, up2, vp2,                                             // InOut
                  thlp2, rtp2, rtpthlp,                                      // InOut
                  sclrp2, sclrprtp, sclrpthlp,                               // InOut
                  err_info);                                                 // InOut

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Error calling calc_sfc_varnce in advance_clubb_core\n";
      return;
    }
  }

  // #######################################################################
  // Advance prognostic variables one timestep.
  // #######################################################################
  // Cx is irrelevant unless a Richardson-number closure uses it. Set it to
  // zero otherwise, to avoid NaN problems.
  if (clubb_config_flags.l_use_C7_Richardson || clubb_config_flags.l_use_C11_Richardson) {
    compute_Cx_fnc_Richardson(workspace, nzm, nzt, ngrdcol, gr,                   // In
                              Lscale_zm, ddzt_umvm_sqd, rho_ds_zm,                // In
                              brunt_vaisala_freq_sqd,                             // In
                              brunt_vaisala_freq_sqd_mixed,                       // In
                              clubb_params,                                       // In
                              clubb_config_flags.l_use_shear_Richardson,          // In
                              clubb_config_flags.l_modify_limiters_for_cnvg_test, // In
                              Cx_fnc_Richardson);                                 // Out
  } else {
    Kokkos::deep_copy(Cx_fnc_Richardson, real_t(0.0));
  }

  // Advance the four main solvers in the native compile-time order.
  for (advance_order_loop_iter = 1; advance_order_loop_iter <= 4; ++advance_order_loop_iter) {
    if (advance_order_loop_iter == order_xm_wpxp) {
      // ----------------------------------------------------------------
      // Advance scalar means/fluxes (rtm/wprtp, thlm/wpthlp, sclrm/wpsclrp).
      // ----------------------------------------------------------------
      // advance_xm_wpxp_bad_wp2 ! Test error comment, DO NOT modify or move
      advance_xm_wpxp(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt_advance, // In
                      sigma_sqd_w, wm_zm, wm_zt, wp2, Lscale_zm,                        // In
                      wp3, Kh_zt, Kh_zm,                                                // In
                      stability_correction,                                             // In
                      invrs_tau_C6_zm, tau_max_zm, wp2rtp, rtpthvp,                     // In
                      rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp,                     // In
                      thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref,                 // In
                      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                            // In
                      invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,                          // In
                      w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm,                     // In
                      mixt_frac_zm, l_implemented, em, wp2sclrp,                        // In
                      sclrpthvp, sclrm_forcing, sclrp2, Cx_fnc_Richardson,              // In
                      pdf_implicit_coefs_terms,                                         // In
                      um_forcing, vm_forcing, ug, vg, wpthvp,                           // In
                      fcor, fcor_y, um_ref, vm_ref, up2, vp2,                           // In
                      uprcp, vprcp, rc_coef_zm,                                         // In
                      clubb_params, nu_vert_res_dep, ts_nudge,                          // In
                      clubb_config_flags.iiPDF_type,                                    // In
                      clubb_config_flags.penta_solve_method,                            // In
                      clubb_config_flags.tridiag_solve_method,                          // In
                      clubb_config_flags.fill_holes_type,                               // In
                      clubb_config_flags.l_predict_upwp_vpwp,                           // In
                      clubb_config_flags.l_ho_nontrad_coriolis,                         // In
                      clubb_config_flags.l_ho_trad_coriolis,                            // In
                      clubb_config_flags.l_diffuse_rtm_and_thlm,                        // In
                      clubb_config_flags.l_stability_correct_Kh_N2_zm,                  // In
                      clubb_config_flags.l_godunov_upwind_wpxp_ta,                      // In
                      clubb_config_flags.l_upwind_xm_ma,                                // In
                      clubb_config_flags.l_uv_nudge,                                    // In
                      clubb_config_flags.Lscale_method == Lscale_tau,                   // In
                      clubb_config_flags.l_use_C7_Richardson,                           // In
                      clubb_config_flags.l_lmm_stepping,                                // In
                      clubb_config_flags.l_enable_relaxed_clipping,                     // In
                      clubb_config_flags.l_linearize_pbl_winds,                         // In
                      clubb_config_flags.l_mono_flux_lim_thlm,                          // In
                      clubb_config_flags.l_mono_flux_lim_rtm,                           // In
                      clubb_config_flags.l_mono_flux_lim_um,                            // In
                      clubb_config_flags.l_mono_flux_lim_vm,                            // In
                      clubb_config_flags.l_mono_flux_lim_spikefix,                      // In
                      wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num,            // InOut
                      stats,                                                            // InOut
                      rtm, wprtp, thlm, wpthlp,                                         // InOut
                      sclrm, wpsclrp, um, upwp, vm, vpwp,                               // InOut
                      um_pert, vm_pert, upwp_pert, vpwp_pert, err_info);                // InOut

      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << "Error calling advance_xm_wpxp in advance_clubb_core\n";
          return;
        }
      }

      // Vince Larson clipped rcm to prevent rvm < 0, 5 Apr 2008.
      // This requires rtm >= 0. rcm_in_layer affects only radiation and is
      // deliberately not clipped/recomputed here (6 Aug 2009).
      clip_rcm(nzt, ngrdcol, rtm,              // In
               "rtm < rcm in advance_xm_wpxp", // In
               rcm);                           // InOut
    } else if (advance_order_loop_iter == order_xp2_xpyp) {
      // ----------------------------------------------------------------
      // Advance total-water/temperature variances and covariance, and
      // horizontal-wind variances. Vertical-wind moments advance separately
      // according to the compile-time order. Certain cases need these
      // tendencies at shorter steps; advancing this stage first can allow a
      // longer timestep.
      // ----------------------------------------------------------------
      advance_xp2_xpyp(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx, // In
                       invrs_tau_xp2_zm, invrs_tau_C4_zm,                              // In
                       invrs_tau_C14_zm, wm_zm,                                        // In
                       rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,                       // In
                       wp2, wp3, upwp, vpwp,                                           // In
                       sigma_sqd_w, wprtp2, wpthlp2,                                   // In
                       wprtpthlp, Kh_zt, rtp2_forcing,                                 // In
                       thlp2_forcing, rtpthlp_forcing,                                 // In
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                          // In
                       thv_ds_zm, cloud_frac,                                          // In
                       pdf_implicit_coefs_terms,                                       // In
                       dt_advance, fcor_y,                                             // In
                       sclrm, wpsclrp,                                                 // In
                       wpsclrp2, wpsclrprtp, wpsclrpthlp,                              // In
                       lhs_splat_wp2,                                                  // In
                       clubb_params, nu_vert_res_dep,                                  // In
                       clubb_config_flags.iiPDF_type,                                  // In
                       clubb_config_flags.tridiag_solve_method,                        // In
                       clubb_config_flags.fill_holes_type,                             // In
                       clubb_config_flags.l_ho_nontrad_coriolis,                       // In
                       clubb_config_flags.l_min_xp2_from_corr_wx,                      // In
                       clubb_config_flags.l_C2_cloud_frac,                             // In
                       clubb_config_flags.l_upwind_xpyp_ta,                            // In
                       clubb_config_flags.l_godunov_upwind_xpyp_ta,                    // In
                       clubb_config_flags.l_lmm_stepping,                              // In
                       l_implemented,                                                  // In
                       stats,                                                          // InOut
                       rtp2, thlp2, rtpthlp, up2, vp2,                                 // InOut
                       sclrp2, sclrprtp, sclrpthlp, err_info);                         // InOut

      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << "Error calling advance_xp2_xpyp in advance_clubb_core\n";
          return;
        }
      }

      // Clip flux covariances after updating their scalar denominators xp2.
      clip_covars_denom(workspace, nzm, ngrdcol, sclr_dim,        // In
                        dt,                                       // In
                        rtp2, thlp2, up2, vp2, wp2,               // In
                        sclrp2,                                   // In
                        clubb_config_flags.l_linearize_pbl_winds, // In
                        clubb_config_flags.l_predict_upwp_vpwp,   // In
                        stats,                                    // InOut
                        wprtp_cl_num, wpthlp_cl_num,              // InOut
                        upwp_cl_num, vpwp_cl_num,                 // InOut
                        wprtp, wpthlp, upwp, vpwp, wpsclrp,       // InOut
                        upwp_pert, vpwp_pert);                    // InOut
    } else if (advance_order_loop_iter == order_wp2_wp3) {
      // ----------------------------------------------------------------
      // Advance the second and third moments of vertical velocity.
      // ----------------------------------------------------------------
      // advance_wp2_wp3_bad_wp2 ! Test error comment, DO NOT modify or move
      advance_wp2_wp3(workspace, nzm, nzt, ngrdcol, gr, dt_advance,         // In
                      sfc_elevation, fcor_y, sigma_sqd_w, wm_zm,            // In
                      wm_zt,                                                // In
                      wpup2, wpvp2, wp2up2, wp2vp2, wp4,                    // In
                      wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,           // In
                      em, Kh_zm, Kh_zt, invrs_tau_C4_zm,                    // In
                      invrs_tau_wp3_zt, invrs_tau_C1_zm,                    // In
                      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                // In
                      invrs_rho_ds_zt, thv_ds_zm,                           // In
                      thv_ds_zt, pdf_params.mixt_frac, Cx_fnc_Richardson,   // In
                      lhs_splat_wp2, lhs_splat_wp3,                         // In
                      pdf_implicit_coefs_terms,                             // In
                      wprtp, wpthlp, rtp2, thlp2,                           // In
                      clubb_params, nu_vert_res_dep,                        // In
                      clubb_config_flags.iiPDF_type,                        // In
                      clubb_config_flags.penta_solve_method,                // In
                      clubb_config_flags.fill_holes_type,                   // In
                      clubb_config_flags.l_min_wp2_from_corr_wx,            // In
                      clubb_config_flags.l_upwind_xm_ma,                    // In
                      clubb_config_flags.l_standard_term_ta,                // In
                      clubb_config_flags.l_partial_upwind_wp3,              // In
                      clubb_config_flags.l_damp_wp2_using_em,               // In
                      clubb_config_flags.l_use_C11_Richardson,              // In
                      clubb_config_flags.l_damp_wp3_Skw_squared,            // In
                      clubb_config_flags.l_lmm_stepping,                    // In
                      clubb_config_flags.l_use_tke_in_wp3_pr_turb_term,     // In
                      clubb_config_flags.l_use_tke_in_wp2_wp3_K_dfsn,       // In
                      clubb_config_flags.l_use_wp3_lim_with_smth_Heaviside, // In
                      clubb_config_flags.l_wp2_fill_holes_tke,              // In
                      clubb_config_flags.l_ho_nontrad_coriolis,             // In
                      l_implemented,                                        // In
                      stats,                                                // InOut
                      up2, vp2, wp2, wp3, err_info);                        // InOut

      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << "Error calling advance_wp2_wp3 in advance_clubb_core\n";
          return;
        }
      }

      // Clip flux covariances after updating the vertical-wind variance wp2.
      clip_covars_denom(workspace, nzm, ngrdcol, sclr_dim,        // In
                        dt,                                       // In
                        rtp2, thlp2, up2, vp2, wp2,               // In
                        sclrp2,                                   // In
                        clubb_config_flags.l_linearize_pbl_winds, // In
                        clubb_config_flags.l_predict_upwp_vpwp,   // In
                        stats,                                    // InOut
                        wprtp_cl_num, wpthlp_cl_num,              // InOut
                        upwp_cl_num, vpwp_cl_num,                 // InOut
                        wprtp, wpthlp, upwp, vpwp, wpsclrp,       // InOut
                        upwp_pert, vpwp_pert);                    // InOut
    } else if (advance_order_loop_iter == order_windm) {
      // ----------------------------------------------------------------
      // Advance mean winds, eddy-diffused scalar means and their fluxes.
      // ----------------------------------------------------------------
      advance_windm_edsclrm(workspace, nzm, nzt, ngrdcol, edsclr_dim, gr, sclr_idx, dt, // In
                            wm_zt, Kh_zm, clubb_params,                                 // In
                            ug, vg, um_ref, vm_ref,                                     // In
                            wp2, up2, vp2, um_forcing, vm_forcing,                      // In
                            edsclrm_forcing, p_in_Pa,                                   // In
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,                      // In
                            fcor, l_implemented,                                        // In
                            nu_vert_res_dep, ts_nudge,                                  // In
                            clubb_config_flags.tridiag_solve_method,                    // In
                            clubb_config_flags.l_predict_upwp_vpwp,                     // In
                            clubb_config_flags.l_upwind_xm_ma,                          // In
                            clubb_config_flags.l_uv_nudge,                              // In
                            clubb_config_flags.l_lmm_stepping,                          // In
                            clubb_config_flags.l_linearize_pbl_winds,                   // In
                            clubb_config_flags.l_do_expldiff_rtm_thlm,                  // In
                            clubb_config_flags.fill_holes_type,                         // In
                            upwp_cl_num, vpwp_cl_num,                                   // InOut
                            stats,                                                      // InOut
                            um, vm, thlm, rtm, edsclrm,                                 // InOut
                            upwp, vpwp, wpedsclrp,                                      // InOut
                            um_pert, vm_pert, upwp_pert,                                // InOut
                            vpwp_pert, err_info);                                       // InOut

      if (clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << "Error calling advance_windm_edsclrm in advance_clubb_core\n";
          return;
        }
      }
    }
  }

  // ----------------------------------------------------------------
  // Advance or diagnose <thl'^3>, <rt'^3>, and <sclr'^3>.
  // ----------------------------------------------------------------
  if (l_advance_xp3 && clubb_config_flags.iiPDF_type != iiPDF_ADG1) {
    // Use the simplified third-moment equation, either prognostic or steady-state.
    advance_xp3(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt, // In
                rtm, thlm, rtp2, thlp2, wprtp,                            // In
                wpthlp, wprtp2, wpthlp2, rho_ds_zm,                       // In
                invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,                // In
                sclrm, sclrp2, wpsclrp, wpsclrp2,                         // In
                wp2, wp3, upwp, vpwp, up2, vp2,                           // In
                thvm, clubb_params,                                       // In
                clubb_config_flags.l_lmm_stepping,                        // In
                stats,                                                    // InOut
                rtp3, thlp3, sclrp3, up3, vp3);                           // InOut
  } else {
    diagnose_xp3(workspace, nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, // In
                 clubb_config_flags.iiPDF_type, clubb_params,          // In
                 wp2, wp3, thvm,                                       // In
                 wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2,     // In
                 sigma_sqd_w, wpsclrp, sclrp2,                         // In
                 stats,                                                // InOut
                 rtp3, thlp3, up3, vp3,                                // InOut
                 sclrp3);                                              // InOut
  }

  if (clubb_config_flags.ipdf_call_placement == ipdf_post_advance_fields ||
      clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields) {
    // Sample post-only PDF calls. For pre/post placement the first call
    // already sampled stats, so this second call must not sample again.
    if (clubb_config_flags.ipdf_call_placement == ipdf_post_advance_fields) {
      l_samp_stats_in_pdf_call = true;
    } else if (clubb_config_flags.ipdf_call_placement == ipdf_pre_post_advance_fields) {
      l_samp_stats_in_pdf_call = false;
    }

    // ########################################################################
    // Diagnose PDF parameters and integrated quantities from the updated
    // prognostic moments: clouds, buoyancy and turbulent advection.
    // ########################################################################
    compute_sigma_sqd_w(workspace, nzm, nzt, ngrdcol, gr,           // In
                        wp3, wp2, thlp2, rtp2,                      // In
                        up2, vp2, wpthlp, wprtp, upwp, vpwp,        // In
                        clubb_params,                               // In
                        clubb_config_flags.l_predict_upwp_vpwp,     // In
                        sigma_sqd_w);                               // Out
    pdf_closure_driver(workspace, gr, nzm, nzt, ngrdcol,            // In
                       dt, hydromet_dim, sclr_dim, sclr_tol,        // In
                       wprtp, thlm, wpthlp, rtp2, rtp3,             // In
                       thlp2, thlp3, rtpthlp, wp2,                  // In
                       wp3, wm_zm, wm_zt,                           // In
                       um, up2, upwp, up3,                          // In
                       vm, vp2, vpwp, vp3,                          // In
                       p_in_Pa, exner,                              // In
                       thv_ds_zm, thv_ds_zt, rtm_ref,               // In
                       wphydrometp,                                 // In
                       wp2hmp, rtphmp_zt, thlphmp_zt,               // In
                       sclrm, wpsclrp, sclrp2,                      // In
                       sclrprtp, sclrpthlp, sclrp3,                 // In
                       p_sfc, l_samp_stats_in_pdf_call,             // In
                       mixt_frac_max_mag, ts_nudge,                 // In
                       rtm_min, rtm_nudge_max_altitude,             // In
                       clubb_params,                                // In
                       clubb_config_flags.iiPDF_type,               // In
                       clubb_config_flags.saturation_formula,       // In
                       clubb_config_flags.l_rtm_nudge,              // In
                       clubb_config_flags.l_trapezoidal_rule_zt,    // In
                       clubb_config_flags.l_trapezoidal_rule_zm,    // In
                       clubb_config_flags.l_call_pdf_closure_twice, // In
                       clubb_config_flags.l_use_cloud_cover,        // In
                       clubb_config_flags.l_rcm_supersat_adj,       // In
                       l_mix_rat_hm,                                // In
                       stats,                                       // InOut
                       rtm, sigma_sqd_w,                            // InOut
                       pdf_implicit_coefs_terms,                    // InOut
                       pdf_params, pdf_params_zm, err_info,         // InOut
                       rcm, cloud_frac,                             // Out
                       ice_supersat_frac, wprcp,                    // Out
                       wpthvp, wp2thvp, wp2up,                      // Out
                       rtpthvp, thlpthvp,                           // Out
                       rcm_in_layer, cloud_cover,                   // Out
                       rcp2_zt, thlprcp,                            // Out
                       rc_coef_zm, sclrpthvp,                       // Out
                       wpup2, wpvp2,                                // Out
                       wp2up2, wp2vp2, wp4,                         // Out
                       wp2rtp, wprtp2, wp2thlp,                     // Out
                       wpthlp2, wprtpthlp,                          // Out
                       uprcp, vprcp,                                // Out
                       w_up_in_cloud, w_down_in_cloud,              // Out
                       cloudy_updraft_frac,                         // Out
                       cloudy_downdraft_frac,                       // Out
                       wp2sclrp, wpsclrp2,                          // Out
                       wpsclrprtp, wpsclrpthlp);                    // Out

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << '\n';
        std::cerr << "Error calling pdf_closure_driver in advance_clubb_core\n";
        return;
      }
    }
  }

  // CLUBB_CAM qclvar is outside the standalone/JAX feature scope. Its source
  // branch only copies rcp2_zt; no host-model-only argument is exposed here.

  // #######################################################################
  // Accumulate statistics and close the same budget windows.
  // #######################################################################
  if (stats.l_sample) {
    auto budget_scope = workspace.scope();
    stats_accumulate(workspace, nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, gr, dt,    // In
                     l_implemented, clubb_config_flags.l_host_applies_sfc_fluxes,   // In
                     clubb_config_flags.l_stability_correct_tau_zm,                 // In
                     clubb_params,                                                  // In
                     um, vm, upwp, vpwp, up2, vp2,                                  // In
                     thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing, // In
                     wpthlp_sfc, wprtp_sfc, wprtp, wpthlp,                          // In
                     wp2, wp3, rtp2, rtp3, thlp2, thlp3,                            // In
                     rtpthlp,                                                       // In
                     p_in_Pa, exner, rho, rho_zm,                                   // In
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,                    // In
                     wm_zt, wm_zm, rcm,                                             // In
                     cloud_frac,                                                    // In
                     thvm, ug, vg,                                                  // In
                     ddzt_umvm_sqd, stability_correction,                           // In
                     Kh_zt,                                                         // In
                     rsat,                                                          // In
                     Kh_zm,                                                         // In
                     em,                                                            // In
                     sclrm, sclrp2,                                                 // In
                     sclrprtp, sclrpthlp, sclrm_forcing,                            // In
                     wpsclrp, wpedsclrp, edsclrm,                                   // In
                     edsclrm_forcing,                                               // In
                     clubb_config_flags.saturation_formula,                         // In
                     stats);                                                        // InOut
    stats_finalize_budget("wp2_bt", per_dt(wp2), stats);
    stats_finalize_budget("vp2_bt", per_dt(vp2), stats);
    stats_finalize_budget("up2_bt", per_dt(up2), stats);
    stats_finalize_budget("wprtp_bt", per_dt(wprtp), stats);
    stats_finalize_budget("wpthlp_bt", per_dt(wpthlp), stats);
    stats_finalize_budget("rtp2_bt", per_dt(rtp2), stats);
    stats_finalize_budget("thlp2_bt", per_dt(thlp2), stats);
    stats_finalize_budget("rtpthlp_bt", per_dt(rtpthlp), stats);
    stats_finalize_budget("rtm_bt", per_dt(rtm), stats);
    stats_finalize_budget("thlm_bt", per_dt(thlm), stats);
    stats_finalize_budget("um_bt", per_dt(um), stats);
    stats_finalize_budget("vm_bt", per_dt(vm), stats);
    stats_finalize_budget("wp3_bt", per_dt(wp3), stats);
    if (clubb_config_flags.l_predict_upwp_vpwp) {
      stats_finalize_budget("upwp_bt", per_dt(upwp), stats);
      stats_finalize_budget("vpwp_bt", per_dt(vpwp), stats);
    }
  }

  if (clubb_at_least_debug_level_api(2)) {
    parameterization_check(nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,               // In
                           thlm_forcing, rtm_forcing, um_forcing,                 // In
                           vm_forcing, wm_zm, wm_zt, p_in_Pa,                     // In
                           rho_zm, rho, exner, rho_ds_zm,                         // In
                           rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,           // In
                           thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, // In
                           vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,         // In
                           rtm, wprtp, thlm, wpthlp, wp2, wp3,                    // In
                           rtp2, thlp2, rtpthlp,                                  // In
                           "end of ",                                             // In
                           wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,            // In
                           sclrp2,                                                // In
                           sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,           // In
                           edsclrm_forcing,                                       // In
                           err_info);                                             // InOut
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << "Error occurred during parameterization_check at end of advance_clubb_core\n";
      return;
    }
  }
}

namespace {
void set_sfc_value_of_flux_profiles(int nzm, int ngrdcol, int sclr_dim, int edsclr_dim,
                                    const grid& gr, // In
                                    bool l_host_applies_sfc_fluxes,
                                    bool l_linearize_pbl_winds, // In
                                    const Real1D& wpthlp_sfc, const Real1D& wprtp_sfc,
                                    const Real1D& upwp_sfc,                                   // In
                                    const Real1D& vpwp_sfc,                                   // In
                                    const Real1D& upwp_sfc_pert, const Real1D& vpwp_sfc_pert, // In
                                    const Real2D& wpsclrp_sfc, const Real2D& wpedsclrp_sfc,   // In
                                    Real2D wpthlp, Real2D wprtp, Real2D upwp, Real2D vpwp, // InOut
                                    Real2D upwp_pert, Real2D vpwp_pert,                    // InOut
                                    Real3D wpsclrp, Real3D wpedsclrp) {                    // InOut
  // Set or clear surface flux profiles according to host surface-flux ownership.
  // Native momentum-level indices are one-based; convert once for device access.
  const int k_lb_zm = gr.k_lb_zm;

  // SET SURFACE VALUES OF FLUXES (BROUGHT IN).
  // Only do this when the host does not apply these fluxes elsewhere (e.g. WRF).
  // Otherwise the _sfc inputs only set surface variances. -dschanen 8 Sept 2009
  if (!l_host_applies_sfc_fluxes) {
    Kokkos::parallel_for(
        "set_sfc_value_of_flux_profiles::fluxes", Kokkos::RangePolicy<>(0, ngrdcol),
        KOKKOS_LAMBDA(const int i) {
          wpthlp(i, k_lb_zm) = wpthlp_sfc(i);
          wprtp(i, k_lb_zm)  = wprtp_sfc(i);
          upwp(i, k_lb_zm)   = upwp_sfc(i);
          vpwp(i, k_lb_zm)   = vpwp_sfc(i);
        });

    if (l_linearize_pbl_winds) {
      Kokkos::parallel_for(
          "set_sfc_value_of_flux_profiles::perturbed_fluxes", Kokkos::RangePolicy<>(0, ngrdcol),
          KOKKOS_LAMBDA(const int i) {
            upwp_pert(i, k_lb_zm) = upwp_sfc_pert(i);
            vpwp_pert(i, k_lb_zm) = vpwp_sfc_pert(i);
          });
    }

    // Set fluxes for passive scalars, if enabled.
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "set_sfc_value_of_flux_profiles::scalar_fluxes",
          Kokkos::RangePolicy<>(0, ngrdcol * sclr_dim), KOKKOS_LAMBDA(const int cell) {
            const int i               = cell % ngrdcol;
            const int sclr            = cell / ngrdcol;
            wpsclrp(i, k_lb_zm, sclr) = wpsclrp_sfc(i, sclr);
          });
    }

    if (edsclr_dim > 0) {
      // This local field is not saved between timesteps. Clear it before
      // overwriting its surface, preserving the native all-level definition.
      Kokkos::deep_copy(wpedsclrp, zero);
      Kokkos::parallel_for(
          "set_sfc_value_of_flux_profiles::eddy_scalar_fluxes",
          Kokkos::RangePolicy<>(0, ngrdcol * edsclr_dim), KOKKOS_LAMBDA(const int cell) {
            const int i                   = cell % ngrdcol;
            const int edsclr              = cell / ngrdcol;
            wpedsclrp(i, k_lb_zm, edsclr) = wpedsclrp_sfc(i, edsclr);
          });
    }
  } else {
    Kokkos::parallel_for(
        "set_sfc_value_of_flux_profiles::host_fluxes", Kokkos::RangePolicy<>(0, ngrdcol),
        KOKKOS_LAMBDA(const int i) {
          wpthlp(i, k_lb_zm) = real_t(0.0);
          wprtp(i, k_lb_zm)  = real_t(0.0);
          upwp(i, k_lb_zm)   = real_t(0.0);
          vpwp(i, k_lb_zm)   = real_t(0.0);
        });

    // Set fluxes for passive scalars, if enabled.
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "set_sfc_value_of_flux_profiles::host_scalar_fluxes",
          Kokkos::RangePolicy<>(0, ngrdcol * sclr_dim), KOKKOS_LAMBDA(const int cell) {
            const int i               = cell % ngrdcol;
            const int sclr            = cell / ngrdcol;
            wpsclrp(i, k_lb_zm, sclr) = real_t(0.0);
          });
    }

    if (edsclr_dim > 0) {
      // This local field is not saved between timesteps; clear every level.
      Kokkos::deep_copy(wpedsclrp, zero);
      Kokkos::parallel_for(
          "set_sfc_value_of_flux_profiles::host_eddy_scalar_fluxes",
          Kokkos::RangePolicy<>(0, ngrdcol * edsclr_dim), KOKKOS_LAMBDA(const int cell) {
            const int i                   = cell % ngrdcol;
            const int edsclr              = cell / ngrdcol;
            wpedsclrp(i, k_lb_zm, edsclr) = real_t(0.0);
          });
    }
  }
}

} // namespace
} // namespace clubb_kokkos
