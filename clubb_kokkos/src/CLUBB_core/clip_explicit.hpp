// Limit moments and cloud water with the Native clipping contracts.
//
// The six public routines preserve Native field/counter ownership and sampled
// budgets. All physical fields use column-first LayoutLeft execution-space Views;
// optional threshold_hi remains a scalar. Generic covariance/variance overloads
// adapt actual caller subviews through the existing array owner once per call.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
#include "stats_netcdf.hpp"
#include <string>
namespace clubb_kokkos {
// Native public selectors and clipping-event counts.
inline constexpr int clip_rtp2        = 1;  // Named constant for rtp2 clipping
inline constexpr int clip_thlp2       = 2;  // Named constant for thlp2 clipping
inline constexpr int clip_rtpthlp     = 3;  // Named constant for rtpthlp clipping
inline constexpr int clip_up2         = 5;  // Named constant for up2 clipping
inline constexpr int clip_vp2         = 6;  // Named constant for vp2 clipping
inline constexpr int clip_wprtp       = 8;  // Named constant for wprtp clipping
inline constexpr int clip_wpthlp      = 9;  // Named constant for wpthlp clipping
inline constexpr int clip_upwp        = 10; // Named constant for upwp clipping
inline constexpr int clip_vpwp        = 11; // Named constant for vpwp clipping
inline constexpr int clip_wp2         = 12; // Named constant for wp2 clipping
inline constexpr int clip_wpsclrp     = 13; // Named constant for wp scalar clipping
inline constexpr int clip_sclrp2      = 14; // Named constant for sclrp2 clipping
inline constexpr int clip_sclrprtp    = 15; // Named constant for sclrprtp clipping
inline constexpr int clip_sclrpthlp   = 16; // Named constant for sclrpthlp clipping
inline constexpr int clip_wphydrometp = 17; // Named constant for wphydrometp clipping
inline constexpr int wprtp_cl_max     = 3;  // Number of wprtp clipping events per timestep
inline constexpr int wpthlp_cl_max    = 3;  // Number of wpthlp clipping events per timestep
inline constexpr int upwp_cl_max      = 3;  // Number of upwp clipping events per timestep
inline constexpr int vpwp_cl_max      = 3;  // Number of vpwp clipping events per timestep
// Native clip_scalar = 7 remains commented out.
// rtp2: r_t'^2         [(kg/kg)^2]
// thlp2: theta_l'^2     [K^2]
// up2: u'^2           [m^2/s^2]
// vp2: v'^2           [m^2/s^2]
// wp2: w'^2           [m^2/s^2]
// sclrp2: sclr'^2  [{units vary}^2]
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_predict_upwp_vpwp: Prognose u'w' and v'w'
// wprtp: w'r_t'        [(kg/kg) m/s]
// wpthlp: w'theta_l'    [K m/s]
// upwp: u'w'          [m^2/s^2]
// vpwp: v'w'          [m^2/s^2]
// wpsclrp: w'sclr'         [units m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
void clip_covars_denom(ScratchWorkspace& workspace,                                   // InOut
                       int nzm, int ngrdcol, int sclr_dim,                            // In
                       real_t dt,                                                     // In
                       Real2D rtp2, Real2D thlp2, Real2D up2, Real2D vp2, Real2D wp2, // In
                       Real3D sclrp2,                                                 // In
                       bool l_linearize_pbl_winds,                                    // In
                       bool l_predict_upwp_vpwp,                                      // In
                       stats_type& stats,                                             // InOut
                       int& wprtp_cl_num, int& wpthlp_cl_num, int& upwp_cl_num,
                       int& vpwp_cl_num, // InOut
                       Real2D wprtp, Real2D wpthlp, Real2D upwp, Real2D vpwp,
                       Real3D wpsclrp,                     // InOut
                       Real2D upwp_pert, Real2D vpwp_pert) // InOut
    ;

// solve_type: Variable being solved; used for STATS.
// xp2: Variance of x, x'^2 (momentum levels)         [{x units}^2]
// yp2: Variance of y, y'^2 (momentum levels)         [{y units}^2]
// xpyp: Covariance of x and y, x'y' (momentum levels) [{x units}*{y units}]
// xpyp_chnge: Net change in x'y' due to clipping [{x units}*{y units}]
void clip_covar(int nzm, int ngrdcol, int solve_type, Real2D xp2, Real2D yp2, // In
                Real2D xpyp,                                                  // InOut
                Real2D xpyp_chnge)                                            // Out
    ;

// solve_type: Variable being solved; used for STATS
// dt: Model timestep; used here for STATS     [s]
// threshold_lo: Minimum value of x'^2                   [{x units}^2]
// xp2: Variance of x, x'^2 (momentum levels)   [{x units}^2]
// threshold_hi: Maximum value of x'^2                   [{x units}^2]
void clip_variance(ScratchWorkspace& workspace,                                     // InOut
                   int nzm, int ngrdcol, const grid& gr, int solve_type, real_t dt, // In
                   Real2D threshold_lo,                                             // In
                   stats_type& stats,                                               // InOut
                   Real2D xp2,                                                      // InOut
                   const real_t* threshold_hi = nullptr)                            // Optional In
    ;

// dt: Model timestep; used here for STATS        [s]
// sfc_elevation: Elevation of ground level                  [m AMSL]
// Skw_max_mag: Maximum allowable magnitude of Skewness    [-]
// wp2_zt: w'^2 interpolated to thermodyamic levels   [m^2/s^2]
// wp3: w'^3 (thermodynamic levels)                [m^3/s^3]
void clip_skewness(ScratchWorkspace& workspace,                                           // InOut
                   int nzt, int ngrdcol, const grid& gr, real_t dt, Real1D sfc_elevation, // In
                   Real1D Skw_max_mag, Real2D wp2_zt,                                     // In
                   bool l_use_wp3_lim_with_smth_Heaviside,                                // In
                   stats_type& stats,                                                     // InOut
                   Real2D wp3)                                                            // InOut
    ;

// sfc_elevation: Elevation of ground level                  [m AMSL]
// Skw_max_mag: Maximum allowable magnitude of Skewness    [-]
// wp2_zt: w'^2 interpolated to thermodyamic levels   [m^2/s^2]
// wp3: w'^3 (thermodynamic levels)                [m^3/s^3]
void clip_skewness_core(int nzt, int ngrdcol, const grid& gr, Real1D sfc_elevation, // In
                        Real1D Skw_max_mag, Real2D wp2_zt,                          // In
                        bool l_use_wp3_lim_with_smth_Heaviside,                     // In
                        Real2D wp3)                                                 // InOut
    ;

// rtm: Total water mixing ratio             [kg/kg]
// rcm: Cloud water mixing ratio  [kg/kg]
void clip_rcm(int nzt, int ngrdcol, Real2D rtm, // In
              const std::string& message,       // In
              Real2D rcm)                       // InOut
    ;

// Actual moment callers pass scalar sections; adapt at this public API boundary.
template <class X, class Y, class XY, class Delta>
void clip_covar(int nzm, int ngrdcol, int solve_type, X xp2, Y yp2, XY xpyp, Delta xpyp_chnge) {
  auto xp2_exec    = to_execution_space(xp2);
  auto yp2_exec    = to_execution_space(yp2);
  auto xpyp_exec   = to_execution_space(xpyp);
  auto change_exec = to_execution_space(xpyp_chnge, false);
  clip_covar(nzm, ngrdcol, solve_type, xp2_exec, yp2_exec, xpyp_exec, change_exec);
  from_execution_space(xpyp, xpyp_exec);
  from_execution_space(xpyp_chnge, change_exec);
}
template <class Threshold, class Field>
void clip_variance(ScratchWorkspace& workspace, int nzm, int ngrdcol, const grid& gr,
                   int solve_type, real_t dt, Threshold threshold_lo, stats_type& stats, Field xp2,
                   const real_t* threshold_hi = nullptr) {
  auto workspace_scope = workspace.scope();
  auto threshold_exec  = to_execution_space(workspace, threshold_lo);
  auto xp2_exec        = to_execution_space(workspace, xp2);
  clip_variance(workspace, nzm, ngrdcol, gr, solve_type, dt, threshold_exec, stats, xp2_exec,
                threshold_hi);
  from_execution_space(workspace, xp2, xp2_exec);
}
} // namespace clubb_kokkos
