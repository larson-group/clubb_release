// Public Native-order mean/flux timestep interface; see the source description.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "parameters_tunable.hpp"
#include "pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// dt: Timestep                                 [s]
// sigma_sqd_w: sigma_sqd_w on momentum levels           [-]
// wm_zm: w wind component on momentum levels      [m/s]
// wm_zt: w wind component on thermodynamic levels [m/s]
// wp2: w'^2 (momentum levels)                   [m^2/s^2]
// Lscale_zm: Turbulent mixing length interp to m levs [m]
// wp3: <w'^3> on thermodynamic levels           [m^3/s^3]
// Kh_zt: Eddy diffusivity on thermodynamic levels [m^2/s]
// Kh_zm: Eddy diffusivity on momentum levels
// stability_correction: Stability correction factor
// invrs_tau_C6_zm: Inverse time-scale on mom. levels applied to C6 term [1/s]
// tau_max_zm: Max. allowable eddy dissipation time scale on m-levs  [s]
// wp2rtp: <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
// rtpthvp: r_t'th_v' (momentum levels)              [(kg/kg) K]
// rtm_forcing: r_t forcing (thermodynamic levels)       [(kg/kg)/s]
// wprtp_forcing: <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
// rtm_ref: rtm for nudging                          [kg/kg]
// wp2thlp: <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
// thlpthvp: th_l'th_v' (momentum levels)             [K^2]
// thlm_forcing: th_l forcing (thermodynamic levels)      [K/s]
// wpthlp_forcing: <w'th_l'> forcing (momentum levels)      [K/s^2]
// thlm_ref: thlm for nudging                         [K]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on moment. levs. [K]
// rtp2: r_t'^2 (momentum levels)                 [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                [K^2]
// w_1_zm: Mean w (1st PDF component)              [m/s]
// w_2_zm: Mean w (2nd PDF component)              [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)       [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)       [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// l_implemented: True if CLUBB is being implemented and run in a host model
// em: Turbulent Kinetic Energy (TKE)           [m^2/s^2]
// wp2sclrp: <w'^2 sclr'> (thermodynamic levels)   [Units vary]
// sclrpthvp: <sclr' th_v'> (momentum levels)       [Units vary]
// sclrm_forcing: sclrm forcing (thermodynamic levels)  [Units vary]
// sclrp2: For clipping Vince Larson             [Units vary]
// Cx_fnc_Richardson: Cx_fnc computed from Richardson_num       [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// um_forcing: <u> forcing term (thermodynamic levels)      [m/s^2]
// vm_forcing: <v> forcing term (thermodynamic levels)      [m/s^2]
// ug: <u> geostrophic wind (thermodynamic levels)  [m/s]
// vg: <v> geostrophic wind (thermodynamic levels)  [m/s]
// wpthvp: <w'thv'> (momentum levels)                   [m/s K]
// fcor: Traditional Coriolis parameter               [s^-1] Vertical planetary vorticity.
// Proportional to sin(latitude)
// fcor_y: Nontraditional Coriolis parameter            [s^-1] Meridional planetary vorticity.
// Proportional to cos(latitude)
// um_ref: Reference u wind component for nudging       [m/s]
// vm_ref: Reference v wind component for nudging       [m/s]
// up2: Variance of the u wind component             [m^2/s^2]
// vp2: Variance of the v wind component             [m^2/s^2]
// uprcp: < u' r_c' >                               [(m kg)/(s kg)]
// vprcp: < v' r_c' >                               [(m kg)/(s kg)]
// rc_coef_zm: Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// nu_vert_res_dep: Vertical resolution dependent nu values
// ts_nudge: Timescale of u/v nudging             [s]
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use
// for the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
// two-component PDF.
// penta_solve_method: Method to solve then penta-diagonal system
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity
// when <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
// equations of <w'w'>, <u'w'>, and <u'u'>.
// l_ho_trad_coriolis: Flag to implement the traditional Coriolis terms in the prognostic
// equations of <v'w'> and <u'w'>.
// l_diffuse_rtm_and_thlm: This flag determines whether or not we want CLUBB to do diffusion on
// rtm and thlm
// l_stability_correct_Kh_N2_zm: This flag determines whether or not we want CLUBB to apply a
// stability correction
// l_godunov_upwind_wpxp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent advection terms. It affects
// wpxp only.
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_uv_nudge: For wind speed nudging
// l_diag_Lscale_from_tau: First diagnose dissipation time tau, and then diagnose the mixing
// length scale as Lscale = tau * tke
// l_use_C7_Richardson: Parameterize C7 based on Richardson number
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm
// l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm
// l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um
// l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// rtm: r_t  (total water mixing ratio)           [kg/kg]
// wprtp: w'r_t'                                    [(kg/kg) m/s]
// thlm: th_l (liquid water potential temperature) [K]
// wpthlp: w'th_l'                                   [K m/s]
// sclrm: [Units vary]
// wpsclrp: [Units vary]
// um: <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
// upwp: <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// vm: <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
// vpwp: <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
// um_pert: perturbed <u>    [m/s]
// vm_pert: perturbed <v>    [m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_xm_wpxp(
    int nzm, int nzt, int ngrdcol, int sclr_dim, Real1D sclr_tol, const grid& gr, // In
    real_t dt,                                                                    // In
    Real2D sigma_sqd_w, Real2D wm_zm, Real2D wm_zt, Real2D wp2, Real2D Lscale_zm, // In
    Real2D wp3, Real2D Kh_zt, Real2D Kh_zm,                                       // In
    Real2D stability_correction,                                                  // In
    Real2D invrs_tau_C6_zm, Real2D tau_max_zm, Real2D wp2rtp, Real2D rtpthvp,     // In
    Real2D rtm_forcing, Real2D wprtp_forcing, Real2D rtm_ref, Real2D wp2thlp,     // In
    Real2D thlpthvp, Real2D thlm_forcing, Real2D wpthlp_forcing, Real2D thlm_ref, // In
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,                   // In
    Real2D invrs_rho_ds_zt, Real2D thv_ds_zm, Real2D rtp2, Real2D thlp2,          // In
    Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm, Real2D varnce_w_2_zm,     // In
    Real2D mixt_frac_zm, bool l_implemented, Real2D em, Real3D wp2sclrp,          // In
    Real3D sclrpthvp, Real3D sclrm_forcing, Real3D sclrp2,                        // In
    Real2D Cx_fnc_Richardson,                                                     // In
    const implicit_coefs_terms& pdf_implicit_coefs_terms,                         // In
    Real2D um_forcing, Real2D vm_forcing, Real2D ug, Real2D vg, Real2D wpthvp,    // In
    Real1D fcor, Real1D fcor_y, Real2D um_ref, Real2D vm_ref, Real2D up2,         // In
    Real2D vp2,                                                                   // In
    Real2D uprcp, Real2D vprcp, Real2D rc_coef_zm,                                // In
    Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,              // In
    real_t ts_nudge,                                                              // In
    int iiPDF_type,                                                               // In
    int penta_solve_method,                                                       // In
    int tridiag_solve_method,                                                     // In
    int fill_holes_type,                                                          // In
    bool l_predict_upwp_vpwp,                                                     // In
    bool l_ho_nontrad_coriolis,                                                   // In
    bool l_ho_trad_coriolis,                                                      // In
    bool l_diffuse_rtm_and_thlm,                                                  // In
    bool l_stability_correct_Kh_N2_zm,                                            // In
    bool l_godunov_upwind_wpxp_ta,                                                // In
    bool l_upwind_xm_ma,                                                          // In
    bool l_uv_nudge,                                                              // In
    bool l_diag_Lscale_from_tau,                                                  // In
    bool l_use_C7_Richardson,                                                     // In
    bool l_lmm_stepping,                                                          // In
    bool l_enable_relaxed_clipping,                                               // In
    bool l_linearize_pbl_winds,                                                   // In
    bool l_mono_flux_lim_thlm,                                                    // In
    bool l_mono_flux_lim_rtm,                                                     // In
    bool l_mono_flux_lim_um,                                                      // In
    bool l_mono_flux_lim_vm,                                                      // In
    bool l_mono_flux_lim_spikefix,                                                // In
    int& wprtp_cl_num, int& wpthlp_cl_num, int& upwp_cl_num, int& vpwp_cl_num,    // InOut
    f2k::stats_type& stats,                                                       // InOut
    Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp,                         // InOut
    Real3D sclrm, Real3D wpsclrp, Real2D um, Real2D upwp, Real2D vm, Real2D vpwp, // InOut
    Real2D um_pert, Real2D vm_pert, Real2D upwp_pert, Real2D vpwp_pert,           // InOut
    err_info_type& err_info                                                       // InOut
);
} // namespace clubb_kokkos
