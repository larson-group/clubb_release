// Prognose scalar and horizontal variances, retaining Native budgets/errors.
//
// advance_xp2_xpyp updates variances/covariances, clipping and passive scalars.
// update_xp2_mc updates their rain-evaporation tendency fields with the original
// ordered precipitation-fraction recurrence. Fields use column-first LayoutLeft;
// grid and scalar category metadata remain Native one based. Stats stays Native.
// All arguments are required and output Views remain owned by the caller.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "array_index.hpp"
#include "grid_class.hpp"
#include "parameters_tunable.hpp"
#include "pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"
namespace clubb_kokkos {
// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// invrs_tau_xp2_zm: Inverse time-scale for xp2 on momentum levels [1/s]
// invrs_tau_C4_zm: Inverse time-scale for C4 terms on mom. levels [1/s]
// invrs_tau_C14_zm: Inverse time-scale for C14 terms on mom. levels [1/s]
// wm_zm: w-wind component on momentum levels   [m/s]
// rtm: Total water mixing ratio (t-levs)     [kg/kg]
// wprtp: <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
// thlm: Liquid potential temp. (t-levs)       [K]
// wpthlp: <w'th_l'> (momentum levels)           [(m K)/s]
// wpthvp: <w'th_v'> (momentum levels)           [(m K)/s]
// um: u wind (thermodynamic levels)         [m/s]
// vm: v wind (thermodynamic levels)         [m/s]
// wp2: <w'^2> (momentum levels)              [m^2/s^2]
// wp3: <w'^3> (thermodynamic levels)         [m^3/s^3]
// upwp: <u'w'> (momentum levels)              [m^2/s^2]
// vpwp: <v'w'> (momentum levels)              [m^2/s^2]
// sigma_sqd_w: sigma_sqd_w (momentum levels)         [-]
// wprtp2: <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
// wpthlp2: <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
// wprtpthlp: <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
// Kh_zt: Eddy diffusivity on thermo. levels    [m^2/s]
// rtp2_forcing: <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
// thlp2_forcing: <th_l'^2> forcing (momentum levels)   [K^2/s]
// rtpthlp_forcing: <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
// rho_ds_zm: Dry, static density on momentum levs. [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ mom. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on mom. levs. [K]
// cloud_frac: Cloud fraction (thermodynamic levels) [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// dt: Model timestep                                [s]
// fcor_y: Nontraditional Coriolis parameter         [s^-1] Meridional planetary vorticity.
// Proportional to cos(latitude)
// sclrm: Mean value; pass. scalar (t-levs.) [{sclr units}]
// wpsclrp: <w'sclr'> (momentum levels)        [m/s{sclr units}]
// wpsclrp2: <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
// wpsclrprtp: <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
// wpsclrpthlp: <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]
// lhs_splat_wp2: LHS coefficient of wp2 splatting term  [1/s]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// nu_vert_res_dep: Vertical resolution dependent nu values
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
// PDF.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
// equations of <w'w'>, <u'w'>, and <u'u'>.
// l_min_xp2_from_corr_wx: Flag to base the threshold minimum value of xp2 (rtp2 and thlp2) on
// keeping the overall correlation of w and x within the limits of -max_mag_correlation_flux to
// max_mag_correlation_flux.
// l_C2_cloud_frac: Flag to use cloud fraction to adjust the value of the turbulent dissipation
// coefficient, C2.
// l_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// l_godunov_upwind_xpyp_ta: This flag determines whether we want to use a Godunov-like upwind
// differencing approximation rather than a centered differencing for turbulent or mean advection
// terms. It affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_implemented: True if CLUBB is being implemented and run in a host model
// rtp2: <r_t'^2>                      [(kg/kg)^2]
// thlp2: <th_l'^2>                     [K^2]
// rtpthlp: <r_t'th_l'>                   [(kg K)/kg]
// up2: <u'^2>                        [m^2/s^2]
// vp2: <v'^2>                        [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_xp2_xpyp(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, int sclr_dim,
                      Real1D sclr_tol,
                      const grid& gr,                                  // In
                      const sclr_idx_type& sclr_idx,                   // In
                      Real2D invrs_tau_xp2_zm, Real2D invrs_tau_C4_zm, // In
                      Real2D invrs_tau_C14_zm, Real2D wm_zm,           // In
                      Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp, Real2D wpthvp,
                      Real2D um,                                                         // In
                      Real2D vm,                                                         // In
                      Real2D wp2, Real2D wp3, Real2D upwp, Real2D vpwp,                  // In
                      Real2D sigma_sqd_w, Real2D wprtp2, Real2D wpthlp2,                 // In
                      Real2D wprtpthlp, Real2D Kh_zt, Real2D rtp2_forcing,               // In
                      Real2D thlp2_forcing, Real2D rtpthlp_forcing,                      // In
                      Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,        // In
                      Real2D thv_ds_zm, Real2D cloud_frac,                               // In
                      const implicit_coefs_terms& pdf_implicit_coefs_terms,              // In
                      real_t dt, Real1D fcor_y,                                          // In
                      Real3D sclrm, Real3D wpsclrp,                                      // In
                      Real3D wpsclrp2, Real3D wpsclrprtp, Real3D wpsclrpthlp,            // In
                      Real2D lhs_splat_wp2,                                              // In
                      Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,   // In
                      int iiPDF_type,                                                    // In
                      int tridiag_solve_method,                                          // In
                      int fill_holes_type,                                               // In
                      bool l_ho_nontrad_coriolis,                                        // In
                      bool l_min_xp2_from_corr_wx,                                       // In
                      bool l_C2_cloud_frac,                                              // In
                      bool l_upwind_xpyp_ta,                                             // In
                      bool l_godunov_upwind_xpyp_ta,                                     // In
                      bool l_lmm_stepping,                                               // In
                      bool l_implemented,                                                // In
                      stats_type& stats,                                                 // InOut
                      Real2D rtp2, Real2D thlp2, Real2D rtpthlp, Real2D up2, Real2D vp2, // InOut
                      Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp,
                      err_info_type& err_info) // InOut
    ;

// nzm: Momentum points in the Vertical        [-]
// nzt: Thermodynamic points in the Vertical        [-]
// dt: Model timestep        [s]
// cloud_frac: Cloud fraction                        [-]
// rcm: Cloud water mixing ratio              [kg/kg]
// rvm: Vapor water mixing ratio              [kg/kg]
// thlm: Liquid potential temperature          [K]
// wm: Mean vertical velocity                [m/s]
// exner: Exner function                        [-]
// rrm_evap: Evaporation of rain                   [kg/kg/s] It is expected that this variable is
// negative, as that is the convention in Morrison microphysics
// pdf_params: PDF parameters
// rtp2_mc: Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
// thlp2_mc: Tendency of <thl'^2> due to evaporation  [K^2/s]
// wprtp_mc: Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
// wpthlp_mc: Tendency of <w'thl'> due to evaporation  [m*K/s^2]
// rtpthlp_mc: Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]
void update_xp2_mc(ScratchWorkspace& workspace, const grid& gr, int nzm, int nzt, int ngrdcol,
                   real_t dt,
                   Real2D cloud_frac,                                                         // In
                   Real2D rcm, Real2D rvm, Real2D thlm,                                       // In
                   Real2D wm, Real2D exner, Real2D rrm_evap, const pdf_parameter& pdf_params, // In
                   Real2D rtp2_mc, Real2D thlp2_mc, Real2D wprtp_mc, Real2D wpthlp_mc, // InOut
                   Real2D rtpthlp_mc)                                                  // InOut
    ;
} // namespace clubb_kokkos
