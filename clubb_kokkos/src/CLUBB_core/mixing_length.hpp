// Mixing distances and turbulence dissipation times for CLUBB.
//
// calc_Lscale coordinates one timestep: select a mixing-length method, compute
// lengths and inverse damping times, then sample registered diagnostics.
// set_Lscale_max returns a cap from host horizontal spacing or 100000 m in
// standalone mode. calc_Lscale_directly selects the legacy or coherent parcel
// method, and diagnose_Lscale_from_tau diagnoses lengths from damping times.
//
// Fields are column-first LayoutLeft execution-space Views. Thermodynamic and
// momentum fields have nzt and nzm levels; lengths are in metres, time scales
// in seconds, and inverse times in 1/s. Grid metadata retains Native one-based
// values. Statistics and errors update the caller's supplied state, and all
// parcel/scratch storage belongs to the current call.
//
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
#include "pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"

namespace clubb_kokkos {

// nzm: Number of momentum levels
// nzt: Number of thermodynamic levels
// ngrdcol: Number of grid columns
// gr: Grid structure
// l_implemented: True when current configuration is implemented
// host_dx: Host-model grid spacing in x direction [m]
// host_dy: Host-model grid spacing in y direction [m]
// p_in_Pa: Pressure [Pa]
// exner: Exner function [-]
// rtm: Total water mixing ratio [kg/kg]
// thlm: Liquid water potential temperature [K]
// thvm: Virtual potential temperature [K]
// thlp2: Variance of thl [K^2]
// rtp2: Variance of rt [(kg/kg)^2]
// rtpthlp: Covariance of rt and thl [K kg/kg]
// pdf_params: PDF parameters
// em: Turbulent kinetic energy [m^2/s^2]
// thv_ds_zt: Reference-state virtual potential temperature [K]
// lmin: Minimum mixing length [m]
// upwp_sfc: Surface value of u'w' [m^2/s^2]
// vpwp_sfc: Surface value of v'w' [m^2/s^2]
// ddzt_umvm_sqd: Squared vertical shear of horizontal wind [s^-2]
// ice_supersat_frac: Ice supersaturation fraction [-]
// ufmin: Minimum velocity scale [m/s]
// tau_const: Constant eddy dissipation time scale [s]
// sfc_elevation: Surface elevation [m]
// clubb_params: Tunable CLUBB parameters
// saturation_formula: Saturation vapor pressure formula selector
// l_Lscale_plume_centered: Whether to use plume-centered Lscale calculation
// Lscale_method: Mixing-length formulation
// l_e3sm_config: Whether E3SM configuration options are active
// l_smooth_Heaviside_tau_wpxp: Whether to smooth Heaviside functions in wpxp tau
// l_modify_limiters_for_cnvg_test: Whether to modify limiters for convergence tests
// l_use_invrs_tau_N2_iso: Whether to use isotropic N2 inverse tau
// brunt_vaisala_freq_sqd_smth: Smoothed Brunt-Vaisala frequency squared [s^-2]
// stats: Statistics accumulator
// err_info: Error information
// invrs_tau_zt: Inverse time scale on thermodynamic levels [1/s]
// invrs_tau_zm: Inverse time scale on momentum levels [1/s]
// invrs_tau_xp2_zm: Inverse time scale for xp2 terms on momentum levels [1/s]
// invrs_tau_wp3_zt: Inverse time scale for wp3 terms on thermodynamic levels [1/s]
// invrs_tau_C1_zm: Inverse time scale used for C1 term on momentum levels [1/s]
// invrs_tau_C4_zm: Inverse time scale used for C4 terms on momentum levels [1/s]
// invrs_tau_C6_zm: Inverse time scale used for C6 term on momentum levels [1/s]
// invrs_tau_C14_zm: Inverse time scale used for C14 terms on momentum levels [1/s]
// tau_max_zm: Maximum allowable eddy dissipation time scale on momentum levels [s]
// tau_max_zt: Maximum allowable eddy dissipation time scale on thermodynamic levels [s]
// tau_zm: Eddy dissipation time scale on momentum levels [s]
// Lscale: Turbulent mixing length on thermodynamic levels [m]
// Lscale_zm: Turbulent mixing length on momentum levels [m]
// Lscale_up: Upward component of turbulent mixing length [m]
// Lscale_down: Downward component of turbulent mixing length [m]
void calc_Lscale(ScratchWorkspace& workspace,                                               // InOut
                 int nzm, int nzt, int ngrdcol, const grid& gr, bool l_implemented,         // In
                 Real1D host_dx, Real1D host_dy,                                            // In
                 Real2D p_in_Pa, Real2D exner, Real2D rtm, Real2D thlm, Real2D thvm,        // In
                 Real2D thlp2, Real2D rtp2, Real2D rtpthlp,                                 // In
                 const pdf_parameter& pdf_params, Real2D em, Real2D thv_ds_zt, Real1D lmin, // In
                 Real1D upwp_sfc, Real1D vpwp_sfc, Real2D ddzt_umvm_sqd,                    // In
                 Real2D ice_supersat_frac,                                                  // In
                 real_t ufmin, real_t tau_const, Real1D sfc_elevation, Real2D clubb_params, // In
                 int saturation_formula, bool l_Lscale_plume_centered,                      // In
                 int Lscale_method,                                                         // In
                 bool l_e3sm_config,                                                        // In
                 bool l_smooth_Heaviside_tau_wpxp,                                          // In
                 bool l_modify_limiters_for_cnvg_test,                                      // In
                 bool l_use_invrs_tau_N2_iso,                                               // In
                 Real2D brunt_vaisala_freq_sqd_smth,                                        // In
                 stats_type& stats, err_info_type& err_info,                                // InOut
                 Real2D invrs_tau_zt, Real2D invrs_tau_zm, Real2D invrs_tau_xp2_zm,         // Out
                 Real2D invrs_tau_wp3_zt,                                                   // Out
                 Real2D invrs_tau_C1_zm, Real2D invrs_tau_C4_zm,                            // Out
                 Real2D invrs_tau_C6_zm, Real2D invrs_tau_C14_zm,                           // Out
                 Real2D tau_max_zm, Real2D tau_max_zt, Real2D tau_zm,                       // Out
                 Real2D Lscale, Real2D Lscale_zm, Real2D Lscale_up, Real2D Lscale_down      // Out
);

// ngrdcol: Number of grid columns
// l_implemented: True if CLUBB is being implemented and run in a host model
// host_dx: Host model's east-west horizontal grid spacing     [m]
// host_dy: Host model's north-south horizontal grid spacing   [m]
// Lscale_max: Maximum allowable value for Lscale   [m]
void set_Lscale_max(int ngrdcol, bool l_implemented, Real1D host_dx, Real1D host_dy, // In
                    Real1D Lscale_max                                                // Out
);

// ngrdcol: Number of grid columns
// nzm: Number of momentum levels
// nzt: Number of thermodynamic levels
// gr: Grid structure
// l_implemented: True if CLUBB is being implemented and run in a host model
// p_in_Pa: Air pressure (thermodynamic levels)       [Pa]
// exner: Exner function [-]
// rtm: Mean total-water mixing ratio [kg/kg]
// thlm: Mean liquid-water potential temperature [K]
// thvm: Mean virtual potential temperature [K]
// rtp2_zt: Total-water variance [(kg/kg)^2]
// thlp2_zt: Liquid-water potential-temperature variance [K^2]
// rtpthlp_zt: Total-water/temperature covariance [K kg/kg]
// pdf_params: PDF Parameters  [units vary]
// em: Turbulent kinetic energy [m^2/s^2]
// thv_ds_zt: Reference-state virtual potential temperature [K]
// Lscale_max: Per-column cap on combined mixing length [m]
// lmin: Surface value of the height-dependent length floor [m]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// saturation_formula: Integer that stores the saturation formula to be used
// l_Lscale_plume_centered: Alternate that uses the PDF to compute the perturbed values
// Lscale_method: Mixing-length formulation
// stats: Statistics accumulator
// err_info: err_info struct containing err_code and err_header
// Lscale: Mixing length      [m]
// Lscale_up: Mixing length up   [m]
// Lscale_down: Mixing length down [m]
void calc_Lscale_directly(ScratchWorkspace& workspace,                                  // InOut
                          int ngrdcol, int nzm, int nzt, const grid& gr,                // In
                          bool l_implemented, Real2D p_in_Pa, Real2D exner, Real2D rtm, // In
                          Real2D thlm, Real2D thvm, Real2D rtp2_zt, Real2D thlp2_zt,
                          Real2D rtpthlp_zt,                                            // In
                          const pdf_parameter& pdf_params, Real2D em, Real2D thv_ds_zt, // In
                          Real1D Lscale_max, Real1D lmin,                               // In
                          Real2D clubb_params,                                          // In
                          int saturation_formula,                                       // In
                          bool l_Lscale_plume_centered,                                 // In
                          int Lscale_method,                                            // In
                          stats_type& stats, err_info_type& err_info,                   // InOut
                          Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down           // Out
);

// nzm: Number of momentum levels
// nzt: Number of thermodynamic levels
// ngrdcol: Number of grid columns
// gr: Grid structure
// upwp_sfc: Surface value of u'w' [m^2/s^2]
// vpwp_sfc: Surface value of v'w' [m^2/s^2]
// ddzt_umvm_sqd: Squared vertical shear of horizontal wind [s^-2]
// ice_supersat_frac: Ice supersaturation fraction [-]
// em: Turbulent kinetic energy [m^2/s^2]
// ufmin: Minimum velocity scale [m/s]
// tau_const: Constant eddy dissipation time scale [s]
// sfc_elevation: Surface elevation [m]
// Lscale_max: Maximum allowable mixing length [m]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// stats: Statistics accumulator
// l_e3sm_config: Whether E3SM configuration options are active
// l_smooth_Heaviside_tau_wpxp: Use the smoothed Heaviside 'Peskin' function to compute
// invrs_tau_wpxp_zm
// brunt_vaisala_freq_sqd_smth: smoothed Buoyancy frequency squared, N^2     [s^-2]
// Ri_zm: Richardson number
// err_info: err_info struct containing err_code and err_header
// invrs_tau_zt: Inverse time scale on thermodynamic levels [1/s]
// invrs_tau_zm: Inverse time scale on momentum levels [1/s]
// invrs_tau_sfc: Surface contribution to inverse tau [1/s]
// invrs_tau_no_N2_zm: Inverse tau without N2 contribution on momentum levels [1/s]
// invrs_tau_bkgnd: Background inverse tau [1/s]
// invrs_tau_shear: Shear contribution to inverse tau [1/s]
// invrs_tau_N2_iso: Isotropic N2 contribution to inverse tau [1/s]
// invrs_tau_wp2_zm: Inverse time scale for wp2 terms on momentum levels [1/s]
// invrs_tau_xp2_zm: Inverse time scale for xp2 terms on momentum levels [1/s]
// invrs_tau_wp3_zm: Inverse time scale for wp3 terms on momentum levels [1/s]
// invrs_tau_wp3_zt: Inverse time scale for wp3 terms on thermodynamic levels [1/s]
// invrs_tau_wpxp_zm: Inverse time scale for wpxp terms on momentum levels [1/s]
// tau_max_zm: Maximum allowable eddy dissipation time scale on momentum levels [s]
// tau_max_zt: Maximum allowable eddy dissipation time scale on thermodynamic levels [s]
// tau_zm: Eddy dissipation time scale on momentum levels [s]
// tau_zt: Eddy dissipation time scale on thermodynamic levels [s]
// Lscale: Turbulent mixing length on thermodynamic levels [m]
// Lscale_up: Upward component of turbulent mixing length [m]
// Lscale_down: Downward component of turbulent mixing length [m]
void diagnose_Lscale_from_tau(ScratchWorkspace& workspace,                            // InOut
                              int nzm, int nzt, int ngrdcol, const grid& gr,          // In
                              Real1D upwp_sfc, Real1D vpwp_sfc, Real2D ddzt_umvm_sqd, // In
                              Real2D ice_supersat_frac,                               // In
                              Real2D em,                                              // In
                              real_t ufmin, real_t tau_const,                         // In
                              Real1D sfc_elevation, Real1D Lscale_max,                // In
                              Real2D clubb_params,                                    // In
                              stats_type& stats,                                      // InOut
                              bool l_e3sm_config,                                     // In
                              bool l_smooth_Heaviside_tau_wpxp,                       // In
                              Real2D brunt_vaisala_freq_sqd_smth, Real2D Ri_zm,       // In
                              err_info_type& err_info,                                // InOut
                              Real2D invrs_tau_zt, Real2D invrs_tau_zm,               // Out
                              Real2D invrs_tau_sfc, Real2D invrs_tau_no_N2_zm,
                              Real2D invrs_tau_bkgnd,                           // Out
                              Real2D invrs_tau_shear, Real2D invrs_tau_N2_iso,  // Out
                              Real2D invrs_tau_wp2_zm, Real2D invrs_tau_xp2_zm, // Out
                              Real2D invrs_tau_wp3_zm, Real2D invrs_tau_wp3_zt,
                              Real2D invrs_tau_wpxp_zm, // Out
                              Real2D tau_max_zm, Real2D tau_max_zt, Real2D tau_zm,
                              Real2D tau_zt,                                      // Out
                              Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down // Out
);

} // namespace clubb_kokkos
