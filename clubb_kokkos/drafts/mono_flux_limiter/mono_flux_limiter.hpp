// Bound turbulent transport using the fields that can reach each vertical level.
//
// monotonic_turbulent_flux_limit limits wpxp and re-solves/adjusts xm, preserving
// Native budgets, spike removal and fatal ordering. calc_turb_adv_range computes
// Native1based lower/upper influencing-level indices from the PDF velocities.
// Fields and integer ranges are column-first LayoutLeft execution Views. Matrix
// scratch stays (band,column,level) until the actual solver boundary. No mutable
// scientific state survives a call. See INTEGRATION.md: draft only, uncompiled.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// solve_type: Variables being solved for.
// dt: Model timestep length                           [s]
// xm_old: xm at previous time step (thermo. levs.) [units vary]
// xp2: x'^2 (momentum levels)                   [units vary]
// wm_zt: w wind component on thermodynamic levels [m/s]
// xm_forcing: xm forcings (thermodynamic levels)       [units vary]
// rho_ds_zm: Dry, static density on momentum levels   [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ moment. levs. [m^3/kg]
// invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs. [m^3/kg]
// xp2_threshold: Lower limit of x'^2                      [units vary]
// xm_tol: Lower limit of maxdev                    [units vary]
// l_implemented: True if CLUBB is being implemented and run in a host model
// low_lev_effect: Index of lowest level that has an effect (for lev. k)
// high_lev_effect: Index of highest level that has an effect (for lev. k)
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
// spurious drying tendencies at model top
// stats: Native statistics context, mutated by sampled budget calls [-]
// xm: xm at current time step (thermodynamic levels)  [units vary]
// wpxp: w'x' (momentum levels)                          [units vary]
// err_info: err_info struct containing err_code and err_header
void monotonic_turbulent_flux_limit(int nzm, int nzt, int ngrdcol, const grid &gr, int solve_type,
                                    real_t dt, // In
                                    Real2D xm_old, Real2D xp2, Real2D wm_zt,
                                    Real2D xm_forcing,                                        // In
                                    Real2D rho_ds_zm, Real2D rho_ds_zt,                       // In
                                    Real2D invrs_rho_ds_zm, Real2D invrs_rho_ds_zt,           // In
                                    real_t xp2_threshold, real_t xm_tol, bool l_implemented,  // In
                                    Kokkos::View<int **, Kokkos::LayoutLeft> low_lev_effect,  // In
                                    Kokkos::View<int **, Kokkos::LayoutLeft> high_lev_effect, // In
                                    int tridiag_solve_method,                                 // In
                                    bool l_upwind_xm_ma,                                      // In
                                    bool l_mono_flux_lim_spikefix,                            // In
                                    f2k::stats_type &stats,                          // InOut
                                    Real2D xm, Real2D wpxp, err_info_type &err_info) // InOut
    ;

// nzm: Number of momentum levels [-]
// nzt: Number of thermodynamic levels [-]
// ngrdcol: Number of columns [-]
// gr: Native grid coordinates and spacings [-,m]
// dt: Model timestep length                       [s]
// w_1_zm: Mean w (1st PDF component)                   [m/s]
// w_2_zm: Mean w (2nd PDF component)                   [m/s]
// varnce_w_1_zm: Variance of w (1st PDF component)            [m^2/s^2]
// varnce_w_2_zm: Variance of w (2nd PDF component)            [m^2/s^2]
// mixt_frac_zm: Weight of 1st PDF component (Sk_w dependent) [-]
// stats: Native statistics context, mutated by sampled budget calls [-]
// low_lev_effect: Index of lowest level that has an effect (for lev. k)
// high_lev_effect: Index of highest level that has an effect (for lev. k)
void calc_turb_adv_range(int nzm, int nzt, int ngrdcol, const grid &gr, real_t dt, // In
                         Real2D w_1_zm, Real2D w_2_zm, Real2D varnce_w_1_zm,
                         Real2D varnce_w_2_zm,                                     // In
                         Real2D mixt_frac_zm,                                      // In
                         f2k::stats_type &stats,                                   // InOut
                         Kokkos::View<int **, Kokkos::LayoutLeft> low_lev_effect,  // Out
                         Kokkos::View<int **, Kokkos::LayoutLeft> high_lev_effect) // Out
    ;
} // namespace clubb_kokkos
