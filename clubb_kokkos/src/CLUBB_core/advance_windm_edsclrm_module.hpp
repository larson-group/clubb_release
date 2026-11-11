// Advance mean horizontal winds and eddy scalars using Native diagnostics.
//
// advance_windm_edsclrm applies diffusion, momentum surface fluxes and nudging,
// with optional perturbed winds and explicit thlm/rtm copies. All 52 Native
// arguments are required; output Views and scalar counters belong to callers.
// Physical fields use column-first LayoutLeft; grid/category metadata retains
// Native one-based values. Statistics remains the Native backend. Standalone
// initialization rejects UV sponge damping before allocating model fields.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "array_index.hpp"
#include "grid_class.hpp"
#include "parameters_tunable.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"
namespace clubb_kokkos {
// sclr_idx: Indices identifying the eddy scalars for thlm and rtm [-]
// dt: Model timestep                             [s]
// wm_zt: w wind component on thermodynamic levels    [m/s]
// Kh_zm: Eddy diffusivity for scalars on momentum levs. [m^2/s]
// clubb_params: Array of CLUBB's tunable parameters            [units vary]
// ug: u (west-to-east) geostrophic wind comp.     [m/s]
// vg: v (south-to-north) geostrophic wind comp.   [m/s]
// um_ref: Reference u wind component for nudging      [m/s]
// vm_ref: Reference v wind component for nudging      [m/s]
// wp2: w'^2 (momentum levels)                      [m^2/s^2]
// up2: u'^2 (momentum levels)                      [m^2/s^2]
// vp2: v'^2 (momentum levels)                      [m^2/s^2]
// um_forcing: u forcing                                   [m/s/s]
// vm_forcing: v forcing                                   [m/s/s]
// edsclrm_forcing: Eddy scalar large-scale forcing        [{units vary}/s]
// p_in_Pa: Air pressure (thermodynamic levels)          [Pa]
// rho_ds_zm: Dry, static density on momentum levels      [kg/m^3]
// rho_ds_zt: Dry, static density on thermodynamic levels  [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels  [m^3/kg]
// fcor: Coriolis parameter                            [s^-1]
// l_implemented: True if CLUBB is being implemented and run in a host model
// nu_vert_res_dep: Vertical resolution dependent nu values
// ts_nudge: Timescale of u/v nudging             [s]
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity when
// <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_uv_nudge: For wind speed nudging
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_do_expldiff_rtm_thlm: Explicitly diffuse rtm and thlm with eddy scalars
// um: Mean u (west-to-east) wind component         [m/s]
// vm: Mean v (south-to-north) wind component       [m/s]
// thlm: Liquid potential temperature                 [K]
// rtm: Total water mixing ratio                     [kg/kg]
// edsclrm: Mean eddy scalar quantity             [units vary]
// upwp: <u'w'> (momentum levels)                     [m^2/s^2]
// vpwp: <v'w'> (momentum levels)                     [m^2/s^2]
// wpedsclrp: w'edsclr' (momentum levels)           [m/s {units vary}]
// um_pert: perturbed <u>    [m/s]
// vm_pert: perturbed <v>    [m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_windm_edsclrm(ScratchWorkspace& workspace,                                   // InOut
                           int nzm, int nzt, int ngrdcol, int edsclr_dim, const grid& gr, // In
                           const sclr_idx_type& sclr_idx, real_t dt,                      // In
                           Real2D wm_zt, Real2D Kh_zm, Real2D clubb_params,               // In
                           Real2D ug, Real2D vg, Real2D um_ref, Real2D vm_ref,            // In
                           Real2D wp2, Real2D up2, Real2D vp2, Real2D um_forcing,
                           Real2D vm_forcing,                                             // In
                           Real3D edsclrm_forcing, Real2D p_in_Pa,                        // In
                           Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt,    // In
                           Real1D fcor, bool l_implemented,                               // In
                           const nu_vertical_res_dep& nu_vert_res_dep, real_t ts_nudge,   // In
                           int tridiag_solve_method,                                      // In
                           bool l_predict_upwp_vpwp,                                      // In
                           bool l_upwind_xm_ma,                                           // In
                           bool l_uv_nudge,                                               // In
                           bool l_lmm_stepping,                                           // In
                           bool l_linearize_pbl_winds,                                    // In
                           bool l_do_expldiff_rtm_thlm,                                   // In
                           int fill_holes_type,                                           // In
                           int& upwp_cl_num, int& vpwp_cl_num,                            // InOut
                           stats_type& stats,                                             // InOut
                           Real2D um, Real2D vm, Real2D thlm, Real2D rtm, Real3D edsclrm, // InOut
                           Real2D upwp, Real2D vpwp, Real3D wpedsclrp,                    // InOut
                           Real2D um_pert, Real2D vm_pert, Real2D upwp_pert,
                           Real2D vpwp_pert,        // InOut
                           err_info_type& err_info) // InOut
    ;
} // namespace clubb_kokkos
