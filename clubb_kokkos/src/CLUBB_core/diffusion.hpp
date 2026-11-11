// Fill Native band operators while retaining caller-owned output Views.
//
// diffusion_zt_lhs and diffusion_zm_lhs fill caller-owned lhs arrays in Native
// (band,column,level) order, with super/main/sub bands. Physical diffusivity,
// density and spacing fields are column first. All arguments are required. Generic
// entry overloads preserve existing caller sections through the shared array
// owner, without any band transpose.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
namespace clubb_kokkos {
// K_zm: Coef. of eddy diffusivity at momentum levels   [m^2/s]
// K_zt: Coef. of eddy diffusivity at thermo. levels    [m^2/s]
// nu: Background constant coef. of eddy diffusivity  [m^2/s]
// invrs_rho_ds_ztzxt: Inverse dry statis density on thermo. levels   [m^3/kg]
// rho_ds_zm: Dry statis density on momentum levels          [kg/m^3]
// lhs: LHS coefficient of diffusion term    [1/s]
void diffusion_zt_lhs(ScratchWorkspace& workspace, // InOut
                      int nzm, int nzt, int ngrdcol, const grid& gr, Real2D K_zm, Real2D K_zt, // In
                      Real1D nu,                                                               // In
                      Real2D invrs_rho_ds_ztzxt, Real2D rho_ds_zm,                             // In
                      Real3D lhs) // Out
    ;

// K_zt: Coef. of eddy diffusivity at thermo. levels   [m^2/s]
// K_zm: Coef. of eddy diffusivity at momentum levels [m^2/s]
// nu: Background constant coef. of eddy diffusivity [m^2/s]
// invrs_rho_ds_zm: Inverse dry, static density on momentum levels [m^3/kg]
// rho_ds_zt: Dry, static density on thermodynamic levels [kg/m^3]
// lhs: LHS coefficient of diffusion term    [1/s]
void diffusion_zm_lhs(ScratchWorkspace& workspace, // InOut
                      int nzm, int nzt, int ngrdcol, const grid& gr, Real2D K_zt, Real2D K_zm, // In
                      Real1D nu,                                                               // In
                      Real2D invrs_rho_ds_zm, Real2D rho_ds_zt,                                // In
                      Real3D lhs) // Out
    ;

// Preserve actual operator-call sections at the public array API boundary.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void diffusion_zt_lhs(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                      Array0 K_zm, Array1 K_zt, Array2 nu, Array3 invrs_rho_ds_ztzxt,
                      Array4 rho_ds_zm, Array5 lhs) {
  auto workspace_scope         = workspace.scope();
  auto K_zm_exec               = to_execution_space(workspace, K_zm);
  auto K_zt_exec               = to_execution_space(workspace, K_zt);
  auto nu_exec                 = to_execution_space(workspace, nu);
  auto invrs_rho_ds_ztzxt_exec = to_execution_space(workspace, invrs_rho_ds_ztzxt);
  auto rho_ds_zm_exec          = to_execution_space(workspace, rho_ds_zm);
  auto lhs_exec                = to_execution_space(workspace, lhs, false);
  diffusion_zt_lhs(workspace, nzm, nzt, ngrdcol, gr, K_zm_exec, K_zt_exec, nu_exec,
                   invrs_rho_ds_ztzxt_exec, rho_ds_zm_exec, lhs_exec);
  from_execution_space(workspace, lhs, lhs_exec);
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void diffusion_zm_lhs(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                      Array0 K_zt, Array1 K_zm, Array2 nu, Array3 invrs_rho_ds_zm, Array4 rho_ds_zt,
                      Array5 lhs) {
  auto workspace_scope      = workspace.scope();
  auto K_zt_exec            = to_execution_space(workspace, K_zt);
  auto K_zm_exec            = to_execution_space(workspace, K_zm);
  auto nu_exec              = to_execution_space(workspace, nu);
  auto invrs_rho_ds_zm_exec = to_execution_space(workspace, invrs_rho_ds_zm);
  auto rho_ds_zt_exec       = to_execution_space(workspace, rho_ds_zt);
  auto lhs_exec             = to_execution_space(workspace, lhs, false);
  diffusion_zm_lhs(workspace, nzm, nzt, ngrdcol, gr, K_zt_exec, K_zm_exec, nu_exec,
                   invrs_rho_ds_zm_exec, rho_ds_zt_exec, lhs_exec);
  from_execution_space(workspace, lhs, lhs_exec);
}
} // namespace clubb_kokkos
