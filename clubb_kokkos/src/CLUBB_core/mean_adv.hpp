// Fill Native band operators while retaining caller-owned output Views.
//
// term_ma_zt_lhs and term_ma_zm_lhs fill caller-owned lhs_ma arrays in Native
// (band,column,level) order, with super/main/sub bands. Velocities, spacings and
// interpolation weights are column first. All arguments are required. Generic
// entry overloads preserve existing caller sections through the shared array
// owner, without any band transpose.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
namespace clubb_kokkos {
// wm_zt: wm_zt                        [m/s]
// weights_zt2zm: Thermodynamic-to-momentum interpolation weights (column,momentum
// level,upper/lower category) [-]; Native categories1:2
// invrs_dzt: Inverse of grid spacing      [1/m]
// invrs_dzm: Inverse of grid spacing      [1/m]
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// grid_dir: Grid direction multiplier (+1 = ascending; -1 = descending)
// lhs_ma: Mean advection contributions to lhs    [1/s]
void term_ma_zt_lhs(int nzm, int nzt, int ngrdcol, Real2D wm_zt, Real3D weights_zt2zm, // In
                    Real2D invrs_dzt, Real2D invrs_dzm,                                // In
                    bool l_upwind_xm_ma, real_t grid_dir,                              // In
                    Real3D lhs_ma)                                                     // Out
    ;

// wm_zm: wm_zm                        [m/s]
// invrs_dzm: Inverse of grid spacing      [1/m]
// weights_zm2zt: Momentum-to-thermodynamic interpolation weights (column,thermodynamic
// level,upper/lower category) [-]; Native categories1:2
// lhs_ma: Mean advection contributions to lhs  [1/s]
void term_ma_zm_lhs(int nzm, int nzt, int ngrdcol, Real2D wm_zm, // In
                    Real2D invrs_dzm, Real3D weights_zm2zt,      // In
                    Real3D lhs_ma)                               // Out
    ;

// Preserve actual operator-call sections at the public array API boundary.
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void term_ma_zt_lhs(int nzm, int nzt, int ngrdcol, Array0 wm_zt, Array1 weights_zt2zm,
                    Array2 invrs_dzt, Array3 invrs_dzm, bool l_upwind_xm_ma, real_t grid_dir,
                    Array4 lhs_ma) {
  auto wm_zt_exec         = to_execution_space(wm_zt);
  auto weights_zt2zm_exec = to_execution_space(weights_zt2zm);
  auto invrs_dzt_exec     = to_execution_space(invrs_dzt);
  auto invrs_dzm_exec     = to_execution_space(invrs_dzm);
  auto lhs_ma_exec        = to_execution_space(lhs_ma, false);
  term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt_exec, weights_zt2zm_exec, invrs_dzt_exec, invrs_dzm_exec,
                 l_upwind_xm_ma, grid_dir, lhs_ma_exec);
  from_execution_space(lhs_ma, lhs_ma_exec);
}
template <class Array0, class Array1, class Array2, class Array3>
void term_ma_zm_lhs(int nzm, int nzt, int ngrdcol, Array0 wm_zm, Array1 invrs_dzm,
                    Array2 weights_zm2zt, Array3 lhs_ma) {
  auto wm_zm_exec         = to_execution_space(wm_zm);
  auto invrs_dzm_exec     = to_execution_space(invrs_dzm);
  auto weights_zm2zt_exec = to_execution_space(weights_zm2zt);
  auto lhs_ma_exec        = to_execution_space(lhs_ma, false);
  term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm_exec, invrs_dzm_exec, weights_zm2zt_exec, lhs_ma_exec);
  from_execution_space(lhs_ma, lhs_ma_exec);
}
} // namespace clubb_kokkos
