// Compute PDF-based turbulent advection for variances, covariances and fluxes.
//
// xpyp_term_ta_pdf_lhs and its Godunov variant fill three-band implicit
// operators; xpyp_term_ta_pdf_rhs and its Godunov variant fill explicit
// tendencies. Physical fields use column-first LayoutLeft Views. LHS output
// retains Native (band,column,momentum level) order, with super/main/sub bands.
// Both storage endpoints are zero; grid direction selects upwind differences.
// Generic entry overloads retain caller-owned Views through the existing array
// owner. All arguments are required; see INTEGRATION.md for pending validation.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"

namespace clubb_kokkos {

// coef_wpxpyp_implicit: Coef. of <x'y'> in <w'x'y'>; t-levs  [m/s]
// rho_ds_zt: Dry, static density at t-levels      [kg/m^3]
// rho_ds_zm: Dry, static density at m-levs       [kg/m^3]
// invrs_rho_ds_zm: Inv dry, static density at m-levels  [m^3/kg]
// l_upwind_xpyp_turbulent_adv: Flag to use "upwind" discretization
// sgn_turbulent_vel: Sign of the turbulent velocity      [-]
// coef_wpxpyp_implicit_zm: coef_wpxpyp_implicit interp m-levs  [m/s]
// lhs_ta: LHS coefficient of xpyp turbulent advection  [1/s]
// lhs_ta storage: (band,column,momentum level), super/main/sub bands [1/s].
void xpyp_term_ta_pdf_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Real2D coef_wpxpyp_implicit,        // In
                          Real2D rho_ds_zt, Real2D rho_ds_zm, // In
                          Real2D invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Real2D sgn_turbulent_vel,           // In
                          Real2D coef_wpxpyp_implicit_zm,     // In
                          Real3D lhs_ta);                     // Out

// Preserve Native array sections at the existing shared ownership boundary.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void xpyp_term_ta_pdf_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Array0 coef_wpxpyp_implicit,        // In
                          Array1 rho_ds_zt, Array2 rho_ds_zm, // In
                          Array3 invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Array4 sgn_turbulent_vel,           // In
                          Array5 coef_wpxpyp_implicit_zm,     // In
                          Array6 lhs_ta)                      // Out
{
  // Retain caller allocations across the shared array boundary.
  auto coef_wpxpyp_implicit_exec = to_execution_space(coef_wpxpyp_implicit);
  auto rho_ds_zt_exec = to_execution_space(rho_ds_zt);
  auto rho_ds_zm_exec = to_execution_space(rho_ds_zm);
  auto invrs_rho_ds_zm_exec = to_execution_space(invrs_rho_ds_zm);
  auto sgn_turbulent_vel_exec = to_execution_space(sgn_turbulent_vel);
  auto coef_wpxpyp_implicit_zm_exec = to_execution_space(coef_wpxpyp_implicit_zm);
  auto lhs_ta_exec = to_execution_space(lhs_ta, false);
  xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wpxpyp_implicit_exec, rho_ds_zt_exec,
                       rho_ds_zm_exec, invrs_rho_ds_zm_exec, l_upwind_xpyp_turbulent_adv,
                       sgn_turbulent_vel_exec, coef_wpxpyp_implicit_zm_exec, lhs_ta_exec);
  from_execution_space(lhs_ta, lhs_ta_exec);
}

// coef_wpxpyp_implicit: Coef. of <x'y'> in <w'x'y'>; t-lev [m/s]
// invrs_rho_ds_zm: Inv dry, static density @ m-level [m^3/kg]
// rho_ds_zm: Dry, static density at m-lev [kg/m^3]
// lhs_ta storage: (band,column,momentum level), super/main/sub bands [1/s].
void xpyp_term_ta_pdf_lhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Real2D coef_wpxpyp_implicit,                   // In
                                  Real2D invrs_rho_ds_zm, Real2D rho_ds_zm,      // In
                                  Real3D lhs_ta);                                // Out

// Preserve Native array sections at the existing shared ownership boundary.
template <class Array0, class Array1, class Array2, class Array3>
void xpyp_term_ta_pdf_lhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Array0 coef_wpxpyp_implicit,                   // In
                                  Array1 invrs_rho_ds_zm, Array2 rho_ds_zm,      // In
                                  Array3 lhs_ta)                                 // Out
{
  // Retain caller allocations across the shared array boundary.
  auto coef_wpxpyp_implicit_exec = to_execution_space(coef_wpxpyp_implicit);
  auto invrs_rho_ds_zm_exec = to_execution_space(invrs_rho_ds_zm);
  auto rho_ds_zm_exec = to_execution_space(rho_ds_zm);
  auto lhs_ta_exec = to_execution_space(lhs_ta, false);
  xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ngrdcol, gr, coef_wpxpyp_implicit_exec,
                               invrs_rho_ds_zm_exec, rho_ds_zm_exec, lhs_ta_exec);
  from_execution_space(lhs_ta, lhs_ta_exec);
}

// term_wpxpyp_explicit: RHS: <w'x'y'> eq; t-levs      [m/s(x un)(y un)]
// rho_ds_zt: Dry, static density at t-levs          [kg/m^3]
// rho_ds_zm: Dry, static density at m-levs       [kg/m^3]
// invrs_rho_ds_zm: Inv dry, static density at m-levs   [m^3/kg]
// l_upwind_xpyp_turbulent_adv: Flag to use "upwind" discretization
// sgn_turbulent_vel: Sign of the turbulent velocity           [-]
// term_wpxpyp_explicit_zm: term_wpxpyp_expl. zm       [m/s(x un)(y un)]
// rhs_ta: RHS portion of xpyp turbulent advection      [(x un)(y un)/s]
void xpyp_term_ta_pdf_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Real2D term_wpxpyp_explicit,        // In
                          Real2D rho_ds_zt, Real2D rho_ds_zm, // In
                          Real2D invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Real2D sgn_turbulent_vel,           // In
                          Real2D term_wpxpyp_explicit_zm,     // In
                          Real2D rhs_ta);                     // Out

// Preserve Native array sections at the existing shared ownership boundary.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void xpyp_term_ta_pdf_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                          Array0 term_wpxpyp_explicit,        // In
                          Array1 rho_ds_zt, Array2 rho_ds_zm, // In
                          Array3 invrs_rho_ds_zm,             // In
                          bool l_upwind_xpyp_turbulent_adv,   // In
                          Array4 sgn_turbulent_vel,           // In
                          Array5 term_wpxpyp_explicit_zm,     // In
                          Array6 rhs_ta)                      // Out
{
  // Retain caller allocations across the shared array boundary.
  auto term_wpxpyp_explicit_exec = to_execution_space(term_wpxpyp_explicit);
  auto rho_ds_zt_exec = to_execution_space(rho_ds_zt);
  auto rho_ds_zm_exec = to_execution_space(rho_ds_zm);
  auto invrs_rho_ds_zm_exec = to_execution_space(invrs_rho_ds_zm);
  auto sgn_turbulent_vel_exec = to_execution_space(sgn_turbulent_vel);
  auto term_wpxpyp_explicit_zm_exec = to_execution_space(term_wpxpyp_explicit_zm);
  auto rhs_ta_exec = to_execution_space(rhs_ta, false);
  xpyp_term_ta_pdf_rhs(nzm, nzt, ngrdcol, gr, term_wpxpyp_explicit_exec, rho_ds_zt_exec,
                       rho_ds_zm_exec, invrs_rho_ds_zm_exec, l_upwind_xpyp_turbulent_adv,
                       sgn_turbulent_vel_exec, term_wpxpyp_explicit_zm_exec, rhs_ta_exec);
  from_execution_space(rhs_ta, rhs_ta_exec);
}

// term_wpxpyp_explicit_zm: RHS: <w'x'y'> eq; m-lev(k)   [m/s(x un)(y un)]
// invrs_rho_ds_zm: Inv dry, static density at m-lev (k) [m^3/kg]
// sgn_turbulent_vel: Sign of turbulent velocity on thermodynamic levels here [-].
// rho_ds_zm: Dry, static density at m-lev (k) [kg/m^3]
// rhs_ta: Explicit turbulent-advection tendency [(x units)(y units)/s].
void xpyp_term_ta_pdf_rhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Real2D term_wpxpyp_explicit_zm,                // In
                                  Real2D invrs_rho_ds_zm,                        // In
                                  Real2D sgn_turbulent_vel,                      // In
                                  Real2D rho_ds_zm,                              // In
                                  Real2D rhs_ta);                                // Out

// Preserve Native array sections at the existing shared ownership boundary.
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void xpyp_term_ta_pdf_rhs_godunov(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                  Array0 term_wpxpyp_explicit_zm,                // In
                                  Array1 invrs_rho_ds_zm,                        // In
                                  Array2 sgn_turbulent_vel,                      // In
                                  Array3 rho_ds_zm,                              // In
                                  Array4 rhs_ta)                                 // Out
{
  // Retain caller allocations across the shared array boundary.
  auto term_wpxpyp_explicit_zm_exec = to_execution_space(term_wpxpyp_explicit_zm);
  auto invrs_rho_ds_zm_exec = to_execution_space(invrs_rho_ds_zm);
  auto sgn_turbulent_vel_exec = to_execution_space(sgn_turbulent_vel);
  auto rho_ds_zm_exec = to_execution_space(rho_ds_zm);
  auto rhs_ta_exec = to_execution_space(rhs_ta, false);
  xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, gr, term_wpxpyp_explicit_zm_exec,
                               invrs_rho_ds_zm_exec, sgn_turbulent_vel_exec, rho_ds_zm_exec,
                               rhs_ta_exec);
  from_execution_space(rhs_ta, rhs_ta_exec);
}

} // namespace clubb_kokkos
