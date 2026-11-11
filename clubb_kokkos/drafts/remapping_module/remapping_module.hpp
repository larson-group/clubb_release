// Remap vertical fields while preserving their column mass.
//
// remap_vals_to_target maps cell means between two grids using the Ullrich
// operator or piecewise parabolic method (PPM). remap_iv_from_var_name selects
// wind/positive/ordinary PPM behavior from the Native case-sensitive name list.
// Fields are (column,level), LayoutLeft; returned arrays own their storage in
// the input field's memory space. Density splines/edges ascend in altitude;
// field units are retained. All scratch is call-owned. The private cached
// operator and full diagnostics retain Native behavior. The isolated draft
// has source-audit and Serial compilation evidence; runtime/BFB remain pending.
#pragma once

#include "array_types.hpp"
#include "grid_class.hpp"
#include <string>

namespace clubb_kokkos {

// var_name: Case-sensitive field name; ADJUSTL/TRIM and Native LEN64 apply.
int remap_iv_from_var_name(
    const std::string &var_name                                                     // In
);

// ngrdcol: Number of grid columns []
// gr_source: Source grid, where source_values are given. Read-only in the complete Native
// body, despite omitted INTENT.
// gr_target: Target grid, where source_values are remapped. Read-only in the complete
// Native body, despite omitted INTENT.
// source_values_idx: dimension of source_values []
// source_values: given values on the source grid that should be remapped to the target grid
// target_values_idx: dimension of remapped source_values on target grid []
// total_idx_rho_lin_spline: number of indices for the linear spline definition arrays []
// rho_lin_spline_vals: rho values at the given altitudes [kg/m^3]
// rho_lin_spline_levels: altitudes for the given rho values [m]
// iv: setting for E3SM's PPM method, that sets what type of value is remapped, the following are
// the three valid settings: iv =-1: winds iv = 0: positive definite scalars iv = 1: others
// p_sfc: prescribed pressure at the surface [Pa]
// grid_remap_method: specifies which remapping method should be used options: 1. Ullrich linear
// remapping scheme (first order) 2. Piecewise Parabolic Method (PPM) (high order)
// l_zt_variable: if .true., the variable that was input is given on the zt levels of the grid, if
// it is .false., it is given on the zm levels
Real2D remap_vals_to_target(
    int ngrdcol,                                                                    // In
    const grid &gr_source, const grid &gr_target,                                   // In
    int source_values_idx,                                                          // In
    Real2D source_values,                                                           // In
    int target_values_idx,                                                          // In
    int total_idx_rho_lin_spline,                                                   // In
    Real2D rho_lin_spline_vals,                                                     // In
    Real2D rho_lin_spline_levels,                                                   // In
    int iv, Real1D p_sfc,                                                           // In
    int grid_remap_method,                                                          // In
    bool l_zt_variable                                                              // In
);

template <class Source, class Density, class Altitudes, class Pressure>
Kokkos::View<real_t **, Kokkos::LayoutLeft, typename Source::memory_space>
remap_vals_to_target(
    int ngrdcol,                                                                    // In
    const grid &gr_source, const grid &gr_target,                                   // In
    int source_values_idx,                                                          // In
    Source source_values,                                                           // In
    int target_values_idx,                                                          // In
    int total_idx_rho_lin_spline,                                                   // In
    Density rho_lin_spline_vals,                                                    // In
    Altitudes rho_lin_spline_levels,                                                // In
    int iv, Pressure p_sfc,                                                         // In
    int grid_remap_method,                                                          // In
    bool l_zt_variable                                                              // In
) {
  // Stage the existing public array boundary once; retain caller memory space.
  // Arguments (Native order):
  // ngrdcol: Number of grid columns []
  // gr_source: Source grid, where source_values are given. Read-only in the complete Native
  // body, despite omitted INTENT.
  // gr_target: Target grid, where source_values are remapped. Read-only in the complete
  // Native body, despite omitted INTENT.
  // source_values_idx: dimension of source_values []
  // source_values: given values on the source grid that should be remapped to the target grid
  // target_values_idx: dimension of remapped source_values on target grid []
  // total_idx_rho_lin_spline: number of indices for the linear spline definition arrays []
  // rho_lin_spline_vals: rho values at the given altitudes [kg/m^3]
  // rho_lin_spline_levels: altitudes for the given rho values [m]
  // iv: setting for E3SM's PPM method, that sets what type of value is remapped, the following are
  // the three valid settings: iv =-1: winds iv = 0: positive definite scalars iv = 1: others
  // p_sfc: prescribed pressure at the surface [Pa]
  // grid_remap_method: specifies which remapping method should be used options: 1. Ullrich linear
  // remapping scheme (first order) 2. Piecewise Parabolic Method (PPM) (high order)
  // l_zt_variable: if .true., the variable that was input is given on the zt levels of the grid,
  // if it is .false., it is given on the zm levels

  const auto remapped = remap_vals_to_target(
      ngrdcol, // In
      gr_source, gr_target, // In
      source_values_idx, // In
      to_execution_space(source_values), // In
      target_values_idx, // In
      total_idx_rho_lin_spline, // In
      to_execution_space(rho_lin_spline_vals), // In
      to_execution_space(rho_lin_spline_levels), // In
      iv, to_execution_space(p_sfc), // In
      grid_remap_method, // In
      l_zt_variable // In
  );
  using Result = Kokkos::View<real_t **, Kokkos::LayoutLeft, typename Source::memory_space>;
  Result result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "remap_vals_to_target"), ngrdcol,
                target_values_idx);
  from_execution_space(result, remapped);
  return result;
}

} // namespace clubb_kokkos
