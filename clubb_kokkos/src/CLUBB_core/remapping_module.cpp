// Conservative vertical remapping; see remapping_module.hpp for public APIs.
//
// Native profile integration, pressure recurrences, limiter stages and PPM cell
// sums remain ordered. Independent cells/columns use Kokkos kernels; host code
// retains diagnostic and early-return order. Private PPM coefficients use
// (column,level,coefficient), adapting Native (coefficient,column,level).
// The public PPM wrapper's sequential profile calls retain its reused-output
// failure behavior. No cached/global storage or JAX cumulative sums are used.
// References: Ullrich/Taylor remapping Parts I/II, formulae9/12/14/30;
// E3SM mapz_module.F90, the PPM source cited in the canonical Fortran file.
#include "error_code.hpp"
#include "remapping_module.hpp"
#include "constants_clubb.hpp"
#include "interpolation.hpp"
#include "model_flags.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {
namespace {
constexpr real_t tol = real_t(1.0e-6); // Tolerance for real-number comparisons.

template <class SourceLevels, class TargetLevels, class Output>
void remapping_matrix(ScratchWorkspace& workspace,                            // InOut
                      int levels_source_idx, int levels_target_idx,           // In
                      SourceLevels levels_source, TargetLevels levels_target, // In
                      Output remapping_matrix                                 // Out
);
void remap_vals_ppm(ScratchWorkspace& workspace,                  // InOut
                    int ngrdcol,                                  // In
                    int levels_source_idx, int levels_target_idx, // In
                    Real2D levels_source, Real2D levels_target,   // In
                    Real2D source_values, int iv,                 // In
                    Real2D remap_vals_ppm                         // Out
);
} // namespace

//-----------------------------------------------------------------------------
// GENERAL REMAPPING FUNCTIONS AND SUBROUTINES
//-----------------------------------------------------------------------------
int remap_iv_from_var_name(
    const std::string &var_name                                                     // In
) {
  // Return the Native PPM selector: -1winds,0positive-definite,1other fields.
  // ADJUSTL/TRIM remove blanks only; assignment to LEN64 then comparisons
  // ignore trailing blanks, including any introduced by truncating a long name.
  // Arguments (Native order):
  // var_name: Case-sensitive field name; ADJUSTL/TRIM and Native LEN64 apply.

  std::string name = var_name;
  name.erase(0, name.find_first_not_of(' '));
  name.erase(name.find_last_not_of(' ') + 1);
  name.resize(std::min(name.size(), std::size_t(64)));
  name.erase(name.find_last_not_of(' ') + 1);
  int iv = 1;
  if (name == "rho" || name == "rho_ds_zt" || name == "invrs_rho_ds_zt" || name == "thv_ds_zt" ||
      name == "rfrzm" || name == "rtm_ref" || name == "thlm_ref" || name == "rtm" ||
      name == "thlm" || name == "thvm" || name == "rcm" || name == "cloud_frac" ||
      name == "ice_supersat_frac" || name == "rcm_in_layer" || name == "cloud_cover" ||
      name == "cloudy_updraft_frac" || name == "cloudy_downdraft_frac" || name == "Kh_zt" ||
      name == "Lscale" || name == "rho_zm" || name == "rho_ds_zm" || name == "invrs_rho_ds_zm" ||
      name == "thv_ds_zm" || name == "up2" || name == "vp2" || name == "rtp2" || name == "thlp2" ||
      name == "wp2" || name == "rc_coef_zm" || name == "wp4" || name == "wp2up2" ||
      name == "wp2vp2" || name == "invrs_tau_zm" || name == "Kh_zm")
    iv = 0;
  if (name == "wm_zt" || name == "um_ref" || name == "vm_ref" || name == "ug" || name == "vg" ||
      name == "um" || name == "vm" || name == "um_pert" || name == "vm_pert" || name == "wm_zm")
    iv = -1;
  return iv;
}

namespace {
//=============================================================================
template <class Density, class Altitudes, class GridLevels, class Output>
void calc_mass_over_grid_intervals(ScratchWorkspace& workspace,                 // InOut
                                   int total_idx_lin_spline,                    // In
                                   Density lin_spline_rho_vals,                 // In
                                   Altitudes lin_spline_rho_levels,             // In
                                   int grid_levels_idx, GridLevels grid_levels, // In
                                   Output mass_per_interval                     // Out
) {
  // Density[kg/m^3] and ascending altitudes[m] yield mass per unit area[kg/m^2].
  // The moving spline index/segment sum retain Native profile order.
  // Arguments (Native order):
  // total_idx_lin_spline: The total number of indices of the spline points []
  // lin_spline_rho_vals: Dry, static density                   [kg/m^3]
  // lin_spline_rho_levels: The altitudes of the given rho values [m]
  // grid_levels_idx: The total number of indices of the grid levels   []
  // grid_levels: levels of grid over which the masses should be calculated [m]

  auto workspace_scope = workspace.scope();
  using Space = typename Density::memory_space;

  Kokkos::View<int*, Kokkos::LayoutLeft, Space> failure =
      workspace.view<Kokkos::View<int*, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "calc_mass_over_grid_intervals::failure", 1);
  using Policy = Kokkos::RangePolicy<typename Density::execution_space>;
  Kokkos::parallel_for(
      "calc_mass_over_grid_intervals::profile", Policy(0, 1), KOKKOS_LAMBDA(int) {
        failure(0) = 0;
        // The spline must cover the highest requested grid edge.
        if (grid_levels(grid_levels_idx - 1) - lin_spline_rho_levels(total_idx_lin_spline - 1) >
            tol) {
          failure(0) = 1;
          return;
        }
        real_t val_below, level_below; // Density[kg/m^3] and altitude[m] below the interval.
        real_t mass_over_interval;     // Current cell mass per unit area[kg/m^2].
        real_t upper_zm_level_rho;     // Density at the upper cell edge[kg/m^3].
        int l              = 0;
        bool l_level_found = false;
        // Find the first level and initialize the lower segment endpoint.
        while (!l_level_found && l < total_idx_lin_spline) {
          if (Kokkos::abs(grid_levels(0) - lin_spline_rho_levels(l)) < tol) {
            val_below     = lin_spline_rho_vals(l);
            level_below = grid_levels(0);
            l_level_found = true;
          } else if (grid_levels(0) < lin_spline_rho_levels(l)) {
            if (l == 0) {
              failure(0) = 2;
              return;
            } else {
              val_below     = lin_interpolate_two_points(grid_levels(0),               // In
                                                         lin_spline_rho_levels(l),     // In
                                                         lin_spline_rho_levels(l - 1), // In
                                                         lin_spline_rho_vals(l),       // In
                                                         lin_spline_rho_vals(l - 1)    // In
              );
              level_below = grid_levels(0);
              l_level_found = true;
              l = l - 1;
            }
          }
          l = l + 1;
        }
        // Add all spline segments inside each cell, then its final partial segment.
        for (int k = 0; k < grid_levels_idx - 1; ++k) {
          mass_over_interval = real_t(0.0);
          while (lin_spline_rho_levels(l) < grid_levels(k + 1) && l < total_idx_lin_spline - 1) {
            mass_over_interval = mass_over_interval + (lin_spline_rho_levels(l) - level_below) *
                                                          (val_below + lin_spline_rho_vals(l)) / 2;
            level_below        = lin_spline_rho_levels(l);
            val_below          = lin_spline_rho_vals(l);
            l = l + 1;
          }
          if (Kokkos::abs(grid_levels(k + 1) - lin_spline_rho_levels(l)) < tol) {
            // Approximately equal altitudes: copy the spline value.
            upper_zm_level_rho = lin_spline_rho_vals(l);
          } else {
            // Spline node l lies above the upper cell edge.
            upper_zm_level_rho = lin_interpolate_two_points(grid_levels(k + 1),           // In
                                                            lin_spline_rho_levels(l),     // In
                                                            lin_spline_rho_levels(l - 1), // In
                                                            lin_spline_rho_vals(l),       // In
                                                            lin_spline_rho_vals(l - 1)    // In
            );
          }
          mass_over_interval   = mass_over_interval + (grid_levels(k + 1) - level_below) *
                                                          (val_below + upper_zm_level_rho) / 2;
          mass_per_interval(k) = mass_over_interval;
          val_below = upper_zm_level_rho;
          level_below          = grid_levels(k + 1);
        }
      });
  const auto failure_host = to_host_space(workspace, failure);
  if (failure_host(0)) {
    if (failure_host(0) == 1)
      std::cerr << "Wrong input: highest level of lin_spline_rho_levels must"
                   " be higher or equal to highest zm level of grid_levels\n";
    else
      std::cerr << "Wrong input: lowest level of lin_spline_rho_levels must"
                   " be lower or equal to lowest zm level of grid_levels\n";
    std::cerr << "Cannot compute the mass to get the remapping operator if the"
                 "density function is not defined over the whole grid region\n";
    throw std::runtime_error("Cannot compute the remapping operator");
  }
}

//=============================================================================
template <class Vectors, class Matrices, class Output>
void matrix_vector_mult(int ngrdcol,                   // In
                        int dim_input, int dim_output, // In
                        Vectors x_vectors,             // In
                        Matrices A_matrices,           // In
                        Output y_vectors               // Out
) {
  // Compute y=Ax for every column; each dot product retains ascending j order.
  // Arguments (Native order):
  // ngrdcol: Number of grid columns []
  // dim_input: dimension of the input vector []
  // dim_output: dimension of the output vector []
  // x_vectors: values of the x vectors (ngrdcol many vectors)
  // A_matrices: values of the matrices (ngrdcol many matrices)

  using Policy = Kokkos::RangePolicy<typename Vectors::execution_space>;
  Kokkos::parallel_for(
      "matrix_vector_mult", Policy(0, ngrdcol * dim_output), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t value = 0;
        for (int j = 0; j < dim_input; ++j)
          value = value + A_matrices(i, k, j) * x_vectors(i, j);
        y_vectors(i, k) = value;
      });
}

//=============================================================================
template <class Mass, class Field>
real_t vertical_integral_conserve_mass(ScratchWorkspace& workspace,
                                       int grid_levels_idx,                // InOut, In
                                       Field field, Mass mass_per_interval // In
) {
  // Integrate field times the exact density-spline cell masses [kg/m^2]. The
  // remapping helper already computed these masses with the same profile fold.
  auto workspace_scope = workspace.scope();
  using Space          = typename Mass::memory_space;
  auto integral        = workspace.view<Kokkos::View<real_t*, Kokkos::LayoutLeft, Space>>(
      Kokkos::WithoutInitializing, "vertical_integral_conserve_mass", 1);
  using Policy = Kokkos::RangePolicy<typename Mass::execution_space>;
  Kokkos::parallel_for(
      "vertical_integral_conserve_mass::ordered_sum", Policy(0, 1), KOKKOS_LAMBDA(int) {
        real_t vertical_integral_conserve_mass = zero;
        // Keep the ascending Native SUM fold.
        for (int k = 0; k < grid_levels_idx - 1; ++k)
          vertical_integral_conserve_mass =
              vertical_integral_conserve_mass + field(k) * mass_per_interval(k);
        integral(0) = vertical_integral_conserve_mass;
      });
  return to_host_space(workspace, integral)(0);
}

//=============================================================================
template <class Source, class Target, class SourceMass, class TargetMass>
void check_conservation(ScratchWorkspace& workspace,                                     // InOut
                        int levels_source_idx, int levels_target_idx,                    // In
                        Source field_source, Target field_target,                        // In
                        SourceMass mass_on_source_cells, TargetMass mass_on_target_cells // In
) {
  // Compare spline-weighted vertical integrals before/after remapping. Reuse
  // the exact cell masses already calculated by the remapping helper.
  // Both relative and absolute tolerances are retained in Native order.
  auto workspace_scope         = workspace.scope();
  const real_t integral_source = vertical_integral_conserve_mass(
      workspace, levels_source_idx, field_source, mass_on_source_cells);
  const real_t integral_target = vertical_integral_conserve_mass(
      workspace, levels_target_idx, field_target, mass_on_target_cells);
  const real_t tol_cons = real_t(1.e-3);
  const real_t err_percentage = real_t(1.e-7);
  if (Kokkos::abs(integral_target - integral_source) >
          Kokkos::abs(integral_source * err_percentage) &&
      Kokkos::abs(integral_target - integral_source) > tol_cons) {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10)
            << "WARNING! Integral for field was not conserved.\n"
            << "Integral was " << integral_target << " instead of " << integral_source << '\n';
    std::cerr << message.str();
  }
}

//=============================================================================
template <class Source, class Target>
void check_monotonicity(ScratchWorkspace& workspace,                // InOut
                        int field_source_idx, int field_target_idx, // In
                        Source field_source, Target field_target,   // In
                        int grid_remap_method                       // In
) {
  // Warn for target values outside source extrema. Native PPM omits the first
  // and last three target cells from this global-monotonicity diagnostic.
  // Arguments (Native order):
  // field_source_idx: dimension of field_source []
  // field_target_idx: dimension of field_target []
  // field_source: values for a variable on the source grid
  // field_target: remapped values for a variable on the target grid
  // grid_remap_method: integer specifying which remapping method was used

  auto workspace_scope = workspace.scope();
  int start_index = 0, end_index = field_target_idx - 1;
  if (grid_remap_method == ppm_remap) {
    start_index = 3;
    end_index   = field_target_idx - 4;
  }
  const auto source_host = to_host_space(workspace, field_source);
  const auto target_host = to_host_space(workspace, field_target);
  real_t source_maxval = -std::numeric_limits<real_t>::max();
  real_t source_minval = std::numeric_limits<real_t>::max();
  for (int k = 0; k < field_source_idx; ++k) {
    if (source_host(k + 1) > source_maxval)
      source_maxval = source_host(k + 1);
    if (source_host(k + 1) < source_minval)
      source_minval = source_host(k + 1);
  }
  for (int k = start_index; k <= end_index; ++k) {
    const real_t value = target_host(k);
    if ((value < source_minval && Kokkos::abs(value - source_minval) > tol) ||
        (value > source_maxval && Kokkos::abs(value - source_maxval) > tol)) {
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10)
              << "WARNING! The remapped values are not monotone.\n"
              << "The extrema on the source field were: " << source_minval << " and "
              << source_maxval << ", but the remapped value is: " << value << '\n';
      std::cerr << message.str();
    }
  }
}

//=============================================================================
template <class Matrix>
void check_consistency(ScratchWorkspace& workspace,   // InOut
                       int dim_input, int dim_output, // In
                       Matrix R_ij                    // In
) {
  // Apply the actual matrix-vector owner to a vector of ones and inspect every
  // output in Native order. Singleton dimensions replace sequence association.
  // Arguments (Native order):
  // dim_input: dimension of the input values on the source grid []
  // dim_output: dimension of the remapped values on the target grid []
  // R_ij: matrix to apply to input values for remapping (R_ij) []

  auto workspace_scope = workspace.scope();
  using Space = typename Matrix::memory_space;
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> ones =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "check_consistency::ones", 1, dim_input);
  Kokkos::deep_copy(ones, one);

  using Matrix3 =
      Kokkos::View<real_t***, Kokkos::LayoutStride, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  const Matrix3 R_matrices(R_ij.data(), Kokkos::LayoutStride(1, 1, dim_output, R_ij.stride(0),
                                                             dim_input, R_ij.stride(1)));
  auto output_ones = workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(
      Kokkos::WithoutInitializing, "matrix_vector_mult", 1, dim_output);
  matrix_vector_mult(1, dim_input, dim_output, ones, R_matrices, output_ones);
  const auto output_host = to_host_space(workspace, output_ones);
  bool l_consistent = true;
  int k                  = 0;
  while (l_consistent && k < dim_output) {
    if (Kokkos::abs(output_host(0, k) - one) > tol)
      l_consistent = false;
    k = k + 1;
  }
  if (!l_consistent) {
    std::ostringstream message;
    message << std::setprecision(std::numeric_limits<real_t>::max_digits10)
            << "WARNING! The matrix_vector_mult function is not consistent\n"
            << "The output should be all ones, but was: ";
    for (int n = 0; n < dim_output; ++n)
      message << output_host(0, n) << ' ';
    std::cerr << message.str() << '\n';
  }
}

//=============================================================================
void remap_vals_to_target_helper(ScratchWorkspace& workspace,                  // InOut
                                 int ngrdcol,                                  // In
                                 int levels_source_idx, int levels_target_idx, // In
                                 Real2D levels_source, Real2D levels_target,   // In
                                 int total_idx_rho_lin_spline,                 // In
                                 Real2D rho_lin_spline_vals,                   // In
                                 Real2D rho_lin_spline_levels,                 // In
                                 Real2D source_values,                         // In
                                 int iv, Real1D p_sfc,                         // In
                                 int grid_remap_method,                        // In
                                 Real2D target_values,                         // Out
                                 const Real3D* R_ij = nullptr                  // In Optional
) {
  // Integrate the density spline, build pressures and remap the source field.
  // An optional Ullrich operator is reused; debug2 checks still run afterward.
  // Arguments (Native order):
  // ngrdcol: Number of grid columns []
  // levels_source_idx: Number of levels in the source grid [].
  // levels_target_idx: Number of levels in the target grid [].
  // levels_source: altitudes for the source grid [m]
  // levels_target: altitudes for the target grid [m]
  // total_idx_rho_lin_spline: dimension of the linear spline definition arrays for rho []
  // rho_lin_spline_vals: rho values at the given altitudes  [kg/m^3]
  // rho_lin_spline_levels: altitudes for the given rho values [m]
  // source_values: given values on the source grid that should be remapped to the target grid
  // iv: setting for E3SM's PPM method, that sets what type of value is remapped, the following are
  // the three valid settings: iv =-1: winds iv = 0: positive definite scalars iv = 1: others
  // p_sfc: prescribed pressure at the surface [Pa]
  // grid_remap_method: specifies which remapping method should be used options: 1. Ullrich linear
  // remapping scheme (first order) 2. Piecewise Parabolic Method (PPM) (high order)
  // R_ij: Optional read-only (column,target cell,source cell) remapping operator [].

  auto workspace_scope        = workspace.scope();
  Real2D mass_on_source_cells = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_to_target_helper::mass_on_source_cells", ngrdcol,
      levels_source_idx - 1);
  Real2D mass_on_target_cells = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_to_target_helper::mass_on_target_cells", ngrdcol,
      levels_target_idx - 1);
  Real2D pressure_levels_source = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_to_target_helper::pressure_levels_source", ngrdcol,
      levels_source_idx);
  Real2D pressure_levels_target = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_to_target_helper::pressure_levels_target", ngrdcol,
      levels_target_idx);
  Real3D R_ij_internal;
  if (grid_remap_method == cons_ullrich_remap && R_ij == nullptr)
    R_ij_internal = workspace.view<Real3D>(Kokkos::WithoutInitializing,
                                           "remap_vals_to_target_helper::R_ij_internal", ngrdcol,
                                           levels_target_idx - 1, levels_source_idx - 1);

  for (int i = 0; i < ngrdcol; ++i) {
    // Calculate pressure levels from altitudes and densities. These operations
    // are also required when the optional Ullrich operator has been supplied.
    calc_mass_over_grid_intervals(workspace, total_idx_rho_lin_spline,
                                  Kokkos::subview(rho_lin_spline_vals, i, Kokkos::ALL),
                                  Kokkos::subview(rho_lin_spline_levels, i, Kokkos::ALL),
                                  levels_source_idx, Kokkos::subview(levels_source, i, Kokkos::ALL),
                                  Kokkos::subview(mass_on_source_cells, i, Kokkos::ALL));
    calc_mass_over_grid_intervals(workspace, total_idx_rho_lin_spline,
                                  Kokkos::subview(rho_lin_spline_vals, i, Kokkos::ALL),
                                  Kokkos::subview(rho_lin_spline_levels, i, Kokkos::ALL),
                                  levels_target_idx, Kokkos::subview(levels_target, i, Kokkos::ALL),
                                  Kokkos::subview(mass_on_target_cells, i, Kokkos::ALL));
    Kokkos::parallel_for(
        "remap_vals_to_target_helper::pressure", 1, KOKKOS_LAMBDA(int) {
          pressure_levels_source(i, 0) = p_sfc(i);
          for (int k = 1; k < levels_source_idx; ++k)
            pressure_levels_source(i, k) =
                pressure_levels_source(i, k - 1) - mass_on_source_cells(i, k - 1) * grav;
          pressure_levels_target(i, 0) = p_sfc(i);
          for (int k = 1; k < levels_target_idx; ++k)
            pressure_levels_target(i, k) =
                pressure_levels_target(i, k - 1) - mass_on_target_cells(i, k - 1) * grav;
        });
    // Construct an Ullrich operator only when none was supplied.
    if (grid_remap_method == cons_ullrich_remap && R_ij == nullptr)
      remapping_matrix(workspace, levels_source_idx, levels_target_idx,
                       Kokkos::subview(pressure_levels_source, i, Kokkos::ALL),
                       Kokkos::subview(pressure_levels_target, i, Kokkos::ALL),
                       Kokkos::subview(R_ij_internal, i, Kokkos::ALL, Kokkos::ALL));
  }

  // Remap values with the selected remapping method.
  if (grid_remap_method == cons_ullrich_remap) {
    if (R_ij != nullptr)
      matrix_vector_mult(ngrdcol, levels_source_idx - 1, levels_target_idx - 1, source_values,
                         *R_ij, target_values);
    else
      matrix_vector_mult(ngrdcol, levels_source_idx - 1, levels_target_idx - 1, source_values,
                         R_ij_internal, target_values);
  } else if (grid_remap_method == ppm_remap) {
    remap_vals_ppm(workspace, ngrdcol, levels_source_idx, levels_target_idx, pressure_levels_source,
                   pressure_levels_target, source_values, iv, target_values);
  } else {
    std::cerr << "There is no remapping method implemented for grid_remap_method="
              << grid_remap_method << '\n';
    throw std::runtime_error("Please adjust the configurable model flags.");
  }

  // Check remapped values, including a supplied operator, in Native column order.
  if (clubb_at_least_debug_level_api(2)) {
    for (int i = 0; i < ngrdcol; ++i) {
      check_conservation(workspace, levels_source_idx, levels_target_idx,
                         Kokkos::subview(source_values, i, Kokkos::ALL),
                         Kokkos::subview(target_values, i, Kokkos::ALL),
                         Kokkos::subview(mass_on_source_cells, i, Kokkos::ALL),
                         Kokkos::subview(mass_on_target_cells, i, Kokkos::ALL));
      check_monotonicity(workspace, levels_source_idx - 1, levels_target_idx - 1, // In
                         Kokkos::subview(source_values, i, Kokkos::ALL),          // In
                         Kokkos::subview(target_values, i, Kokkos::ALL),          // In
                         grid_remap_method                                        // In
      );
      if (grid_remap_method == cons_ullrich_remap) {
        if (R_ij != nullptr)
          check_consistency(workspace, levels_source_idx - 1, levels_target_idx - 1, // In
                            Kokkos::subview(*R_ij, i, Kokkos::ALL, Kokkos::ALL)      // In
          );
        else
          check_consistency(workspace, levels_source_idx - 1, levels_target_idx - 1,    // In
                            Kokkos::subview(R_ij_internal, i, Kokkos::ALL, Kokkos::ALL) // In
          );
      }
    }
  }
}
} // namespace

//=============================================================================
static void remap_vals_to_target_values(ScratchWorkspace& workspace,                  // InOut
                                        int ngrdcol,                                  // In
                                        const grid& gr_source, const grid& gr_target, // In
                                        int source_values_idx,                        // In
                                        Real2D source_values,                         // In
                                        int target_values_idx,                        // In
                                        int total_idx_rho_lin_spline,                 // In
                                        Real2D rho_lin_spline_vals,                   // In
                                        Real2D rho_lin_spline_levels,                 // In
                                        int iv, Real1D p_sfc,                         // In
                                        int grid_remap_method,                        // In
                                        bool l_zt_variable,                           // In
                                        Real2D target_values                          // Out
) {
  // Remap cell means between grids. Density splines define exact cell masses;
  // momentum-level fields first receive surrounding thermodynamic grid edges.
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

  auto workspace_scope = workspace.scope();
  if (l_zt_variable) {
    // Surround thermodynamic-level values with the momentum levels.
    remap_vals_to_target_helper(workspace, ngrdcol, gr_source.nzm, gr_target.nzm, gr_source.zm,
                                gr_target.zm, total_idx_rho_lin_spline, rho_lin_spline_vals,
                                rho_lin_spline_levels, source_values, iv, p_sfc, grid_remap_method,
                                target_values);
  } else if (!l_zt_variable) {
    Real2D levels_source =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "remap_vals_to_target::levels_source",
                               ngrdcol, source_values_idx + 1);
    Real2D levels_target =
        workspace.view<Real2D>(Kokkos::WithoutInitializing, "remap_vals_to_target::levels_target",
                               ngrdcol, target_values_idx + 1);
    // For momentum-level values use thermodynamic levels inside the grid;
    // at either exterior boundary use the momentum and thermodynamic levels.
    const auto source_zm = gr_source.zm, source_zt = gr_source.zt;
    const auto target_zm = gr_target.zm, target_zt = gr_target.zt;
    const int source_nzt = gr_source.nzt, source_nzm = gr_source.nzm;
    const int target_nzt = gr_target.nzt, target_nzm = gr_target.nzm;
    Kokkos::parallel_for(
        "remap_vals_to_target::levels", ngrdcol, KOKKOS_LAMBDA(int i) {
          levels_source(i, 0) = source_zm(i, 0);
          levels_target(i, 0) = target_zm(i, 0);
          for (int k = 0; k < source_nzt; ++k)
            levels_source(i, k + 1) = source_zt(i, k);
          for (int k = 0; k < target_nzt; ++k)
            levels_target(i, k + 1) = target_zt(i, k);
          levels_source(i, source_nzt + 1) = source_zm(i, source_nzm - 1);
          levels_target(i, target_nzt + 1) = target_zm(i, target_nzm - 1);
        });
    remap_vals_to_target_helper(workspace, ngrdcol, gr_source.nzt + 2, gr_target.nzt + 2,
                                levels_source, levels_target, total_idx_rho_lin_spline,
                                rho_lin_spline_vals, rho_lin_spline_levels, source_values, iv,
                                p_sfc, grid_remap_method, target_values);
  }
}

namespace {
//-----------------------------------------------------------------------------
// LINEAR REMAPPING METHOD BY ULLRICH ET AL.
// Formula (30), Arbitrary-Order Conservative and Consistent Remapping and a
// Theory of Linear Maps: Part II.
//-----------------------------------------------------------------------------
template <class Matrix, class SourceLevels, class TargetLevels>
void check_remap_matrix_conservation(ScratchWorkspace& workspace,                           // InOut
                                     Matrix R_ij, int dim1, int dim2,                       // In
                                     SourceLevels levels_source, TargetLevels levels_target // In
) {
  // Conservation condition (9), Ullrich/Taylor Part I: maintain global mass.
  // Arguments (Native order):
  // R_ij: matrix to apply to input values for remapping (R_ij) []
  // dim1: Number of target cells (matrix first axis) []
  // dim2: dimensions of the remapping matrix []
  // levels_source: the pressure levels in the source grid [kg/(ms^2)]
  // levels_target: the pressure levels in the target grid [kg/(ms^2)]

  auto workspace_scope = workspace.scope();
  const auto R         = to_host_space(workspace, R_ij);
  const auto source    = to_host_space(workspace, levels_source);
  const auto target    = to_host_space(workspace, levels_target);
  bool l_conservative = true;
  int l                = 0;
  while (l_conservative && l < dim2) {
    real_t weighted_col_sum = 0;
    for (int k = 0; k < dim1; ++k)
      weighted_col_sum = weighted_col_sum + R(k, l) * Kokkos::abs(target(k + 1) - target(k));
    if (Kokkos::abs(weighted_col_sum - Kokkos::abs(source(l + 1) - source(l))) > tol)
      l_conservative = false;
    l = l + 1;
  }
  if (!l_conservative) {
    std::cerr << "WARNING! The remapping operator is not conservative\n";
    throw std::runtime_error("Operator should be conservative, something went wrong...");
  }
}

template <class Matrix>
void check_remap_matrix_consistency(ScratchWorkspace& workspace,    // InOut
                                    Matrix R_ij, int dim1, int dim2 // In
) {
  // Consistency condition (12), Ullrich/Taylor Part I: a constant field stays
  // constant. In particular, applying the operator to ones must return ones.
  // Arguments (Native order):
  // R_ij: matrix to apply to input values for remapping (R_ij) []
  // dim1: Number of target cells (matrix first axis) []
  // dim2: dimension of the remapping matrix []

  auto workspace_scope = workspace.scope();
  const auto R         = to_host_space(workspace, R_ij);
  bool l_consistent = true;
  int k                = 0;
  while (l_consistent && k < dim1) {
    real_t row_sum = 0;
    for (int l = 0; l < dim2; ++l)
      row_sum = row_sum + R(k, l);
    if (Kokkos::abs(row_sum - 1) > tol)
      l_consistent = false;
    k = k + 1;
  }
  if (!l_consistent) {
    std::cerr << "WARNING! The remapping operator is not consistent\n";
    throw std::runtime_error("Operator should be consistent, something went wrong...");
  }
}

template <class Matrix>
void check_remap_matrix_monotonicity(ScratchWorkspace& workspace,    // InOut
                                     Matrix R_ij, int dim1, int dim2 // In
) {
  // Monotonicity condition (14), Ullrich/Taylor Part I: introduce no additional
  // global extrema. A remapping operator is monotone when all entries are >=0.
  // Arguments (Native order):
  // R_ij: matrix to apply to input values for remapping (R_ij) []
  // dim1: Number of target cells (matrix first axis) []
  // dim2: dimensions of the remapping matrix []

  auto workspace_scope = workspace.scope();
  const auto R         = to_host_space(workspace, R_ij);
  bool l_monotone = true;
  int k                = 0;
  while (l_monotone && k < dim1) {
    for (int l = 0; l < dim2; ++l)
      if (R(k, l) < 0)
        l_monotone = false;
    k = k + 1;
  }
  if (!l_monotone) {
    std::cerr << "WARNING! The remapping operator is not monotone\n";
    throw std::runtime_error("Operator should be monotone, something went wrong...");
  }
}

template <class SourceLevels, class TargetLevels, class Output>
void remapping_matrix(ScratchWorkspace& workspace,                            // InOut
                      int levels_source_idx, int levels_target_idx,           // In
                      SourceLevels levels_source, TargetLevels levels_target, // In
                      Output remapping_matrix                                 // Out
) {
  // Ullrich et al. Part II, formula (30). Pressure edges are surface-to-top.
  // Arguments (Native order):
  // levels_source_idx: Number of levels in the source grid [].
  // levels_target_idx: Number of levels in the target grid [].
  // levels_source: the pressure levels in the source grid [kg/(ms^2)]
  // levels_target: the pressure levels in the target grid [kg/(ms^2)]

  auto workspace_scope = workspace.scope();
  using Exec = typename SourceLevels::execution_space;

  Kokkos::parallel_for(
      "remapping_matrix",
      Kokkos::RangePolicy<Exec>(0, (levels_target_idx - 1) * (levels_source_idx - 1)),
      KOKKOS_LAMBDA(int index) {
        const int k = index % (levels_target_idx - 1), l = index / (levels_target_idx - 1);
        const real_t omega_ov_upper = Kokkos::max(levels_target(k + 1), levels_source(l + 1));
        const real_t omega_ov_lower = Kokkos::min(levels_target(k), levels_source(l));
        real_t omega_ov;
        if (omega_ov_upper < omega_ov_lower)
          omega_ov = Kokkos::abs(omega_ov_upper - omega_ov_lower);
        else
          omega_ov = 0;
        remapping_matrix(k, l) = omega_ov / Kokkos::abs(levels_target(k + 1) - levels_target(k));
      });
  if (clubb_at_least_debug_level_api(2)) {
    check_remap_matrix_conservation(workspace,
                                    remapping_matrix,                             // In
                                    levels_target_idx - 1, levels_source_idx - 1, // In
                                    levels_source, levels_target                  // In
    );
    check_remap_matrix_consistency(workspace,
                                   remapping_matrix,                            // In
                                   levels_target_idx - 1, levels_source_idx - 1 // In
    );
    check_remap_matrix_monotonicity(workspace,
                                    remapping_matrix,                            // In
                                    levels_target_idx - 1, levels_source_idx - 1 // In
    );
  }
}

//-----------------------------------------------------------------------------
// PIECEWISE PARABOLIC METHOD BY ULLRICH AND WOODWARD
// The Piecewise Parabolic Method (PPM) for gas-dynamical simulations.
//
// Obtained from E3SM, sponsored by the U.S. Department of Energy, Office of
// Science, Biological/Environmental Research, Earth Systems Model Development.
// https://github.com/E3SM-Project/E3SM/blob/master/components/eam/src/dynamics/fv/mapz_module.F90
// kmppm, steepz, ppm2m and map1_ppm below retain that implementation.
//-----------------------------------------------------------------------------
// Scalar-cell overload is the coefficient-shape/device adaptation of kmppm;
// the Native-signature array overload below dispatches independent longitudes.
template <class Coefficients>
KOKKOS_INLINE_FUNCTION void kmppm(real_t dm, Coefficients a4, int lmt) {
  constexpr real_t r12 = one / twelve;
  // Developer: S.-J. Lin, NASA-GSFC.
  if (lmt == 3)
    return;
  if (lmt == 0) {
    // Standard PPM constraint.
    if (dm == zero) {
      a4(1) = a4(0);
      a4(2) = a4(0);
      a4(3) = zero;
    } else {
      const real_t da1 = a4(2) - a4(1);
      const real_t da2 = da1 * da1;
      const real_t a6da = a4(3) * da1;
      if (a6da < -da2) {
        a4(3) = three * (a4(1) - a4(0));
        a4(2) = a4(1) - a4(3);
      } else if (a6da > da2) {
        a4(3) = three * (a4(2) - a4(0));
        a4(1) = a4(2) - a4(3);
      }
    }
  } else if (lmt == 1) {
    // Improved full monotonicity constraint (Lin). No first guess of A6 is read.
    const real_t qmp = two * dm;
    a4(1) =
        a4(0) - Kokkos::copysign(Kokkos::min(Kokkos::abs(qmp), Kokkos::abs(a4(1) - a4(0))), qmp);
    a4(2) =
        a4(0) + Kokkos::copysign(Kokkos::min(Kokkos::abs(qmp), Kokkos::abs(a4(2) - a4(0))), qmp);
    a4(3) = three * (two * a4(0) - (a4(1) + a4(2)));
  } else if (lmt == 2) {
    // Positive definite constraint.
    if (Kokkos::abs(a4(2) - a4(1)) < -a4(3)) {
      const real_t fmin =
          a4(0) + one_fourth * ((a4(2) - a4(1)) * (a4(2) - a4(1))) / a4(3) + a4(3) * r12;
      if (fmin < zero) {
        if (a4(0) < a4(2) && a4(0) < a4(1)) {
          a4(2) = a4(0);
          a4(1) = a4(0);
          a4(3) = zero;
        } else if (a4(2) > a4(1)) {
          a4(3) = three * (a4(1) - a4(0));
          a4(2) = a4(1) - a4(3);
        } else {
          a4(3) = three * (a4(2) - a4(0));
          a4(1) = a4(2) - a4(3);
        }
      }
    }
  }
}

template <class Difference, class Coefficients>
void kmppm(
    Difference dm,                                                                  // In
    Coefficients a4,                                                                // InOut
    int itot, int lmt                                                               // In
) {
  // Apply the PPM monotonicity constraint independently to each longitude.
  // Arguments (Native order):
  // dm: Monotonic mismatch; same field units as a4.
  // a4: PPM coefficients AA, AL, AR, A6, stored as (column,coefficient).
  // itot: Total Longitudes
  // lmt: 0: Standard PPM constraint 1: Improved full monotonicity constraint (Lin) 2: Positive
  // definite constraint 3: do nothing (return immediately)

  // a4 is (column,coefficient); lmt3 returns without reading either input.
  if (lmt == 3)
    return;
  using Exec = typename Coefficients::execution_space;
  Kokkos::parallel_for(
      "kmppm", Kokkos::RangePolicy<Exec>(0, itot),
      KOKKOS_LAMBDA(int i) { kmppm(dm(i), Kokkos::subview(a4, i, Kokkos::ALL), lmt); });
}

//=============================================================================
template <class Coefficients, class Mismatch, class Monotonic, class Difference, class Pressure,
          class Sum>
void steepz(ScratchWorkspace& workspace,                                   // InOut
            int i1, int i2, int km,                                        // In
            Coefficients a4,                                               // InOut
            Mismatch df2, Monotonic dm, Difference dq, Pressure dp, Sum d4 // In
) {
  // Discontinuity adjustment defined in the PPM paper. Private arrays map
  // the inclusive longitude range onto local columns; levels are zero based.
  // Arguments (Native order):
  // i1: Starting longitude
  // i2: Finishing longitude
  // km: Total levels
  // a4: PPM coefficients AA, AL, AR, A6, stored as (column,level,coefficient).
  // df2: first guess mismatch
  // dm: monotonic mismatch
  // dq: backward diff of q
  // dp: grid size
  // d4: backward sum:  dp(k)+ dp(k-1)

  auto workspace_scope = workspace.scope();
  using Space = typename Coefficients::memory_space;
  using Exec = typename Coefficients::execution_space;
  const int columns = i2 - i1 + 1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> alfa =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "steepz::alfa", columns, km);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> f =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "steepz::f", columns, km);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> rat =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "steepz::rat", columns, km);
  // Compute ratio of dq/dp.
  Kokkos::parallel_for(
      "steepz::rat", Kokkos::RangePolicy<Exec>(0, columns * (km - 1)), KOKKOS_LAMBDA(int index) {
        const int i = index % columns, k = index / columns + 1;
        rat(i, k) = dq(i, k - 1) / d4(i, k);
      });
  // Compute F.
  Kokkos::parallel_for(
      "steepz::f", Kokkos::RangePolicy<Exec>(0, columns * (km - 2)), KOKKOS_LAMBDA(int index) {
        const int i = index % columns, k = index / columns + 1;
        f(i, k) = (rat(i, k + 1) - rat(i, k)) / (dp(i, k - 1) + dp(i, k) + dp(i, k + 1));
      });
  Kokkos::parallel_for(
      "steepz::alfa", Kokkos::RangePolicy<Exec>(0, columns * Kokkos::max(0, km - 4)),
      KOKKOS_LAMBDA(int index) {
        const int i = index % columns, k = index / columns + 2;
        if (f(i, k + 1) * f(i, k - 1) < zero && df2(i, k) != zero) {
          const real_t dg2 = (f(i, k + 1) - f(i, k - 1)) *
                             ((dp(i, k + 1) - dp(i, k - 1)) * (dp(i, k + 1) - dp(i, k - 1)) +
                              d4(i, k) * d4(i, k + 1));
          // Prvalues keep device min/max from binding host constants by reference.
          alfa(i, k) = Kokkos::max(
              real_t(zero), Kokkos::min(real_t(one_half), -three_sixteenth * dg2 / df2(i, k)));
        } else
          alfa(i, k) = zero;
      });
  Kokkos::parallel_for(
      "steepz::a4", Kokkos::RangePolicy<Exec>(0, columns * Kokkos::max(0, km - 5)),
      KOKKOS_LAMBDA(int index) {
        const int i = index % columns, k = index / columns + 3;
        a4(i, k, 1) = (one - alfa(i, k - 1) - alfa(i, k)) * a4(i, k, 1) +
                      alfa(i, k - 1) * (a4(i, k, 0) - dm(i, k)) +
                      alfa(i, k) * (a4(i, k - 1, 0) + dm(i, k - 1));
      });
}

//=============================================================================
template <class Coefficients, class Pressure>
void ppm2m(ScratchWorkspace& workspace,                            // InOut
           Coefficients a4,                                        // InOut
           Pressure delp, int km, int i1, int i2, int iv, int kord // In
) {
  // Reconstruct PPM cell edges and constrain their monotonicity.
  // a4 is (column,level,coefficient): AA,AL,AR,A6 retain Native order 1:4.
  // Arguments (Native order):
  // a4: PPM coefficients AA, AL, AR, A6, stored as (column,level,coefficient).
  // delp: layer pressure thickness
  // km: vertical dimension
  // i1: Starting longitude
  // i2: Finishing longitude
  // iv: iv =-1: winds iv = 0: positive definite scalars iv = 1: others
  // kord: Order (or more accurately method no.):

  auto workspace_scope = workspace.scope();
  using Space = typename Coefficients::memory_space;
  using Exec = typename Coefficients::execution_space;
  const int km1 = km - 1, it = i2 - i1 + 1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> dc =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "ppm2m::dc", it, km);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> h2;
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> delq =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "ppm2m::delq", it, km);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> df2 =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "ppm2m::df2", it, km);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> d4 =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                        "ppm2m::d4", it, km);

  Kokkos::parallel_for(
      "ppm2m::delq_d4", Kokkos::RangePolicy<Exec>(0, it * (km - 1)), KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it + 1;
        delq(i, k - 1) = a4(i, k, 0) - a4(i, k - 1, 0);
        d4(i, k)       = delp(i, k - 1) + delp(i, k);
      });
  Kokkos::parallel_for(
      "ppm2m::dc_df2", Kokkos::RangePolicy<Exec>(0, it * (km - 2)), KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it + 1;
        const real_t c1 = (delp(i, k - 1) + one_half * delp(i, k)) / d4(i, k + 1);
        const real_t c2 = (delp(i, k + 1) + one_half * delp(i, k)) / d4(i, k);
        const real_t tmp =
            delp(i, k) * (c1 * delq(i, k) + c2 * delq(i, k - 1)) / (d4(i, k) + delp(i, k + 1));
        const real_t qmax =
            Kokkos::max(Kokkos::max(a4(i, k - 1, 0), a4(i, k, 0)), a4(i, k + 1, 0)) - a4(i, k, 0);
        const real_t qmin =
            a4(i, k, 0) - Kokkos::min(Kokkos::min(a4(i, k - 1, 0), a4(i, k, 0)), a4(i, k + 1, 0));
        dc(i, k)  = Kokkos::copysign(Kokkos::min(Kokkos::min(Kokkos::abs(tmp), qmax), qmin), tmp);
        df2(i, k) = tmp;
      });

  // 4th order interpolation of the provisional cell edge value.
  Kokkos::parallel_for(
      "ppm2m::provisional_edges", Kokkos::RangePolicy<Exec>(0, it * (km - 3)),
      KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it + 2;
        const real_t c1 = delq(i, k - 1) * delp(i, k - 1) / d4(i, k);
        const real_t a1 = d4(i, k - 1) / (d4(i, k) + delp(i, k - 1));
        const real_t a2 = d4(i, k + 1) / (d4(i, k) + delp(i, k));
        a4(i, k, 1)     = a4(i, k - 1, 0) + c1 +
                          two / (d4(i, k - 1) + d4(i, k + 1)) *
                              (delp(i, k) * (c1 * (a1 - a2) + a2 * dc(i, k - 1)) -
                               delp(i, k - 1) * a1 * dc(i, k));
      });
  steepz(workspace, i1, i2, km,  // In
         a4,                     // InOut
         df2, dc, delq, delp, d4 // In
  );

  // Area-preserving cubic with second derivative zero at the top boundary.
  Kokkos::parallel_for(
      "ppm2m::top", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
        const real_t d1 = delp(i, 0), d2 = delp(i, 1);
        const real_t qm = (d2 * a4(i, 0, 0) + d1 * a4(i, 1, 0)) / (d1 + d2);
        const real_t dq = two * (a4(i, 1, 0) - a4(i, 0, 0)) / (d1 + d2);
        const real_t c1 =
            four * (a4(i, 2, 1) - qm - d2 * dq) / (d2 * (two * d2 * d2 + d1 * (d2 + three * d1)));
        const real_t c3 = dq - one_half * c1 * (d2 * (five * d1 + d2) - three * (d1 * d1));
        a4(i, 1, 1) = qm - one_fourth * c1 * d1 * d2 * (d2 + three * d1);
        a4(i, 0, 1) = d1 * (two * c1 * (d1 * d1) - c3) + a4(i, 1, 1);
        dc(i, 0) = a4(i, 0, 0) - a4(i, 0, 1);
        // No over- and undershoot condition.
        const real_t cmax = Kokkos::max(a4(i, 0, 0), a4(i, 1, 0));
        const real_t cmin = Kokkos::min(a4(i, 0, 0), a4(i, 1, 0));
        a4(i, 1, 1) = Kokkos::max(cmin, a4(i, 1, 1));
        a4(i, 1, 1) = Kokkos::min(cmax, a4(i, 1, 1));
      });
  if (iv == 0) {
    Kokkos::parallel_for(
        "ppm2m::top_positive", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
          // WS 02.04.23: FVGCM sets both edges to AA; CAM instead clips AL.
          a4(i, 0, 1) = Kokkos::max(real_t(zero), a4(i, 0, 1));
          a4(i, 1, 1) = Kokkos::max(real_t(zero), a4(i, 1, 1));
        });
  } else if (iv == -1) {
    // Winds: retain Native km>32 sponge, which the JAX implementation omits.
    if (km > 32) {
      Kokkos::parallel_for(
          "ppm2m::top_wind_sponge", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
            a4(i, 0, 1) = a4(i, 0, 0);
            a4(i, 0, 2) = a4(i, 0, 0);
          });
    } else {
      Kokkos::parallel_for(
          "ppm2m::top_wind", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
            if (a4(i, 0, 0) * a4(i, 0, 1) <= zero)
              a4(i, 0, 1) = zero;
            else
              a4(i, 0, 1) = Kokkos::copysign(
                  Kokkos::min(Kokkos::abs(a4(i, 0, 0)), Kokkos::abs(a4(i, 0, 1))), a4(i, 0, 0));
          });
    }
  }

  // Area-preserving cubic with second derivative zero at the surface.
  Kokkos::parallel_for(
      "ppm2m::bottom", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
        const real_t d1 = delp(i, km - 1), d2 = delp(i, km1 - 1);
        const real_t qm = (d2 * a4(i, km - 1, 0) + d1 * a4(i, km1 - 1, 0)) / (d1 + d2);
        const real_t dq = two * (a4(i, km1 - 1, 0) - a4(i, km - 1, 0)) / (d1 + d2);
        const real_t c1 =
            (a4(i, km1 - 1, 1) - qm - d2 * dq) / (d2 * (two * d2 * d2 + d1 * (d2 + three * d1)));
        const real_t c3 = dq - two * c1 * (d2 * (five * d1 + d2) - three * (d1 * d1));
        a4(i, km - 1, 1) = qm - c1 * d1 * d2 * (d2 + three * d1);
        a4(i, km - 1, 2) = d1 * (eight * c1 * (d1 * d1) - c3) + a4(i, km - 1, 1);
        dc(i, km - 1) = a4(i, km - 1, 2) - a4(i, km - 1, 0);
        // No over- and undershoot condition.
        const real_t cmax = Kokkos::max(a4(i, km - 1, 0), a4(i, km1 - 1, 0));
        const real_t cmin = Kokkos::min(a4(i, km - 1, 0), a4(i, km1 - 1, 0));
        a4(i, km - 1, 1) = Kokkos::max(cmin, a4(i, km - 1, 1));
        a4(i, km - 1, 1) = Kokkos::min(cmax, a4(i, km - 1, 1));
      });
  // Enforce the surface constraint.
  if (iv == 0) {
    Kokkos::parallel_for(
        "ppm2m::bottom_positive", Kokkos::RangePolicy<Exec>(0, it),
        KOKKOS_LAMBDA(int i) {
          a4(i, km - 1, 2) = Kokkos::max(real_t(zero), a4(i, km - 1, 2));
        });
  } else if (iv == -1) {
    Kokkos::parallel_for(
        "ppm2m::bottom_wind", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
          if (a4(i, km - 1, 0) * a4(i, km - 1, 2) <= zero)
            a4(i, km - 1, 2) = zero;
          else
            a4(i, km - 1, 2) = Kokkos::copysign(
                Kokkos::min(Kokkos::abs(a4(i, km - 1, 0)), Kokkos::abs(a4(i, km - 1, 2))),
                a4(i, km - 1, 0));
        });
  }
  Kokkos::parallel_for(
      "ppm2m::right_edges", Kokkos::RangePolicy<Exec>(0, it * km1), KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it;
        a4(i, k, 2) = a4(i, k + 1, 1);
      });

  // f(s) = AL+s*[(AR-AL)+A6*(1-s)] (0<=s<=1).
  // The top two and bottom two layers always use monotonic mapping.
  // The two boundary cells are independent. Compute A6 and apply the same
  // scalar limiter in one kernel, avoiding a separate launch for each stage/cell.
  Kokkos::parallel_for(
      "ppm2m::top_limiter", Kokkos::RangePolicy<Exec>(0, it * 2), KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it;
        a4(i, k, 3) = three * (two * a4(i, k, 0) - (a4(i, k, 1) + a4(i, k, 2)));
        kmppm(dc(i, k), Kokkos::subview(a4, i, k, Kokkos::ALL), 0);
      });
  if (kord >= 7) {
    h2 = workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(
        Kokkos::WithoutInitializing, "ppm2m::h2", it, km);
    // Huynh's second constraint. The public wrapper fixes kord=4, but this
    // Native private branch is retained; no public higher-order option is added.
    Kokkos::parallel_for(
        "ppm2m::Huynh_h2", Kokkos::RangePolicy<Exec>(0, it * (km - 2)), KOKKOS_LAMBDA(int index) {
          const int i = index % it, k = index / it + 1;
          // Method 1: h2=delq(k)-delq(k-1).
          // Method 2: two*(dc(k+1)/delp(k+1)-dc(k-1)/delp(k-1))
          //           /(delp(k)+one_half*(delp(k-1)+delp(k+1)))*delp(k)**2.
          // Method 3:
          h2(i, k) = dc(i, k + 1) - dc(i, k - 1);
        });
    real_t fac;
    if (kord == 7)
      fac = three_halves; // Original quasi-monotone.
    else
      fac = one_eighth; // Full monotone.
    for (int k = 2; k < km - 2; ++k) {
      Kokkos::parallel_for(
          "ppm2m::Huynh_edges", Kokkos::RangePolicy<Exec>(0, it), KOKKOS_LAMBDA(int i) {
            // Right edges. Alternative guesses in E3SM use two*delq(k-1) and
            // fac*h2(k-1)+one_half*delq(k-1), respectively.
            const real_t pmp = two * dc(i, k);
            real_t qmp       = a4(i, k, 0) + pmp;
            real_t lac       = a4(i, k, 0) + fac * h2(i, k - 1) + dc(i, k);
            real_t qmin      = Kokkos::min(Kokkos::min(a4(i, k, 0), qmp), lac);
            real_t qmax      = Kokkos::max(Kokkos::max(a4(i, k, 0), qmp), lac);
            a4(i, k, 2)      = Kokkos::min(Kokkos::max(a4(i, k, 2), qmin), qmax);
            // Left edges. Alternative guesses use two*delq(k) and
            // fac*h2(k+1)-one_half*delq(k), respectively.
            qmp         = a4(i, k, 0) - pmp;
            lac         = a4(i, k, 0) + fac * h2(i, k + 1) - dc(i, k);
            qmin        = Kokkos::min(Kokkos::min(a4(i, k, 0), qmp), lac);
            qmax        = Kokkos::max(Kokkos::max(a4(i, k, 0), qmp), lac);
            a4(i, k, 1) = Kokkos::min(Kokkos::max(a4(i, k, 1), qmin), qmax);
            // Recompute A6.
            a4(i, k, 3) = three * (two * a4(i, k, 0) - (a4(i, k, 1) + a4(i, k, 2)));
          });
      // Additional constraint to prevent negatives when kord=7.
      if (iv == 0 && kord == 7)
        kmppm(Kokkos::subview(dc, Kokkos::ALL, k),              // In
              Kokkos::subview(a4, Kokkos::ALL, k, Kokkos::ALL), // InOut
              it, 2                                             // In
        );
    }
  } else {
    int lmt = kord - 3;
    lmt = Kokkos::max(0, lmt);
    if (iv == 0)
      lmt = Kokkos::min(2, lmt);
    // Fuse independent interior-layer limiter calls. Keep A6 unread until the
    // branch writes/uses it: kord4's Lin limiter initializes A6 itself.
    Kokkos::parallel_for(
        "ppm2m::interior_limiter", Kokkos::RangePolicy<Exec>(0, it * Kokkos::max(0, km - 4)),
        KOKKOS_LAMBDA(int index) {
          const int i = index % it, k = index / it + 2;
          if (kord != 4)
            a4(i, k, 3) = three * (two * a4(i, k, 0) - (a4(i, k, 1) + a4(i, k, 2)));
          kmppm(dc(i, k),                               // In
                Kokkos::subview(a4, i, k, Kokkos::ALL), // InOut
                lmt                                     // In
          );
        });
  }
  Kokkos::parallel_for(
      "ppm2m::bottom_limiter", Kokkos::RangePolicy<Exec>(0, it * 2), KOKKOS_LAMBDA(int index) {
        const int i = index % it, k = index / it + km1 - 1;
        a4(i, k, 3) = three * (two * a4(i, k, 0) - (a4(i, k, 1) + a4(i, k, 2)));
        kmppm(dc(i, k), Kokkos::subview(a4, i, k, Kokkos::ALL), 0);
      });
}

//=============================================================================
template <class SourcePressure, class Source, class TargetPressure, class Target>
void map1_ppm(ScratchWorkspace& workspace,                                       // InOut
              int km, SourcePressure pe1, Source q1, int kn, TargetPressure pe2, // In
              Target q2,                                                         // InOut
              int ng_s, int ng_n, int itot, int i1, int i2,                      // In
              int j, int jfirst, int jlast, int iv, int kord                     // In
) {
  // Perform PPM on one latitude. Pressure edges increase from top to surface.
  // q1/q2 retain (longitude,latitude,level), with latitude's lower bound
  // jfirst-ng_s mapped to zero. Work arrays map i1:i2 onto local columns.
  // Arguments (Native order):
  // km: Original vertical dimension
  // pe1: pressure at layer edges (from model top to bottom surface) in the original vertical
  // coordinate
  // q1: Field input
  // kn: Target vertical dimension
  // pe2: pressure at layer edges (from model top to bottom surface) in the new vertical coordinate
  // q2: Field output
  // ng_s: Ghosted latitudes south
  // ng_n: Ghosted latitudes north
  // itot: Total latitudes
  // i1: Starting longitude
  // i2: Finishing longitude
  // j: Current latitude
  // jfirst: Starting latitude
  // jlast: Finishing latitude
  // iv: Mode: 0 ==  constituents  1 == ???
  // kord: Method order

  auto workspace_scope = workspace.scope();
  using Space = typename Source::memory_space;
  using Exec = typename Source::execution_space;
  constexpr real_t r3 = one / three, r23 = two / three;
  const int columns = i2 - i1 + 1, latitude = j - (jfirst - ng_s);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Space> dp1 =
      workspace.view<Kokkos::View<real_t**, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "map1_ppm::dp1", columns, km);
  Kokkos::View<real_t***, Kokkos::LayoutLeft, Space> q4 =
      workspace.view<Kokkos::View<real_t***, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "map1_ppm::q4", columns, km, 4);
  // k0 stores zero-based source intervals; slot zero seeds the preceding search.
  Kokkos::View<int**, Kokkos::LayoutLeft, Space> k0 =
      workspace.view<Kokkos::View<int**, Kokkos::LayoutLeft, Space>>(
          Kokkos::WithoutInitializing, "map1_ppm::k0", columns, kn + 2);
  Kokkos::View<int*, Kokkos::LayoutLeft, Space> failure_k =
      workspace.view<Kokkos::View<int*, Kokkos::LayoutLeft, Space>>(Kokkos::WithoutInitializing,
                                                                    "map1_ppm::failure_k", columns);
  Kokkos::parallel_for(
      "map1_ppm::input", Kokkos::RangePolicy<Exec>(0, columns * km), KOKKOS_LAMBDA(int index) {
        const int column = index % columns, i = column + i1, k = index / columns;
        dp1(column, k) = pe1(i, k + 1) - pe1(i, k);
        q4(column, k, 0) = q1(i, latitude, k);
      });
  // Compute vertical subgrid distribution.
  ppm2m(workspace,
        q4,                       // InOut
        dp1, km, i1, i2, iv, kord // In
  );

  // For each pe2(i,k), find the smallest k0 such that
  // pe1(i,k0)<=pe2(i,k)<=pe1(i,k0+1). Edges are assumed monotone and share
  // endpoints. Retain the ascending search from the preceding k0.
  Kokkos::parallel_for(
      "map1_ppm::search", Kokkos::RangePolicy<Exec>(0, columns), KOKKOS_LAMBDA(int column) {
        const int i       = column + i1;
        failure_k(column) = -1;
        k0(column, 0)     = 0;
        for (int k = 0; k <= kn; ++k) {
          int k0found = -1;
          for (int kk = k0(column, k); kk < km; ++kk) {
            if (pe2(i, k) <= pe1(i, kk + 1)) {
              k0(column, k + 1) = kk;
              k0found = kk;
              break;
            }
          }
          if (k0found < 0) {
            failure_k(column) = k;
            return;
          }
        }
      });
  const auto failure_host = to_host_space(workspace, failure_k);
  for (int column = 0; column < columns; ++column) {
    if (failure_host(column) >= 0) {
      const int i = column + i1, k = failure_host(column);
      const auto pe1_host = to_host_space(workspace, pe1), pe2_host = to_host_space(workspace, pe2);
      std::ostringstream message;
      message << std::setprecision(std::numeric_limits<real_t>::max_digits10)
              << "mapz error - k0found i j k (kk,pe1,pe2) = " << -1 << ' ' << i + 1 << ' ' << j + 1
              << ' ' << k + 1 << ' ';
      // Keep Native diagnostic's km+1 range, including its undefined pe2 read
      // if kn<km on this failed-search path. Do not clamp or invent a fallback.
      for (int kk = 0; kk <= km; ++kk)
        message << kk + 1 << ' ' << pe1_host(i, kk) << ' ' << pe2_host(i, kk) << ' ';
      std::cerr << message.str() << '\n';
      // Native returns before any output interpolation. q2 is unchanged.
      return;
    }
  }

  // Interpolate. Each target cell's whole-source-cell sum remains ascending;
  // no parallel reduction or cumulative-sum reassociation is substituted.
  Kokkos::parallel_for(
      "map1_ppm::interpolate", Kokkos::RangePolicy<Exec>(0, columns * kn),
      KOKKOS_LAMBDA(int index) {
        const int column = index % columns, i = column + i1, k = index / columns;
        // Contribution from whole source cells strictly between the target edges.
        real_t qsumk = zero;
        for (int kl = k0(column, k + 1) + 1; kl < k0(column, k + 2); ++kl)
          qsumk = qsumk + dp1(column, kl) * q4(column, kl, 0);
        const int kk = k0(column, k + 1);
        // Contribution between pe1(i,kk) and pe2(i,k).
        const real_t pl = (pe2(i, k) - pe1(i, kk)) / dp1(column, kk);
        // Are both target edges in the same source interval?
        if (k0(column, k + 2) == k0(column, k + 1)) {
          const real_t pr = (pe2(i, k + 1) - pe1(i, kk)) / dp1(column, kk);
          q2(i, latitude, k) =
              q4(column, kk, 1) +
              one_half * (q4(column, kk, 3) + q4(column, kk, 2) - q4(column, kk, 1)) * (pr + pl) -
              q4(column, kk, 3) * r3 * (pr * (pr + pl) + pl * pl);
        } else {
          // Contribution between pe2(i,k) and pe1(i,kk+1).
          real_t qsum =
              (pe1(i, kk + 1) - pe2(i, k)) *
              (q4(column, kk, 1) +
               one_half * (q4(column, kk, 3) + q4(column, kk, 2) - q4(column, kk, 1)) * (one + pl) -
               q4(column, kk, 3) * (r3 * (one + pl * (one + pl))));
          // Contribution across the intervening whole source cells.
          qsum = qsum + qsumk;
          // Contribution from the last source interval to the upper target edge.
          const int kl      = k0(column, k + 2);
          const real_t delp = pe2(i, k + 1) - pe1(i, kl);
          const real_t esl  = delp / dp1(column, kl);
          qsum = qsum + delp * (q4(column, kl, 1) + one_half * esl *
                                                        (q4(column, kl, 2) - q4(column, kl, 1) +
                                                         q4(column, kl, 3) * (one - r23 * esl)));
          q2(i, latitude, k) = qsum / (pe2(i, k + 1) - pe2(i, k));
        }
      });
}

//=============================================================================
void remap_vals_ppm(ScratchWorkspace& workspace,                  // InOut
                    int ngrdcol,                                  // In
                    int levels_source_idx, int levels_target_idx, // In
                    Real2D levels_source, Real2D levels_target,   // In
                    Real2D source_values, int iv,                 // In
                    Real2D remap_vals_ppm                         // Out
) {
  // E3SM PPM wrapper accepting CLUBB pressure levels from surface to top.
  // Arguments (Native order):
  // ngrdcol: Number of grid columns []
  // levels_source_idx: Number of levels in the source grid [].
  // levels_target_idx: Number of levels in the target grid [].
  // levels_source: the pressure levels of the source grid [kg/(ms^2)]
  // levels_target: the pressure levels of the target grid [kg/(ms^2)]
  // source_values: given values on the source grid that should be remapped to the target grid
  // iv: setting for E3SM's PPM method, that sets what type of value is remapped, the following are
  // the three valid settings: iv =-1: winds iv = 0: positive definite scalars iv = 1: others

  auto workspace_scope = workspace.scope();

  // Explicit singleton dimensions replace Native scalar sequence association
  // with map1_ppm's rank2 edges/rank3 latitude fields. Output scratch is reused
  // between columns, preserving the Native failed-search return semantics.
  Real2D source_flipped = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_ppm::source_flipped", 1, levels_source_idx);
  Real2D target_flipped = workspace.view<Real2D>(
      Kokkos::WithoutInitializing, "remap_vals_ppm::target_flipped", 1, levels_target_idx);
  Real3D source_values_flipped =
      workspace.view<Real3D>(Kokkos::WithoutInitializing, "remap_vals_ppm::source_values_flipped",
                             1, 1, levels_source_idx - 1);
  Real3D target_values_flipped =
      workspace.view<Real3D>(Kokkos::WithoutInitializing, "remap_vals_ppm::target_values_flipped",
                             1, 1, levels_target_idx - 1);
  const int ghost_lat_south = 0, ghost_lat_north = 0, total_lat = 1;
  const int start_long = 0, finish_long = 0, current_lat = 0, start_lat = 0, finish_lat = 0,
            kord = 4;
  for (int i = 0; i < ngrdcol; ++i) {
    // E3SM uses top-to-surface edges and values, opposite to CLUBB's grids.
    Kokkos::parallel_for(
        "remap_vals_ppm::input", Kokkos::max(levels_source_idx, levels_target_idx),
        KOKKOS_LAMBDA(int k) {
          if (k < levels_source_idx) {
            source_flipped(0, k) = levels_source(i, levels_source_idx - k - 1);
            if (k < levels_source_idx - 1)
              source_values_flipped(0, 0, k) = source_values(i, levels_source_idx - k - 2);
          }
          if (k < levels_target_idx)
            target_flipped(0, k) = levels_target(i, levels_target_idx - k - 1);
        });
    map1_ppm(workspace, levels_source_idx - 1, source_flipped, source_values_flipped, // In
             levels_target_idx - 1, target_flipped,                                   // In
             target_values_flipped,                                                   // InOut
             ghost_lat_south, ghost_lat_north, total_lat, start_long, finish_long,    // In
             current_lat, start_lat, finish_lat, iv, kord                             // In
    );
    // Restore surface-to-top ordering.
    Kokkos::parallel_for(
        "remap_vals_ppm::output", levels_target_idx - 1, KOKKOS_LAMBDA(int k) {
          remap_vals_ppm(i, k) = target_values_flipped(0, 0, levels_target_idx - k - 2);
        });
  }
}
} // namespace
namespace {
template <class Source> bool remap_overlaps(const Source& source, Real2D output) {
  // Compare complete same-space spans before replacing an owned-return snapshot.
  if (!source.size() || !output.size())
    return false;
  const auto begin       = reinterpret_cast<std::uintptr_t>(source.data());
  const auto destination = reinterpret_cast<std::uintptr_t>(output.data());
  return begin <= destination ? destination - begin < source.span() * sizeof(real_t)
                              : begin - destination < output.span() * sizeof(real_t);
}
struct OwnedRemapCompletion {
  // Legacy APIs destroy their local workspace on return, so complete its queued
  // use first, including exception exits. Case-owned workspace paths need no fence.
  ~OwnedRemapCompletion() { Kokkos::fence("remap_vals_to_target::owned_workspace"); }
};
} // namespace

void remap_vals_to_target(ScratchWorkspace& workspace,                  // InOut
                          int ngrdcol,                                  // In
                          const grid& gr_source, const grid& gr_target, // In
                          int source_values_idx,                        // In
                          Real2D source_values,                         // In
                          int target_values_idx,                        // In
                          int total_idx_rho_lin_spline,                 // In
                          Real2D rho_lin_spline_vals,                   // In
                          Real2D rho_lin_spline_levels,                 // In
                          int iv, Real1D p_sfc,                         // In
                          int grid_remap_method,                        // In
                          bool l_zt_variable,                           // In
                          Real2D target_values                          // Out
) {
  // Fill an existing target field; preserve source snapshots when any read input overlaps it.
  auto workspace_scope = workspace.scope();
  if (remap_overlaps(source_values, target_values) ||
      remap_overlaps(rho_lin_spline_vals, target_values) ||
      remap_overlaps(rho_lin_spline_levels, target_values) ||
      remap_overlaps(p_sfc, target_values) || remap_overlaps(gr_source.zm, target_values) ||
      remap_overlaps(gr_source.zt, target_values) || remap_overlaps(gr_target.zm, target_values) ||
      remap_overlaps(gr_target.zt, target_values)) {
    auto result = workspace.view<Real2D>(
        Kokkos::WithoutInitializing, "remap_vals_to_target::overlap", ngrdcol, target_values_idx);
    remap_vals_to_target_values(workspace, ngrdcol, gr_source, gr_target, source_values_idx,
                                source_values, target_values_idx, total_idx_rho_lin_spline,
                                rho_lin_spline_vals, rho_lin_spline_levels, iv, p_sfc,
                                grid_remap_method, l_zt_variable, result);
    Kokkos::deep_copy(target_values, result);
  } else {
    remap_vals_to_target_values(workspace, ngrdcol, gr_source, gr_target, source_values_idx,
                                source_values, target_values_idx, total_idx_rho_lin_spline,
                                rho_lin_spline_vals, rho_lin_spline_levels, iv, p_sfc,
                                grid_remap_method, l_zt_variable, target_values);
  }
}

Real2D remap_vals_to_target(ScratchWorkspace& workspace,                  // InOut
                            int ngrdcol,                                  // In
                            const grid& gr_source, const grid& gr_target, // In
                            int source_values_idx,                        // In
                            Real2D source_values,                         // In
                            int target_values_idx,                        // In
                            int total_idx_rho_lin_spline,                 // In
                            Real2D rho_lin_spline_vals,                   // In
                            Real2D rho_lin_spline_levels,                 // In
                            int iv, Real1D p_sfc,                         // In
                            int grid_remap_method,                        // In
                            bool l_zt_variable                            // In
) {
  // Reserve the returned target in the active caller frame, beyond this helper's temporary scopes.
  auto result = workspace.view<Real2D>(Kokkos::WithoutInitializing, "remap_vals_to_target", ngrdcol,
                                       target_values_idx);
  remap_vals_to_target(workspace, ngrdcol, gr_source, gr_target, source_values_idx, source_values,
                       target_values_idx, total_idx_rho_lin_spline, rho_lin_spline_vals,
                       rho_lin_spline_levels, iv, p_sfc, grid_remap_method, l_zt_variable, result);
  return result;
}

Real2D remap_vals_to_target(int ngrdcol,                                  // In
                            const grid& gr_source, const grid& gr_target, // In
                            int source_values_idx,                        // In
                            Real2D source_values,                         // In
                            int target_values_idx,                        // In
                            int total_idx_rho_lin_spline,                 // In
                            Real2D rho_lin_spline_vals,                   // In
                            Real2D rho_lin_spline_levels,                 // In
                            int iv, Real1D p_sfc,                         // In
                            int grid_remap_method,                        // In
                            bool l_zt_variable                            // In
) {
  // Retain the legacy independently owned target; temporary storage ends at this return boundary.
  Real2D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "remap_vals_to_target"), ngrdcol,
                target_values_idx);
  ScratchWorkspace workspace;
  auto workspace_scope = workspace.scope();
  OwnedRemapCompletion completion;
  remap_vals_to_target(workspace, ngrdcol, gr_source, gr_target, source_values_idx, source_values,
                       target_values_idx, total_idx_rho_lin_spline, rho_lin_spline_vals,
                       rho_lin_spline_levels, iv, p_sfc, grid_remap_method, l_zt_variable, result);
  return result;
}

} // namespace clubb_kokkos
