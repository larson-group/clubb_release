// Native-order hole filling and hydrometeor clipping.
// See fill_holes.hpp for public entries, layouts and supported vertical methods.
#include "scratch_workspace.hpp"
#include "error_code.hpp"
#include "fill_holes.hpp"
#include "model_flags.hpp"
#include "stats_netcdf.hpp"
#include "index_mapping.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
namespace clubb_kokkos {
namespace {

// Shared host rendering of the Native one-level warnings. It changes no fields.
template <class Input, class Output>
void print_hm_fill_diagnostics(int num_hm_fill, Input hm_one_lev, Output hm_one_lev_filled) {
  std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  real_t total_hole = zero, total_mass = zero;
  for (int i = 0; i < num_hm_fill; ++i) {
    if (hm_one_lev(i) < zero)
      total_hole = total_hole + hm_one_lev(i);
    else
      total_mass = total_mass + hm_one_lev(i);
  }
  if (Kokkos::abs(total_mass) < eps) {
    std::cout << "Warning: One-level hole filling was not successful! total_mass ~= 0\n";
    return;
  }
  if (Kokkos::abs(total_hole) > total_mass) {
    for (int i = 0; i < num_hm_fill; ++i)
      std::cout << "Warning: One-level hole filling was not able to fill holes completely!"
                   " The holes were filled partially. |total_hole| > total_mass\n";
  }
  real_t sum_before = zero, sum_after = zero;
  bool l_negative = false;
  for (int i = 0; i < num_hm_fill; ++i) {
    sum_before = sum_before + hm_one_lev(i);
    sum_after  = sum_after + hm_one_lev_filled(i);
    if (hm_one_lev_filled(i) < zero)
      l_negative = true;
  }
  if (Kokkos::abs(sum_before - sum_after) > Kokkos::abs(sum_before + sum_after) * eps / 2)
    std::cout << "Warning: Hole filling was not conservative!\n";
  if (l_negative)
    std::cout << "Warning: Hole filling failed! A hole could not be filled.\n";
}
constexpr bool l_debug = false;
constexpr bool l_print_before_after = false;
// Native ANY visits logical elements, including column windows whose
// LayoutLeft allocation retains the parent array's padded leading dimension.
template <class View> bool any_below(View field, real_t threshold) {
  static_assert(View::rank >= 1 && View::rank <= 3);
  int below = 0;
  Kokkos::parallel_reduce(
      "fill_holes::any_below", field.size(),
      KOKKOS_LAMBDA(int cell, int &count) {
        // CUDA extended lambdas must capture the View before a constexpr branch.
        const int first_extent = field.extent_int(0);
        real_t value;
        if constexpr (View::rank == 1)
          value = field(cell);
        else if constexpr (View::rank == 2)
          value = field(((cell % first_extent)), ((cell / first_extent)));
        else
          value = field(((cell % first_extent)), (((cell / first_extent))) % field.extent_int(1),
                        cell / (first_extent * field.extent_int(1)));
        if (value < threshold)
          ++count;
      },
      below);
  return below != 0;
}
template <class View> void print_field(View field) {
  auto values = to_host_space(field);
  std::cerr << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t k = 0; k < values.extent(1); ++k)
    for (std::size_t i = 0; i < values.extent(0); ++i)
      std::cerr << values(i, k) << " ";
  std::cerr << "\n";
}
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_global(ScratchWorkspace& workspace, int nz, int ngrdcol, real_t threshold, // In
                       int lower_hf_level, int upper_hf_level,                             // In
                       Real2D dz, Real2D rho_ds, int grid_dir_indx,                        // In
                       Real2D field)                                                       // InOut
    ;

// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_sliding_window(ScratchWorkspace& workspace, int nz, int ngrdcol,
                               real_t threshold,                            // In
                               int lower_hf_level, int upper_hf_level,      // In
                               Real2D dz, Real2D rho_ds, int grid_dir_indx, // In
                               Real2D field)                                // InOut
    ;
} // namespace

//=============================================================================
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// fill_holes_type: Option for which type of hole filler to use
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_vertical_api(ScratchWorkspace& workspace, int nz, int ngrdcol,
                             real_t threshold,                            // In
                             int lower_hf_level, int upper_hf_level,      // In
                             Real2D dz, Real2D rho_ds, int grid_dir_indx, // In
                             int fill_holes_type,                         // In
                             Real2D field)                                // InOut
{
  // Description:
  //   This subroutine calls a hole filling method, specified by fill_holes_type.

  //   The lowest storage level (k=0) should not be included, as the hole-filling scheme
  //   should not alter the set value of 'field' at the surface (for momentum
  //   level variables), or consider the value of 'field' at a level below the
  //   surface (for thermodynamic level variables).
  //
  //   For momentum level variables only, the hole-filling scheme should not
  //   alter the set value of 'field' at the upper boundary level (k=nz-1).
  //   So for momemtum level variables, call with upper_hf_level=nz-2, and
  //   for thermodynamic level variables, call with upper_hf_level=nz-1.
  //
  //-----------------------------------------------------------------------
  auto workspace_scope = workspace.scope();
  Real2D field_initial;
  if (l_debug || l_print_before_after) {
    field_initial =
        workspace_scope.view<Real2D>("fill_holes_vertical_api::field_initial", ngrdcol, nz);
    Kokkos::deep_copy(field_initial, field);
    std::cerr << "field threshold: " << threshold << "\nfield: ";
    print_field(field);
    std::cerr << " -- over\n";
  }
  // Only bother with a fill call if there are values below threshold.
  bool l_field_below_threshold = any_below(field, threshold);
  if (l_field_below_threshold) {
    if (fill_holes_type == global_fill) {
      // Modify the entire range: maximally effective, computationally cheap, minimally local.
      fill_holes_global(workspace, nz, ngrdcol, threshold, lower_hf_level, upper_hf_level, dz,
                        rho_ds, grid_dir_indx, field);
    } else if (fill_holes_type == sliding_window) {
      // Consecutive overlapping windows are serial; global fill is the Native fallback.
      fill_holes_sliding_window(workspace, nz, ngrdcol, threshold, lower_hf_level, upper_hf_level,
                                dz, rho_ds, grid_dir_indx, field);
    } else {
      // Same supported scope as JAX: the four other Native methods are excluded.
      // Preserve the Native no-hole return before method dispatch/termination.
      std::cerr << "ERROR: unknown fill_holes_type: " << fill_holes_type << "\n";
      std::exit(EXIT_SUCCESS);
    }
  }
  if (l_debug || l_print_before_after) {
    auto values  = to_host_space(workspace, field);
    auto initial = to_host_space(workspace, field_initial);
    auto density = to_host_space(workspace, rho_ds);
    auto spacing = to_host_space(workspace, dz);
    if (any_below(field, threshold)) {
      std::cerr << "THERE ARE BELOW THRESHOLD VALUES IN THE FIELD AFTER THE FILL\n"
                   "EVEN THOUGH THE FIELD AVERAGE IS BELOW THRESHOLD. This is not\n"
                   "a gauranteed error, since smoothing contraints may prevent \n"
                   "some above threshold levels from having all their mass exhausted.\n";
      for (int i = 0; i < ngrdcol; ++i) {
        real_t numer = zero, denom = zero;
        for (int k = 0; k < nz; ++k) {
          numer += values(i, k) * density(i, k) * spacing(i, k);
          denom += density(i, k) * spacing(i, k);
        }
        if (numer / denom > threshold) {
          for (int k = 0; k < nz; ++k) {
            if (values(i, k) < threshold)
              std::cout << "field(" << (((i + 1))) << " " << (((k + 1))) << ") = " << values(i, k)
                        << " < " << threshold << "\n";
          }
          std::cerr << "column " << (((i + 1))) << " field average = " << numer / denom
                    << " -- threshold = " << threshold << "\n";
        } else {
          std::cerr << "FIELD AVERAGE BELOW THRESHOLD - FILLING INCOMPLETE AS EXPECTED\n";
        }
      }
    }
    real_t before = zero, after = zero;
    for (int k = 0; k < nz; ++k)
      for (int i = 0; i < ngrdcol; ++i) {
        before += initial(i, k) * density(i, k) * spacing(i, k);
        after += values(i, k) * density(i, k) * spacing(i, k);
      }
    std::cerr << "difference in mass (before-after): " << before - after << "\nfield after: ";
    print_field(field);
    std::cerr << " -- over\n";
  }
}

namespace {

//=============================================================================
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_global(ScratchWorkspace& workspace, int nz, int ngrdcol, real_t threshold, // In
                       int lower_hf_level, int upper_hf_level,                             // In
                       Real2D dz, Real2D rho_ds, int grid_dir_indx,                        // In
                       Real2D field)                                                       // InOut
{
  // Description:
  //   This subroutine clips values of 'field' that are below 'threshold' using
  //   the whole range [lower_hf_level:upper_hf_level] as the fill window. This
  //   maximized effectiveness, but minimized locality.
  //
  //   Mass is conserved by reducing the clipped field everywhere by a constant
  //   multiplicative coefficient.
  //
  //   This subroutine does not guarantee that the clipped field will exceed
  //   threshold everywhere; blunt clipping is needed for that.
  //-----------------------------------------------------------------------
  auto workspace_scope = workspace.scope();
  if (l_debug)
    std::cout << "fill_holes_type: GLOBAL FILL\n";
  Real2D rho_ds_dz     = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "fill_holes_global::rho_ds_dz", ngrdcol, nz);
  Real2D field_clipped = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "fill_holes_global::field_clipped", ngrdcol, nz);
  Real1D denom_integral_global = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "fill_holes_global::denom_integral_global", ngrdcol);
  Real1D numer_integral_global = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "fill_holes_global::numer_integral_global", ngrdcol);
  Real1D field_avg_global = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "fill_holes_global::field_avg_global", ngrdcol);
  Real1D mass_fraction_global = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "fill_holes_global::mass_fraction_global", ngrdcol);
  // Independent columns share a kernel; each vertical integral stays sequential in Native order.
  Kokkos::parallel_for(
      "fill_holes_global::columns", ngrdcol, KOKKOS_LAMBDA(int column) {
        const int i = column;
        for (int k = 0; k < nz; ++k)
          rho_ds_dz(i, k) = rho_ds(i, k) * dz(i, k);
        // Compute the numerator and denominator integrals.
        numer_integral_global(i) = real_t(0.0);
        denom_integral_global(i) = real_t(0.0);
        for (int k = lower_hf_level; grid_dir_indx * k - grid_dir_indx * upper_hf_level < 1;
             k += grid_dir_indx) {
          numer_integral_global(i) = numer_integral_global(i) + rho_ds_dz(i, k) * field(i, k);
          denom_integral_global(i) = denom_integral_global(i) + rho_ds_dz(i, k);
        }
        // Vertical average: see advance_helper_module::vertical_avg.
        field_avg_global(i) = numer_integral_global(i) / denom_integral_global(i);
        // Clip upward if all holes can be filled; otherwise remove all above-threshold mass.
        for (int k = lower_hf_level; grid_dir_indx * k - grid_dir_indx * upper_hf_level < 1;
             k += grid_dir_indx) {
          if (field_avg_global(i) >= threshold)
            field_clipped(i, k) = threshold > field(i, k) ? threshold : field(i, k);
          else
            field_clipped(i, k) = threshold < field(i, k) ? threshold : field(i, k);
        }
        // Only the numerator changes for the clipped vertical integral.
        numer_integral_global(i) = real_t(0.0);
        for (int k = lower_hf_level; grid_dir_indx * k - grid_dir_indx * upper_hf_level < 1;
             k += grid_dir_indx)
          numer_integral_global(i) =
              numer_integral_global(i) + rho_ds_dz(i, k) * field_clipped(i, k);
        bool l_column_below_threshold = false;
        for (int k = lower_hf_level; grid_dir_indx * k - grid_dir_indx * upper_hf_level < 1;
             k += grid_dir_indx)
          if (field(i, k) < threshold)
            l_column_below_threshold = true;
        // Do not update a column without holes inside the filling range.
        if (l_column_below_threshold) {
          const real_t field_clipped_avg = numer_integral_global(i) / denom_integral_global(i);
          // Avoid division by zero when the clipped average is approximately threshold.
          if (Kokkos::abs(field_clipped_avg - threshold) >
              Kokkos::abs(field_clipped_avg + threshold) * eps / 2) {
            // Multiplicative normalization conserves the original integrated mass.
            mass_fraction_global(i) =
                (field_avg_global(i) - threshold) / (field_clipped_avg - threshold);
            for (int k = lower_hf_level; grid_dir_indx * k - grid_dir_indx * upper_hf_level < 1;
                 k += grid_dir_indx)
              field(i, k) = threshold + mass_fraction_global(i) * (field_clipped(i, k) - threshold);
          }
        }
      });
}
} // namespace

namespace {

//=============================================================================
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// dz: Spacing between thermodynamic grid levels; centered over momentum grid levels OR Spcaing
// between momentum grid levels; centered over thermodynamic grid levels
// rho_ds: Dry, static density on thermodynamic or momentum levels    [kg/m^3]
// grid_dir_indx: Grid direction index (+1 ascending; -1 descending)
// field: The field (e.g. wp2) that contains holes [Units same as threshold]
void fill_holes_sliding_window(ScratchWorkspace& workspace, int nz, int ngrdcol,
                               real_t threshold,                            // In
                               int lower_hf_level, int upper_hf_level,      // In
                               Real2D dz, Real2D rho_ds, int grid_dir_indx, // In
                               Real2D field)                                // InOut
{
  // Description:
  //   This subroutine clips values of 'field' that are below 'threshold' as much
  //   as possible (i.e. "fills holes"), but conserves the total integrated mass
  //   of 'field'.  This prevents clipping from acting as a spurious source.
  //
  //   This performs a sliding window technique, modifying consecutive ranges of
  //   vertical levels in serial, this is computationally expensive, but highly local.
  //   This high locally has a tradeoff with effectiveness, and can often fail to fill
  //   all the holes, especially if there is more than ~5, as a result, this relies on
  //   the global fill if the first pass of the sliding window fails to fill all holes
  //
  //   Mass is conserved by reducing the clipped field everywhere by a constant
  //   multiplicative coefficient.
  //
  // References:
  //   ``Numerical Methods for Wave Equations in Geophysical Fluid
  //     Dynamics'', Durran (1999), p. 292.
  //-----------------------------------------------------------------------
  auto workspace_scope = workspace.scope();
  if (l_debug)
    std::cout << "fill_holes_type: SLIDING WINDOW FILL WITH GLOBAL FALLBACK\n";
  bool l_field_below_threshold = any_below(field, threshold);
  if (l_field_below_threshold) {
    Real2D rho_ds_dz = workspace_scope.view<Real2D>(
        Kokkos::WithoutInitializing, "fill_holes_sliding_window::rho_ds_dz", ngrdcol, nz);
    Real2D invrs_denom_integral = workspace_scope.view<Real2D>(
        Kokkos::WithoutInitializing, "fill_holes_sliding_window::invrs_denom_integral", ngrdcol,
        nz);
    Real2D field_clipped = workspace_scope.view<Real2D>(
        Kokkos::WithoutInitializing, "fill_holes_sliding_window::field_clipped", ngrdcol, nz);
    Kokkos::parallel_for(
        "fill_holes_sliding_window::rho_ds_dz", ((ngrdcol * nz)), KOKKOS_LAMBDA(int cell) {
          const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
          rho_ds_dz(i, k) = rho_ds(i, k) * dz(i, k);
        });
    const int start_indx = ((grid_dir_indx * num_hf_draw_points + lower_hf_level));
    const int stop_indx  = ((-grid_dir_indx * num_hf_draw_points + upper_hf_level));
    // The denominator does not change. Explicit source-order sums retain the Native
    // OpenMP-offload improvement:
    // https://github.com/larson-group/clubb/issues/1138#issuecomment-1974918151
    Kokkos::parallel_for(
        "fill_holes_sliding_window::denominator", ngrdcol, KOKKOS_LAMBDA(int column) {
          const int i = column;
          for (int k = start_indx; (grid_dir_indx * k - grid_dir_indx * stop_indx) < 1;
               k += grid_dir_indx) {
            real_t rho_k_sum = real_t(0.0);
            const int start_indx_j = ((-grid_dir_indx * num_hf_draw_points + k));
            const int stop_indx_j  = ((grid_dir_indx * num_hf_draw_points + k));
            for (int j = start_indx_j; (grid_dir_indx * j - grid_dir_indx * stop_indx_j) < 1;
                 j += grid_dir_indx)
              rho_k_sum = rho_k_sum + rho_ds_dz(i, j);
            invrs_denom_integral(i, k) = one / rho_k_sum;
          }
        });
    // Make one pass through each profile. Overlapping windows read prior updates;
    // the k loop must stay serial, including for the descending grid.
    Kokkos::parallel_for(
        "fill_holes_sliding_window::columns", ngrdcol, KOKKOS_LAMBDA(int column) {
          const int i = column;
          for (int k = start_indx; (grid_dir_indx * k - grid_dir_indx * stop_indx) < 1;
               k += grid_dir_indx) {
            const int k_start             = ((-grid_dir_indx * num_hf_draw_points + k));
            const int k_end               = ((grid_dir_indx * num_hf_draw_points + k));
            bool l_window_below_threshold = false;
            for (int j = k_start; (grid_dir_indx * j - grid_dir_indx * k_end) < 1;
                 j += grid_dir_indx)
              if (field(i, j) < threshold)
                l_window_below_threshold = true;
            if (l_window_below_threshold) {
              // Compute the original field average, which must be conserved.
              real_t field_avg = real_t(0.0);
              for (int j = k_start; (grid_dir_indx * j - grid_dir_indx * k_end) < 1;
                   j += grid_dir_indx)
                field_avg = field_avg + rho_ds_dz(i, j) * field(i, j);
              field_avg = field_avg * invrs_denom_integral(i, k);
              // Remove either the holes or the available above-threshold material.
              for (int j = k_start; (grid_dir_indx * j - grid_dir_indx * k_end) < 1;
                   j += grid_dir_indx) {
                if (field_avg >= threshold)
                  field_clipped(i, j) = threshold > field(i, j) ? threshold : field(i, j);
                else
                  field_clipped(i, j) = threshold < field(i, j) ? threshold : field(i, j);
              }
              // The clipped total mass is no smaller when filling all holes is possible.
              real_t field_clipped_avg = real_t(0.0);
              for (int j = k_start; (grid_dir_indx * j - grid_dir_indx * k_end) < 1;
                   j += grid_dir_indx)
                field_clipped_avg = field_clipped_avg + rho_ds_dz(i, j) * field_clipped(i, j);
              field_clipped_avg = field_clipped_avg * invrs_denom_integral(i, k);
              if (Kokkos::abs(field_clipped_avg - threshold) >
                  Kokkos::abs(field_clipped_avg + threshold) * eps / 2) {
                const real_t mass_fraction =
                    (field_avg - threshold) / (field_clipped_avg - threshold);
                for (int j = k_start; (grid_dir_indx * j - grid_dir_indx * k_end) < 1;
                     j += grid_dir_indx)
                  field(i, j) = threshold + mass_fraction * (field_clipped(i, j) - threshold);
              }
            }
          }
        });
  }
  // Check the entire field, as Native does, before the global fallback.
  l_field_below_threshold = any_below(field, threshold);
  if (l_field_below_threshold) {
    if (l_debug)
      std::cout << "first pass insufficient -- doing global fill\n";
    fill_holes_global(workspace, nz, ngrdcol, threshold, lower_hf_level, upper_hf_level, dz, rho_ds,
                      grid_dir_indx, field);
  }
}
} // namespace

//=============================================================================
// nz: Number of grid levels [-]
// ngrdcol: Number of columns [-]
// threshold: A threshold (e.g. w_tol*w_tol) below which field must not fall
// [Units vary; same as field]
// lower_hf_level: Lower zero-based slot of the hole-filling range      []
// upper_hf_level: Upper zero-based slot of the hole-filling range      []
// wp2: The field (here wp2) that contains holes [m^2/s^2]
// up2: Zonal TKE component supplying the hole [m^2/s^2]
// vp2: The fields (here up2, vp2) that will be used to fill holes [m^2/s^2]
void fill_holes_wp2_from_horz_tke(ScratchWorkspace& workspace, int nz, int ngrdcol,
                                  real_t threshold,                       // In
                                  int lower_hf_level, int upper_hf_level, // In
                                  Real2D wp2, Real2D up2, Real2D vp2)     // InOut
{
  // Description:
  //   This subroutine clips values of wp2 that are below 'threshold' as much
  //   as possible (i.e. "fills holes"), but conserves the turbulent kinetic energy
  //   (up2+vp2+wp2). This prevents clipping from acting as a spurious source.
  //
  //   Turbulent kinetic energy at each height level is conserved by reducing up2 and vp2
  //   by a multiplicative coefficient.
  //
  //   This subroutine does not guarantee that the clipped field will exceed
  //   threshold everywhere; blunt clipping is needed for that.
  //
  //   The lowest storage level (k=0) should not be included, as the hole-filling scheme
  //   should not alter the set value of 'field' at the surface (for momentum
  //   level variables), or consider the value of 'field' at a level below the
  //   surface (for thermodynamic level variables).
  //
  //   For momentum level variables only, the hole-filling scheme should not
  //   alter the set value of 'field' at the upper boundary level (k=nz-1).
  //   So for momemtum level variables, call with upper_hf_level=nz-2, and
  //   for thermodynamic level variables, call with upper_hf_level=nz-1.
  //-----------------------------------------------------------------------
  auto workspace_scope         = workspace.scope();
  bool l_field_below_threshold = any_below(wp2, threshold);
  if (!l_field_below_threshold)
    return;
  const bool l_clubb_at_least_debug_level_3 = clubb_at_least_debug_level_api(3);
  Real2D wp2_initial, up2_initial, vp2_initial;
  if (l_clubb_at_least_debug_level_3) {
    wp2_initial =
        workspace_scope.view<Real2D>("fill_holes_wp2_from_horz_tke::wp2_initial", ngrdcol, nz);
    up2_initial =
        workspace_scope.view<Real2D>("fill_holes_wp2_from_horz_tke::up2_initial", ngrdcol, nz);
    vp2_initial =
        workspace_scope.view<Real2D>("fill_holes_wp2_from_horz_tke::vp2_initial", ngrdcol, nz);
    Kokkos::deep_copy(wp2_initial, wp2);
    Kokkos::deep_copy(up2_initial, up2);
    Kokkos::deep_copy(vp2_initial, vp2);
  }
  // Each level transfers TKE independently; numeric lower..upper bounds are unchanged.
  Kokkos::parallel_for(
      "fill_holes_wp2_from_horz_tke::cells",
      ngrdcol * std::max(0, ((-lower_hf_level + upper_hf_level + 1))), KOKKOS_LAMBDA(int cell) {
        const int i = ((cell % ngrdcol)), k = ((lower_hf_level + (cell / ngrdcol)));
        real_t up2_old, vp2_old, wp2_old, tke_x2_old, missing_wp2, up2_avail, vp2_avail,
            up2_vp2_avail, ratio;
        tke_x2_old = up2(i, k) + vp2(i, k) + wp2(i, k);
        up2_old    = up2(i, k);
        vp2_old    = vp2(i, k);
        wp2_old    = wp2(i, k);
        // Check if we have a hole to fill at level k and
        // there is buffer TKE in up2 and/or vp2 available
        if (wp2(i, k) < threshold && (up2(i, k) > threshold || vp2(i, k) > threshold)) {
          // Size of the hole to fill
          missing_wp2 = threshold - wp2(i, k);
          // Calculate available TKE that can be moved to wp2
          // At least one of these is strictly > 0
          up2_avail     = ((up2(i, k) - threshold) > (zero) ? (up2(i, k) - threshold) : (zero));
          vp2_avail     = ((vp2(i, k) - threshold) > (zero) ? (vp2(i, k) - threshold) : (zero));
          up2_vp2_avail = up2_avail + vp2_avail;
          if (missing_wp2 >= up2_vp2_avail) {
            // Not enough TKE available to fill hole
            // Take everything up to threshold from up2 and vp2, if possible
            up2(i, k) = ((up2(i, k)) < (threshold) ? (up2(i, k)) : (threshold));
            vp2(i, k) = ((vp2(i, k)) < (threshold) ? (vp2(i, k)) : (threshold));
            // Add available TKE to wp2
            wp2(i, k) = wp2(i, k) + up2_vp2_avail;
          } else {
            // We have enough available TKE to fill hole
            // Fill
            wp2(i, k) = threshold;
            // Remove TKE from up2 and vp2
            if (Kokkos::abs(up2_avail) < std::numeric_limits<real_t>::epsilon() * 1000) {
              // No TKE available from up2
              // up2 stays the same
              // Take all from vp2
              vp2(i, k) = vp2(i, k) - missing_wp2;
            } else {
              if (Kokkos::abs(vp2_avail) < std::numeric_limits<real_t>::epsilon() * 1000) {
                // No TKE available from vp2
                // Take all from up2
                up2(i, k) = up2(i, k) - missing_wp2;
                // vp2 statys the same
              } else {
                // Calculate portion of up2/vp2 that we want to take away
                ratio = missing_wp2 / up2_vp2_avail;
                // TKE available in both, take according to ratio
                // Since we know that there is enough available TKE,
                // this reduction will not reduce either variables below threshold
                up2(i, k) = threshold + up2_avail * (one - ratio);
                vp2(i, k) = threshold + vp2_avail * (one - ratio);
              }
            }
          }
        }
      });
  if (l_clubb_at_least_debug_level_3) {
    std::cerr << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    // Mirror old/new values for Native ordered diagnostics; no I/O inside kernels.
    auto wp2_new = to_host_space(workspace, wp2), up2_new = to_host_space(workspace, up2),
         vp2_new    = to_host_space(workspace, vp2);
    auto wp2_before = to_host_space(workspace, wp2_initial),
         up2_before = to_host_space(workspace, up2_initial),
         vp2_before = to_host_space(workspace, vp2_initial);
    for (int i = 0; i < ngrdcol; ++i)
      for (int k = lower_hf_level; k < upper_hf_level + 1; ++k) {
        const real_t wp2_old = wp2_before(i, k), up2_old = up2_before(i, k),
                     vp2_old    = vp2_before(i, k);
        const real_t tke_x2_old = up2_old + vp2_old + wp2_old;
        real_t tke_diff;
        tke_diff = Kokkos::abs(tke_x2_old - (up2_new(i, k) + vp2_new(i, k) + wp2_new(i, k)));
        if (tke_diff > std::numeric_limits<real_t>::epsilon() * 1000) {
          std::cerr << "Warning! fill_holes_wp2_from_horz_tke is not conservative " << " "
                    << "at level " << " " << (((i + 1))) << " " << (((k + 1))) << " "
                    << ". Difference: " << " " << tke_diff << "\n";
        }
        if (wp2_new(i, k) > threshold) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! " << " "
                    << "Too much TKE added to wp2. old: " << " " << wp2_old << " "
                    << ", new: " << " " << wp2_new(i, k) << "\n";
        }
        if (wp2_new(i, k) < wp2_old) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! wp2 was reduced. " << " "
                    << "old: " << " " << wp2_old << " " << ", new: " << " " << wp2_new(i, k)
                    << "\n";
        }
        if (up2_new(i, k) > up2_old) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! up2 increased. old: " << " "
                    << up2_old << " " << ", new: " << " " << up2_new(i, k) << "\n";
        }
        if (up2_old <= threshold && up2_new(i, k) < up2_old) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! TKE was falsely " << " "
                    << "taken from up2. Threshold: " << " " << threshold << " " << ", old: " << " "
                    << up2_old << " " << ", new: " << " " << up2_new(i, k) << "\n";
        }
        if (vp2_new(i, k) > vp2_old) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! vp2 increased. old: " << " "
                    << vp2_old << " " << ", new: " << " " << vp2_new(i, k) << "\n";
        }
        if (vp2_old <= threshold && vp2_new(i, k) < vp2_old) {
          std::cerr << "Warning in fill_holes_wp2_from_horz_tke! TKE was falsely " << " "
                    << "taken from vp2. Threshold: " << " " << threshold << " " << ", old: " << " "
                    << vp2_old << " " << ", new: " << " " << vp2_new(i, k) << "\n";
        }
      }
  }
}

//=============================================================================
// num_hm_fill: number of hydrometeors involved
// hm_one_lev: Same-phase species mixing ratios at one level [kg/kg]
// hm_one_lev_filled: Same-phase mixing ratios after redistribution [kg/kg]
void hole_filling_hm_one_lev(int num_hm_fill, Real1D hm_one_lev, // In
                             Real1D hm_one_lev_filled)           // Out
{
  // Description:
  // Fills holes between same-phase (i.e. either liquid or frozen) hydrometeors for
  // one height level.
  //
  // Warning: Do not input hydrometeors of different phases, e.g. liquid and frozen.
  // Otherwise heat will not be conserved.
  //
  // References:
  //
  // None
  //-----------------------------------------------------------------------
  // The concrete entry keeps diagnostics on the host. A batched caller invokes
  // the identical device overload on its per-level sections without kernel launches.
  Kokkos::parallel_for(
      "hole_filling_hm_one_lev::level", 1, KOKKOS_LAMBDA(int) {
        auto input = Kokkos::subview(hm_one_lev, Kokkos::ALL);
        auto output = Kokkos::subview(hm_one_lev_filled, Kokkos::ALL);
        hole_filling_hm_one_lev<decltype(input), decltype(output)>(num_hm_fill, input, output);
      });
  if (clubb_at_least_debug_level_api(2)) {
    auto input = to_host_space(hm_one_lev);
    auto output = to_host_space(hm_one_lev_filled);
    print_hm_fill_diagnostics(num_hm_fill, input, output);
  }
}

//=============================================================================
// nzt: Number of thermodynamic levels [-]
// hydromet_dim: Number of hydrometeor species [-]
// hydromet: Hydrometeor fields (level,species) [units vary]
// l_frozen_hm: if true, then the hydrometeor is frozen; otherwise liquid
// l_mix_rat_hm: if true, then the quantity is a hydrometeor mixing ratio
// hydromet_filled: Hydrometeor fields after frozen-species redistribution [units vary]
void fill_holes_hydromet_api(ScratchWorkspace& workspace, int nzt, int hydromet_dim,
                             Real2D hydromet,                         // In
                             Bool1D l_frozen_hm, Bool1D l_mix_rat_hm, // In
                             Real2D hydromet_filled)                  // Out
{
  // Description:
  // Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  // The hole filling conserves water substance between all same-phase (frozen or liquid)
  // hydrometeors at each height level.
  //
  // Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  //
  // Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  //
  // References:
  //
  // None
  //-----------------------------------------------------------------------
  auto workspace_scope = workspace.scope();
  auto frozen_flags    = to_host_space(workspace, l_frozen_hm);
  auto mixing_flags    = to_host_space(workspace, l_mix_rat_hm);
  int num_frozen_hm = 0;
  // Determine the number of frozen mixing ratios, preserving species order.
  for (int i = 0; i < hydromet_dim; ++i)
    if (frozen_flags(i) && mixing_flags(i))
      num_frozen_hm = ((num_frozen_hm + 1));
  Real2D hydromet_frozen = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "fill_holes_hydromet_api::hydromet_frozen", nzt, num_frozen_hm);
  Real2D hydromet_frozen_filled = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "fill_holes_hydromet_api::hydromet_frozen_filled", nzt,
      num_frozen_hm);
  // Pack the exact Native frozen subset; zero padding would alter sum association.
  int j = 0;
  for (int i = 0; i < hydromet_dim; ++i) {
    if (frozen_flags(i) && mixing_flags(i)) {
      Kokkos::deep_copy(Kokkos::subview(hydromet_frozen, Kokkos::ALL, j),
                        Kokkos::subview(hydromet, Kokkos::ALL, i));
      j = ((j + 1));
    }
  }
  // Same-phase redistribution is sequential in species, independent in height.
  Kokkos::parallel_for(
      "fill_holes_hydromet_api::levels", nzt, KOKKOS_LAMBDA(int level) {
        bool l_holes = false;
        for (int j = 0; j < num_frozen_hm; ++j)
          if (hydromet_frozen(level, j) < zero)
            l_holes = true;
        auto input = Kokkos::subview(hydromet_frozen, level, Kokkos::ALL);
        auto output = Kokkos::subview(hydromet_frozen_filled, level, Kokkos::ALL);
        if (l_holes) {
          hole_filling_hm_one_lev<decltype(input), decltype(output)>(num_frozen_hm, input, output);
        } else {
          for (int j = 0; j < num_frozen_hm; ++j)
            output(j) = input(j);
        }
      });
  if (clubb_at_least_debug_level_api(2)) {
    auto input  = to_host_space(workspace, hydromet_frozen);
    auto output = to_host_space(workspace, hydromet_frozen_filled);
    for (int i = 0; i < nzt; ++i) {
      bool l_holes = false;
      for (int j = 0; j < num_frozen_hm; ++j)
        if (input(i, j) < zero)
          l_holes = true;
      if (l_holes)
        print_hm_fill_diagnostics(num_frozen_hm, Kokkos::subview(input, i, Kokkos::ALL),
                                  Kokkos::subview(output, i, Kokkos::ALL));
    }
  }
  // Reassemble the filled array; liquid species and concentrations are unchanged.
  j = 0;
  for (int i = 0; i < hydromet_dim; ++i) {
    if (frozen_flags(i) && mixing_flags(i)) {
      Kokkos::deep_copy(Kokkos::subview(hydromet_filled, Kokkos::ALL, i),
                        Kokkos::subview(hydromet_frozen_filled, Kokkos::ALL, j));
      j = ((j + 1));
    } else {
      Kokkos::deep_copy(Kokkos::subview(hydromet_filled, Kokkos::ALL, i),
                        Kokkos::subview(hydromet, Kokkos::ALL, i));
    }
  }
  // Native liquid-phase same-level hole filling remains unimplemented.
}

//=============================================================================
// nzt: Number of thermodynamic levels [-]
// dt: Timestep         [s]
// exner: Exner function                            [-]
// hydromet_name: Fixed-width species name [-]
// rvm_mc: Water-vapor tendency [kg/kg/s]
// thlm_mc: Liquid-water potential-temperature tendency [K/s]
// hydromet: Hydrometeor array                         [units vary]
void fill_holes_wv(int nzt, real_t dt, Real1D exner, const std::string &hydromet_name, // In
                   Real1D rvm_mc, Real1D thlm_mc, Real1D hydromet)                     // InOut
{
  // Description:
  // Fills holes using the cloud water mixing ratio from the current height level.
  //
  // References:
  //
  // None
  //-----------------------------------------------------------------------
  // An unknown species terminates only upon reaching an actual hole, after the
  // Native vapor-tendency update. Keep this exceptional path ordered on the host.
  std::string name = hydromet_name;
  while (!name.empty() && name.back() == ' ')
    name.pop_back();
  const bool l_rain = name == "rrm";
  const bool l_frozen = name == "rim" || name == "rsm" || name == "rgm";
  if (!l_rain && !l_frozen) {
    auto hm = to_host_space(hydromet);
    for (int k = 0; k < nzt; ++k) {
      if (hm(k) < zero_threshold) {
        auto vapor = to_host_space(rvm_mc);
        const real_t rvm_clip_tndcy = hm(k) / dt;
        vapor(k)                    = vapor(k) + rvm_clip_tndcy;
        from_execution_space(rvm_mc, to_execution_space(vapor));
        std::cerr << "Fatal error in microphys_driver\n";
        std::exit(EXIT_FAILURE); // Native ERROR STOP, after the vapor update.
      }
    }
    return;
  }
  Kokkos::parallel_for(
      "fill_holes_wv::levels", nzt, KOKKOS_LAMBDA(int level) {
        const int k = level;
        if (hydromet(k) < zero_threshold) {
          // Tendency transferred from the hydrometeor to vapor.
          const real_t rvm_clip_tndcy = hydromet(k) / dt;
          rvm_mc(k)                   = rvm_mc(k) + rvm_clip_tndcy;
          // Evaporation and sublimation have different latent heat compensation.
          if (l_rain)
            thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * (Lv / (Cp * exner(k)));
          else
            thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * (Ls / (Cp * exner(k)));
          hydromet(k) = zero_threshold;
        }
      });
}

//=============================================================================
// gr: Native grid coordinates and spacings [-,m]
// ngrdcol: Number of model columns
// nzt: Number of vertical levels
// dt: Timestep         [s]
// hydromet_dim: Number of hydrometeor fields
// hm_metadata: Hydrometeor species indices, flags and tolerances [-,units vary]
// l_fill_holes_hm: Enable same-level frozen-species redistribution [-]
// rho_ds_zt: Dry, static density on thermo. levels    [kg/m^3]
// exner: Exner function                           [-]
// stats: Native statistics context, mutated by sampled budget calls [-]
// thlm_mc: Microphysics contributions to liquid potential temp [K/s]
// rvm_mc: Microphysics contributions to vapor water           [kg/kg/s]
// hydromet: Mean of hydrometeor fields    [units vary]
void fill_holes_driver_api(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, int nzt,
                           real_t dt, int hydromet_dim,                               // In
                           const hm_metadata_type& hm_metadata, bool l_fill_holes_hm, // In
                           Real2D rho_ds_zt, Real2D exner, int fill_holes_type,       // In
                           stats_type& stats, Real2D thlm_mc, Real2D rvm_mc,
                           Real3D hydromet) // InOut
{
  // Description:
  // Apply hydrometeor hole filling and clipping to all model columns.
  // Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  // The hole filling conserves water substance between all same-phase (frozen or liquid)
  // hydrometeors at each height level.
  //
  // Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  //
  // Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  //
  // References:
  //
  // None
  //-----------------------------------------------------------------------
  auto workspace_scope    = workspace.scope();
  Real3D hydromet_filled  = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing,
                                                         "fill_holes_driver_api::hydromet_filled",
                                                         ngrdcol, nzt, hydromet_dim);
  Real3D hydromet_clipped = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing,
                                                         "fill_holes_driver_api::hydromet_clipped",
                                                         ngrdcol, nzt, hydromet_dim);
  real_t max_velocity;
  std::string name_bt, name_hf, name_wvhf, name_cl, name_mc, hydromet_name;
  constexpr bool l_hole_fill = true;
  auto hydromet_list_host    = to_host_space(workspace, hm_metadata.hydromet_list);
  auto l_mix_rat_hm_host     = to_host_space(workspace, hm_metadata.l_mix_rat_hm);
  auto l_frozen_hm_host      = to_host_space(workspace, hm_metadata.l_frozen_hm);
  auto hydromet_tol_host     = to_host_space(workspace, hm_metadata.hydromet_tol);
  // Stats payloads are staged only at the existing host stats API boundary.
  auto budget = [&](const std::string& name, int ihm, bool begin) -> void {
    auto budget_scope = workspace.scope();
    auto values = budget_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                            "fill_holes_driver_api::stats_values", ngrdcol, nzt);
    Kokkos::parallel_for(
        "fill_holes_driver_api::stats_values", ((ngrdcol * nzt)), KOKKOS_LAMBDA(int cell) {
          const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
          values(i, k) = hydromet(i, k, ihm) / dt;
        });
    if (begin)
      stats_begin_budget(name, values, stats);
    else
      stats_finalize_budget(name, values, stats);
  };
  // Start _hf budgets before any same-phase or vertical filling.
  if (stats.l_sample) {
    for (int ihm = 0; ihm < hydromet_dim; ++ihm) {
      setup_stats_names(ihm, hydromet_dim, hydromet_list_host, max_velocity, &name_bt, &name_hf,
                        &name_wvhf, &name_cl, &name_mc);
      budget(name_hf, ihm, true);
    }
  }
  if (any_below(hydromet, zero_threshold) && l_fill_holes_hm) {
    // Preserve single-column public APIs at their actual boundaries. Packing a
    // strided row reuses the established array owner; no Fortran scalar call.
    for (int i = 0; i < ngrdcol; ++i) {
      auto input =
          to_execution_space(workspace, Kokkos::subview(hydromet, i, Kokkos::ALL, Kokkos::ALL));
      auto output = to_execution_space(
          workspace, Kokkos::subview(hydromet_filled, i, Kokkos::ALL, Kokkos::ALL), false);
      fill_holes_hydromet_api(workspace, nzt, hydromet_dim, input, hm_metadata.l_frozen_hm,
                              hm_metadata.l_mix_rat_hm, output);
      from_execution_space(workspace, Kokkos::subview(hydromet_filled, i, Kokkos::ALL, Kokkos::ALL),
                           output);
    }
    Kokkos::deep_copy(hydromet, hydromet_filled);
  }
  for (int ihm = 0; ihm < hydromet_dim; ++ihm) {
    setup_stats_names(ihm, hydromet_dim, hydromet_list_host, max_velocity, &name_bt, &name_hf,
                      &name_wvhf, &name_cl, &name_mc);
    hydromet_name.assign(hydromet_list_host.data() + ((hydromet_list_host.extent(0) * ihm)),
                         hydromet_list_host.extent(0));
    while (!hydromet_name.empty() && hydromet_name.back() == ' ')
      hydromet_name.pop_back();
    if (clubb_at_least_debug_level_api(1)) {
      auto values = to_host_space(workspace, hydromet);
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          if (values(i, k, ihm) < zero_threshold)
            std::cerr << hydromet_name << " < " << zero_threshold
                      << " in fill_holes_driver_api at i, k= " << (((i + 1))) << " " << (((k + 1)))
                      << "\n";
    }
    auto species =
        to_execution_space(workspace, Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm));
    if (any_below(species, zero_threshold) && !hydromet_name.empty() && hydromet_name[0] == 'r' &&
        l_hole_fill) {
      for (int i = 0; i < ngrdcol; ++i) {
        auto column = to_execution_space(
            workspace, Kokkos::subview(hydromet, std::make_pair(i, (((i + 1)))), Kokkos::ALL, ihm));
        if (any_below(column, zero_threshold)) {
          auto dz = to_execution_space(
              workspace, Kokkos::subview(gr.dzt, std::make_pair(i, (((i + 1)))), Kokkos::ALL));
          auto rho = to_execution_space(
              workspace, Kokkos::subview(rho_ds_zt, std::make_pair(i, (((i + 1)))), Kokkos::ALL));
          // Native deliberately passes1..nzt and positive direction for these rows.
          fill_holes_vertical_api(workspace, nzt, 1, zero_threshold, 0, (nzt - 1), dz, rho, 1,
                                  fill_holes_type, column);
          from_execution_space(
              workspace,
              Kokkos::subview(hydromet, std::make_pair(i, (((i + 1)))), Kokkos::ALL, ihm), column);
        }
      }
    }
    if (stats.l_sample) {
      budget(name_hf, ihm, false);
      budget(name_wvhf, ihm, true);
    }
    if (!hydromet_name.empty() && hydromet_name[0] == 'r' && l_hole_fill) {
      for (int i = 0; i < ngrdcol; ++i) {
        auto hm = to_execution_space(workspace, Kokkos::subview(hydromet, i, Kokkos::ALL, ihm));
        if (any_below(hm, zero_threshold)) {
          // Residual holes can require vapor, notably ASEX A209 with LHS sampling.
          auto exner_column = to_execution_space(workspace, Kokkos::subview(exner, i, Kokkos::ALL));
          auto vapor = to_execution_space(workspace, Kokkos::subview(rvm_mc, i, Kokkos::ALL));
          auto temperature =
              to_execution_space(workspace, Kokkos::subview(thlm_mc, i, Kokkos::ALL));
          fill_holes_wv(nzt, dt, exner_column, hydromet_name, vapor, temperature, hm);
          from_execution_space(workspace, Kokkos::subview(rvm_mc, i, Kokkos::ALL), vapor);
          from_execution_space(workspace, Kokkos::subview(thlm_mc, i, Kokkos::ALL), temperature);
          from_execution_space(workspace, Kokkos::subview(hydromet, i, Kokkos::ALL, ihm), hm);
        }
      }
    }
    if (stats.l_sample)
      budget(name_wvhf, ihm, false);
    if (l_mix_rat_hm_host(ihm)) {
      if (stats.l_sample)
        budget(name_cl, ihm, true);
      const bool l_frozen       = l_frozen_hm_host(ihm);
      const real_t hydromet_tol = hydromet_tol_host(ihm);
      Kokkos::parallel_for(
          "fill_holes_driver_api::mixing_ratio_clip", ((ngrdcol * nzt)), KOKKOS_LAMBDA(int cell) {
            const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
            // Clip negative values, then transfer very small positive mass to vapor.
            if (hydromet(i, k, ihm) < zero_threshold)
              hydromet(i, k, ihm) = zero_threshold;
            if (hydromet(i, k, ihm) <= hydromet_tol) {
              rvm_mc(i, k) = rvm_mc(i, k) + hydromet(i, k, ihm) / dt;
              if (!l_frozen)
                thlm_mc(i, k) =
                    thlm_mc(i, k) - (Lv / (Cp * exner(i, k))) * (hydromet(i, k, ihm) / dt);
              else
                thlm_mc(i, k) =
                    thlm_mc(i, k) - (Ls / (Cp * exner(i, k))) * (hydromet(i, k, ihm) / dt);
              hydromet(i, k, ihm) = zero;
            }
          });
      if (stats.l_sample)
        budget(name_cl, ihm, false);
    }
  }
  // Concentration clipping does not change the mixing-ratio fields.
  for (int i = 0; i < ngrdcol; ++i) {
    auto input =
        to_execution_space(workspace, Kokkos::subview(hydromet, i, Kokkos::ALL, Kokkos::ALL));
    auto output = to_execution_space(
        workspace, Kokkos::subview(hydromet_clipped, i, Kokkos::ALL, Kokkos::ALL), false);
    clip_hydromet_conc_mvr(workspace, nzt, hydromet_dim, hm_metadata, input, output);
    from_execution_space(workspace, Kokkos::subview(hydromet_clipped, i, Kokkos::ALL, Kokkos::ALL),
                         output);
  }
  for (int ihm = 0; ihm < hydromet_dim; ++ihm) {
    if (!l_mix_rat_hm_host(ihm)) {
      setup_stats_names(ihm, hydromet_dim, hydromet_list_host, max_velocity, &name_bt, &name_hf,
                        &name_wvhf, &name_cl, &name_mc);
      if (stats.l_sample)
        budget(name_cl, ihm, true);
      Kokkos::deep_copy(Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm),
                        Kokkos::subview(hydromet_clipped, Kokkos::ALL, Kokkos::ALL, ihm));
      if (stats.l_sample)
        budget(name_cl, ihm, false);
    }
  }
}

//=============================================================================
// nzt: Number of thermodynamic levels [-]
// hydromet_dim: Number of hydrometeor fields
// hm_metadata: Hydrometeor species indices, flags and tolerances [-,units vary]
// hydromet: Mean of hydrometeor fields    [units vary]
// hydromet_clipped: Clipped mean of hydrometeor fields    [units vary]
void clip_hydromet_conc_mvr(ScratchWorkspace& workspace, int nzt, int hydromet_dim,
                            const hm_metadata_type& hm_metadata, // In
                            Real2D hydromet,                     // In
                            Real2D hydromet_clipped)             // Out
{
  auto workspace_scope = workspace.scope();
  // Description:
  // Increases the value of a hydrometeor concentration when it is too small,
  // according to the hydrometeor mixing ratio (which remains unchanged) and
  // the maximum acceptable drop or particle mean volume radius for that
  // hydrometeor species.

  // References:
  //-----------------------------------------------------------------------
  auto mix    = to_host_space(workspace, hm_metadata.l_mix_rat_hm);
  auto frozen = to_host_space(workspace, hm_metadata.l_frozen_hm);
  for (int idx = 0; idx < hydromet_dim; ++idx) {
    if (!mix(idx)) {
      // Precompute public Native metadata functions on the host once per species.
      // No scalar bridge enters the independent level kernel.
      real_t Nxm_min_coef;
      const real_t mvr_hydromet_max = mvr_hm_max(idx, hm_metadata);
      if (!frozen(idx)) {
        // <Nr>=<rr>/((4/3)*pi*rho_lw*mvr_rain_max^3).
        Nxm_min_coef = one / (four_thirds * pi * rho_lw *
                              ((mvr_hydromet_max * mvr_hydromet_max) * mvr_hydromet_max));
      } else {
        // Frozen species use rho_ice and their species-specific maximum radius.
        Nxm_min_coef = one / (four_thirds * pi * rho_ice *
                              ((mvr_hydromet_max * mvr_hydromet_max) * mvr_hydromet_max));
      }
      const int rx_idx = ((Nx2rx_hm_idx(idx, hm_metadata) + 1)) - 1;
      Kokkos::parallel_for(
          "clip_hydromet_conc_mvr::levels", nzt, KOKKOS_LAMBDA(int level) {
            const int k = level;
            if (hydromet(k, rx_idx) > zero) {
              const real_t minimum     = Nxm_min_coef * hydromet(k, rx_idx);
              hydromet_clipped(k, idx) = hydromet(k, idx) > minimum ? hydromet(k, idx) : minimum;
            } else {
              hydromet_clipped(k, idx) = zero;
            }
          });
    } else {
      Kokkos::deep_copy(Kokkos::subview(hydromet_clipped, Kokkos::ALL, idx),
                        Kokkos::subview(hydromet, Kokkos::ALL, idx));
    }
  }
}

//=============================================================================
// ihm: Zero-based species slot [-]
// hydromet_dim: Number of hydrometeor species [-]
// hydromet_list: Fixed-width species names (character,species) [-]
// max_velocity: Maximum sedimentation velocity [m/s]
void setup_stats_names(
    int ihm, int hydromet_dim,                                                 // In
    Kokkos::View<char**, Kokkos::LayoutLeft, Kokkos::HostSpace> hydromet_list, // In
    real_t& max_velocity,                                                      // InOut
    std::string* name_bt, std::string* name_hf,
    std::string* name_wvhf,                     // Optional Out
    std::string* name_cl, std::string* name_mc) // Optional Out
{
  // Description:
  //
  // Determines the stats output names depending on the hydrometeor.
  // Attention: hydromet_list needs to be set up before this routine is called.
  //
  // Bogus example
  std::string name_bt_local, name_hf_local, name_wvhf_local, name_cl_local, name_mc_local;
  max_velocity = zero;
  name_bt_local = "";
  name_hf_local = "";
  name_wvhf_local = "";
  name_cl_local = "";
  name_mc_local = "";
  std::string species(hydromet_list.data() + (hydromet_list.extent(0) * ihm),
                      hydromet_list.extent(0));
  while (!species.empty() && species.back() == ' ')
    species.pop_back();
  if (species == "rrm") {
    name_bt_local = "rrm_bt";
    name_hf_local = "rrm_hf";
    name_wvhf_local = "rrm_wvhf";
    name_cl_local = "rrm_cl";
    name_mc_local = "rrm_mc";
    max_velocity = -real_t(9.1);
  } else if (species == "rim") {
    name_bt_local = "rim_bt";
    name_hf_local = "rim_hf";
    name_wvhf_local = "rim_wvhf";
    name_cl_local = "rim_cl";
    name_mc_local = "rim_mc";
    // Morrison limit
    //         max_velocity = -1.2_core_rknd ! m/s
    // Made up limit.  The literature suggests that it is quite possible
    // that snow flake might achieve a terminal velocity of 2 m/s, and this
    // happens in the COAMPS microphysics -dschanen 29 Sept 2009
    max_velocity = -real_t(1.2);
  } else if (species == "rsm") {
    name_bt_local = "rsm_bt";
    name_hf_local = "rsm_hf";
    name_wvhf_local = "rsm_wvhf";
    name_cl_local = "rsm_cl";
    name_mc_local = "rsm_mc";
    max_velocity = -real_t(2.0);
  } else if (species == "rgm") {
    name_bt_local = "rgm_bt";
    name_hf_local = "rgm_hf";
    name_wvhf_local = "rgm_wvhf";
    name_cl_local = "rgm_cl";
    name_mc_local = "rgm_mc";
    max_velocity = -real_t(20.);
  } else if (species == "Nrm") {
    name_bt_local = "Nrm_bt";
    name_cl_local = "Nrm_cl";
    name_mc_local = "Nrm_mc";
    max_velocity = -real_t(9.1);
  } else if (species == "Nim") {
    name_bt_local = "Nim_bt";
    name_cl_local = "Nim_cl";
    name_mc_local = "Nim_mc";
    // Morrison limit
    //         max_velocity = -1.2_core_rknd ! m/s
    // Made up limit.  The literature suggests that it is quite possible
    // that snow flake might achieve a terminal velocity of 2 m/s, and this
    // happens in the COAMPS microphysics -dschanen 29 Sept 2009
    max_velocity = -real_t(1.2);
  } else if (species == "Nsm") {
    name_bt_local = "Nsm_bt";
    name_cl_local = "Nsm_cl";
    name_mc_local = "Nsm_mc";
    max_velocity = -real_t(2.0);
  } else if (species == "Ngm") {
    name_bt_local = "Ngm_bt";
    name_cl_local = "Ngm_cl";
    name_mc_local = "Ngm_mc";
    max_velocity = -real_t(20.);
  } else if (species == "Ncm") {
    name_bt_local = "Ncm_bt";
    name_cl_local = "Ncm_cl";
    name_mc_local = "Ncm_mc";
    // Use the rain water limit, since Morrison has no explicit limit on
    // cloud water.  Presumably these numbers are never large.
    // -dschanen 28 Sept 2009
    max_velocity = -real_t(9.1);
  } else {
    max_velocity = zero;
  }
  if (name_bt != nullptr)
    *name_bt = name_bt_local;
  if (name_hf != nullptr)
    *name_hf = name_hf_local;
  if (name_wvhf != nullptr)
    *name_wvhf = name_wvhf_local;
  if (name_cl != nullptr)
    *name_cl = name_cl_local;
  if (name_mc != nullptr)
    *name_mc = name_mc_local;
}
} // namespace clubb_kokkos
