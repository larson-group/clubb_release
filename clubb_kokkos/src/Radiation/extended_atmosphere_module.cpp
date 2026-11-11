// Atmosphere extension bounds and grids; see the matching public header.
#include "extended_atmosphere_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

void determine_extended_atmos_bounds(
    ScratchWorkspace& workspace,                               // InOut
    int grid_size, Real1D zt_grid,                             // In
    Real1D zm_grid, Real1D zm_grid_spacing, Real1D p_in_Pa_zm, // In
    int& extended_atmos_bottom_level,                          // Out
    int& extended_atmos_top_level,                             // Out
    int& extended_atmos_range_size,                            // Out
    int& lin_int_buffer_size,                                  // Out
    const parameters_radiation_state& parameters_radiation,    // Native module input
    extended_atmosphere_state& extended_atmosphere             // Native module InOut
) {
  // Determine the bottom/top extension levels and linear-interpolation buffer.
  // The buffer counters cooling spikes from a large temperature difference
  // between the extension profile and the computational profile.
  // Reference: McClatchey et al. (1972), Environmental Research Papers411,p.94.

  // Source bounds are scalar initialization control. Host mirrors preserve
  // ordered diagnostics; complete-grid construction below uses the array owner.
  auto workspace_scope             = workspace.scope();
  const auto extended_alt_host     = to_host_space(workspace, extended_atmosphere.extended_alt);
  const auto extended_p_in_mb_host = to_host_space(workspace, extended_atmosphere.extended_p_in_mb);
  const auto zm_grid_host          = to_host_space(workspace, zm_grid);
  const auto spacing_host          = to_host_space(workspace, zm_grid_spacing);
  const auto pressure_host         = to_host_space(workspace, p_in_Pa_zm);
  const real_t radiation_top       = parameters_radiation.radiation_top;
  const int extended_atmos_dim     = extended_atmosphere.extended_atmos_dim;

  real_t dz10, dz_model, dz, zm_grid_top, extended_bottom, buffer_size;
  int j, k;

  // ---- Begin Code ----
  if (radiation_top < zm_grid_host(grid_size - 1)) {
    std::cerr << "In subroutine determine_extended_atmos_bounds\n";
    throw std::runtime_error(
        "top of the radiation grid is below the top of the computational grid");
  }

  // Determine the bounds to use for the extended atmosphere.
  j = 0;

  // Altitudes increase with increasing grid level.
  while (extended_alt_host(j) < zm_grid_host(grid_size - 1) && j < extended_atmos_dim - 1) {
    ++j;
  }

  // Pressure decreases with increasing grid level. Native valid input must
  // contain a bracketing extension level; no extra clipping is introduced.
  while (pressure_host(grid_size - 1) < extended_p_in_mb_host(j) * pascal_per_mb) {
    ++j;
  }

  if (extended_alt_host(j) < zm_grid_host(grid_size - 1)) {
    std::cerr << "In subroutine determine_extended_atmos_bounds\n";
    throw std::runtime_error("Extended atmosphere is below the top of the computational grid");
  }
  if (extended_alt_host(extended_atmos_dim - 1) < radiation_top) {
    std::cerr
        << "In subroutine determine_extended_atmos_bounds\n"
        << "Atmosphere cannot be extended because extension data does not reach radiation_top\n";
    throw std::runtime_error("ERROR STOP");
  }
  if (pressure_host(grid_size - 1) < extended_p_in_mb_host(j) * pascal_per_mb) {
    std::cerr << "In subroutine determine_extended_atmos_bounds\n";
    throw std::runtime_error(
        "pressure at top of computational grid less than pressure at base of radiative grid");
  }

  k = 0;
  if (j < extended_atmos_dim) {
    while (extended_alt_host(k) < radiation_top && k < extended_atmos_dim - 1) {
      ++k;
    }
    // Roll back when the selected level is above the specified radiation top.
    if (extended_alt_host(k) > radiation_top)
      --k;
  } else {
    k = j;
  }

  extended_atmos_bottom_level = j;
  extended_atmos_top_level    = k;
  extended_atmos_range_size   = k - j + 1;
  if (extended_atmos_range_size < 1) {
    std::cerr << "In subroutine determine_extended_atmos_bounds\n";
    throw std::runtime_error("radiation top below computational grid");
  }

  // Determine a buffer with no more than ten levels.
  zm_grid_top         = zm_grid_host(grid_size - 1);
  extended_bottom     = extended_alt_host(extended_atmos_bottom_level);
  dz10                = (extended_bottom - zm_grid_top) / real_t(10);
  dz_model            = spacing_host(grid_size - 2);
  dz                  = std::max(dz10, dz_model);
  buffer_size         = (extended_bottom - zm_grid_top) / dz;
  lin_int_buffer_size = static_cast<int>(buffer_size);

  extended_atmosphere.total_atmos_dim = grid_size + lin_int_buffer_size + extended_atmos_range_size;
  const int total_atmos_dim           = extended_atmosphere.total_atmos_dim;
  extended_atmosphere.complete_momentum = Real1D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "extended_atmosphere::complete_momentum"),
      total_atmos_dim + 1);
  const auto complete_momentum = extended_atmosphere.complete_momentum;
  const auto extended_alt      = extended_atmosphere.extended_alt;
  const int bottom_level       = extended_atmos_bottom_level;

  // Build the complete momentum grid with zero-based storage, retaining
  // its post-extension overwrite: the final extrapolation can start at a level
  // assigned by the preceding loop, depending on where that loop exits.
  Kokkos::parallel_for(
      "determine_extended_atmos_bounds::momentum", Kokkos::RangePolicy<>(0, 1),
      KOKKOS_LAMBDA(const int) {
        for (int j = 0; j < grid_size; ++j)
          complete_momentum(j) = zm_grid(j);

        int i = 1;
        for (int j = grid_size; j < grid_size + lin_int_buffer_size; ++j) {
          complete_momentum(j) = zm_grid_top + dz * real_t(i);
          ++i;
        }

        i = 0;
        for (int j = grid_size + lin_int_buffer_size; j < total_atmos_dim; ++j) {
          complete_momentum(j) = extended_alt(bottom_level + i);
          ++i;
          if (i + bottom_level == extended_atmos_dim)
            break;
        }

        // Use a linear extension for the topmost points (usually one or two).
        for (int j = grid_size + lin_int_buffer_size + i - 1; j < total_atmos_dim + 1; ++j) {
          const real_t dz_extension = complete_momentum(j - 1) - complete_momentum(j - 2);
          complete_momentum(j)      = complete_momentum(j - 1) + dz_extension;
        }
      });

  extended_atmosphere.complete_alt =
      Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "extended_atmosphere::complete_alt"),
             total_atmos_dim);
  const auto complete_alt = extended_atmosphere.complete_alt;
  Kokkos::parallel_for(
      "determine_extended_atmos_bounds::thermodynamic", Kokkos::RangePolicy<>(0, 1),
      KOKKOS_LAMBDA(const int) {
        for (int j = 0; j < grid_size - 1; ++j)
          complete_alt(j) = zt_grid(j);
        for (int j = grid_size - 1; j < total_atmos_dim; ++j) {
          complete_alt(j) = (complete_momentum(j + 1) + complete_momentum(j)) / real_t(2);
        }
      });
}

void finalize_extended_atm(extended_atmosphere_state& extended_atmosphere) {
  // Free the extension and complete-grid arrays. Native leaves its dimension
  // integers unchanged; a new case obtains a new state owner.
  extended_atmosphere.extended_alt      = Real1D{};
  extended_atmosphere.extended_T_in_K   = Real1D{};
  extended_atmosphere.extended_sp_hmdty = Real1D{};
  extended_atmosphere.extended_p_in_mb  = Real1D{};
  extended_atmosphere.extended_o3l      = Real1D{};
  extended_atmosphere.complete_alt      = Real1D{};
  extended_atmosphere.complete_momentum = Real1D{};
}

} // namespace clubb_kokkos
