// Construct the optional fixed dycore/output grid from its Native height file.
//
// setup_gr_dycore reads the existing 73-level ../input/grid/dycore.grd file,
// trims its endpoints to the requested column surface/top heights [m], and
// calls the C++ grid constructor. File parsing uses the C++ height-file reader.
// Heights and grid storage belong to the caller; no state is retained here.
// Only this fixed-grid setup is provided. Adaptive/generalized-grid algorithms
// remain outside the supported standalone phase.
#pragma once

#include "array_types.hpp"
#include "grid_class.hpp"
#include "error_code.hpp"
#include "err_info_type_module.hpp"
#include <iostream>
#include <string>

namespace clubb_kokkos {

// =============================================================================
template <class Surface, class Top>
void setup_gr_dycore(ScratchWorkspace& workspace,    // InOut
                     int iunit, int ngrdcol,         // In
                     Surface grid_sfc, Top grid_top, // In
                     bool l_ascending_grid,          // In
                     grid& gr,                       // Out
                     err_info_type& err_info) {      // InOut
  // Description:
  //   Set up the fixed dycore grid, replacing the caller's previous grid.
  // iunit: Native logical unit used by the height-file reader [-].
  // ngrdcol: Number of grid columns [-].
  // grid_sfc: Requested lower momentum height by column [m].
  // grid_top: Requested upper momentum height by column [m].
  // l_ascending_grid: True for ascending storage; false for descending storage [-].
  // gr: Caller-owned grid coordinates, spacings and interpolation weights [units vary].
  // err_info: Caller-owned status and diagnostic headers [-].

  // Native IntentOut releases existing allocatable grid components at entry.
  auto workspace_scope = workspace.scope();
  gr                   = {};

  int nlevel, level_min; // Selected zero-based upper/lower file levels [-].
  using HostProfile = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  using HostFields  = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  // TODO: Retain Native's fixed file/extent until its input contract is generalized.
  constexpr int grid_type         = 3;
  constexpr int nzmax             = 73;
  constexpr bool l_implemented    = false;
  const std::string zm_grid_fname = "../input/grid/dycore.grd";
  const std::string zt_grid_fname = "";

  // These call-local buffers are the explicit file-I/O boundary. Stage the
  // input height profiles once; the grid constructor owns execution-space data.
  auto momentum_heights = workspace_scope.view<HostFields>(
      Kokkos::WithoutInitializing, "setup_gr_dycore::momentum_heights", ngrdcol, nzmax);
  auto thermodynamic_heights = workspace_scope.view<HostFields>(
      Kokkos::WithoutInitializing, "setup_gr_dycore::thermodynamic_heights", ngrdcol, nzmax - 1);
  auto sfc_elevation = workspace_scope.view<HostProfile>(Kokkos::WithoutInitializing,
                                                         "setup_gr_dycore::sfc_elevation", ngrdcol);
  auto deltaz        = workspace_scope.view<HostProfile>(Kokkos::WithoutInitializing,
                                                         "setup_gr_dycore::deltaz", ngrdcol);
  const auto grid_sfc_host = to_host_space(workspace, grid_sfc);
  const auto grid_top_host = to_host_space(workspace, grid_top);

  for (int i = 0; i < ngrdcol; ++i) {
    // Named mutable subviews preserve the reader's two Out arguments.
    auto momentum_column      = Kokkos::subview(momentum_heights, i, Kokkos::ALL);
    auto thermodynamic_column = Kokkos::subview(thermodynamic_heights, i, Kokkos::ALL);
    read_grid_heights(workspace, nzmax, grid_type,  // In
                      zm_grid_fname, zt_grid_fname, // In
                      iunit,                        // In
                      momentum_column,              // Out
                      thermodynamic_column,         // Out
                      err_info);                    // InOut
    sfc_elevation(i) = real_t(0.0);
    deltaz(i)        = real_t(0.0);
  }

  // Native reports this error without an early return; preserve that order.
  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << "Fatal error calling read_grid_heights in setup_gr_dycore\n";
  }

  nlevel = nzmax - 1;
  while (momentum_heights(0, nlevel) > grid_top_host(0) && nlevel > 0)
    nlevel = nlevel - 1;

  for (int i = 0; i < ngrdcol; ++i)
    momentum_heights(i, nlevel) = grid_top_host(i);

  level_min = nlevel;
  while (momentum_heights(0, level_min) > grid_sfc_host(0) && level_min > 0)
    level_min = level_min - 1;

  for (int i = 0; i < ngrdcol; ++i)
    momentum_heights(i, level_min) = grid_sfc_host(i);

  setup_grid(nlevel + 1, ngrdcol, sfc_elevation, l_implemented,                 // In
             l_ascending_grid, grid_type, deltaz, grid_sfc_host, grid_top_host, // In
             momentum_heights, thermodynamic_heights,                           // In
             gr, err_info);                                                     // InOut

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Fatal error calling setup_grid in setup_gr_dycore\n";
    }
  }
}

} // namespace clubb_kokkos
