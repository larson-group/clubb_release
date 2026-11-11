// Construct the optional fixed dycore/output grid from its Native height file.
//
// setup_gr_dycore reads the existing 73-level ../input/grid/dycore.grd file,
// trims its endpoints to the requested column surface/top heights [m], and
// calls the C++ grid constructor. File parsing stays in the public f2k reader.
// Heights and grid storage belong to the caller; no state is retained here.
// Only this fixed-grid setup is provided. Adaptive/generalized-grid algorithms
// remain outside the supported standalone phase.
#pragma once

#include "array_types.hpp"
#include "grid_class.hpp"
#include "error_code.hpp"
#include "err_info_type_module.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iostream>
#include <string>

namespace clubb_kokkos {

// =============================================================================
template <class Surface, class Top>
void setup_gr_dycore(int iunit, int ngrdcol,         // In
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
  gr = {};

  int nlevel, level_min; // Native one-based selected upper/lower file levels [-].
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
  HostFields momentum_heights(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "setup_gr_dycore::momentum_heights"), ngrdcol,
      nzmax);
  HostFields thermodynamic_heights(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "setup_gr_dycore::thermodynamic_heights"),
      ngrdcol, nzmax - 1);
  HostProfile sfc_elevation(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "setup_gr_dycore::sfc_elevation"), ngrdcol);
  HostProfile deltaz(Kokkos::view_alloc(Kokkos::WithoutInitializing, "setup_gr_dycore::deltaz"),
                     ngrdcol);
  const auto grid_sfc_host = to_host_space(grid_sfc);
  const auto grid_top_host = to_host_space(grid_top);

  for (int i = 1; i <= ngrdcol; ++i) {
    // Named mutable subviews preserve the reader's two Out arguments.
    auto momentum_column      = Kokkos::subview(momentum_heights, i - 1, Kokkos::ALL);
    auto thermodynamic_column = Kokkos::subview(thermodynamic_heights, i - 1, Kokkos::ALL);
    f2k::read_grid_heights(nzmax, grid_type,             // In
                           zm_grid_fname, zt_grid_fname, // In
                           iunit,                        // In
                           momentum_column,              // Out
                           thermodynamic_column,         // Out
                           err_info);                    // InOut
    sfc_elevation(i - 1) = real_t(0.0);
    deltaz(i - 1)        = real_t(0.0);
  }

  // Native reports this error without an early return; preserve that order.
  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << "Fatal error calling read_grid_heights in setup_gr_dycore\n";
  }

  nlevel = nzmax;
  while (momentum_heights(0, nlevel - 1) > grid_top_host(0) && nlevel > 1)
    nlevel = nlevel - 1;

  for (int i = 1; i <= ngrdcol; ++i)
    momentum_heights(i - 1, nlevel - 1) = grid_top_host(i - 1);

  level_min = nlevel;
  while (momentum_heights(0, level_min - 1) > grid_sfc_host(0) && level_min > 1)
    level_min = level_min - 1;

  for (int i = 1; i <= ngrdcol; ++i)
    momentum_heights(i - 1, level_min - 1) = grid_sfc_host(i - 1);

  setup_grid(nlevel, ngrdcol, sfc_elevation, l_implemented,                     // In
             l_ascending_grid, grid_type, deltaz, grid_sfc_host, grid_top_host, // In
             momentum_heights, thermodynamic_heights,                           // In
             gr, err_info);                                                     // InOut

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Fatal error calling setup_grid in setup_gr_dycore\n";
    }
  }
}

} // namespace clubb_kokkos
