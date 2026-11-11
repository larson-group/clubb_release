// Interpolate and differentiate fields on CLUBB's staggered vertical grid.
//
// zt2zm_api/zm2zt_api select scalar-level, profile or column/profile operators;
// zt2zm2zt/zm2zt2zm smooth by interpolating to the other grid and back. ddzm/ddzt
// return vertical derivatives on thermodynamic/momentum levels; flip reverses
// a profile. Heights and spacing are in metres; derivatives divide field units
// by metres. Linear weights support stretched grids, with boundary extensions
// and optional minimum values. Native cubic/quintic compile flags are retained.
// Arrays use column-first LayoutLeft in the execution space. The 2D interpolation
// generics also stage host/strided input and return results in the caller's
// memory space. Use these weights,
// rather than half-sums of neighboring fields, on stretched grids.
// Grid endpoints and API level k are zero-based storage coordinates.
// The temporary UDT bridge converts endpoint values at the Fortran boundary.
// setup_grid/cleanup_grid_api construct and release geometry; setup_grid_heights
// refreshes time-varying host heights and interpolation weights. read_grid_heights
// reads stretched-grid files and checks their level counts and increasing heights.
// References: Golaz et al. (2002), JAS59,3540–3551, sec.3c; Equations for CLUBB,
// Grid Configuration; Steffen (1990), A&A239,443–450.
#pragma once

#include <Kokkos_Core.hpp>
#include "scratch_workspace.hpp"
#include <cstdint>
#include <string>
#include "array_types.hpp"
#include "err_info_type_module.hpp"
#include <optional>

namespace clubb_kokkos {

// BEGIN GENERATED GRID TYPE
// Mirrors src/CLUBB_core/grid_class.F90; arrays use Fortran-order LayoutLeft.
struct grid {
  int nzm{};
  int nzt{};
  Kokkos::View<real_t**, Kokkos::LayoutLeft> zm;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> zt;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> invrs_dzm;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> invrs_dzt;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> dzm;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> dzt;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> weights_zm2zt;
  Kokkos::View<real_t***, Kokkos::LayoutLeft> weights_zt2zm;
  int k_lb_zm{};
  int k_lb_zt{};
  int k_ub_zm{};
  int k_ub_zt{};
  int grid_dir_indx{};
  real_t grid_dir{};
};

// END GENERATED GRID TYPE

void setup_grid(int nzmax, int ngrdcol, Real1D sfc_elevation, bool l_implemented, // In
                bool l_ascending_grid, int grid_type, Real1D deltaz, Real1D zm_init,
                Real1D zm_top,                                                 // In
                Real2D momentum_heights, Real2D thermodynamic_heights,         // In
                grid& gr, err_info_type& err_info);                            // InOut
void cleanup_grid_api(grid& gr);                                               // InOut
void setup_grid_heights(int nzm, int nzt, int ngrdcol,                         // In
                        bool l_implemented, bool l_ascending_grid,             // In
                        int grid_type, Real1D deltaz, Real1D zm_init,          // In
                        Real2D momentum_heights, Real2D thermodynamic_heights, // In
                        grid& gr, err_info_type& err_info);                    // InOut

void read_grid_heights(
    int nzmax, int grid_type,                                                           // In
    const std::string& zm_grid_fname,                                                   // In
    const std::string& zt_grid_fname, int file_unit,                                    // In
    Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> momentum_heights,      // Out
    Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> thermodynamic_heights, // Out
    err_info_type& err_info);                                                           // InOut

template <class Momentum, class Thermodynamic>
void read_grid_heights(ScratchWorkspace& workspace,                     // InOut
                       int nzmax, int grid_type,                        // In
                       const std::string& zm_grid_fname,                // In
                       const std::string& zt_grid_fname, int file_unit, // In
                       Momentum momentum_heights,                       // Out
                       Thermodynamic thermodynamic_heights,             // Out
                       err_info_type& err_info) {                       // InOut
  // File I/O uses independent pooled host buffers; restore caller-owned storage.
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> momentum_host =
      to_host_space(workspace, momentum_heights);
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> thermodynamic_host =
      to_host_space(workspace, thermodynamic_heights);
  read_grid_heights(nzmax, grid_type, zm_grid_fname, zt_grid_fname, file_unit, // In
                    momentum_host, thermodynamic_host,                         // Out
                    err_info);                                                 // InOut
  from_execution_space(workspace, momentum_heights, momentum_host);
  from_execution_space(workspace, thermodynamic_heights, thermodynamic_host);
}

// Shape/memory-space adapters retain the caller-owned grid and error state.
// Compatible execution Views alias; host and strided inputs stage once.
template <class Surface, class Spacing, class Start, class Top, class Momentum, class Thermodynamic>
void setup_grid(int nzmax, int ngrdcol, Surface sfc_elevation, bool l_implemented, // In
                bool l_ascending_grid, int grid_type, Spacing deltaz, Start zm_init,
                Top zm_top,                                                     // In
                Momentum momentum_heights, Thermodynamic thermodynamic_heights, // In
                grid& gr, err_info_type& err_info) {                            // InOut
  setup_grid(nzmax, ngrdcol, to_execution_space(sfc_elevation), l_implemented,  // In
             l_ascending_grid, grid_type, to_execution_space(deltaz),           // In
             to_execution_space(zm_init), to_execution_space(zm_top),           // In
             to_execution_space(momentum_heights), to_execution_space(thermodynamic_heights), // In
             gr, err_info); // InOut
}
template <class Spacing, class Start, class Momentum, class Thermodynamic>
void setup_grid_heights(int nzm, int nzt, int ngrdcol,                                   // In
                        bool l_implemented, bool l_ascending_grid,                       // In
                        int grid_type, Spacing deltaz, Start zm_init,                    // In
                        Momentum momentum_heights, Thermodynamic thermodynamic_heights,  // In
                        grid& gr, err_info_type& err_info) {                             // InOut
  setup_grid_heights(nzm, nzt, ngrdcol,                                                  // In
                     l_implemented, l_ascending_grid,                                    // In
                     grid_type, to_execution_space(deltaz), to_execution_space(zm_init), // In
                     to_execution_space(momentum_heights),                               // In
                     to_execution_space(thermodynamic_heights),                          // In
                     gr, err_info);                                                      // InOut
}

// C++ overloads implement the native generic interfaces. Native procedure bodies
// retain their definition order in grid_class.cpp; these declarations dispatch
// rank and preserve optional scalar inputs.
Real2D zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                 std::optional<real_t> zm_min = std::nullopt);
Real1D zt2zm_api(const grid& gr, Real1D azt, std::optional<real_t> zm_min = std::nullopt);
Real1D zt2zm_api(const grid& gr, Real1D azt, real_t zm_min);
real_t zt2zm_api(const grid& gr, Real1D azt, int k, std::optional<real_t> zm_min = std::nullopt);
Real2D zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                 std::optional<real_t> zt_min = std::nullopt);
Real1D zm2zt_api(const grid& gr, Real1D azm, std::optional<real_t> zt_min = std::nullopt);
Real1D zm2zt_api(const grid& gr, Real1D azm, real_t zt_min);
real_t zm2zt_api(const grid& gr, Real1D azm, int k, std::optional<real_t> zt_min = std::nullopt);
Real2D zt2zm2zt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
                std::optional<real_t> zt_min = std::nullopt);
Real2D zm2zt2zm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
                std::optional<real_t> zm_min = std::nullopt);
Real2D ddzm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm);
Real1D ddzm(const grid& gr, Real1D azm);
Real2D ddzt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt);
Real1D ddzt(const grid& gr, Real1D azt);
Real1D flip(Real1D x, int xdim);

// Fill existing output storage instead of constructing and copying a returned
// field. Overlapping inputs retain snapshot semantics through a temporary.
void zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt, // In
               Real2D azm,                                                // Out
               std::optional<real_t> zm_min = std::nullopt);              // Optional In
void zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm, // In
               Real2D azt,                                                // Out
               std::optional<real_t> zt_min = std::nullopt);              // Optional In
void ddzm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,      // In
          Real2D derivative);                                             // Out
void ddzt(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,      // In
          Real2D derivative);                                             // Out

// Explicit workspace variants reserve only routine-local results/packing. Existing
// owning interfaces above remain available for persistent grid/case/snapshot fields.
Real2D zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                 Real2D azt, std::optional<real_t> zm_min = std::nullopt);
void zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Real2D azt, Real2D azm, std::optional<real_t> zm_min = std::nullopt);
Real2D zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                 Real2D azm, std::optional<real_t> zt_min = std::nullopt);
void zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Real2D azm, Real2D azt, std::optional<real_t> zt_min = std::nullopt);
Real2D ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm);
void ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azm,
          Real2D derivative);
Real2D ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt);
void ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Real2D azt,
          Real2D derivative);
Real2D zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                Real2D azt, std::optional<real_t> zt_min = std::nullopt);
void zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Real2D azt, Real2D smoothed, std::optional<real_t> zt_min = std::nullopt);
Real2D zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                Real2D azm, std::optional<real_t> zm_min = std::nullopt);
void zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Real2D azm, Real2D smoothed, std::optional<real_t> zm_min = std::nullopt);

template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, Source azt, // In
               Output azm,                                                // Out
               std::optional<real_t> zm_min = std::nullopt) {             // Optional In
  // Stage host/strided callers at the boundary; dense execution Views alias.
  auto source = to_execution_space(azt);
  auto output = to_execution_space(azm, false);
  zt2zm_api(nzm, nzt, ngrdcol, gr, source, output, zm_min);
  from_execution_space(azm, output);
}

template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, Source azm, // In
               Output azt,                                                // Out
               std::optional<real_t> zt_min = std::nullopt) {             // Optional In
  // Fill caller-owned storage and restore a host/strided output when needed.
  auto source = to_execution_space(azm);
  auto output = to_execution_space(azt, false);
  zm2zt_api(nzm, nzt, ngrdcol, gr, source, output, zt_min);
  from_execution_space(azt, output);
}

template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void ddzm(int nzm, int nzt, int ngrdcol, const grid& gr, Source azm, // In
          Output derivative) {                                       // Out
  // Fill a derivative in the caller's memory space without an intermediate result.
  auto source = to_execution_space(azm);
  auto output = to_execution_space(derivative, false);
  ddzm(nzm, nzt, ngrdcol, gr, source, output);
  from_execution_space(derivative, output);
}

template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void ddzt(int nzm, int nzt, int ngrdcol, const grid& gr, Source azt, // In
          Output derivative) {                                       // Out
  // Preserve derivative arithmetic while adapting only host/strided boundaries.
  auto source = to_execution_space(azt);
  auto output = to_execution_space(derivative, false);
  ddzt(nzm, nzt, ngrdcol, gr, source, output);
  from_execution_space(derivative, output);
}

// Driver profiles may be host-resident or strided; core inputs already on the
// execution space alias their allocation through the shared array adapter.
template <class View>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename View::memory_space>
zt2zm_api(int nzm, int nzt, int ngrdcol, const grid& gr, View azt,
          std::optional<real_t> zm_min = std::nullopt) {
  // Stage the generic call at its array boundary. A host result can then be
  // assigned to the driver's strided host profile without a cross-space copy.
  auto result = zt2zm_api(nzm, nzt, ngrdcol, gr, to_execution_space(azt), zm_min);
  if constexpr (std::is_same_v<typename View::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename View::memory_space> output(
        "zt2zm_api::output", ngrdcol, nzm);
    from_execution_space(output, result);
    return output;
  }
}
template <class View>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename View::memory_space>
zm2zt_api(int nzm, int nzt, int ngrdcol, const grid& gr, View azm,
          std::optional<real_t> zt_min = std::nullopt) {
  // Preserve the caller's memory space across the generic array boundary.
  auto result = zm2zt_api(nzm, nzt, ngrdcol, gr, to_execution_space(azm), zt_min);
  if constexpr (std::is_same_v<typename View::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename View::memory_space> output(
        "zm2zt_api::output", ngrdcol, nzt);
    from_execution_space(output, result);
    return output;
  }
}

template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Source azt, Output azm, std::optional<real_t> zm_min = std::nullopt) {
  // Stage only transient host/strided inputs and restore the existing output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, azm);
  zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, source, output, zm_min);
  from_execution_space(workspace, azm, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
zt2zm_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azt,
          std::optional<real_t> zm_min = std::nullopt) {
  // Retain caller memory space and caller-frame lifetime for returned temporaries.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto result = zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, source, zm_min);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}
template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
               Source azm, Output azt, std::optional<real_t> zt_min = std::nullopt) {
  // Stage only transient host/strided inputs and restore the existing output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, azt);
  zm2zt_api(workspace, nzm, nzt, ngrdcol, gr, source, output, zt_min);
  from_execution_space(workspace, azt, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
zm2zt_api(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azm,
          std::optional<real_t> zt_min = std::nullopt) {
  // Retain caller memory space and caller-frame lifetime for returned temporaries.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto result = zm2zt_api(workspace, nzm, nzt, ngrdcol, gr, source, zt_min);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}
template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azm,
          Output derivative) {
  // Stage only transient host/strided inputs and restore the existing output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, derivative);
  ddzm(workspace, nzm, nzt, ngrdcol, gr, source, output);
  from_execution_space(workspace, derivative, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
ddzm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azm) {
  // Retain caller memory space and caller-frame lifetime for returned temporaries.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto result = ddzm(workspace, nzm, nzt, ngrdcol, gr, source);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}
template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azt,
          Output derivative) {
  // Stage only transient host/strided inputs and restore the existing output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, derivative);
  ddzt(workspace, nzm, nzt, ngrdcol, gr, source, output);
  from_execution_space(workspace, derivative, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
ddzt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azt) {
  // Retain caller memory space and caller-frame lifetime for returned temporaries.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto result = ddzt(workspace, nzm, nzt, ngrdcol, gr, source);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}
template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Source azt, Output smoothed, std::optional<real_t> zt_min = std::nullopt) {
  // Keep smoothing intermediates transient and restore host/strided caller output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, smoothed);
  zt2zm2zt(workspace, nzm, nzt, ngrdcol, gr, source, output, zt_min);
  from_execution_space(workspace, smoothed, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
zt2zm2zt(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azt,
         std::optional<real_t> zt_min = std::nullopt) {
  // Return a caller-frame field in the source memory space without an owning allocation.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azt);
  auto result = zt2zm2zt(workspace, nzm, nzt, ngrdcol, gr, source, zt_min);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}
template <class Source, class Output,
          std::enable_if_t<Kokkos::is_view_v<Source> && Kokkos::is_view_v<Output>, int> Enabled = 0>
void zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
              Source azm, Output smoothed, std::optional<real_t> zm_min = std::nullopt) {
  // Keep smoothing intermediates transient and restore host/strided caller output.
  static_assert(Source::rank == 2 && Output::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t> &&
                std::is_same_v<typename Output::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto output = to_execution_space(workspace, Kokkos::WithoutInitializing, smoothed);
  zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, source, output, zm_min);
  from_execution_space(workspace, smoothed, output);
}

template <class Source, std::enable_if_t<Kokkos::is_view_v<Source>, int> Enabled = 0>
Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Source::memory_space>
zm2zt2zm(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr, Source azm,
         std::optional<real_t> zm_min = std::nullopt) {
  // Return a caller-frame field in the source memory space without an owning allocation.
  static_assert(Source::rank == 2);
  static_assert(std::is_same_v<typename Source::non_const_value_type, real_t>);
  auto source = to_execution_space(workspace, azm);
  auto result = zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, source, zm_min);
  if constexpr (std::is_same_v<typename Source::memory_space, typename Real2D::memory_space>)
    return result;
  else
    return to_host_space(workspace, result);
}

} // namespace clubb_kokkos
