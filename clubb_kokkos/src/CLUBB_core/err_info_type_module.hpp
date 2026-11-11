#pragma once

#include "scratch_workspace.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"

#include "array_types.hpp"
#include <optional>
#include <type_traits>

namespace clubb_kokkos {

// BEGIN GENERATED ERR INFO TYPE
// Mirrors src/CLUBB_core/err_info_type_module.F90; arrays use Fortran-order LayoutLeft.
struct err_info_type {
  Kokkos::View<int*, Kokkos::LayoutLeft> err_code;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> host_lat;
  Kokkos::View<real_t*, Kokkos::LayoutLeft> host_lon;
  int chunk_idx{};
  int mpi_rank{};
  Kokkos::View<char**, Kokkos::LayoutLeft> err_header;
  std::string err_header_global{};
};

// END GENERATED ERR INFO TYPE

void init_default_err_info_api(ScratchWorkspace& workspace, // InOut
                               int ngrdcol,                 // In
                               err_info_type& err_info      // InOut
);
void init_err_info_api(ScratchWorkspace& workspace, // InOut
                       int ngrdcol, int chunk_idx_in, int mpi_rank_in, Real1D lat_in,
                       Real1D lon_in,          // In
                       err_info_type& err_info // InOut
);
void set_err_info_values_api(ScratchWorkspace& workspace,                    // InOut
                             int ngrdcol,                                    // In
                             err_info_type& err_info,                        // InOut
                             std::optional<int> chunk_idx_in = std::nullopt, // Optional In
                             std::optional<int> mpi_rank_in  = std::nullopt, // Optional In
                             std::optional<Real1D> lat_in    = std::nullopt, // Optional In
                             std::optional<Real1D> lon_in    = std::nullopt  // Optional In
);
void cleanup_err_info_api(err_info_type& err_info // InOut
);
template <class Lat, class Lon>
inline void init_err_info_api(ScratchWorkspace& workspace, // InOut
                              int ngrdcol, int chunk_idx_in, int mpi_rank_in, Lat lat_in,
                              Lon lon_in,             // In
                              err_info_type& err_info // InOut
) {
  // Stage host/strided inputs once through the existing array owner.
  auto workspace_scope = workspace.scope();
  init_err_info_api(workspace, ngrdcol, chunk_idx_in, mpi_rank_in,
                    to_execution_space(workspace, lat_in), // In
                    to_execution_space(workspace, lon_in), // In
                    err_info                               // InOut
  );
}
template <class Lat = Real1D, class Lon = Real1D>
inline void set_err_info_values_api(ScratchWorkspace& workspace,                    // InOut
                                    int ngrdcol,                                    // In
                                    err_info_type& err_info,                        // InOut
                                    std::optional<int> chunk_idx_in = std::nullopt, // Optional In
                                    std::optional<int> mpi_rank_in  = std::nullopt, // Optional In
                                    std::optional<Lat> lat_in       = std::nullopt, // Optional In
                                    std::optional<Lon> lon_in       = std::nullopt  // Optional In
) {
  // Optional HostSpace/strided Views retain independent presence semantics.
  auto workspace_scope = workspace.scope();
  std::optional<Real1D> lat_exec;
  std::optional<Real1D> lon_exec;

  if (lat_in.has_value())
    lat_exec = to_execution_space(workspace, *lat_in);

  if (lon_in.has_value())
    lon_exec = to_execution_space(workspace, *lon_in);
  set_err_info_values_api(workspace,
                          ngrdcol,                                      // In
                          err_info,                                     // InOut
                          chunk_idx_in, mpi_rank_in, lat_exec, lon_exec // Optional In
  );
}
template <class Lat, class Lon,
          std::enable_if_t<Kokkos::is_view_v<Lat> && Kokkos::is_view_v<Lon>, int> = 0>
inline void set_err_info_values_api(ScratchWorkspace& workspace, // InOut
                                    int ngrdcol,                 // In
                                    err_info_type& err_info,     // InOut
                                    std::optional<int> chunk_idx_in, std::optional<int> mpi_rank_in,
                                    Lat lat_in, // Optional In
                                    Lon lon_in  // Optional In
) {
  // Existing driver supplies both coordinate arrays as ordinary Views.
  auto workspace_scope = workspace.scope();
  set_err_info_values_api(workspace,
                          ngrdcol,                                               // In
                          err_info,                                              // InOut
                          chunk_idx_in, mpi_rank_in, std::optional<Lat>(lat_in), // Optional In
                          std::optional<Lon>(lon_in)                             // Optional In
  );
}

} // namespace clubb_kokkos
