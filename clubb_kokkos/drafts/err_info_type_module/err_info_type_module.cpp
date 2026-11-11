// Initialize, update and release per-case CLUBB error context.
//
// Public routines:
//   init_default_err_info_api: allocate context with default MPI/chunk and coordinates.
//   init_err_info_api: allocate context using caller-provided model information.
//   set_err_info_values_api: update selected context values and regenerate error headers.
//   cleanup_err_info_api: release column arrays while retaining scalar context.
//
// Initialize an unallocated object before use. Change model information through
// set_err_info_values_api so stored headers stay synchronized; err_code may be
// assigned directly. Column headers use (character,column) storage, with 200
// characters per column. Fixed-width header formatting temporarily remains in
// the Native public setter, called with no optional context updates.
#include "err_info_type_module.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iostream>
#include <optional>
#include <stdexcept>

namespace clubb_kokkos {

namespace {
void update_headers(
    int ngrdcol, // In
    err_info_type& err_info // InOut
);
} // namespace

// =============================================================================
void init_default_err_info_api(
    int ngrdcol, // In
    err_info_type& err_info // InOut
) {
  // Description: Allocate error context with default model information.
  //              Model information and stored headers are initialized together.
  //
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // err_info: Caller-owned error status, model context and stored headers [-].

  // Dummy defaults are call-local; Native initializes both integer values to 1.
  const int mpi_rank_dummy = 1;
  const int chunk_idx_dummy = 1;
  Real1D lat_default("init_default_err_info_api::lat", ngrdcol);
  Real1D lon_default("init_default_err_info_api::lon", ngrdcol);
  Kokkos::deep_copy(lat_default, unused_var);
  Kokkos::deep_copy(lon_default, unused_var);

  // err_code is no-error; MPI/chunk are 1; every coordinate is unused_var.
  // Preserve the Native actual argument order of the two equal dummy values.
  init_err_info_api(
      ngrdcol, mpi_rank_dummy, chunk_idx_dummy, // In
      lat_default, // In
      lon_default, // In
      err_info // InOut
  );
}

// =============================================================================
void init_err_info_api(
    int ngrdcol, int chunk_idx_in, int mpi_rank_in, Real1D lat_in, Real1D lon_in, // In
    err_info_type& err_info // InOut
) {
  // Description: Allocate error context with the supplied model information.
  //              Model information and stored headers are initialized together.
  //
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // chunk_idx_in: Index of the chunk of columns that CLUBB is running on [-].
  // mpi_rank_in: MPI rank of the process [-].
  // lat_in: Input column latitudes [degrees].
  // lon_in: Input column longitudes [degrees].
  // err_info: Caller-owned error status, model context and stored headers [-].

  // Native InOut ALLOCATE cannot silently replace already allocated fields.
  if (err_info.host_lat.is_allocated() || err_info.host_lon.is_allocated() ||
      err_info.err_code.is_allocated() || err_info.err_header.is_allocated())
    throw std::logic_error("init_err_info_api: Native ALLOCATE requires unallocated fields");

  // Allocate arrays in Native order. Header storage is (character,column).
  err_info.host_lat = decltype(err_info.host_lat)("err_info::host_lat", ngrdcol);
  err_info.host_lon = decltype(err_info.host_lon)("err_info::host_lon", ngrdcol);
  err_info.err_code = decltype(err_info.err_code)("err_info::err_code", ngrdcol);
  err_info.err_header = decltype(err_info.err_header)("err_info::err_header", 200, ngrdcol);

  // Populate err_code, then the rest of the model context and headers.
  Kokkos::deep_copy(err_info.err_code, clubb_no_error);
  set_err_info_values_api(
      ngrdcol, // In
      err_info, // InOut
      chunk_idx_in, mpi_rank_in, lat_in, lon_in // Optional In
  );
}

// =============================================================================
void set_err_info_values_api(
    int ngrdcol, // In
    err_info_type& err_info, // InOut
    std::optional<int> chunk_idx_in, std::optional<int> mpi_rank_in, // Optional In
    std::optional<Real1D> lat_in, std::optional<Real1D> lon_in // Optional In
) {
  // Description: Set the info variables and rewrite the headers
  //
  // USE THIS SUBROUTINE TO SET INFO VALUES FOR ERR_INFO
  // ERR_CODE CAN BE SET IN THE USUAL WAY BY ASSIGNMENT
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // err_info: Caller-owned error status, model context and stored headers [-].
  // chunk_idx_in: Index of the chunk of columns that CLUBB is running on [-].
  // mpi_rank_in: MPI rank of the process [-].
  // lat_in: Input column latitudes [degrees].
  // lon_in: Input column longitudes [degrees].

  // Preserve Native's actual guard: err_header is repeated; host_lon and
  // err_code are not tested. No repair of partially allocated source states.
  if (!err_info.err_header.is_allocated() || !err_info.host_lat.is_allocated() ||
      !err_info.err_header.is_allocated()) {
    std::cerr << "Warning! Error assigning values to allocatable arrays in CLUBB's "
                 "err_info struct. Arrays were not allocated before use. They will be "
                 "allocated now, but it is advisable to call an init_api "
                 "subroutines first!"
              << '\n';
    // Native assumes the allocation state is shared by all array components.
    init_default_err_info_api(
        ngrdcol, // In
        err_info // InOut
    );
  }

  // Populate parallelization info, when provided.
  if (mpi_rank_in.has_value())
    err_info.mpi_rank = *mpi_rank_in;

  if (chunk_idx_in.has_value())
    err_info.chunk_idx = *chunk_idx_in;

  // Populate and normalize latitude, if provided.
  if (lat_in.has_value()) {
    const auto lat_values = *lat_in;
    const auto host_lat = err_info.host_lat;
    Kokkos::parallel_for(
        "set_err_info_values_api::latitude", ngrdcol, KOKKOS_LAMBDA(const int i) {
          real_t tmp; // Adjusted coordinate for this column [degrees].
          // Values above 90 use the Native assumed [0,180] input range.
          if (lat_values(i) > real_t(90.))
            tmp = lat_values(i) - real_t(90.);
          else
            tmp = lat_values(i);
          host_lat(i) = tmp;
        });
  }

  // Populate and normalize longitude, if provided.
  if (lon_in.has_value()) {
    const auto lon_values = *lon_in;
    const auto host_lon = err_info.host_lon;
    Kokkos::parallel_for(
        "set_err_info_values_api::longitude", ngrdcol, KOKKOS_LAMBDA(const int i) {
          real_t tmp; // Adjusted coordinate for this column [degrees].
          // Values above 180 use the Native assumed [0,360] input range.
          if (lon_values(i) > real_t(180.))
            tmp = lon_values(i) - real_t(180.);
          else
            tmp = lon_values(i);
          host_lon(i) = tmp;
        });
  }

  // Write updated info to headers.
  update_headers(
      ngrdcol, // In
      err_info // InOut
  );
}

// =============================================================================
namespace {
void update_headers(
    int ngrdcol, // In
    err_info_type& err_info // InOut
) {
  // Description: Rewrite the headers from the current info values
  //
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // err_info: Caller-owned error status, model context and stored headers [-].

  // Transient I/O boundary: exact Native I3/F6.2 padding, overflow and character
  // widths remain in Fortran. The alias selects PUBLIC set_err_info_values_api
  // with all four optional inputs ABSENT, so copied coordinates are not
  // normalized a second time. Private Native converter accessors stay transient.
  auto formatted = err_info;
  f2k::set_err_info_values_api(ngrdcol, formatted);

  // Keep original case arrays/metadata and column-header allocation. Only the
  // two formatted header outputs return from this call-local snapshot.
  Kokkos::deep_copy(err_info.err_header, formatted.err_header);
  err_info.err_header_global = formatted.err_header_global;
}
} // namespace

// =============================================================================
void cleanup_err_info_api(
    err_info_type& err_info // InOut
) {
  // Description:
  //   De-allocates the memory for err_info
  //
  // References:
  //   None
  //------------------------------------------------------------------------------
  // err_info: Caller-owned error status, model context and stored headers [-].

  // Kokkos releases reference-counted allocations rather than reporting Native
  // DEALLOCATE STAT values. Keep the source release order and scalar context.
  err_info.err_header = {};
  err_info.err_code = {};
  err_info.host_lat = {};
  err_info.host_lon = {};
}

} // namespace clubb_kokkos
