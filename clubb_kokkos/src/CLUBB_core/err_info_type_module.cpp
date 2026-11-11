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
// characters per column. Fixed-width headers are formatted on the host from
// case-owned context, then copied to the original column-header View.
#include "scratch_workspace.hpp"
#include "err_info_type_module.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace clubb_kokkos {

namespace {
std::string a_field(std::string_view value, std::size_t width) {
  // Implement Native Aw: left-pad short strings or retain the first w characters.
  if (value.size() >= width)
    return std::string(value.substr(0, width));
  return std::string(width - value.size(), ' ') + std::string(value);
}

std::string i_field(int value, std::size_t width) {
  // Implement Native Iw, including all-star overflow instead of widening a field.
  const auto text = std::to_string(value);
  return text.size() > width ? std::string(width, '*') : a_field(text, width);
}

std::string f6_2_field(real_t value) {
  // Implement Native F6.2: signed zero, half-away rounding, nonfinite spelling
  // and all-star overflow. std::fixed rounds exact ties to even on this stack.
  if (std::isnan(value))
    return a_field("NaN", 6);
  if (std::isinf(value))
    return a_field(std::signbit(value) ? "-Inf" : "Inf", 6);
  if (std::abs(value) >= real_t(1000))
    return "******";

  // Scale the exact binary significand by 25 and shift by the remaining power
  // of two (100=25*4). Integer half-away rounding avoids false ties from a
  // rounded floating product and does not depend on extended host precision.
  constexpr int digits = std::numeric_limits<real_t>::digits;
  static_assert(std::numeric_limits<real_t>::radix == 2 && digits <= 53);
  int exponent;
  const auto fraction    = std::frexp(std::abs(value), &exponent);
  const auto significand = static_cast<std::uint64_t>(std::ldexp(fraction, digits));
  const auto numerator   = significand * 25;
  const int shift        = digits - exponent - 2;
  const auto hundredths =
      shift >= 64
          ? 0u
          : static_cast<unsigned int>((numerator + (std::uint64_t(1) << (shift - 1))) >> shift);
  std::string text = (std::signbit(value) ? "-" : "") + std::to_string(hundredths / 100) + '.';
  text += static_cast<char>('0' + (hundredths % 100) / 10);
  text += static_cast<char>('0' + hundredths % 10);
  return text.size() > 6 ? "******" : a_field(text, 6);
}

void update_headers(ScratchWorkspace& workspace, // InOut
                    int ngrdcol,                 // In
                    err_info_type& err_info      // InOut
);
} // namespace

// =============================================================================
void init_default_err_info_api(ScratchWorkspace& workspace, // InOut
                               int ngrdcol,                 // In
                               err_info_type& err_info      // InOut
) {
  // Description: Allocate error context with default model information.
  //              Model information and stored headers are initialized together.
  //
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // err_info: Caller-owned error status, model context and stored headers [-].

  // Dummy defaults are call-local; Native initializes both integer values to 1.
  auto workspace_scope      = workspace.scope();
  const int mpi_rank_dummy  = 1;
  const int chunk_idx_dummy = 1;
  Real1D lat_default = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                    "init_default_err_info_api::lat", ngrdcol);
  Real1D lon_default = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                    "init_default_err_info_api::lon", ngrdcol);
  Kokkos::deep_copy(lat_default, unused_var);
  Kokkos::deep_copy(lon_default, unused_var);

  // err_code is no-error; MPI/chunk are 1; every coordinate is unused_var.
  // Preserve the Native actual argument order of the two equal dummy values.
  init_err_info_api(workspace, ngrdcol, mpi_rank_dummy, chunk_idx_dummy, // In
                    lat_default,                                         // In
                    lon_default,                                         // In
                    err_info                                             // InOut
  );
}

// =============================================================================
void init_err_info_api(ScratchWorkspace& workspace, // InOut
                       int ngrdcol, int chunk_idx_in, int mpi_rank_in, Real1D lat_in,
                       Real1D lon_in,          // In
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
  auto workspace_scope = workspace.scope();
  if (err_info.host_lat.is_allocated() || err_info.host_lon.is_allocated() ||
      err_info.err_code.is_allocated() || err_info.err_header.is_allocated())
    throw std::logic_error("init_err_info_api: Native ALLOCATE requires unallocated fields");

  // Allocate arrays in Native order. Header storage is (character,column).
  err_info.host_lat = decltype(err_info.host_lat)(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "err_info::host_lat"), ngrdcol);
  err_info.host_lon = decltype(err_info.host_lon)(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "err_info::host_lon"), ngrdcol);
  err_info.err_code   = decltype(err_info.err_code)("err_info::err_code", ngrdcol);
  err_info.err_header = decltype(err_info.err_header)(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "err_info::err_header"), 200, ngrdcol);

  // Populate err_code, then the rest of the model context and headers.
  Kokkos::deep_copy(err_info.err_code, clubb_no_error);
  set_err_info_values_api(workspace,
                          ngrdcol,                                  // In
                          err_info,                                 // InOut
                          chunk_idx_in, mpi_rank_in, lat_in, lon_in // Optional In
  );
}

// =============================================================================
void set_err_info_values_api(ScratchWorkspace& workspace, // InOut
                             int ngrdcol,                 // In
                             err_info_type& err_info,     // InOut
                             std::optional<int> chunk_idx_in,
                             std::optional<int> mpi_rank_in, // Optional In
                             std::optional<Real1D> lat_in,
                             std::optional<Real1D> lon_in // Optional In
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
  auto workspace_scope = workspace.scope();
  if (!err_info.err_header.is_allocated() || !err_info.host_lat.is_allocated() ||
      !err_info.err_header.is_allocated()) {
    std::cerr << "Warning! Error assigning values to allocatable arrays in CLUBB's "
                 "err_info struct. Arrays were not allocated before use. They will be "
                 "allocated now, but it is advisable to call an init_api "
                 "subroutines first!"
              << '\n';
    // Native assumes the allocation state is shared by all array components.
    init_default_err_info_api(workspace,
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
    const auto host_lat   = err_info.host_lat;
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
    const auto host_lon   = err_info.host_lon;
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
  update_headers(workspace, ngrdcol, // In
                 err_info            // InOut
  );
}

// =============================================================================
namespace {
void update_headers(ScratchWorkspace& workspace, // InOut
                    int ngrdcol,                 // In
                    err_info_type& err_info      // InOut
) {
  // Description: Rewrite the headers from the current info values
  //
  //-----------------------------------------------------------------------
  // ngrdcol: Number of grid columns [-].
  // err_info: Caller-owned error status, model context and stored headers [-].

  auto workspace_scope = workspace.scope();
  const auto host_lat  = to_host_space(workspace, err_info.host_lat);
  const auto host_lon  = to_host_space(workspace, err_info.host_lon);
  using HostHeaders    = Kokkos::View<char**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto err_header = workspace_scope.view<HostHeaders>(Kokkos::WithoutInitializing,
                                                      "update_headers::err_header", 200, ngrdcol);

  // Populate error header character string that we want to append
  // to the error messages throughout CLUBB.

  // Global (column-independent) header.
  err_info.err_header_global =
      a_field("Fatal error in CLUBB - MPI Process ", 38) + i_field(err_info.mpi_rank, 3) +
      a_field(" / Chunk index ", 15) + i_field(err_info.chunk_idx, 3) + '\n' +
      a_field("Column index not available. ", 28) +
      a_field("Latitude/Longitude for column 1: ", 33) + f6_2_field(host_lat(0)) + '/' +
      f6_2_field(host_lon(0)) + '\n' + a_field("Latitude range from -90 to 90 degrees, ", 39) +
      a_field("Longitude range from -180 to 180 degrees", 40);
  err_info.err_header_global.resize(214, ' ');
  // Preserve the C++ scalar-character contract established by the old converter.
  const auto last = err_info.err_header_global.find_last_not_of(' ');
  err_info.err_header_global.resize(last == std::string::npos ? 0 : last + 1);

  // Populate column-specific headers.
  for (int i = 0; i < ngrdcol; ++i) {
    std::string header =
        a_field("Fatal error in CLUBB in - MPI Process ", 38) + i_field(err_info.mpi_rank, 3) +
        a_field(" / Chunk index ", 15) + i_field(err_info.chunk_idx, 3) + '\n' +
        a_field("Grid column index i = ", 22) + i_field(i + 1, 2) + a_field(", Latitude ", 11) +
        f6_2_field(host_lat(i)) + a_field(" / Longitude ", 13) + f6_2_field(host_lon(i)) + '\n' +
        a_field("Latitude range from -90 to 90 degrees, ", 39) +
        a_field("Longitude range from -180 to 180 degrees", 40);
    header.resize(200, ' ');
    for (int c = 0; c < 200; ++c)
      err_header(c, i) = header[c];
  }
  auto updated_headers =
      Kokkos::subview(err_info.err_header, std::make_pair(0, 200), std::make_pair(0, ngrdcol));
  from_execution_space(workspace, updated_headers, err_header);
}
} // namespace

// =============================================================================
void cleanup_err_info_api(err_info_type& err_info // InOut
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
  err_info.err_code   = {};
  err_info.host_lat   = {};
  err_info.host_lon   = {};
}

} // namespace clubb_kokkos
