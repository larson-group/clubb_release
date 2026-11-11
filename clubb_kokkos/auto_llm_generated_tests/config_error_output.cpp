// Check diagnostic text edges that successful whole-model comparisons cannot reach.
// Compare public err_info headers with Native I2/I3/F6.2 formatting, including
// rounding, signed zero, nonfinite values, field overflow and metadata-only updates.
// Configuration printers must retain Native row mappings and caller stream state.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "SILHS/parameters_silhs.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" void native_error_headers_c(int ngrdcol, int chunk_idx, int mpi_rank,
                                       const clubb_kokkos::real_t* lat,
                                       const clubb_kokkos::real_t* lon, char* global_header,
                                       char* column_headers);

namespace {
void require(bool condition, const std::string& message) {
  // Identify the output contract that failed without altering the production text.
  if (!condition)
    throw std::runtime_error(message);
}

template <class Lat, class Lon>
void compare(const clubb_kokkos::err_info_type& err_info, int chunk_idx, int mpi_rank,
             const Lat& lat, const Lon& lon) {
  // Native constructs its own context from the same raw coordinates; C++ must
  // retain all padding, line breaks and overflow fields without changing storage.
  const int ngrdcol = lat.extent(0);
  std::array<char, 214> global;
  std::vector<char> columns(200 * ngrdcol);
  native_error_headers_c(ngrdcol, chunk_idx, mpi_rank, lat.data(), lon.data(), global.data(),
                         columns.data());
  std::string expected(global.data(), global.size());
  const auto last = expected.find_last_not_of(' ');
  expected.resize(last == std::string::npos ? 0 : last + 1);
  require(err_info.err_header_global == expected, "Global header differs from Native");
  const auto headers = clubb_kokkos::to_host_space(err_info.err_header);
  for (int i = 0; i < ngrdcol; ++i)
    for (int c = 0; c < 200; ++c)
      require(headers(c, i) == columns[c + 200 * i],
              "Column header differs from Native at column " + std::to_string(i + 1));
}
} // namespace

int main(int argc, char** argv) {
  // Cover field edges with 104 columns (I2 overflow), then compare a metadata-only
  // update and default context. Printer assertions use deliberately altered streams.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    using Host1D          = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    constexpr int ngrdcol = 104;
    const std::array<real_t, 24> values{
        real_t(.125),
        real_t(-.125),
        real_t(.015),
        real_t(-.015),
        real_t(.005),
        real_t(-.005),
        real_t(0),
        -real_t(0),
        real_t(99.995),
        real_t(-99.995),
        real_t(999.99),
        real_t(-999),
        real_t(90),
        real_t(-180),
        real_t(1e20),
        real_t(-1e20),
        std::numeric_limits<real_t>::infinity(),
        -std::numeric_limits<real_t>::infinity(),
        std::numeric_limits<real_t>::quiet_NaN(),
        std::nextafter(real_t(.125), real_t(0)),
        std::numeric_limits<real_t>::denorm_min(),
        -std::numeric_limits<real_t>::min(),
        std::nextafter(real_t(.005), real_t(0)),
        std::nextafter(real_t(.125), std::numeric_limits<real_t>::infinity())};
    Host1D lat("latitude", ngrdcol), lon("longitude", ngrdcol);
    for (int i = 0; i < ngrdcol; ++i) {
      lat(i) = values[i % values.size()];
      lon(i) = values[(i + 1) % values.size()];
    }
    err_info_type err_info;
    init_err_info_api(workspace, ngrdcol, 7, 42, lat, lon, err_info);
    compare(err_info, 7, 42, lat, lon);
    const auto* allocation = err_info.err_header.data();
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    set_err_info_values_api(workspace, ngrdcol, err_info, -999, 1000);
    compare(err_info, -999, 1000, lat, lon);
    require(err_info.err_header.data() == allocation, "Header update replaced caller storage");
    const auto codes = to_host_space(err_info.err_code);
    for (int i = 0; i < ngrdcol; ++i)
      require(codes(i) == clubb_fatal_error, "Header update changed an error code");
    cleanup_err_info_api(err_info);
    init_default_err_info_api(workspace, 1, err_info);
    Host1D unused("unused coordinate", 1);
    unused(0) = real_t(-999);
    compare(err_info, 1, 1, unused, unused);
    cleanup_err_info_api(err_info);

    clubb_config_flags_type clubb_flags{};
    clubb_flags.penta_solve_method = 26;
    clubb_flags.l_mono_flux_lim_um = true;
    silhs_config_flags_type silhs_flags{};
    silhs_flags.cluster_allocation_strategy = 26;
    silhs_flags.l_lh_importance_sampling    = true;
    std::ostringstream stream;
    stream << std::hex << std::showpos << std::boolalpha << std::setfill('#') << std::setw(20000);
    const auto flags = stream.flags();
    print_clubb_config_flags_api(stream, clubb_flags);
    print_silhs_config_flags_api(stream, silhs_flags);
    const auto text = stream.str();
    require(stream.flags() == flags && stream.width() == 20000 && stream.fill() == '#',
            "Configuration printer changed caller formatting");
    require(text.find("penta_solve_method = 26\n") != std::string::npos &&
                text.find("cluster_allocation_strategy = 26\n") != std::string::npos &&
                text.find("l_lh_importance_sampling = T\n") != std::string::npos,
            "Configuration printer inherited caller hex/sign/logical formatting");
    require(text.find("l_mono_flux_lim_vm = T\n") != std::string::npos,
            "Native um/vm diagnostic label mapping changed");
    std::cout << "Configuration streams and Native fixed-width header edges pass\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
