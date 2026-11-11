// Check NetCDF-C input conventions on a tiny file with unequal axis lengths.
//
// The fixture exercises Native IDs/dimension order, one-based hyperslabs,
// float/double/int conversion, untouched output tails, character padding and
// missing-data errors. Whole-model IOP/restart comparisons remain the model gate.
#include "netcdf.hpp"
#include "Input_fields/stat_file_utils.hpp"
#include "Input_fields/inputfields_state.hpp"
#include <cmath>
#include <array>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const char* message) {
  // Fail with the specific file-boundary contract that changed.
  if (!condition)
    throw std::runtime_error(message);
}

void success(int status) {
  // Preserve the library's explanation for any unexpected NetCDF failure.
  if (status != NC_NOERR)
    throw std::runtime_error(nc_strerror(status));
}

void check_model_readers(const std::filesystem::path& path) {
  // Reuse the same file to check zero-based level metadata, strided prefix writes,
  // Native full-output unit conversion and batched exact/interpolated brackets.
  using namespace clubb_kokkos;
  ScratchWorkspace workspace;
  auto scope = workspace.scope();
  inputfields_state state;
  state.clubb_day   = 1;
  state.clubb_month = 1;
  state.clubb_year  = 2001;
  stat_file descriptor;
  bool error;
  open_netcdf_read(workspace, "field", path.string(), descriptor, error, state);
  check(!error && descriptor.ia == 0 && descriptor.iz == 2 && descriptor.ntimes == 2,
        "Reader level endpoints must be zero-based; time records remain counts");
  check(descriptor.dtwrite == real_t(60) && descriptor.time == 0.,
        "Reader date origin or record interval changed");
  Real2D backing("netcdf_input::strided", 2, 5);
  Kokkos::deep_copy(backing, real_t(-99));
  auto row = Kokkos::subview(backing, 1, Kokkos::ALL);
  get_netcdf_var(workspace, descriptor, "field", 2, false, row, error);
  auto host = to_host_space(workspace, backing);
  check(!error && host(1, 0) == real_t(6.1) && host(1, 1) == real_t(8.1) &&
            host(1, 2) == real_t(10.1) && host(1, 3) == real_t(-99) && host(1, 4) == real_t(-99),
        "Strided reader must write only the selected vertical prefix");
  for (int k = 0; k < 5; ++k)
    check(host(0, k) == real_t(-99), "Reader wrote a neighboring column");
  get_netcdf_var(workspace, descriptor, "field", 2, true, row, error);
  host = to_host_space(workspace, backing);
  check(!error && host(1, 0) == real_t(6.1) / g_per_kg && host(1, 3) == real_t(-99) / g_per_kg,
        "Native unit conversion must include the full output tail");
  get_netcdf_var(workspace, descriptor, "field", 3, false, row, error);
  check(error, "Reader must reject records beyond ntimes");
  const auto after_error = to_host_space(workspace, backing);
  check(after_error(1, 0) == host(1, 0) && after_error(1, 3) == host(1, 3),
        "Early reader errors must preserve caller storage");

  Real1D heights("netcdf_input::heights", 6);
  auto heights_host = Kokkos::create_mirror_view(heights);
  const real_t values[]{50, 100, 150, 250, 300, 350};
  for (int k = 0; k < 6; ++k)
    heights_host(k) = values[k];
  Kokkos::deep_copy(heights, heights_host);
  int first, last;
  CLUBB_levels_within_LES_domain(workspace, descriptor, heights, first, last);
  check(first == 1 && last == 4, "LES overlap bounds must be inclusive zero-based slots");
  auto exact  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "test::exact", 4);
  auto lower  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "test::lower", 4);
  auto upper  = scope.view<InputIndices>(Kokkos::WithoutInitializing, "test::upper", 4);
  auto linear = scope.view<Bool1D>(Kokkos::WithoutInitializing, "test::linear", 4);
  LES_grid_to_CLUBB_grid(descriptor, ConstInputProfile(heights), first, last, exact, lower, upper,
                         linear);
  const auto exact_host  = to_host_space(workspace, exact);
  const auto lower_host  = to_host_space(workspace, lower);
  const auto upper_host  = to_host_space(workspace, upper);
  const auto linear_host = to_host_space(workspace, linear);
  check(exact_host(0) == 0 && exact_host(3) == 2 && !linear_host(0) && !linear_host(3),
        "Batched exact LES matches changed");
  check(exact_host(1) == -1 && lower_host(1) == 0 && upper_host(1) == 1 && lower_host(2) == 1 &&
            upper_host(2) == 2 && linear_host(1) && linear_host(2),
        "Batched neighboring LES slots or unset sentinel changed");
  Real1D average("netcdf_input::average", 6);
  stat_file_average(workspace, path.string(), 6, 1, 2, heights, "field", 1, false, error, state,
                    average);
  auto average_host           = to_host_space(workspace, average);
  const real_t lowest_average = (real_t(0.1) / g_per_kg + real_t(6.1) / g_per_kg) / real_t(2);
  check(!error && average_host(1) == lowest_average && std::isinf(average_host(0)) &&
            average_host(0) > 0 && std::isinf(average_host(5)) && average_host(5) > 0,
        "Record average or disabled lower/upper extension changed");
  const auto ghost =
      stat_file_average(workspace, path.string(), 6, 1, 2, heights, "field", 1, true, error, state);
  const auto ghost_host = to_host_space(workspace, ghost);
  check(!error && ghost_host(0) == lowest_average,
        "Special lower ghost value must retain the lowest file average");
  check(stat_file_num_vertical_levels(workspace, "field", path.string(), state) == 3,
        "A last-slot endpoint must not replace the vertical count");
  close_netcdf_read(descriptor);
  open_netcdf_read(workspace, "no_height", path.string(), descriptor, error, state);
  check(error && descriptor.iz == -2,
        "Missing Native -1 height endpoint must retain a distinct zero-coordinate sentinel");
  close_netcdf_read(descriptor);
}

} // namespace

int main(int argc, char* argv[]) {
  Kokkos::ScopeGuard kokkos(argc, argv);
  // Write a small independent C-order fixture, then read through the Native-order
  // public service. Unequal dimensions expose an accidental missing reversal.
  namespace nc = clubb_kokkos::netcdf;
  const auto path =
      std::filesystem::temp_directory_path() /
      ("clubb_netcdf_input_" +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".nc");
  int file = -1;
  try {
    int time, lev, col, field, profile, matrix, dates, scalar, time_var, x, y, no_height;
    success(nc_create(path.c_str(), NC_CLOBBER, &file));
    success(nc_def_dim(file, "time", 2, &time));
    success(nc_def_dim(file, "lev", 3, &lev));
    success(nc_def_dim(file, "col", 2, &col));
    const int dimensions[]{time, lev, col};
    success(nc_def_var(file, "field", NC_DOUBLE, 3, dimensions, &field));
    success(nc_def_var(file, "lev", NC_FLOAT, 1, &lev, &profile));
    success(nc_def_var(file, "matrix", NC_DOUBLE, 2, dimensions, &matrix));
    success(nc_def_var(file, "dates", NC_INT, 1, &time, &dates));
    success(nc_def_var(file, "scalar", NC_INT, 0, nullptr, &scalar));
    success(nc_put_att_text(file, field, "units", 4, "g/kg"));
    success(nc_put_att_text(file, field, "long_text", 11, "abcdefghijk"));
    success(nc_put_att_text(file, NC_GLOBAL, "title", 7, "fixture"));
    success(nc_put_att_int(file, field, "numeric", NC_INT, 1, &col));
    success(nc_def_var(file, "time", NC_DOUBLE, 1, &time, &time_var));
    const std::string time_units = "seconds since 2001-01-01 00:00:00.0";
    success(nc_put_att_text(file, time_var, "units", time_units.size(), time_units.data()));
    success(nc_def_dim(file, "x", 1, &x));
    success(nc_def_dim(file, "y", 1, &y));
    const int no_height_dims[]{time, y, x};
    success(nc_def_var(file, "no_height", NC_FLOAT, 3, no_height_dims, &no_height));
    success(nc_enddef(file));
    const double field_values[]{0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1, 11.1};
    const float profile_values[]{100.f, 200.f, 300.f};
    const double matrix_values[]{1., 2., 3., 4., 5., 6.};
    const int date_values[]{20010101, 20010102}, scalar_value = 19991231;
    success(nc_put_var_double(file, field, field_values));
    success(nc_put_var_float(file, profile, profile_values));
    success(nc_put_var_double(file, matrix, matrix_values));
    success(nc_put_var_int(file, dates, date_values));
    success(nc_put_var_int(file, scalar, &scalar_value));
    const double time_values[]{0., 60.};
    success(nc_put_var_double(file, time_var, time_values));
    success(nc_close(file));
    file = -1;

    success(nc::open(path.string() + "  ", file));
    int varid, dimid, length, type, rank;
    success(nc::inq_varid(file, "field  ", varid));
    check(varid == field + 1, "Variable IDs must remain one-based");
    success(nc::inq_dimid(file, "lev", dimid));
    check(dimid == lev + 1, "Dimension IDs must remain one-based");
    std::array<int, NC_MAX_VAR_DIMS> dimids;
    success(nc::inquire_variable(file, varid, type, rank, dimids.data()));
    check(type == NC_DOUBLE && rank == 3 && dimids[0] == col + 1 && dimids[1] == lev + 1 &&
              dimids[2] == time + 1,
          "Variable inquiry must reverse C dimension order");
    std::array<char, 8> name;
    success(nc::inquire_dimension(file, dimid, name.data(), name.size(), length));
    check(length == 3 && std::string(name.data(), name.size()) == "lev     ",
          "Dimension names must retain Fortran blank padding");
    const int start[]{2, 1, 2}, count[]{1, 3, 1};
    std::array<double, 5> values{-99., -99., -99., -99., -99.};
    success(nc::get_var(file, varid, values.size(), 3, start, count, values.data()));
    check(values == std::array<double, 5>{7.1, 9.1, 11.1, -99., -99.},
          "Native hyperslab selection or unread output tail changed");
    std::array<float, 3> single;
    success(nc::get_var(file, varid, single.size(), 3, start, count, single.data()));
    check(single == std::array<float, 3>{float(7.1), float(9.1), float(11.1)},
          "NetCDF double-to-default-real conversion changed");
    check(nc::get_var(file, varid, 2, 3, start, count, single.data()) == NC_EEDGE,
          "An undersized host destination must be rejected before writing");
    const int invalid_start[]{0, 1, 2};
    check(nc::get_var(file, varid, 3, 3, invalid_start, count, single.data()) == NC_EINVALCOORDS,
          "Native start=0 must retain the NetCDF invalid-coordinate error");
    std::array<char, 8> units;
    success(nc::get_att(file, varid, "units", units.data(), units.size()));
    check(std::string(units.data(), units.size()) == "g/kg    ",
          "Text attributes must retain Fortran blank padding");
    success(nc::get_att(file, varid, "long_text", units.data(), units.size()));
    check(std::string(units.data(), units.size()) == "abcdefgh",
          "Overlong character metadata must stay within its caller's buffer");
    check(nc::get_att(file, varid, "missing", units.data(), units.size()) == NC_ENOTATT,
          "Missing-attribute status changed");
    check(nc::get_att(file, varid, "numeric", units.data(), units.size()) == NC_ECHAR,
          "Numeric-to-character attribute rejection changed");
    success(nc::get_att(file, 0, "title", units.data(), units.size()));
    check(std::string(units.data(), units.size()) == "fixture ",
          "Native global attribute ID must remain zero");

    success(nc::inq_varid(file, "matrix", varid));
    const int all_start[]{1, 1}, all_count[]{3, 2};
    std::array<float, 6> raw;
    success(nc::get_var(file, varid, raw.size(), 2, all_start, all_count, raw.data()));
    check(raw == std::array<float, 6>{1.f, 2.f, 3.f, 4.f, 5.f, 6.f},
          "The service must leave the raw IOP sequence for its existing reshape");
    success(nc::inq_varid(file, "lev", varid));
    std::array<double, 4> levels{-1., -1., -1., -1.};
    const int first = 1, three = 3;
    success(nc::get_var(file, varid, levels.size(), 1, &first, &three, levels.data()));
    check(levels == std::array<double, 4>{100., 200., 300., -1.},
          "Default-real-to-core-real conversion or surface extension changed");
    success(nc::inq_varid(file, "dates", varid));
    std::array<int, 2> date;
    success(nc::get_var(file, varid, date.size(), date.data()));
    check(date == std::array<int, 2>{20010101, 20010102}, "Integer vector values changed");
    int value;
    success(nc::get_var(file, varid, value));
    check(value == 20010101, "Scalar reading a vector must select its first element");
    success(nc::inq_varid(file, "scalar", varid));
    success(nc::get_var(file, varid, value));
    check(value == 19991231, "Rank-zero integer date read changed");
    check(nc::inq_varid(file, "absent", varid) == NC_ENOTVAR,
          "Optional missing-variable status changed");
    success(nc::close(file));
    file = -1;
    check_model_readers(path);
    std::filesystem::remove(path);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (file >= 0)
      nc_close(file);
    std::filesystem::remove(path);
    return 1;
  }
}
