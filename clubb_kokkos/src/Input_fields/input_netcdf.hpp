// Read input/restart NetCDF columns with Native dimensions, dates and units.
//
// open_netcdf_read validates the test field and creates a caller-owned descriptor;
// get_netcdf_var reads the first horizontal column at a Native one-based record,
// optionally converting units. close_netcdf_read closes the ID explicitly.
// Raw reads use the NetCDF-C service; date/shape interpretation and math stay here.
#pragma once
#include "stat_file_module.hpp"
#include <string>
#include <array>
#include <iostream>
#include "netcdf.hpp"
#include "CLUBB_core/constants_clubb.hpp"
namespace clubb_kokkos {
struct inputfields_state;
void open_netcdf_read(ScratchWorkspace& workspace, const std::string& test_variable,
                      const std::string& path,               // In
                      stat_file& ncf, bool& l_error,         // Out
                      const inputfields_state& inputfields); // Native module In
void close_netcdf_read(const stat_file& ncf);                // In

template <class Profile>
void get_netcdf_var(ScratchWorkspace& workspace, const stat_file& ncf, const std::string& varname,
                    int itime, bool l_convert_to_MKS, // In
                    Profile x, bool& l_error)         // Out
{
  // Description:
  //   Read in values from a netCDF file.
  auto scope = workspace.scope();
  static_assert(Profile::rank == 1);
  static_assert(std::is_same_v<typename Profile::non_const_value_type, real_t>);
  // Assumptions: NF90_FLOAT is Native real4; variables follow COARDS and have
  // three or four dimensions. The first horizontal column is selected.
  int varid, ierr, xtype, ndims;
  std::array<int, netcdf::nc_max_var_dims> dimIds, start, count;
  l_error = false;
  ierr    = netcdf::inq_varid(ncf.iounit, varname, varid);
  if (ierr != netcdf::nc_noerr) {
    std::cerr << netcdf::strerror(ierr) << '\n';
    l_error = true;
    return;
  }
  ierr = netcdf::inquire_variable(ncf.iounit, varid, xtype, ndims, dimIds.data());
  if (ierr != netcdf::nc_noerr || (ndims != 3 && ndims != 4) ||
      !(xtype == netcdf::nc_double || xtype == netcdf::nc_float)) {
    std::cerr << "input_netcdf.get_var : The netCDF data doesn't "
                 "conform to expected precision, shape, or dimensions\n";
    l_error = true;
    return;
  }
  if (itime > ncf.ntimes) {
    std::cerr << "input_netcdf.get_var: The time specified exceeds the netCDF data\n";
    l_error = true;
    return;
  }
  if (int(x.extent(0)) < ncf.iz + 1) {
    std::cerr << "input_netcdf.get_var: Output array is too small for vertical levels\n";
    l_error = true;
    return;
  }
  int idx_x = -1, idx_y = -1, idx_z = -1, idx_t = -1, idx_col = -1;
  int xdim = -1, ydim = -1, cdim = -1;
  for (int i = 0; i < ndims; ++i) {
    int dim_len;
    std::array<char, netcdf::nc_max_name> name;
    ierr = netcdf::inquire_dimension(ncf.iounit, dimIds[i], name.data(), name.size(), dim_len);
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    std::string dim_name(name.data(), name.size());
    dim_name.erase(dim_name.find_last_not_of(' ') + 1);
    if (dim_name == "Z" || dim_name == "z" || dim_name == "zt" || dim_name == "zm" ||
        dim_name == "altitude" || dim_name == "height" || dim_name == "lev" ||
        dim_name == "lh_zt" || dim_name == "rad_zt" || dim_name == "rad_zm")
      idx_z = i;
    else if (dim_name == "T" || dim_name == "t" || dim_name == "time")
      idx_t = i;
    else if (dim_name == "X" || dim_name == "x" || dim_name == "longitude" || dim_name == "lon") {
      idx_x = i;
      xdim  = dim_len;
    } else if (dim_name == "Y" || dim_name == "y" || dim_name == "latitude" || dim_name == "lat") {
      idx_y = i;
      ydim  = dim_len;
    } else if (dim_name == "col" || dim_name == "column") {
      idx_col = i;
      cdim    = dim_len;
    } else {
      l_error = true;
      return;
    }
  }
  if (idx_z < 0 || idx_t < 0) {
    l_error = true;
    return;
  }
  if (idx_col < 0) {
    if (idx_x < 0 || idx_y < 0 || xdim != 1 || ydim != 1) {
      l_error = true;
      return;
    }
  } else if (cdim < 1) {
    l_error = true;
    return;
  }
  for (int i = 0; i < ndims; ++i) {
    start[i] = 1;
    count[i] = 1;
  }
  start[idx_t] = itime;
  count[idx_z] = ncf.iz + 1;
  // Obtain a vertical column at time itime. Preserve the undefined output tail.
  auto x_host = scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
      Kokkos::WithoutInitializing, "get_netcdf_var::x_host", ncf.iz + 1);
  ierr          = netcdf::get_var(ncf.iounit, varid, ncf.iz + 1, ndims, start.data(), count.data(),
                                  x_host.data());
  auto x_prefix = Kokkos::subview(x, Kokkos::pair<int, int>(0, ncf.iz + 1));
  array_detail::copy_transient(workspace, x_prefix, x_host);
  if (ierr != netcdf::nc_noerr) {
    std::cerr << netcdf::strerror(ierr) << '\n';
    l_error = true;
    return;
  }
  if (l_convert_to_MKS) {
    std::array<char, 10> unit_chars;
    ierr = netcdf::get_att(ncf.iounit, varid, "units", unit_chars.data(), unit_chars.size());
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    } else {
      std::string units(unit_chars.data(), unit_chars.size());
      units.erase(units.find_last_not_of(' ') + 1);
      if (units == "g/kg")
        Kokkos::parallel_for(
            "get_netcdf_var::g_per_kg",
            Kokkos::RangePolicy<typename Profile::execution_space>(0, x.extent(0)),
            KOKKOS_LAMBDA(int i) { x(i) = x(i) / g_per_kg; });
      else if (units == "W/m2") {
        l_error = true;
        std::cerr << "get_netcdf_var: Unable to convert variables of this type to MKS units\n";
        return;
      } else if (units == "K/day")
        Kokkos::parallel_for(
            "get_netcdf_var::per_day",
            Kokkos::RangePolicy<typename Profile::execution_space>(0, x.extent(0)),
            KOKKOS_LAMBDA(int i) { x(i) = x(i) / sec_per_day; });
      // Default: do nothing, exactly as Native.
    }
  }
}
} // namespace clubb_kokkos
