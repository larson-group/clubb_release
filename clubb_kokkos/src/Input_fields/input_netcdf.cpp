// Read Native-compatible NetCDF profiles and interpret their dimensions and time.
#include "input_netcdf.hpp"
#include "inputfields_state.hpp"
#include "CLUBB_core/calendar.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "netcdf.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

void open_netcdf_read(ScratchWorkspace& workspace, const std::string& test_variable,
                      const std::string& path,              // In
                      stat_file& ncf, bool& l_error,        // Out
                      const inputfields_state& inputfields) // Native module In
{
  // Description:
  //   Open a netCDF data set in read-only mode.
  // The appended case context supplies only Native's three module date fields.
  auto scope = workspace.scope();
  // Intent(out) releases allocatable components without defining other scalars.
  ncf.z        = Real1D{};
  ncf.lat_vals = Real1D{};
  ncf.lon_vals = Real1D{};
  ncf.samp_idx = Real1D{};
  ncf.grid_avg_var.clear();
  ncf.samples_of_var.clear();
  ncf.nsamp = 0;

  int ierr, itmp, varid, xtype;
  int xdim, ydim, cdim, ndims;
  double hours, minutes, seconds, multiplier, delta_d;
  int netcdf_day, netcdf_month, netcdf_year, clubb_idate, netcdf_idate;
  std::array<int, netcdf::nc_max_var_dims> dimIds;
  std::string zname, time_name;
  l_error = false;

  ierr = netcdf::open(path, ncf.iounit);
  if (ierr != netcdf::nc_noerr) {
    std::cerr << "Error opening netCDF file: " << path << '\n';
    std::cerr << netcdf::strerror(ierr) << '\n';
    l_error = true;
    return;
  }
  ierr = netcdf::inq_varid(ncf.iounit, test_variable, varid);
  if (ierr != netcdf::nc_noerr) {
    std::cerr << netcdf::strerror(ierr) << '\n';
    l_error = true;
    return;
  }
  ierr = netcdf::inquire_variable(ncf.iounit, varid, xtype, ndims, dimIds.data());
  if (ierr != netcdf::nc_noerr || (ndims != 3 && ndims != 4) ||
      !(xtype == netcdf::nc_double || xtype == netcdf::nc_float)) {
    std::cerr << "input_netcdf.open_netcdf_read : The netCDF data doesn't "
                 "conform to expected precision, shape, or dimensions\n";
    l_error = true;
    return;
  }
  xdim       = -1;
  ydim       = -1;
  cdim       = -1;
  ncf.iz     = -2;
  ncf.ntimes = -1;
  for (int i = 0; i < ndims; ++i) {
    std::array<char, netcdf::nc_max_name> name;
    ierr = netcdf::inquire_dimension(ncf.iounit, dimIds[i], name.data(), name.size(), itmp);
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    std::string dim_name(name.data(), name.size());
    dim_name.erase(dim_name.find_last_not_of(' ') + 1); // Native TRIM.
    if (dim_name == "Z" || dim_name == "z" || dim_name == "zt" || dim_name == "zm" ||
        dim_name == "altitude" || dim_name == "height" || dim_name == "lev" ||
        dim_name == "lh_zt" || dim_name == "rad_zt" || dim_name == "rad_zm") {
      ncf.ia       = 0;
      ncf.iz       = itmp - 1;
      ncf.AltDimId = dimIds[i];
      zname        = dim_name;
    } else if (dim_name == "X" || dim_name == "x" || dim_name == "longitude" || dim_name == "lon") {
      xdim          = itmp;
      ncf.LongDimId = dimIds[i];
    } else if (dim_name == "Y" || dim_name == "y" || dim_name == "latitude" || dim_name == "lat") {
      ydim         = itmp;
      ncf.LatDimId = dimIds[i];
    } else if (dim_name == "T" || dim_name == "t" || dim_name == "time") {
      ncf.ntimes    = itmp;
      ncf.TimeDimId = dimIds[i];
      time_name     = dim_name;
    } else if (dim_name == "col" || dim_name == "column") {
      cdim = itmp;
    } else {
      l_error = true;
      return;
    }
  }
  if (ncf.iz < 0 || ncf.ntimes < 1) {
    l_error = true;
    return;
  }
  if (cdim < 0) {
    if (ydim != 1 || xdim != 1) {
      std::cerr << "input_netcdf.open_netcdf_read : The netCDF data doesn't "
                   "conform to the expected X or Y dimension\n";
      l_error = true;
      return;
    }
  } else if (cdim < 1) {
    std::cerr << "input_netcdf.open_netcdf_read : The netCDF data doesn't "
                 "conform to the expected col dimension\n";
    l_error = true;
    return;
  }
  if (!l_error) {
    // Allocate altitudes and read them through a call-local primitive buffer.
    ncf.z =
        Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "open_netcdf_read::z"), ncf.iz + 1);
    ierr = netcdf::inq_varid(ncf.iounit, zname, ncf.AltVarId);
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    auto z_host = scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
        Kokkos::WithoutInitializing, "open_netcdf_read::z_host", ncf.iz + 1);
    const int first = 1, height_count = ncf.iz + 1;
    ierr = netcdf::get_var(ncf.iounit, ncf.AltVarId, ncf.iz + 1, 1, &first, &height_count,
                           z_host.data());
    Kokkos::deep_copy(ncf.z, z_host);
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    ierr = netcdf::inq_varid(ncf.iounit, time_name, ncf.TimeVarId);
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    std::array<char, 80> time_chars;
    ierr =
        netcdf::get_att(ncf.iounit, ncf.TimeVarId, "units", time_chars.data(), time_chars.size());
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    // Read starting time from the "units" attribute of the netcdf file.
    std::string time(time_chars.data(), time_chars.size());
    time.erase(time.find_last_not_of(' ') + 1);
    const int length = time.size();
    if (length < 29) {
      std::cerr << "The NetCDF file does not have a proper time unit specification. "
                   "The \"units\" attribute for the time variable must be in the form:\n";
      std::cerr << "TIMEUNITS since YYYY-MM-DD HH:MM:SS.S\n";
      l_error = true;
      return;
    }
    // CLUBB files use TIMEUNITS since YYYY-MM-DD HH:MM:SS.S. Some COAMPS
    // files use single-digit month/day without .S; retain Native index rules.
    const int first_hyphen_idx  = time.find('-');
    const int second_hyphen_idx = time.rfind('-');
    const int first_colon_idx   = time.find(':');
    const int second_colon_idx  = time.rfind(':');
    netcdf_year                 = std::stoi(time.substr(first_hyphen_idx - 4, 4));
    std::string dummy_string;
    if (second_hyphen_idx - first_hyphen_idx == 2)
      dummy_string = time.substr(first_hyphen_idx + 1, 1);
    else
      dummy_string = time.substr(first_hyphen_idx + 1, 2);
    netcdf_month = std::stoi(dummy_string);
    if (first_colon_idx - second_hyphen_idx == 5)
      dummy_string = time.substr(second_hyphen_idx + 1, 1);
    else
      dummy_string = time.substr(second_hyphen_idx + 1, 2);
    netcdf_day = std::stoi(dummy_string);
    hours      = std::stod(time.substr(first_colon_idx - 2, 2));
    minutes    = std::stod(time.substr(first_colon_idx + 1, 2));
    seconds    = std::stod(time.substr(second_colon_idx + 1, 2));
    ncf.year   = netcdf_year;
    ncf.month  = netcdf_month;
    ncf.day    = netcdf_day;
    // Difference in days (file minus simulation), fixing different date origins.
    // July 2013 Eric Raut.
    clubb_idate  = gregorian2julian_date(inputfields.clubb_day, inputfields.clubb_month,
                                         inputfields.clubb_year);
    netcdf_idate = gregorian2julian_date(netcdf_day, netcdf_month, netcdf_year);
    delta_d      = double(netcdf_idate - clubb_idate);
    ncf.time     = (hours * double(sec_per_hr)) + (minutes * double(sec_per_min)) + seconds +
                   (delta_d * double(sec_per_day));
    // Get rid of compiler warning.
    multiplier           = 0.0;
    const auto time_unit = time.substr(0, time.find(' '));
    if (time_unit == "hours")
      multiplier = double(sec_per_hr);
    else if (time_unit == "minutes")
      multiplier = double(sec_per_min);
    else if (time_unit == "seconds")
      multiplier = 1.;
    else {
      multiplier = 0.; // Native initialization before unsupported-unit return.
      l_error    = true;
      return;
    }
    std::array<real_t, 2> write_times;
    const int time_count = 2;
    ierr =
        netcdf::get_var(ncf.iounit, ncf.TimeVarId, 2, 1, &first, &time_count, write_times.data());
    if (ierr != netcdf::nc_noerr) {
      std::cerr << netcdf::strerror(ierr) << '\n';
      l_error = true;
      return;
    }
    ncf.dtwrite = (write_times[1] - write_times[0]) * real_t(multiplier);
  }
}

void close_netcdf_read(const stat_file& ncf) // In
{
  // Description:
  //   Close a previously opened netCDF file.
  const int ierr = netcdf::close(ncf.iounit);
  std::cout << "closing " << ncf.iounit << '\n';
  if (ierr != netcdf::nc_noerr) {
    std::cerr << "Error closing a netCDF file\n";
    std::cerr << netcdf::strerror(ierr) << '\n';
    throw std::runtime_error("Fatal error.");
  }
}
} // namespace clubb_kokkos
