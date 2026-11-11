// Convert Native NetCDF identifiers and packed slabs at the C boundary.
#include "netcdf.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>
#include <vector>

namespace clubb_kokkos::netcdf {
namespace {

std::string c_name(std::string_view name) {
  // Fortran character arguments omit trailing blanks before entering NetCDF-C.
  const auto end = name.find_last_not_of(' ');
  return std::string(name.substr(0, end == std::string_view::npos ? 0 : end + 1));
}

template <class Value>
int get_values(int ncid, int varid, int nvalues, int ndims,
               const int* start, const int* count, Value* values) {
  // Reverse a Native slab once and let NetCDF perform the requested kind conversion.
  int file_ndims;
  const int status = nc_inq_varndims(ncid, varid - 1, &file_ndims);
  if (status != NC_NOERR)
    return status;
  if (ndims != file_ndims || nvalues < 0)
    return NC_EINVAL;
  std::vector<std::size_t> c_start(ndims), c_count(ndims);
  std::size_t size = 1;
  for (int dim = 0; dim < ndims; ++dim) {
    if (start[dim] < 1)
      return NC_EINVALCOORDS;
    if (count[dim] < 0)
      return NC_EEDGE;
    c_start[ndims - dim - 1] = std::size_t(start[dim] - 1);
    c_count[ndims - dim - 1] = std::size_t(count[dim]);
    if (count[dim] &&
        size > std::numeric_limits<std::size_t>::max() / std::size_t(count[dim]))
      return NC_EEDGE;
    size *= std::size_t(count[dim]);
  }
  // Native's raw pointer interface cannot report a short destination safely.
  // Reject only that invalid boundary rather than writing past caller storage.
  if (size > std::size_t(nvalues))
    return NC_EEDGE;
  if constexpr (std::is_same_v<Value, float>)
    return nc_get_vara_float(ncid, varid - 1, c_start.data(), c_count.data(), values);
  else if constexpr (std::is_same_v<Value, double>)
    return nc_get_vara_double(ncid, varid - 1, c_start.data(), c_count.data(), values);
  else
    return nc_get_vara_int(ncid, varid - 1, c_start.data(), c_count.data(), values);
}

template <class Value>
int put_values(int ncid, int varid, const Value* values, int nvalues, const std::vector<int>& start,
               const std::vector<int>& count) {
  // Convert the packed source slab once, preserving caller precision and order.
  int ndims;
  int status = nc_inq_varndims(ncid, varid - 1, &ndims);
  if (status != NC_NOERR)
    return status;
  if (nvalues < 0 || start.size() != count.size() ||
      (!start.empty() && start.size() != std::size_t(ndims)))
    return NC_EINVAL;
  std::vector<std::size_t> c_start(ndims, 0), c_count(ndims);
  if (start.empty()) {
    // A Fortran source shape is required for records and partial writes. Only
    // fixed full variables can infer their counts from the file dimensions.
    std::vector<int> dimids(ndims);
    status = nc_inq_vardimid(ncid, varid - 1, dimids.data());
    if (status != NC_NOERR)
      return status;
    int unlimited;
    status = nc_inq_unlimdim(ncid, &unlimited);
    if (status != NC_NOERR)
      return status;
    for (int dim = 0; dim < ndims; ++dim) {
      if (dimids[dim] == unlimited)
        return NC_EINVAL;
      status = nc_inq_dimlen(ncid, dimids[dim], &c_count[dim]);
      if (status != NC_NOERR)
        return status;
    }
  } else {
    for (int dim = 0; dim < ndims; ++dim) {
      if (start[dim] < 1)
        return NC_EINVALCOORDS;
      if (count[dim] < 0)
        return NC_EEDGE;
      c_start[ndims - dim - 1] = std::size_t(start[dim] - 1);
      c_count[ndims - dim - 1] = std::size_t(count[dim]);
    }
  }
  std::size_t size = 1;
  for (auto extent : c_count) {
    if (extent && size > std::numeric_limits<std::size_t>::max() / extent)
      return NC_EEDGE;
    size *= extent;
  }
  if (size > std::size_t(nvalues))
    return NC_EEDGE;
  if constexpr (std::is_same_v<Value, float>)
    return nc_put_vara_float(ncid, varid - 1, c_start.data(), c_count.data(), values);
  else if constexpr (std::is_same_v<Value, double>)
    return nc_put_vara_double(ncid, varid - 1, c_start.data(), c_count.data(), values);
  else if constexpr (std::is_same_v<Value, int>)
    return nc_put_vara_int(ncid, varid - 1, c_start.data(), c_count.data(), values);
  else
    return nc_put_vara_text(ncid, varid - 1, c_start.data(), c_count.data(), values);
}

} // namespace

int open(std::string_view path, int& ncid) {
  // Open a read-only file without adding a C++ handle registry or ownership cache.
  return nc_open(c_name(path).c_str(), NC_NOWRITE, &ncid);
}

int close(int ncid) {
  // Release the caller's explicit handle and return the original library status.
  return nc_close(ncid);
}

int inq_varid(int ncid, std::string_view name, int& varid) {
  // Resolve a name and retain Native's one-based variable identifier.
  int c_varid;
  const int status = nc_inq_varid(ncid, c_name(name).c_str(), &c_varid);
  if (status == NC_NOERR)
    varid = c_varid + 1;
  return status;
}

int inq_dimid(int ncid, std::string_view name, int& dimid) {
  // Resolve a name and retain Native's one-based dimension identifier.
  int c_dimid;
  const int status = nc_inq_dimid(ncid, c_name(name).c_str(), &c_dimid);
  if (status == NC_NOERR)
    dimid = c_dimid + 1;
  return status;
}

int inquire_variable(int ncid, int varid, int& xtype, int& ndims, int* dimids) {
  // Report the stored type and fastest-first dimension IDs used by nf90 callers.
  std::array<int, NC_MAX_VAR_DIMS> c_dimids;
  nc_type c_type;
  const int status = nc_inq_var(ncid, varid - 1, nullptr, &c_type, &ndims,
                               c_dimids.data(), nullptr);
  if (status == NC_NOERR) {
    xtype = c_type;
    for (int dim = 0; dim < ndims; ++dim)
      dimids[dim] = c_dimids[ndims - dim - 1] + 1;
  }
  return status;
}

int inquire_dimension(int ncid, int dimid, char* name, int name_length, int& length) {
  // Return a dimension length and fixed-width, blank-padded Fortran name.
  std::array<char, NC_MAX_NAME + 1> c_dimname;
  std::size_t c_length;
  const int status = nc_inq_dim(ncid, dimid - 1, c_dimname.data(), &c_length);
  if (status == NC_NOERR) {
    std::fill_n(name, name_length, ' ');
    const std::string_view dimname(c_dimname.data());
    std::copy_n(dimname.data(), std::min(dimname.size(), std::size_t(name_length)), name);
    if (c_length > std::size_t(std::numeric_limits<int>::max()))
      return NC_ERANGE;
    length = int(c_length);
  }
  return status;
}

int get_att(int ncid, int varid, std::string_view name, char* text, int text_length) {
  // Read character metadata, retaining a bounded, blank-padded caller result.
  std::fill_n(text, text_length, ' ');
  std::size_t length;
  const auto attribute_name = c_name(name);
  int status = nc_inq_attlen(ncid, varid - 1, attribute_name.c_str(), &length);
  if (status != NC_NOERR)
    return status;
  // Read into the attribute's actual extent, then retain the fixed-length
  // Fortran character result without a terminating NUL or a buffer overrun.
  std::vector<char> value(length);
  status = nc_get_att_text(ncid, varid - 1, attribute_name.c_str(), value.data());
  if (status == NC_NOERR)
    std::copy_n(value.data(), std::min(length, std::size_t(text_length)), text);
  return status;
}

int get_var(int ncid, int varid, int nvalues, int ndims,
            const int* start, const int* count, float* values) {
  // Read Native default-real input without an intermediate double conversion.
  return get_values(ncid, varid, nvalues, ndims, start, count, values);
}

int get_var(int ncid, int varid, int nvalues, int ndims,
            const int* start, const int* count, double* values) {
  // Read double core-real input using the same library conversion as nf90.
  return get_values(ncid, varid, nvalues, ndims, start, count, values);
}

int get_var(int ncid, int varid, int nvalues, int* values) {
  // Match nf90's rank-one default count, selecting one along any trailing axes.
  int ndims;
  const int status = nc_inq_varndims(ncid, varid - 1, &ndims);
  if (status != NC_NOERR)
    return status;
  std::vector<int> start(ndims, 1), count(ndims, 1);
  if (ndims)
    count[0] = nvalues;
  return get_values(ncid, varid, nvalues, ndims, start.data(), count.data(), values);
}

int get_var(int ncid, int varid, int& value) {
  // Match nf90 scalar dispatch: select the first value even in a ranked variable.
  int ndims;
  const int status = nc_inq_varndims(ncid, varid - 1, &ndims);
  if (status != NC_NOERR)
    return status;
  std::vector<std::size_t> start(ndims, 0);
  return nc_get_var1_int(ncid, varid - 1, start.data(), &value);
}

std::string strerror(int status) {
  // Supply the unchanged NetCDF error text to the caller's existing diagnostic.
  return nc_strerror(status);
}

int create(std::string_view path, int mode, int& ncid) {
  // Create the caller's file using exactly the requested NetCDF storage mode.
  return nc_create(c_name(path).c_str(), mode, &ncid);
}

int redef(int ncid) {
  // Reenter definition mode for sample variables added after initialization.
  return nc_redef(ncid);
}

int enddef(int ncid) {
  // Commit definitions before the caller writes coordinate or sampled values.
  return nc_enddef(ncid);
}

int def_dim(int ncid, std::string_view name, int length, int& dimid) {
  // Define a Native-named dimension; zero selects the unlimited record axis.
  if (length < 0)
    return NC_EDIMSIZE;
  int c_dimid;
  const int status = nc_def_dim(ncid, c_name(name).c_str(), std::size_t(length), &c_dimid);
  if (status == NC_NOERR)
    dimid = c_dimid + 1;
  return status;
}

int def_var(int ncid, std::string_view name, int xtype, const std::vector<int>& dimids,
            int& varid) {
  // Reverse Native fastest-first dimension IDs, retaining the library's type.
  std::vector<int> c_dimids;
  c_dimids.reserve(dimids.size());
  for (auto id = dimids.rbegin(); id != dimids.rend(); ++id)
    c_dimids.push_back(*id - 1);
  int c_varid;
  const int status = nc_def_var(ncid, c_name(name).c_str(), static_cast<nc_type>(xtype),
                                int(c_dimids.size()), c_dimids.data(), &c_varid);
  if (status == NC_NOERR)
    varid = c_varid + 1;
  return status;
}

int def_var_fill(int ncid, int varid, int no_fill, double fill_value) {
  // Stats defines double variables and supplies their explicit zero fill value.
  // nc_def_var_fill accepts an untyped pointer, so enforce that storage contract.
  nc_type xtype;
  const int status = nc_inq_vartype(ncid, varid - 1, &xtype);
  if (status != NC_NOERR)
    return status;
  if (xtype != NC_DOUBLE)
    return NC_EBADTYPE;
  return nc_def_var_fill(ncid, varid - 1, no_fill, &fill_value);
}

int put_att(int ncid, int varid, std::string_view name, std::string_view text) {
  // Attribute text length is explicit; only the caller applies Native TRIM.
  return nc_put_att_text(ncid, varid - 1, c_name(name).c_str(), text.size(), text.data());
}

int put_att(int ncid, int varid, std::string_view name, float value) {
  // Preserve the Native single-precision attribute type and value.
  return nc_put_att_float(ncid, varid - 1, c_name(name).c_str(), NC_FLOAT, 1, &value);
}

int put_att(int ncid, int varid, std::string_view name, double value) {
  // Preserve the Native double-precision attribute type and value.
  return nc_put_att_double(ncid, varid - 1, c_name(name).c_str(), NC_DOUBLE, 1, &value);
}

int put_att(int ncid, int varid, std::string_view name, int value) {
  // Write integer metadata without converting it through a floating-point kind.
  return nc_put_att_int(ncid, varid - 1, c_name(name).c_str(), NC_INT, 1, &value);
}

int put_var(int ncid, int varid, const float* values, int nvalues, const std::vector<int>& start,
            const std::vector<int>& count) {
  // Let the same C-library conversion as nf90 write packed core-real values.
  return put_values(ncid, varid, values, nvalues, start, count);
}

int put_var(int ncid, int varid, const double* values, int nvalues, const std::vector<int>& start,
            const std::vector<int>& count) {
  // Write double stats buffers without an intermediate core-real conversion.
  return put_values(ncid, varid, values, nvalues, start, count);
}

int put_var(int ncid, int varid, const int* values, int nvalues, const std::vector<int>& start,
            const std::vector<int>& count) {
  // Preserve SILHS integer mixture indices at the typed C-library boundary.
  return put_values(ncid, varid, values, nvalues, start, count);
}

int put_var(int ncid, int varid, const char* values, int nvalues, const std::vector<int>& start,
            const std::vector<int>& count) {
  // Write already packed fixed-width character fields, including trailing blanks.
  return put_values(ncid, varid, values, nvalues, start, count);
}

} // namespace clubb_kokkos::netcdf
