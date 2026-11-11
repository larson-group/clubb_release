// Read and write NetCDF data with Native Fortran index conventions.
//
// File handles and status codes are NetCDF-C values. Variable/dimension IDs and
// hyperslab starts are one-based; dimension lists put the fastest axis first.
// These call-local services convert that boundary once, leaving file meaning,
// Kokkos staging, missing-value policy and handle ownership with their callers.
//
// open obtains a read-only handle; close releases any caller-owned handle.
// inq_varid/inq_dimid resolve names; inquire_variable/inquire_dimension return Native-ordered metadata.
// get_att reads a blank-padded character result; get_var reads a typed packed
// host slab, integer vector or scalar. strerror returns the library diagnostic.
// create/redef/enddef and def_dim/def_var/def_var_fill define caller-owned files.
// put_att writes explicit text or typed scalar metadata. put_var writes packed
// host values; explicit starts/counts retain Native fastest-first ordering.
// Native global attributes use nc_global (zero); variable/dimension IDs start at
// one. The library owns storage conversion; callers own precision and fill values.
#pragma once

#include <netcdf.h>
#include <string>
#include <string_view>
#include <vector>

namespace clubb_kokkos::netcdf {

inline constexpr int nc_noerr = NC_NOERR;
inline constexpr int nc_int = NC_INT, nc_float = NC_FLOAT, nc_double = NC_DOUBLE;
inline constexpr int nc_char = NC_CHAR, nc_global = 0;
inline constexpr int nc_clobber = NC_CLOBBER, nc_unlimited = NC_UNLIMITED;
inline constexpr int nc_max_var_dims = NC_MAX_VAR_DIMS, nc_max_name = NC_MAX_NAME;

int open(std::string_view path, int& ncid);
int close(int ncid);
int inq_varid(int ncid, std::string_view name, int& varid);
int inq_dimid(int ncid, std::string_view name, int& dimid);
int inquire_variable(int ncid, int varid, int& xtype, int& ndims, int* dimids);
int inquire_dimension(int ncid, int dimid, char* name, int name_length, int& length);
int get_att(int ncid, int varid, std::string_view name, char* text, int text_length);

// Explicit start/count vectors follow Native dimension order. The destination
// is a packed host buffer; nvalues is its capacity, including any unread tail.
int get_var(int ncid, int varid, int nvalues, int ndims,
            const int* start, const int* count, float* values);
int get_var(int ncid, int varid, int nvalues, int ndims,
            const int* start, const int* count, double* values);

// Match the nf90 scalar and rank-one INTEGER overloads without an explicit slab.
int get_var(int ncid, int varid, int nvalues, int* values);
int get_var(int ncid, int varid, int& value);
std::string strerror(int status);

int create(std::string_view path, int mode, int& ncid);
int redef(int ncid);
int enddef(int ncid);
int def_dim(int ncid, std::string_view name, int length, int& dimid);
int def_var(int ncid, std::string_view name, int xtype, const std::vector<int>& dimids, int& varid);
int def_var_fill(int ncid, int varid, int no_fill, double fill_value);
int put_att(int ncid, int varid, std::string_view name, std::string_view text);
int put_att(int ncid, int varid, std::string_view name, float value);
int put_att(int ncid, int varid, std::string_view name, double value);
int put_att(int ncid, int varid, std::string_view name, int value);

// Values use packed Native order. Omitted starts/counts write the full fixed
// variable; pass source-shaped counts for partial arrays or unlimited records.
int put_var(int ncid, int varid, const float* values, int nvalues,
            const std::vector<int>& start = {}, const std::vector<int>& count = {});
int put_var(int ncid, int varid, const double* values, int nvalues,
            const std::vector<int>& start = {}, const std::vector<int>& count = {});
int put_var(int ncid, int varid, const int* values, int nvalues, const std::vector<int>& start = {},
            const std::vector<int>& count = {});
int put_var(int ncid, int varid, const char* values, int nvalues,
            const std::vector<int>& start = {}, const std::vector<int>& count = {});

} // namespace clubb_kokkos::netcdf
