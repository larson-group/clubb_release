// Describe an input NetCDF file and its variable metadata for one CLUBB case.
//
// These Native-named schemas preserve source field order and precision. They
// describe input/restart files; the retained output statistics backend still
// owns its separate stats_type. Views share allocation ownership like the
// Native pointer fields, without copying array elements. File closure remains
// explicit in close_netcdf_read; destruction does not silently close an ncid.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include <string>
#include <vector>

namespace clubb_kokkos {
// Zero-based host-column coordinates; NetCDF IDs retain their file protocol.
inline constexpr int clubb_i = 0, clubb_j = 0;

struct grid_avg_variable {
  Kokkos::View<double***, Kokkos::LayoutLeft> ptr; // stat_rknd array alias.
  std::string name;                                // Native fixed30 variable name.
  std::string description;                         // Native fixed100 description.
  std::string units;                               // Native fixed25 units.
  int indx;                                        // NetCDF variable ID.
  bool l_silhs; // Sample once per SILHS sample, rather than timestep.
};

struct samples_of_variable {
  Kokkos::View<double****, Kokkos::LayoutLeft> ptr; // stat_rknd array alias.
  std::string name;                                 // Native fixed30 variable name.
  std::string description;                          // Native fixed100 description.
  std::string units;                                // Native fixed25 units.
  int indx;                                         // NetCDF variable ID.
  bool l_silhs; // Sample once per SILHS sample, rather than timestep.
};

struct stat_file {
  std::string fname, fdir; // Native fixed200 file name without suffix and directory.
  int iounit;              // Open NetCDF dataset ID; closed explicitly by caller.
  int nrecord, ntimes;     // Number of written records and time coordinates.
  bool l_defined, l_byte_swapped;
  int SampDimId, LatDimId, LongDimId, AltDimId, TimeDimId;
  int SampVarId, LatVarId, LongVarId, AltVarId, TimeVarId;
  int ia, iz;                          // Inclusive zero-based vertical bounds; missing iz is -2.
  int nlat, nlon;                      // Horizontal extents.
  int nsamp = 0;                       // Native initialized SILHS sample count.
  Real1D z;                            // Vertical heights [m].
  int day, month, year;                // Starting date.
  Real1D lat_vals, lon_vals, samp_idx; // Degrees N/E and SILHS sample index.
  real_t dtwrite;                      // Output interval [s].
  double time;                         // time_precision start time [s], always double.
  int nvar;                            // Number of file variables.
  std::vector<grid_avg_variable> grid_avg_var;
  std::vector<samples_of_variable> samples_of_var;
};
} // namespace clubb_kokkos
