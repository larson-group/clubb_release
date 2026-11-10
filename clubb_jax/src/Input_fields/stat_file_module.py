"""Describe an open netCDF statistics file: JAX port of ``stat_file_module.F90``.

Description:
  Contains the derived type that describes the contents and location of a
  netCDF file, such as a CLUBB stats file or an LES (large-eddy simulation)
  benchmark file.

Usage:
  stat_file -- what input_netcdf.open_netcdf_read learns about a file: the open
    dataset, its vertical levels and its output times. Pass it to
    input_netcdf.get_netcdf_var to read profiles and to
    input_netcdf.close_netcdf_read to close the file.

Understanding:
  The Fortran keeps netCDF dimension and variable ids (AltDimId, TimeVarId,
  ...). Python's netCDF4 package looks dimensions and variables up by name, so
  these fields hold names instead.

  Not ported: the grid_avg_variable and samples_of_variable types, and the
  output-only fields of stat_file (fname, fdir, nrecord, l_defined, l_byte_swapped,
  the SILHS sample fields and the variable lists). The Fortran uses them to
  write stats files; the JAX stats output is written by stats_netcdf.StatsWriter.
"""

from dataclasses import dataclass, field

import numpy as np


# These are used in a 2D or 3D host model to output multiple columns
# Set clubb_i and clubb_j according to the column within the host model;
# The indices must not exceed nlon (for i) or nlat (for j).
#
# !!!!!!!!! Unused Warning !!!!!!!!!!!!
# Currently we do not save stats output with lat and lon indexing, and
# because of this, we hardcode these as 1. When implementing the multiple
# column stats outputting in the future, we should collapse the lat and
# lon dimension, and use a generic "column" dimension instead.
clubb_i = 1
clubb_j = 1


# Structure to hold the description of a NetCDF output file
@dataclass
class stat_file:
    """Description of a netCDF file opened by input_netcdf.open_netcdf_read."""

    # File information
    iounit: object = None   # The open netCDF4.Dataset

    ntimes: int = -1        # Number of times written

    # NetCDF datafile dimension and variable names
    LatDimId: str = ""
    LongDimId: str = ""
    AltDimId: str = ""
    TimeDimId: str = ""
    AltVarId: str = ""
    TimeVarId: str = ""

    # Grid information
    ia: int = 1             # Vertical extent (one-based, inclusive)
    iz: int = -1
    nlat: int = 1           # The number of points in the X and Y
    nlon: int = 1

    z: np.ndarray = field(default_factory=lambda: np.zeros(0))  # Height of vertical levels [m]

    # Time information
    day: int = 0            # Date of starting time
    month: int = 0
    year: int = 0

    dtwrite: float = 0.0    # Interval between output    [Seconds]

    time: float = 0.0       # Start time                 [Seconds]
