"""Read profiles from netCDF files: JAX port of ``input_netcdf.F90``.

Description:
  Opens a netCDF file, such as a CLUBB stats file or an LES (large-eddy
  simulation) benchmark file, and reads one vertical profile of a variable at
  a time.

Usage:
  ncf, l_error = open_netcdf_read(test_variable, path) -- open the file and
    learn its vertical levels and output times from the dimensions of
    test_variable.
  x, l_error = get_netcdf_var(ncf, varname, itime, l_convert_to_MKS) -- read
    the profile of varname at the one-based output record itime.
  close_netcdf_read(ncf) -- close the file.

  The module variables clubb_day, clubb_month and clubb_year hold the start
  date of the simulation; the driver sets them before restart files are read.

Understanding:
  A variable must have a time dimension ("T", "t" or "time") and a vertical
  dimension ("Z", "z", "zt", "zm", "altitude", "height", "lev", "lh_zt",
  "rad_zt" or "rad_zm"). It must also have either a column dimension ("col" or
  "column") or X and Y dimensions of length one. Values are read raw, without
  netCDF4's masking or scaling, as the Fortran netCDF library reads them. With
  l_convert_to_MKS, values in g/kg are divided by 1000 and values in K/day by
  the seconds per day.

  JAX adaptations:
  - get_netcdf_var returns every column, with shape (column, level). The
    Fortran reads the first column only.
  - Surface variables without a vertical dimension, or with the length-one
    "sfc" dimension of CLUBB stats files, are read as one level at height 0.
    The Fortran rejects them, so a Fortran restart that reads the soil and
    vegetation temperatures fails.
  - close_netcdf_read does not print the Fortran's "closing" message.
"""

import sys

import numpy as np

from clubb_jax.src.CLUBB_core.calendar import gregorian2julian_date
from clubb_jax.src.CLUBB_core.constants_clubb import (
    g_per_kg,
    sec_per_day,
    sec_per_hr,
    sec_per_min,
)
from clubb_jax.src.Input_fields.stat_file_module import stat_file

from netCDF4 import Dataset


clubb_day = 1    # Start time the of simulation
clubb_month = 1
clubb_year = 2000

_Z_NAMES = ("Z", "z", "zt", "zm", "altitude", "height", "lev", "lh_zt", "rad_zt", "rad_zm")
_X_NAMES = ("X", "x", "longitude", "lon")
_Y_NAMES = ("Y", "y", "latitude", "lat")
_T_NAMES = ("T", "t", "time")
_COL_NAMES = ("col", "column")
# JAX adaptation: the surface dimension of CLUBB stats files.
_SFC_NAMES = ("sfc",)


# -----------------------------------------------------------------------
def open_netcdf_read(test_variable, path):
    """Open a netCDF data set in read-only mode.

    Args:
        test_variable: Variable used to determine the dimensions in the file.
        path: The file name.

    Returns:
        ncf: stat_file describing the netCDF file.
        l_error: Whether the file could not be opened or described.
    """
    # Initialize l_error to false
    l_error = False
    ncf = stat_file()

    try:
        ncf.iounit = Dataset(str(path))
    except OSError as exc:
        print("Error opening netCDF file: " + str(path), file=sys.stderr)
        print(exc, file=sys.stderr)
        l_error = True
        return ncf, l_error

    # nf90_get_var reads raw values; netCDF4 must not mask valid zero fill
    # values or apply scale_factor/add_offset attributes implicitly.
    ncf.iounit.set_auto_mask(False)
    ncf.iounit.set_auto_scale(False)

    if test_variable not in ncf.iounit.variables:
        print(f"Variable not found: {test_variable}", file=sys.stderr)
        l_error = True
        return ncf, l_error

    var = ncf.iounit.variables[test_variable]
    ndims = var.ndim
    dtype = np.dtype(var.dtype)

    # JAX adaptation: surface variables have one dimension fewer.
    l_sfc_var = not any(name in _Z_NAMES for name in var.dimensions)
    if ( ( ndims not in (3, 4) and not ( l_sfc_var and ndims == 2 ) )
         or not ( dtype.kind == "f" and dtype.itemsize in (4, 8) ) ):
        print("input_netcdf.open_netcdf_read : The netCDF data doesn't "
              "conform to expected precision, shape, or dimensions", file=sys.stderr)
        l_error = True
        return ncf, l_error

    xdim = -1
    ydim = -1
    cdim = -1
    ncf.iz = -1
    ncf.ntimes = -1
    zname = ""
    time_name = ""

    for dim_name in var.dimensions:

        itmp = len(ncf.iounit.dimensions[dim_name])

        if dim_name in _Z_NAMES:
            ncf.ia = 1
            ncf.iz = itmp
            ncf.AltDimId = dim_name
            zname = dim_name
        elif dim_name in _X_NAMES:
            xdim = itmp
            ncf.LongDimId = dim_name
        elif dim_name in _Y_NAMES:
            ydim = itmp
            ncf.LatDimId = dim_name
        elif dim_name in _T_NAMES:
            ncf.ntimes = itmp
            ncf.TimeDimId = dim_name
            time_name = dim_name
        elif dim_name in _COL_NAMES:
            cdim = itmp
        elif dim_name in _SFC_NAMES and itmp == 1:
            # JAX adaptation: one surface level.
            pass
        else:
            l_error = True
            return ncf, l_error

    if l_sfc_var:
        # JAX adaptation: one surface level at height 0.
        ncf.ia = 1
        ncf.iz = 1

    if ncf.iz < 1 or ncf.ntimes < 1:
        l_error = True
        return ncf, l_error

    if cdim < 0:
        if ydim != 1 or xdim != 1:
            print("input_netcdf.open_netcdf_read : The netCDF data doesn't "
                  "conform to the expected X or Y dimension", file=sys.stderr)
            l_error = True
            return ncf, l_error
    elif cdim < 1:
        print("input_netcdf.open_netcdf_read : The netCDF data doesn't "
              "conform to the expected col dimension", file=sys.stderr)
        l_error = True
        return ncf, l_error

    if not l_error:

        # Determine the altitudes
        if l_sfc_var:
            ncf.z = np.zeros(1)
        else:
            if zname not in ncf.iounit.variables:
                print(f"Variable not found: {zname}", file=sys.stderr)
                l_error = True
                return ncf, l_error
            ncf.AltVarId = zname
            ncf.z = np.asarray(ncf.iounit.variables[zname][:ncf.iz], dtype=np.float64)

        if time_name not in ncf.iounit.variables:
            print(f"Variable not found: {time_name}", file=sys.stderr)
            l_error = True
            return ncf, l_error
        ncf.TimeVarId = time_name

        time = getattr(ncf.iounit.variables[time_name], "units", None)
        if not isinstance(time, str):
            print("Attribute not found: units", file=sys.stderr)
            l_error = True
            return ncf, l_error

        # Read starting time from the "units" attribute of the netcdf file
        time = time.strip()
        length = len(time)

        if length < 29:
            print("The NetCDF file does not have a proper time unit "
                  "specification. The \"units\" attribute for the "
                  "time variable must be in the form:", file=sys.stderr)
            print("TIMEUNITS since YYYY-MM-DD HH:MM:SS.S", file=sys.stderr)
            l_error = True
            return ncf, l_error

        # The following blocks of code are intended to parse date string from
        # input netcdf files.  The CLUBB output files are in the format
        # TIMEUNITS since YYYY-MM-DD HH:MM:SS.S, while other files, like some
        # COAMPS files, may be formatted TIMEUNITS since YYYY-M-D HH:MM:SS (i.e.
        # with single-digit month or day, and without the ".S" at the end.  Hence
        # we need to parse these string carefully to correctly compare dates.
        # Python indices are zero-based, one less than the Fortran index().

        # Search forward
        first_hyphen_idx = time.find("-")
        # Search backward
        second_hyphen_idx = time.rfind("-")
        # Search forward
        first_colon_idx = time.find(":")
        # Search backward
        second_colon_idx = time.rfind(":")

        try:
            # Year
            netcdf_year = int(time[first_hyphen_idx - 4:first_hyphen_idx])

            # Month
            # If difference between hyphen indices is 2, then the month is a single digit.
            if second_hyphen_idx - first_hyphen_idx == 2:
                dummy_string = time[first_hyphen_idx + 1:first_hyphen_idx + 2]
            else:
                dummy_string = time[first_hyphen_idx + 1:first_hyphen_idx + 3]
            netcdf_month = int(dummy_string)

            # Day
            # If difference between second hyphen and first colon indices is 5, then the
            # day is a single digit.
            if first_colon_idx - second_hyphen_idx == 5:
                dummy_string = time[second_hyphen_idx + 1:second_hyphen_idx + 2]
            else:
                dummy_string = time[second_hyphen_idx + 1:second_hyphen_idx + 3]
            netcdf_day = int(dummy_string)

            # Time
            hours = float(time[first_colon_idx - 2:first_colon_idx])
            minutes = float(time[first_colon_idx + 1:first_colon_idx + 3])
            seconds = float(time[second_colon_idx + 1:second_colon_idx + 3])
        except ValueError:
            print("Unable to read the date from the time units: " + time, file=sys.stderr)
            l_error = True
            return ncf, l_error

        ncf.year = netcdf_year
        ncf.month = netcdf_month
        ncf.day = netcdf_day

        # Compute delta_d, the difference in days between the netcdf file's start
        # date and CLUBB's start date (netcdf_date - clubb_date). This difference
        # should be added to ncf%time to fix an inconsistency that is caused when
        # the netcdf date differs from CLUBB's date.
        # July 2013 Eric Raut
        clubb_idate = int(gregorian2julian_date(clubb_day, clubb_month, clubb_year))
        netcdf_idate = int(gregorian2julian_date(netcdf_day, netcdf_month, netcdf_year))
        delta_d = float(netcdf_idate - clubb_idate)

        ncf.time = (hours * float(sec_per_hr)) + \
                   (minutes * float(sec_per_min)) \
                   + seconds + (delta_d * float(sec_per_day))

        # Figure out what units Time is in so dtwrite can be set correctly
        time_units = time.split(" ", 1)[0]
        if time_units == "hours":
            multiplier = float(sec_per_hr)

        elif time_units == "minutes":
            multiplier = float(sec_per_min)

        elif time_units == "seconds":
            multiplier = 1.0

        else:
            l_error = True
            return ncf, l_error

        if ncf.ntimes < 2:
            print("NetCDF: Start+count exceeds dimension bound", file=sys.stderr)
            l_error = True
            return ncf, l_error
        write_times = np.asarray(ncf.iounit.variables[time_name][:2], dtype=np.float64)

        ncf.dtwrite = (write_times[1] - write_times[0]) * multiplier

    return ncf, l_error


# ----------------------------------------------------------------------
def get_netcdf_var(ncf, varname, itime, l_convert_to_MKS):
    """Read in values from a netCDF file.

    Assumptions:
      We assume a NF90_FLOAT the same as a kind=4 real in Fortran, and
      That the variables will obey COARDS conventions and have 4 dimensions.

    Args:
        ncf: NetCDF file structure from open_netcdf_read.
        varname: The variable name as it occurs in the control file.
        itime: Obtain variable varname at time itime (one-based record).
        l_convert_to_MKS: Convert variables to MKS upon reading.

    Returns:
        x: Result variable, shape (column, ncf.iz); one column for files with
           X and Y dimensions.
        l_error: Whether the variable could not be read.
    """
    # Initialize l_error to false
    l_error = False
    x = None

    if varname not in ncf.iounit.variables:
        print(f"Variable not found: {varname}", file=sys.stderr)
        l_error = True
        return x, l_error

    var = ncf.iounit.variables[varname]
    ndims = var.ndim
    dtype = np.dtype(var.dtype)

    # JAX adaptation: surface variables have one dimension fewer.
    l_sfc_var = not any(name in _Z_NAMES for name in var.dimensions)
    if ( ( ndims not in (3, 4) and not ( l_sfc_var and ndims == 2 ) )
         or not ( dtype.kind == "f" and dtype.itemsize in (4, 8) ) ):
        print("input_netcdf.get_var : The netCDF data doesn't "
              "conform to expected precision, shape, or dimensions", file=sys.stderr)
        l_error = True
        return x, l_error

    if itime > ncf.ntimes:
        print("input_netcdf.get_var: The time specified exceeds the netCDF data", file=sys.stderr)
        l_error = True
        return x, l_error

    idx_x = -1
    idx_y = -1
    idx_z = -1
    idx_t = -1
    idx_col = -1
    xdim = -1
    ydim = -1
    cdim = -1

    for i, dim_name in enumerate(var.dimensions):
        dim_len = len(ncf.iounit.dimensions[dim_name])

        if dim_name in _Z_NAMES:
            idx_z = i
        elif dim_name in _T_NAMES:
            idx_t = i
        elif dim_name in _X_NAMES:
            idx_x = i
            xdim = dim_len
        elif dim_name in _Y_NAMES:
            idx_y = i
            ydim = dim_len
        elif dim_name in _COL_NAMES:
            idx_col = i
            cdim = dim_len
        elif dim_name in _SFC_NAMES and dim_len == 1:
            # JAX adaptation: one surface level.
            idx_z = i
        else:
            l_error = True
            return x, l_error

    if ( idx_z < 0 and not l_sfc_var ) or idx_t < 0:
        l_error = True
        return x, l_error

    if idx_col < 0:
        if idx_x < 0 or idx_y < 0 or xdim != 1 or ydim != 1:
            l_error = True
            return x, l_error
    elif cdim < 1:
        l_error = True
        return x, l_error

    # Obtain a vertical column at time itime
    # JAX adaptation: keep every column; the Fortran count is 1 in the column dimension.
    start = []
    for i in range(ndims):
        if i == idx_t:
            start.append(itime - 1)
        elif i == idx_z:
            start.append(slice(0, ncf.iz))
        elif i == idx_col:
            start.append(slice(None))
        else:
            start.append(0)
    try:
        values = np.asarray(var[tuple(start)], dtype=np.float64)
    except (IndexError, ValueError, RuntimeError) as exc:
        print(exc, file=sys.stderr)
        l_error = True
        return x, l_error

    # Return the kept axes as (column, level).
    if idx_col >= 0 and idx_z >= 0:
        if idx_z < idx_col:
            values = values.T
    elif idx_col >= 0:
        # Surface variable without a vertical dimension: (column,)
        values = values[:, np.newaxis]
    else:
        # File with X and Y dimensions: (level,)
        values = values[np.newaxis, :]
    x = values

    if l_convert_to_MKS:

        units = getattr(var, "units", None)

        if not isinstance(units, str):
            print("Attribute not found: units", file=sys.stderr)
            l_error = True
            return x, l_error

        else:
            units = units.strip()
            if units == "g/kg":
                x = x / g_per_kg

            elif units == "W/m2":
                l_error = True
                print("get_netcdf_var: Unable to convert variables of this type to MKS units",
                      file=sys.stderr)
                return x, l_error

            elif units == "K/day":
                x = x / sec_per_day

            else:
                # Do nothing
                pass

    return x, l_error


# -----------------------------------------------------------------------
def close_netcdf_read(ncf):
    """Close a previously opened netCDF file."""
    # Close file
    try:
        ncf.iounit.close()
    except (OSError, RuntimeError) as exc:
        print("Error closing a netCDF file", file=sys.stderr)
        print(exc, file=sys.stderr)
        raise RuntimeError("Fatal error.") from exc
