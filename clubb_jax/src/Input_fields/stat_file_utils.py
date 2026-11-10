"""Average and interpolate netCDF profiles: JAX port of ``stat_file_utils.F90``.

Description:
  Reads profiles from netCDF statistics files, such as CLUBB stats files or
  LES (large-eddy simulation) output, averages them over output records and
  interpolates them onto the CLUBB grid.

Usage:
  Reading profiles:
    stat_file_average(...) -- average a variable over records t1 to t2 and
      interpolate it to the given heights; the restart reader uses it with
      t1 = t2.
    stat_file_average_interval(...) -- average a variable over several record
      intervals, weighting each by its length.
  File information:
    stat_file_num_vertical_levels(varname, filename) -- number of levels.
    stat_file_vertical_levels(varname, filename, nz) -- heights of the levels.
    l_netcdf_file(filename) -- whether the file name ends in ".nc"; stops
      otherwise.
  Matching grids:
    CLUBB_levels_within_LES_domain(fread_var, nz, ngrdcol, CLUBB_grid) -- each
      column's lowest and highest CLUBB levels inside the height range of the
      file.
    LES_grid_to_CLUBB_grid(fread_var, CLUBB_height) -- the file level equal
      to one CLUBB height, or the two file levels around it.

Understanding:
  Level and record indices are one-based, as in the Fortran. A CLUBB level
  below or above the levels of the file is set to +infinity, except that with
  l_spec_bound_cond the lowest CLUBB level (a ghost level below the ground)
  takes the value of the lowest file level.

  JAX adaptation: file profiles keep a leading column axis, (column, level),
  as returned by input_netcdf.get_netcdf_var. The Fortran reads one file
  column; here a file may also hold one profile per output column.
"""

import sys

import numpy as np

from clubb_jax.src.CLUBB_core.constants_clubb import eps
from clubb_jax.src.CLUBB_core.interpolation import lin_interpolate_two_points
from clubb_jax.src.Input_fields.input_netcdf import (
    close_netcdf_read,
    get_netcdf_var,
    open_netcdf_read,
)


# ----------------------------------------------------------------------
def stat_file_average(filename, out_nz, ngrdcol,
                      t1, t2, out_heights,
                      variable_name,
                      npower, l_spec_bound_cond):
    """Read each shared file profile once and average it on every column grid.
    Average a netCDF file variable over the interval t1 to t2

    References:
      None

    Args:
        filename: Name of the file.
        out_nz: Number of output vertical levels [-].
        ngrdcol: Number of output columns [-].
        t1: Beginning timestep to look at (one-based record).
        t2: Ending timestep to look at (one-based record).
        out_heights: Heights of the output grid, shape (ngrdcol, out_nz) [m].
        variable_name: Name of the variable to read in.
        npower: Exponent operator, must be 1 or 2.
        l_spec_bound_cond: Special boundary condition for when the lowest
            model level is a ghost point comes from a fine vertical resolution
            run, shape (ngrdcol,).  E.g.
                Stat file levels    Model levels
                      50m                50m
                      25m
                     -25m
                                        -50m
            So we need to set the lowest model level to the -25m value.

    Returns:
        stat_file_average: The average, shape (ngrdcol, out_nz).
        l_error: error status for this function.

    JAX adaptation: the file may hold one profile, shared by every column, or
    one profile per column. Any other column count is a read error.
    """
    out_heights = np.asarray(out_heights, dtype=np.float64)
    l_spec_bound_cond = np.broadcast_to(np.asarray(l_spec_bound_cond, dtype=bool), (ngrdcol,))

    # Initialize variables
    l_error = False

    num_timesteps = (t2 - t1) + 1
    stat_file_average = np.zeros((ngrdcol, out_nz))

    # Open netCDF file
    faverage, l_error = open_netcdf_read(variable_name, filename)

    if l_error:
        if faverage.iounit is not None:
            close_netcdf_read(faverage)
        return stat_file_average, l_error

    # Determine variable size
    file_nz = faverage.z.size

    # Determine the interpolation needed on each column grid.
    l_interpolate = np.zeros(ngrdcol, dtype=bool)
    for i in range(1, ngrdcol + 1):
        if out_nz != file_nz:
            l_interpolate[i - 1] = True
        else:
            l_interpolate[i - 1] = np.any(np.abs(faverage.z - out_heights[i - 1, :]) >
                                          np.abs(faverage.z + out_heights[i - 1, :]) * eps / 2)

    k_lowest_input, k_highest_input = CLUBB_levels_within_LES_domain(
        faverage, out_nz, ngrdcol, out_heights,
    )

    exact_lev_idx = np.zeros((ngrdcol, out_nz), dtype=int)
    lower_lev_idx = np.zeros((ngrdcol, out_nz), dtype=int)
    upper_lev_idx = np.zeros((ngrdcol, out_nz), dtype=int)
    l_lin_int = np.zeros((ngrdcol, out_nz), dtype=bool)

    for k in range(1, out_nz + 1):
        for i in range(1, ngrdcol + 1):
            if k < k_lowest_input[i - 1] or k > k_highest_input[i - 1]:
                continue
            (
                exact_lev_idx[i - 1, k - 1], lower_lev_idx[i - 1, k - 1],
                upper_lev_idx[i - 1, k - 1], l_lin_int[i - 1, k - 1],
            ) = LES_grid_to_CLUBB_grid(faverage, out_heights[i - 1, k - 1])

    # Read in variables from netCDF file
    for t in range(t1, t2 + 1):

        l_convert_to_MKS = True
        file_variable, l_error = get_netcdf_var(
            faverage, variable_name, t, l_convert_to_MKS,
        )

        if l_error:
            print("stat_file_average: get_var failed for "
                  + variable_name + " in " + str(filename) + " at time=", t, file=sys.stderr)
            close_netcdf_read(faverage)
            return stat_file_average, l_error

        if file_variable.shape[0] not in (1, ngrdcol):
            print(f"stat_file_average: {variable_name} in {filename} has "
                  f"{file_variable.shape[0]} columns, not 1 or {ngrdcol}", file=sys.stderr)
            l_error = True
            close_netcdf_read(faverage)
            return stat_file_average, l_error

        interp_variable = np.zeros((ngrdcol, out_nz))

        # Interpolate within each column's domain and retain the legacy boundary policy.
        for k in range(1, out_nz + 1):
            for i in range(1, ngrdcol + 1):
                # File profile used by column i
                ifile = min(i, file_variable.shape[0]) - 1
                if not l_interpolate[i - 1]:
                    interp_variable[i - 1, k - 1] = file_variable[ifile, k - 1]
                elif k < k_lowest_input[i - 1]:
                    if k == 1 and l_spec_bound_cond[i - 1]:
                        interp_variable[i - 1, k - 1] = file_variable[ifile, 0]
                    else:
                        interp_variable[i - 1, k - 1] = np.inf
                elif k > k_highest_input[i - 1]:
                    interp_variable[i - 1, k - 1] = np.inf
                elif l_lin_int[i - 1, k - 1]:
                    interp_variable[i - 1, k - 1] = lin_interpolate_two_points(
                        out_heights[i - 1, k - 1],
                        faverage.z[upper_lev_idx[i - 1, k - 1] - 1],
                        faverage.z[lower_lev_idx[i - 1, k - 1] - 1],
                        file_variable[ifile, upper_lev_idx[i - 1, k - 1] - 1],
                        file_variable[ifile, lower_lev_idx[i - 1, k - 1] - 1],
                    )
                else:
                    interp_variable[i - 1, k - 1] = file_variable[ifile, exact_lev_idx[i - 1, k - 1] - 1]

        # Apply an exponent (for the tuner), retaining the per-column time sum order.
        if npower != 1:
            stat_file_average = stat_file_average + interp_variable ** npower
        else:
            stat_file_average = stat_file_average + interp_variable

    # Close the netCDF file
    close_netcdf_read(faverage)

    # Take average over num_timesteps
    stat_file_average = stat_file_average / float(num_timesteps)

    return stat_file_average, l_error


# -------------------------------------------------------------------------
def stat_file_average_interval(filename, nz, t, variable_name,
                               out_heights, npower):
    """Reads in netCDF data from a file and then takes several averages
    over an interval.

    References:
      None

    Notes:
      The record list t holds pairs of start and end records; a 0 ends the
      list early.

    Returns:
        stat_file_average_interval: The weighted average, shape (nz,).
        l_error: status of this function.
    """
    l_error = False
    stat_file_average_interval = np.zeros(nz)

    # Sanity check
    if len(t) < 2:
        print("stat_file_average_interval: Invalid time interval", file=sys.stderr)
        l_error = True
        return stat_file_average_interval, l_error

    tdim = -1
    for i in range(1, len(t) + 1):
        if t[i - 1] == 0:
            break
        tdim = i

    l_spec_bound_cond = np.array([False])
    averaged_profile, l_error = stat_file_average(
        filename, nz, 1,
        t[0], t[1], np.reshape(out_heights, (1, nz)),
        variable_name,
        npower, l_spec_bound_cond,
    )
    stat_file_average_interval = averaged_profile[0, :] * float(t[1] - t[0])

    divisor = t[1] - t[0]

    if l_error:
        return stat_file_average_interval, l_error

    for i in range(3, tdim + 1, 2):
        l_spec_bound_cond = np.array([False])
        averaged_profile, l_error = stat_file_average(
            filename, nz, 1,
            t[i - 1], t[i], np.reshape(out_heights, (1, nz)),
            variable_name,
            npower, l_spec_bound_cond,
        )
        stat_file_temp = averaged_profile[0, :]
        stat_file_average_interval = (
            stat_file_average_interval + stat_file_temp * float(t[i] - t[i - 1])
        )
        divisor = divisor + (t[i] - t[i - 1])

    stat_file_average_interval = stat_file_average_interval / float(divisor)

    return stat_file_average_interval, l_error


# -------------------------------------------------------------------------
def stat_file_num_vertical_levels(varname, filename):
    """Returns a integer for the number of vertical levels in a file.

    References:
      None
    """
    file_nz, l_error = open_netcdf_read(varname, filename)

    if l_error:
        print("Error opening " + str(filename), file=sys.stderr)
        raise RuntimeError("Error opening " + str(filename))

    # Set return variable
    stat_file_num_vertical_levels = file_nz.iz

    # Close file
    close_netcdf_read(file_nz)

    return stat_file_num_vertical_levels


# -------------------------------------------------------------------------
def stat_file_vertical_levels(varname, filename, nz):
    """Return a vector containing the heights of each vertical level [m]."""
    file_levels, l_error = open_netcdf_read(varname, filename)

    if l_error:
        print("Error opening " + str(filename), file=sys.stderr)
        raise RuntimeError("Error opening " + str(filename))

    # Set return variable
    stat_file_vertical_levels = file_levels.z[:nz].copy()

    # Close file
    close_netcdf_read(file_levels)

    return stat_file_vertical_levels


# -------------------------------------------------------------------------
def l_netcdf_file(filename):
    """Return whether filename ends in ".nc"; stop for any other file type."""
    filename = str(filename).strip()

    if filename.endswith(".nc"):
        l_netcdf_file = True
    else:
        print("Unrecognized file type: " + filename, file=sys.stderr)
        raise RuntimeError("Unrecognized file type: " + filename)

    return l_netcdf_file


# ===============================================================================
def CLUBB_levels_within_LES_domain(fread_var, nz, ngrdcol, CLUBB_grid):
    """Find each column's lowest and highest grid levels within the LES domain.
    The input column grids are ascending, as in the profile reader.

    References:
      None

    Args:
        fread_var: LES grid and valid level indices.
        nz: Number of CLUBB vertical levels [-].
        ngrdcol: Number of columns [-].
        CLUBB_grid: Column heights on thermodynamic or momentum levels,
            shape (ngrdcol, nz) [m].

    Returns:
        k_lowest_input: Lowest CLUBB level within the LES domain, shape (ngrdcol,) [-].
        k_highest_input: Highest CLUBB level within the LES domain, shape (ngrdcol,) [-].
    """
    if np.any(CLUBB_grid[:, nz - 1] < fread_var.z[fread_var.ia - 1]):
        print("The lowest LES input level is above the top "
              "of the CLUBB model domain.", file=sys.stderr)
        raise RuntimeError("Error in CLUBB_levels_within_LES_domain")

    if np.any(CLUBB_grid[:, 0] > fread_var.z[fread_var.iz - 1]):
        print("The highest LES input level is below the bottom "
              "of the CLUBB model domain.", file=sys.stderr)
        raise RuntimeError("Error in CLUBB_levels_within_LES_domain")

    k_lowest_input = np.full(ngrdcol, 1, dtype=int)
    k_highest_input = np.full(ngrdcol, nz, dtype=int)

    for k in range(1, nz + 1):
        for i in range(1, ngrdcol + 1):
            if CLUBB_grid[i - 1, k - 1] < fread_var.z[fread_var.ia - 1]:
                k_lowest_input[i - 1] = k + 1
            if CLUBB_grid[i - 1, k - 1] > fread_var.z[fread_var.iz - 1]:
                k_highest_input[i - 1] = min(k_highest_input[i - 1], k - 1)

    return k_lowest_input, k_highest_input


# ===============================================================================
def LES_grid_to_CLUBB_grid(fread_var, CLUBB_height):
    """Finds the level on the LES grid that is exactly even with the CLUBB
    grid level (either thermodynamic or momentum grid) that is input
    height. Else, it finds the two LES levels that sandwich the CLUBB
    height that is input.

    Args:
        fread_var: Information about LES run.
        CLUBB_height: Altitude of one CLUBB grid level
            (either thermodynamic or momentum grid levels)  [m]

    Returns:
        exact_lev_idx: In case of an exact match, index of LES level that is
            exactly even with CLUBB level k.
        lower_lev_idx: In case linear interpolation is needed, index of LES
            level that is immediately below CLUBB level k.
        upper_lev_idx: In case linear interpolation is needed, index of LES
            level that is immediately above CLUBB level k.
        l_lin_int: Flag that is turned on if linear interpolation is needed.
    """
    # Initialize the output quantities.
    exact_lev_idx = 0
    lower_lev_idx = 0
    upper_lev_idx = 0
    l_lin_int = False

    # Initialize LES vertical grid loop index, j, at the lowest LES grid index,
    # which is fread_var%ia.

    for j in range(fread_var.ia, fread_var.iz + 1):

        if abs(fread_var.z[j - 1] - CLUBB_height) < 1e-6:

            # There is an LES level altitude at LES level j that is an exact
            # match to the CLUBB level altitude at CLUBB grid level k.
            exact_lev_idx = j
            l_lin_int = False

        elif fread_var.z[j - 1] < CLUBB_height:

            # The LES level altitude at LES level j is lower than the CLUBB level
            # altitude at CLUBB grid level k.
            lower_lev_idx = j

        else:  # fread_var%z(j) > CLUBB_height

            # The LES level altitude at LES level j is higher than the CLUBB level
            # altitude at CLUBB grid level k.
            upper_lev_idx = j

        if exact_lev_idx > 0:
            break  # An exact answer has been found, exit the loop.

        if upper_lev_idx == lower_lev_idx + 1:

            # CLUBB level k has been found between two successive LES levels.
            # Linear interpolation is needed.  An answer has been found, exit
            # the loop.
            l_lin_int = True
            break

    return exact_lev_idx, lower_lev_idx, upper_lev_idx, l_lin_int
