"""File-reading routines used by the M-PACE A forcing initialization."""

from __future__ import annotations

import numpy as np


# ===============================================================================
def file_read_1d(file_unit, path_and_filename, num_datapts, entries_per_line):
    """Read a data file with a declared number of columns.

    The data are read in the form:
    1 => (row 1, column 1); 2 => (row 1, column 2); etc.

    Example: a data file with 18 total data points (DP1 to DP18), with 4 data
    points per row.

                 i = 1     i = 2     i = 3     i = 4
              ---------------------------------------
        k = 1 |   DP1       DP2       DP3       DP4
        k = 2 |   DP5       DP6       DP7       DP8
        k = 3 |   DP9       DP10      DP11      DP12
        k = 4 |   DP13      DP14      DP15      DP16
        k = 5 |   DP17      DP18

    See Michael Falk's comments below for more information.
    """
    # ---- Begin Code ----
    # Python path-based reads do not use the Fortran unit number. NumPy handles
    # both the full and partial final lines as whitespace-separated input.
    del file_unit, entries_per_line

    # Open data file.
    variable = np.fromfile(path_and_filename, sep=" ", count=num_datapts)
    if variable.size != num_datapts:
        raise ValueError(
            f"Expected {num_datapts} values in {path_and_filename}, "
            f"found {variable.size}"
        )
    return variable


# ===============================================================================
def file_read_2d(
    device,
    file_path,
    file_dimension1,
    file_dimension2,
    file_per_line,
):
    """Read the two-dimensional data-file format used by M-PACE A.

    The files list the file_dimension2 values on a given vertical level, then
    move to the next level. Each line has a specific number of values until the
    last line on a level, which is short. The next level begins with a full line
    again. Michael Falk, 24 September 2007.

    References:
        None.
    """
    # ---- Begin Code ----
    # Python path-based reads do not use Fortran units or per-line bookkeeping.
    del device, file_per_line
    num_datapts = file_dimension1 * file_dimension2
    variable = np.fromfile(file_path, sep=" ", count=num_datapts)
    if variable.size != num_datapts:
        raise ValueError(
            f"Expected {num_datapts} values in {file_path}, found {variable.size}"
        )
    return variable.reshape(file_dimension1, file_dimension2)


__all__ = ["file_read_1d", "file_read_2d"]
