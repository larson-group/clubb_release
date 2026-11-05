"""File-reading routines used by the M-PACE A forcing initialization."""

from __future__ import annotations

import numpy as np
import sys


# ===============================================================================
def file_read_1d(file_unit, path_and_filename, num_datapts, entries_per_line):
    """This subroutine reads in values from a data file with a number of
    rows and a declared number of columns (entries_per_line) of data.
    It reads in the data in the form of:
    1 ==> (row 1, column 1); 2 ==> (row 1, column 2); etc.

    Example:  a diagram of a data file with 18 total data points
              (DP1 to DP18), with 4 data points per row.

             i = 1     i = 2     i = 3     i = 4
          ---------------------------------------
    k = 1 |   DP1       DP2       DP3       DP4
          |
    k = 2 |   DP5       DP6       DP7       DP8
          |
    k = 3 |   DP9       DP10      DP11      DP12
          |
    k = 4 |   DP13      DP14      DP15      DP16
          |
    k = 5 |   DP17      DP18

    See Michael Falk's comments below for more information.

    Returns:
        variable

    Fortran declarations (file_read_1d in src/CLUBB_core/file_functions.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
         file_unit,          & ! Unit number of file being read.
         num_datapts,        & ! Total number of data points being read in.
         entries_per_line   ! Number of data points
        ! on one line of the file being read.

        character(*), intent(in) :: &
         path_and_filename  ! Path to file and filename of file being read.

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(num_datapts), intent(out) :: &
         variable           ! Data values output into variable
    """
    # ---- Begin Code ----
    # Python path-based reads do not use the Fortran unit number. NumPy handles
    # both the full and partial final lines as whitespace-separated input.
    del file_unit, entries_per_line

    # Open data file.
    try:
        variable = np.fromfile(path_and_filename, sep=" ", count=num_datapts)
    except OSError:
        print("CLUBB encountered an error trying to open " + str(path_and_filename), file=sys.stderr)
        raise
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
    """Michael Falk wrote this routine to read data files in a particular format for mpace_a.
    The 2d mpace_a files list the (file_dimension2) values on a given vertical level, then
    moves to the next level to list its values.
    Each line has a specific number of values, until the last line on a level, which
    is short-- it has the last few values and then a line break.  The next line, beginning
    the next level, is full-sized again.  24 September 2007

    Fortran declarations (file_read_2d in src/CLUBB_core/file_functions.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          device, &
          file_dimension1, &
          file_dimension2, &
          file_per_line

        character(*), intent(in) :: &
          file_path

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(file_dimension1,file_dimension2), intent(out) :: &
          variable
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
