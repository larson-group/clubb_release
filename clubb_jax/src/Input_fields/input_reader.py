"""This module is respondsible for the procedures and structures necessary to
read in "SAM-Like" case specific files. Currently only the
<casename>_sounding.in file is formatted to be used by this module.

References:
    None.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from clubb_jax.src.CLUBB_core.interpolation import zlinterp_fnc


@dataclass
class one_dim_read_var:
    """Derived type for representing a rank 1 variable read by a procedure."""

    name: str
    dim_name: str
    values: np.ndarray


@dataclass
class two_dim_read_var:
    """Derived type for representing a rank 2 variable read by a procedure."""

    name: str
    dim1_name: str
    dim2_name: str
    values: np.ndarray


blank_value = -999.9


# -------------------------------------------------------------------------------------------------
def read_two_dim_file(iunit, nCol, filename):
    """Read from a file containing data that varies in two dimensions. These
    dimensions are typically height and time.

    Returns (read_vars, other_dim) in place of the source output arguments.
    """
    del iunit

    # First run through, take names and determine how large the data file is.
    lines = Path(filename).read_text().splitlines()

    # Skip all the comments at the top of the file.
    while lines and (not lines[0].strip() or "!" in lines[0]):
        lines.pop(0)
    if not lines:
        raise ValueError(f"{filename}: empty two-dimensional input file")

    names = [name.strip("'\"") for name in lines.pop(0).split()]
    if len(names) != nCol:
        raise ValueError(f"{filename}: expected {nCol} columns, found {len(names)}")

    other_dim_values = []
    blocks = []
    line_index = 0
    while line_index < len(lines):
        if not lines[line_index].strip():
            line_index += 1
            continue
        header = lines[line_index].split()
        line_index += 1
        if len(header) != 2:
            raise ValueError(f"{filename}: invalid time-block header {header!r}")
        other_dim_values.append(float(header[0]))
        nRowI = int(header[1])
        if nRowI < 1:
            raise ValueError(
                "Number of elements must be an integer and greater than zero "
                "in two-dim input file."
            )
        rows = []
        for _ in range(nRowI):
            if line_index >= len(lines):
                raise ValueError(f"{filename}: block ended before {nRowI} rows were read")
            row = lines[line_index].split()
            line_index += 1
            if len(row) != nCol:
                raise ValueError(f"{filename}: expected {nCol} columns, found {len(row)}")
            rows.append([float(value) for value in row])
        blocks.append(np.asarray(rows))

    if not blocks:
        raise ValueError(f"{filename}: no two-dimensional data blocks found")
    nRowI = blocks[0].shape[0]
    if any(block.shape[0] != nRowI for block in blocks):
        raise ValueError(f"{filename}: all time blocks must have the same row count")

    # Store the names into the structure and allocate accordingly.
    values = np.stack(blocks, axis=2)
    read_vars = [
        two_dim_read_var(
            name=names[k],
            dim1_name="Time[s]",
            dim2_name=names[0],
            values=values[:, k, :],
        )
        for k in range(nCol)
    ]
    other_dim = one_dim_read_var(
        name="Time[s]",
        dim_name="Time[s]",
        values=np.asarray(other_dim_values),
    )
    return read_vars, other_dim


# ------------------------------------------------------------------------------------------------
def read_one_dim_file(iunit, nCol, filename):
    """Read from a file containing data that varies in one dimension. The
    dimension is typically time.

    References:
        None.
    """
    del iunit

    # First run through, take names and determine how large the data file is.
    lines = Path(filename).read_text().splitlines()

    # Skip all the comments at the top of the file.
    while lines and (not lines[0].strip() or "!" in lines[0]):
        lines.pop(0)
    if not lines:
        raise ValueError(f"{filename}: empty one-dimensional input file")

    names = [name.strip("'\"") for name in lines.pop(0).split()]
    if len(names) != nCol:
        raise ValueError(f"{filename}: expected {nCol} columns, found {len(names)}")

    # Count and read the rows into the newly allocated arrays.
    rows = []
    for line_number, line in enumerate(lines, start=2):
        if not line.strip():
            continue
        row = line.split()
        if len(row) != nCol:
            raise ValueError(
                f"{filename}:{line_number}: expected {nCol} columns, found {len(row)}"
            )
        rows.append([float(value) for value in row])
    if not rows:
        raise ValueError(f"{filename}: no one-dimensional data rows found")

    values = np.asarray(rows)
    return [
        one_dim_read_var(name=names[k], dim_name=names[0], values=values[:, k])
        for k in range(nCol)
    ]


# ------------------------------------------------------------------------------------------------
def fill_blanks_one_dim_vars(num_vars, one_dim_vars):
    """Fill blank spots, signified by blank_value, with values linearly
    interpolated using the first array element as a guide.

    References:
        None.
    """
    for i in range(num_vars):
        one_dim_vars[i].values = linear_fill_blanks(
            one_dim_vars[i].values.size,
            one_dim_vars[0].values,
            one_dim_vars[i].values,
            0.0,
        )
    return one_dim_vars


# ------------------------------------------------------------------------------------------------
def fill_blanks_two_dim_vars(num_vars, other_dim, two_dim_vars):
    """Fill blank spots, signified by blank_value, with values linearly
    interpolated using the first array element and other_dim as a guide.

    This is a two-step process. First assume that other_dim has no holes, but
    there are blanks for that variable across that dimension. Then fill holes
    across the dimension whose values are first in two_dim_vars.

    Ex. Time is other_dim and Height in meters is the first element in
    two_dim_vars.

    References:
        None.
    """
    dim_size = two_dim_vars[0].values.shape[0]
    other_dim_size = other_dim.values.size

    for i in range(1, num_vars):
        # Interpolate along main dim.
        for j in range(other_dim_size):
            two_dim_vars[i].values[:, j] = linear_fill_blanks(
                dim_size,
                two_dim_vars[0].values[:, j],
                two_dim_vars[i].values[:, j],
                blank_value,
            )

        # Interpolate along other dim.
        for j in range(dim_size):
            two_dim_vars[i].values[j, :] = linear_fill_blanks(
                other_dim_size,
                other_dim.values,
                two_dim_vars[i].values[j, :],
                blank_value,
            )
    return two_dim_vars


# ------------------------------------------------------------------------------------------------
def linear_fill_blanks(dim_grid, grid, var, default_value):
    """Fill blanks in var using grid as a guide.

    Blank values are signified by being less than or equal to blank_value.

    References:
        None.
    """
    grid = np.asarray(grid)
    var = np.asarray(var)
    if grid.size != dim_grid or var.size != dim_grid:
        raise ValueError("dim_grid must match grid and var")

    reversed_grid = np.any(grid[1:] < grid[:-1])

    # Essentially this code leverages the previously written zlinterp function.
    # A smaller temporary grid and var variable are being created to pass to
    # zlinterp. zlinterp then performs the work of taking the temporary var
    # array and interpolating it to the actual grid array.
    valid = var > blank_value
    amt = int(np.count_nonzero(valid))
    if amt == 0:
        return np.full_like(var, default_value)
    if amt == dim_grid:
        return np.array(var, copy=True)

    temp_grid = grid[valid]
    temp_var = var[valid]
    if reversed_grid:
        return np.asarray(zlinterp_fnc(-grid, -temp_grid, temp_var))
    return np.asarray(zlinterp_fnc(grid, temp_grid, temp_var))


# ----------------------------------------------------------------------------
def deallocate_one_dim_vars(num_vars, one_dim_vars):
    """Deallocate values stored in the whole array of one_dim_vars.

    References:
        None.
    """
    for i in range(num_vars):
        one_dim_vars[i].values = np.empty(0)
    return one_dim_vars


# ------------------------------------------------------------------------------------------------
def deallocate_two_dim_vars(num_vars, two_dim_vars, other_dim):
    """Deallocate values stored in the whole array of two_dim_vars.

    References:
        None.
    """
    for i in range(num_vars):
        two_dim_vars[i].values = np.empty((0, 0))
    other_dim.values = np.empty(0)
    return two_dim_vars, other_dim


# ------------------------------------------------------------------------------------------------
def read_x_table(nvar, xdim, ydim, target_name, retVars):
    """Search for target_name in retVars and return it, or fail when absent.

    References:
        None.
    """
    for i in range(nvar):
        if retVars[i].name == target_name:
            values = retVars[i].values
            if values.shape != (xdim, ydim):
                raise ValueError(
                    f"{target_name}: expected shape {(xdim, ydim)}, found {values.shape}"
                )
            return values
    raise ValueError(f"{target_name} could not be found")


# ------------------------------------------------------------------------------------------------
def read_x_profile(nvar, dim_size, target_name, retVars, input_file=None):
    """Search for target_name in retVars and return it, or fail with a warning.

    Modified by Cavyn, June 2010.
    """
    i = get_target_index(nvar, target_name, retVars)
    if i >= 0:
        values = retVars[i].values
        if values.size != dim_size:
            raise ValueError(
                f"{target_name}: expected {dim_size} values, found {values.size}"
            )
        return values
    if input_file is not None:
        raise ValueError(f"{target_name} could not be found. Check the file {input_file}")
    raise ValueError(f"{target_name} could not be found. Check your sounding.in file")


# ------------------------------------------------------------------------------
def get_target_index(nvar, target_name, retVars):
    """Return the index of target_name in retVars, or -1 when absent.

    The returned index is zero-based in Python.

    References:
        None.

    Created by Cavyn, July 2010.
    """
    for i in range(nvar):
        if retVars[i].name == target_name:
            return i
    return -1


# =============================================================================
def count_columns(iunit, filename):
    """Count columns assuming the first non-comment line contains headers.

    References:
        None.

    Created by Cavyn, July 2010.
    """
    del iunit
    for line in Path(filename).read_text().splitlines():
        if line.strip() and "!" not in line:
            return len(line.split())
    raise ValueError(f"Fatal error reading data in {filename}")


__all__ = [
    "one_dim_read_var",
    "read_one_dim_file",
    "two_dim_read_var",
    "read_two_dim_file",
    "fill_blanks_one_dim_vars",
    "fill_blanks_two_dim_vars",
    "deallocate_one_dim_vars",
    "deallocate_two_dim_vars",
    "read_x_table",
    "read_x_profile",
    "get_target_index",
    "count_columns",
]
