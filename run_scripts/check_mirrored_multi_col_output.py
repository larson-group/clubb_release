#!/usr/bin/python3

# =========================================================================================================
# Description: Checks mirrored multi-column stats output in a directory.
#              By mirrored, we mean column 1 should match column N, column 2 should match column N-1, etc.
#              This script exits non-zero if any mirrored-column comparison exceeds the tolerance.
#
#              A mirrored parameter file can be created by running "./create_multi_col_params.py" with
#              "-mode dup_tweak -mirror true". Then generate output with run_scm.py/run_scm_all.py and run
#              this script on the output directory containing *_stats.nc files.
# =========================================================================================================

import argparse
import os
import sys

import netCDF4
import numpy as np


def _is_numeric_netcdf_var(var):
    return np.issubdtype(var.dtype, np.number)


def _compare_mirrored_columns_vectorized(data_with_col_last):
    """Compare mirrored columns in one vectorized pass.

    data_with_col_last is expected to have the column axis at -1.
    Returns average and max absolute differences for each mirrored pair.
    """
    ncol = data_with_col_last.shape[-1]
    npairs = ncol // 2
    if npairs == 0:
        return np.array([]), np.array([])

    # Compare first half against reversed second half:
    # col 0 vs col n-1, col 1 vs col n-2, ...
    lhs = data_with_col_last[..., :npairs]
    rhs = np.flip(data_with_col_last, axis=-1)[..., :npairs]
    abs_diff = np.ma.abs(lhs - rhs)

    # Reduce over every axis except mirrored-pair axis (last axis).
    if abs_diff.ndim == 1:
        avg_by_pair = abs_diff
        max_by_pair = abs_diff
    else:
        # After moving col to the last axis, all leading axes are sample axes
        # (e.g., time/zt/zm/...), so reduce them and keep one value per col-pair.
        reduce_axes = tuple(range(abs_diff.ndim - 1))
        avg_by_pair = np.ma.mean(abs_diff, axis=reduce_axes)
        max_by_pair = np.ma.max(abs_diff, axis=reduce_axes)

    # Fill masked entries to keep output stable.
    return np.ma.filled(avg_by_pair, 0.0), np.ma.filled(max_by_pair, 0.0)


def check_file(file_path, tolerance, verbose):
    try:
        dset = netCDF4.Dataset(file_path)
    except Exception as err:
        print(f"Error opening file {file_path}: {err}")
        return True

    with dset:
        if "col" not in dset.dimensions:
            print(f"Skipping {file_path}: no 'col' dimension found.")
            return False

        ngrdcol = len(dset.dimensions["col"])
        if ngrdcol < 2:
            print(f"Skipping {file_path}: col dimension size is {ngrdcol}.")
            return False

        print(f"Testing {file_path} with ngrdcol = {ngrdcol}")
        differences_found = False
        checked_var_count = 0
        all_zero_var_count = 0
        differing_var_count = 0

        for var_name, var in dset.variables.items():
            if var_name in dset.dimensions:
                continue
            if "col" not in var.dimensions:
                continue
            if "time" not in var.dimensions:
                # This checker is for time-evolving model output fields only.
                # Static arrays (e.g., param x col) are metadata, not mirrored-evolution checks.
                continue
            if not _is_numeric_netcdf_var(var):
                continue

            checked_var_count += 1
            col_axis = var.dimensions.index("col")
            data = np.ma.array(var[:], copy=False)
            # Normalize to a canonical shape where the mirrored column axis is always last,
            # regardless of whether the file stored it as (time, z, col), (time, col), etc.
            data = np.moveaxis(data, col_axis, -1)

            # Fast all-zero detection for this variable.
            if data.count() > 0 and float(np.ma.max(np.ma.abs(data))) == 0.0:
                all_zero_var_count += 1

            if verbose:
                print(f"checking {var_name} dims={var.dimensions}")

            avg_by_pair, max_by_pair = _compare_mirrored_columns_vectorized(data)
            variable_has_difference = False
            for col in range(avg_by_pair.shape[0]):
                col_begin = col
                col_end = ngrdcol - 1 - col
                avg_abs_diff = float(avg_by_pair[col])
                max_abs_diff = float(max_by_pair[col])

                if avg_abs_diff > tolerance:
                    print(
                        f" - {var_name} col{col_begin} vs col{col_end}: "
                        f"avg_abs_diff = {avg_abs_diff}, max_abs_diff = {max_abs_diff}"
                    )
                    differences_found = True
                    variable_has_difference = True

            if variable_has_difference:
                differing_var_count += 1

        print(
            f"Summary for {file_path}: checked={checked_var_count}, "
            f"all_zero={all_zero_var_count}, differing={differing_var_count}"
        )

        return differences_found


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check mirrored columns in *_stats.nc output files."
    )
    parser.add_argument(
        "directory",
        nargs=1,
        help="Directory containing *_stats.nc files to check",
    )
    parser.add_argument(
        "-t",
        "--tol",
        type=float,
        default=0.0,
        help="Average absolute difference tolerance (default: 0.0)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print each checked variable name.",
    )

    args = parser.parse_args()
    directory = args.directory[0]

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    matching_files = sorted(
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith("_stats.nc")
    )

    if not matching_files:
        print(f"No *_stats.nc files found in directory: {directory}")
        sys.exit(0)

    differences_found_any = False
    for file_path in matching_files:
        if check_file(file_path, args.tol, args.verbose):
            differences_found_any = True

    if differences_found_any:
        print("Differences found in one or more files.")
        sys.exit(1)

    print("All checks passed.")
    sys.exit(0)
