#!/usr/bin/python3

#=========================================================================================================
# Description: Checks the multi_col output files in a directory to see if mirrored columns are identical.
#              By mirrored, we mean that the first columns should be identical to the last one,
#              and the second should be identical to the second to last, and so on. This script 
#              will exit with a non-zero error code if any differences in these comparisons are
#              found.
#
#              A mirrored parameter file can be created by running "./create_multi_col_params.py" with
#              the arguments "-mode dup_tweak" and "-mirror". Then generate the multi_col clubb output
#              with "./run_scm.bash -p clubb_params_multi_col.in CASE". Then run this script on the 
#              clubb/output folder to check all  _multi_col_z{m/t}.nc netcdf files it contains.
#=========================================================================================================

import argparse
import netCDF4
import numpy as np
import os
import sys


def check_file(file_path):
    try:
        dset = netCDF4.Dataset(file_path)
    except Exception as e:
        print(f"Error opening file {file_path}: {e}")
        return False

    ngrdcol, = dset["columns"].shape

    print(f"Testing {file_path} with ngrdcol = {ngrdcol}")

    l_differences_found = False

    # Loop over each variable in dataset
    for var in dset.variables:

        # Only consider variables with 3 dimensions, these are the clubb vars we care about
        if dset[var].ndim == 3:

            print(f"checking {var}")

            for col in range(ngrdcol // 2):

                # This assumes the columns are mirrored, e.g., column 1 is identical to
                # column ngrdcol, and column 2 is identical to column ngrdcol-1, etc
                col_begin = col
                col_end = ngrdcol - 1 - col

                # Calculate the average absolute difference between the columns that should match
                avg_abs_diff = np.average(
                    abs(dset[var][:, :, col_begin] - dset[var][:, :, col_end])
                )

                # If the columns are identical, the average absolute difference should be zero
                if avg_abs_diff > 0.0:
                    print(f" - {var} col{col_begin} vs col{col_end}: avg_abs_diff = {avg_abs_diff}")
                    l_differences_found = True

    return l_differences_found


# Main function to handle directory input
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a test on multi_col files in a directory')

    parser.add_argument("directory", nargs=1,
                        help="Directory containing multi_col files to check")

    args = parser.parse_args()

    directory = args.directory[0]

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    # Find files ending with _multi_col_zm.nc or _multi_col_zt.nc
    patterns = ["_multi_col_zm.nc", "_multi_col_zt.nc"]
    matching_files = [
        os.path.join(directory, f) for f in os.listdir(directory)
        if any(f.endswith(pattern) for pattern in patterns)
    ]

    if not matching_files:
        print(f"No matching files found in directory: {directory}")
        sys.exit(0)

    for file_path in matching_files:
        l_differences_found = check_file(file_path)

    # Exit with error code if any differences were found in any file
    if l_differences_found:
        print("Differences found in one or more files.")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)
