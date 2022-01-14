#!/usr/bin/python3

import subprocess
import argparse
import netCDF4
import sys
import numpy

parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2,
                    help="need two files to diff")
                    
args = parser.parse_args()


# Create diff command, using ncdiff and -O to overwrite file if it exists
ncdiffCommand = "ncdiff -O " + args.files[0] + " " + args.files[1] + " diff.nc"

# Run diff command, printing output (there should be no output unless there's an error)
process = subprocess.Popen(ncdiffCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

# Create dataset from nc file
dset = netCDF4.Dataset("diff.nc")

# Loop through all variables in dataset, and make sure they are all completely 0
all_zero = True
for var in dset.variables:
    
    # The CLUBB netcdf variables we are interested in are all 
    # 4 dimensional (time, altitude, latitude, longitude),
    # but currently we don't actually have latitude or longitude in clubb, so those
    # dimensions are hardcoded to be 1. If in the future we remove those useless
    # dimensions (unlikely), then the variables of interested would be 2D. So 
    # for futureproofing, we will just check all variables with more than 1 dimension.
    if dset[var].ndim > 1:
        if not numpy.all(dset[var][:] == 0):
            print(var + " is NON-ZERO")
            all_zero = False
            
if all_zero:
    sys.exit(0)
else:
    print("FAIL: There were some non-zero fields.")
    sys.exit(1)
