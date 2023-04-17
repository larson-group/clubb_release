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

# Create datasets from nc files
dset0 = netCDF4.Dataset(args.files[0])
dset1 = netCDF4.Dataset(args.files[1])

# Loop through all variables in dataset, and make sure they are all completely 0
all_zero = True
for var in dset0.variables:
    
    # The CLUBB netcdf variables we are interested in are all 
    # 4 dimensional (time, altitude, latitude, longitude),
    # but currently we don't actually have latitude or longitude in clubb, so those
    # dimensions are hardcoded to be 1. If in the future we remove those useless
    # dimensions (unlikely), then the variables of interested would be 2D. So 
    # for futureproofing, we will just check all variables with more than 1 dimension.
    if dset0[var].ndim > 1:
        if numpy.any( dset0[var][:] != dset1[var][:] ):
            print(var + " is NON-ZERO")
            all_zero = False
            
if all_zero:
    sys.exit(0)
else:
    print("FAIL: There were some non-zero fields.")
    sys.exit(1)
