#!/usr/bin/env python

import netCDF4
import sys

def compute_relative_diff(ncfile1, ncfile2, field):
    nc1  = netCDF4.Dataset(ncfile1)
    nc2  = netCDF4.Dataset(ncfile2)
    if field not in nc1.variables or field not in nc2.variables:
        raise Exception("The field " + field + " is not in both NetCDF files.")
    var1 = nc1.variables[field]
    var2 = nc2.variables[field]
    # The shape is assumed to be [num_timesteps,1,1,1]
    var1 = var1[:,0,0,0]
    var2 = var2[:,0,0,0]
    if len(var1) != len(var2):
        raise Exception("The fields don't have the same shape!")
    num_timesteps = len(var1)
    # Sum the field over all timesteps.
    sum1 = sum(var1); sum2 = sum(var2)
    if (sum1 != 0.0):
        relative_diff = abs( (sum2-sum1) / sum1 )
    else:
        relative_diff = -999.
    return relative_diff

###############################################################
usage = "Usage: tolerance_check.py ncfile1 ncfile2 field tolerance"
###############################################################
if len(sys.argv) != 5:
    print(usage)
    sys.exit(2)

ncfile1   = sys.argv[1]
ncfile2   = sys.argv[2]
field     = sys.argv[3]
tolerance = float(sys.argv[4])

relative_diff = compute_relative_diff( ncfile1, ncfile2, field )
print("Relative diff: " + str(relative_diff))
if relative_diff > tolerance:
    print("Relative diff tolerance exceeded!")
    sys.exit(1)
else:
    sys.exit(0)
