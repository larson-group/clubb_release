# $Id$
#
# convert_nc_var.py
#
# Description:
#   Applies a conversion to a variable within a netcdf file.

import netCDF4
import numpy as np
from shutil import copyfile

var_to_edit = 'QTO3'

# Add or multiply for conversion
add_to_convert = False
mult_to_convert = True

# Conversion factor for either addition or multiplication
conversion = 1000.

# File location and name
location = '/home/weberjk/workspace_group_config/CLUBB/input/input_fields/'
filename = 'LBA_250m.nc'
nc_in_file = location+filename

# ---------------------------------------------------------------------------------------
# Begin Code
# --------------------------------------------------------------------------------------

# Create a copy of the file
copyfile(nc_in_file,location+'copy_of_orig.nc')

# Open the input file in 'read/write' mode
nc_in = netCDF4.Dataset(nc_in_file,'r+')
var = nc_in.variables[var_to_edit]

# Copy netcdf object into array
var_convert = var[:]

# Convert array
if(add_to_convert ):
    var_convert = var_convert + conversion
elif(mult_to_convert):
    var_convert = var_convert * conversion
else:
    print "Error: Neither add_to_convert of mult_to_convert are True"

# Copy array into netcdf object.
var[:] = var_convert    

nc_in.close()
