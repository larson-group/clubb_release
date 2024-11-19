#!/usr/bin/python3

#=========================================================================================================
# Description: Compare 2 multicolumn standalone netcdf files.
#              Used to help determine if an error has been introduced in 
#              clubb_standalone. It's reccomended to only use this test
#              for arm/bomex, as they are the most stable cases, but it can
#              be used for any clubb cases.
#   
# Usage: Checkout a clean version of clubb with multicolumns enabled by setting 
#        namelist variables in clubb/input/tunable_parameters/configurable_model_flags.in
#
#               $configurable_multi_column_nl
#               num_standalone_columns = # of columns to run
#
#        Then clubb_standalone will output the netcdf data to ../output/case_multi_col.nc
#        which should be saved for comparison somewhere (say save_case_multi_col.nc). Then make 
#        changes to clubb source code, compile and rerun the case with the same number of 
#        standalone_columns, then compare the output files by calling this script with 
#
#           ./check_multicolumn_error.py save_case_multi_col.nc ../output/case_multi_col.nc
#
#        It's assumed that the largest single difference of any field, over any timestep/level/column,
#        will be less than 1.0e-7. For more stable cases like arm and bomex the largest difference is 
#        usually much smaller than 1e-7, but the tolerance has been relaxed for mpace_b, which is noiser,
#        but seems to be the least noisey case where ice develops.
#
# Author: Gunther Huebler
# Reference: https://github.com/larson-group/clubb/issues/1033
#=========================================================================================================

import argparse
import netCDF4
import numpy as np
import tabulate
import sys

# Threshold used to ignore field values for calculating % diff
field_threshold = 1.0e-7

# Threshold used to determine if the ncfiles are close enough
abs_error_threshold = 1.0e-7

# Var to set if we find that the files differ significantly
files_differ = False


# Parse arguments, we only expect 2 file names
parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2, help="need two files to diff")
parser.add_argument("-s","--scale", action="store_true", help="Scale differences by field value")
                    
args = parser.parse_args()


# Create dataset from nc file
dset0 = netCDF4.Dataset(args.files[0])
dset1 = netCDF4.Dataset(args.files[1])

# Define table header, these are the values we're going to output
table = [[  'Var', 
            'Max Abs Diff', 
            'Max % Diff', 
            'Total Abs Diff', 
            'Avg Abs Diff', 
            'Index of Max Abs Diff (time,nz,ngrdcol)']]

# Loop over each variable, the datasets should contain the same variables
for var in dset0.variables:

    # Only consider variables with 3 dimensions, these are the clubb vars we care about
    if dset0[var].ndim == 3:

        # Calculate absolute differences
        abs_diff = abs( dset0[var][:,:,:] - dset1[var][:,:,:] )

        # Scale differences down by the field average, use ceiling to prevent scaling up fields with small values
        if args.scale:
            field_avg = ( np.average( dset0[var] ) + np.average( dset1[var] ) ) / 2.0
            abs_diff = abs_diff / np.ceil(field_avg)

        # Clip fields to ignore tiny values for the % diff
        field_1_clipped = np.clip( dset0[var][:,:,:], a_min = field_threshold, a_max = 9999999.0  )
        field_2_clipped = np.clip( dset1[var][:,:,:], a_min = field_threshold, a_max = 9999999.0 )

        # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
        percent_diff = 200.0 * ( field_1_clipped-field_2_clipped ) \
                               / ( field_1_clipped+field_2_clipped )

        # Save max abs diff since it is reused a few times
        max_abs_diff = np.max(abs_diff)
        avg_abs_diff = np.average(abs_diff)
                     
        # Append the table array with the values we want to print   
        table.append( [ var, 
                        max_abs_diff, 
                        np.max(percent_diff), 
                        np.sum(abs_diff), 
                        avg_abs_diff, 
                        np.where(abs_diff == max_abs_diff) ] )

        # If the largest absolute difference excedes the threshold, then we consider the files different
        if avg_abs_diff > abs_error_threshold:
            files_differ = True

# Print a very pretty table of the values
print("\n",tabulate.tabulate(table, headers='firstrow'))



if not files_differ:
    print("\nPASSED: The average absolute difference does not excede",abs_error_threshold,"for any field.")
    sys.exit(0)
else:
    print("\n###############################################################################################")
    print("WARNING: The average absolute difference excedes",abs_error_threshold," for some field(s)!")
    print("         It's possible that no error has been introduced, but if the case being compared")
    print("         is arm/bomex/cobra then there is very likely an error.")
    print("###############################################################################################\n")

    #print("Print max absolute difference by timestep?")
    #timestep_jump = int(input("Enter timestep jump to view, or 0 to skip: "))

    #if timestep_jump == 0:
    sys.exit(1)


table = [[ 'Timestep' ]]

# Loop over each variable in the dataset
for var in dset0.variables:
    if dset0[var].ndim == 3:
        table[0].append(var)

timesteps = dset0.dimensions['time'].size


# Loop over each timestep
for time in range(0,timesteps,timestep_jump):

    # Save current timestep to table
    table.append( [time] )

    # Loop over each variable in the dataset
    for var in dset0.variables:

        # Only consider variables with 3 dimensions, these are the clubb vars we care about
        if dset0[var].ndim == 3:

            # Calculate max_abs_diff for the current timestep
            max_abs_diff = np.max( abs( dset0[var][time,:,:] - dset1[var][time,:,:] ) )

            # It may be useful to view the absolute differences scaled by timestep since
            # we expect error to accumulate, by default we will not scale
            scaled_diff = max_abs_diff #/ ( float(time+1) )

            # Append the scaled_diff to the last entry of the current line in the table
            table[-1].append( scaled_diff )


print(tabulate.tabulate(table, headers='firstrow'))
