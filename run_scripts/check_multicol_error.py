#!/usr/bin/python3

import argparse
import netCDF4
import numpy as np
import tabulate
import sys

# Threshold used to ignore field values for calculating % diff
field_threshold = 1.0e-7

# Threshold used to determine if the ncfiles are close enough
total_abs_error_threshold = 1.0e-5

# Var to set if we find that the files differ significantly
files_differ = False


# Parse arguments, we only expect 2 file names
parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2,
                    help="need two files to diff")
                    
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
            'Index of Max Abs Diff']]

# Loop over each variable, the datasets should contain the same variables
for var in dset0.variables:

    # Only consider variables with 3 dimensions, these are the clubb vars we care about
    if dset0[var].ndim == 3:

        # Calculate absolute differences
        abs_diff = abs( dset0[var][:,:,:] - dset1[var][:,:,:] )

        # Clip fields to ignore tiny values for the % diff
        field_1_clipped = np.clip( dset0[var][:,:,:], a_min = field_threshold, a_max = 9999999.0  )
        field_2_clipped = np.clip( dset1[var][:,:,:], a_min = field_threshold, a_max = 9999999.0 )

        # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
        percent_diff = 200.0 * ( field_1_clipped-field_2_clipped ) \
                               / ( field_1_clipped+field_2_clipped )
                     
        # Append the table array with the values we want to print   
        table.append( [ var, 
                        np.max(abs_diff), 
                        np.max(percent_diff), 
                        np.sum(abs_diff), 
                        np.average(abs_diff), 
                        np.where(abs_diff == np.max(abs_diff)) ] )

        # If the total absolute difference excedes the threshold, then we consider the files different
        if np.sum(abs_diff) > total_abs_error_threshold:
            files_differ = True

# Print a very pretty table of the values
print("\n",tabulate.tabulate(table, headers='firstrow'))



if not files_differ:
    print("\nPASSED: Sum of all absolute differences does not excede",total_abs_error_threshold,"for any field.")
    sys.exit()
else:
    print("\n###############################################################################################")
    print("WARNING: Sum of all absolute differences excedes",total_abs_error_threshold," for some fields!")
    print("         It's possible that no error has been introduced, but if the case being compared")
    print("         is arm or bomex, this is not a good sign.")
    print("###############################################################################################\n")

    print("Print max absolute difference by timestep?")
    timestep_jump = int(input("Enter timestep jump to view, or 0 to skip: "))

    if timestep_jump == 0:
        sys.exit()


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
