#!/usr/bin/python3

#=========================================================================================================
# Description: The clubb multi_col output data has the ability to be averaged over a sampling
#              period before being output. This script is intended to compare the time averaged
#              multi_col output files to a "master" output file that was generated with
#              full time resolution. Passing this test confirms that the time averaged netcdf
#              file is being correctly generated.
#   
# Usage: Generate multi_col output data with "calls_per_out = 1" to get a full time resolution 
#        netcdf file and save it somewhere.
#
#               ./create_multi_col_params.py -param_file -n 100 -calls_per_out 1
#               ./run_scm.bash -e -p clubb_params_multi_col.in arm
#               mv ../output ../output_master; mkdir ../output
#
#        Generate multi_col output data with "calls_per_out = n" with n > 1 to get a time averaged
#        netcdf file to compare to master. Use the same parameter generation method for both, 
#        handld by default if you don't input the "mode" argument.
#
#               ./create_multi_col_params.py -param_file -n 100 -calls_per_out 4
#               ./run_scm.bash -e -p clubb_params_multi_col.in arm
#                  
#        Run this script on the gerated multi_col files
#
#               python check_multi_col_avg.py ../output_master/arm_multi_col_zt.nc ../output/arm_multi_col_zt.nc
#
# Author: Gunther Huebler
#=========================================================================================================

import argparse
import netCDF4
import numpy as np
import tabulate
import sys

# Threshold used to determine if the differences are close enough
abs_error_threshold = 1.0e-7

# Parse arguments, we only expect 2 file names
parser = argparse.ArgumentParser(description='Run a test')

parser.add_argument("files", nargs=2, help="need two files to diff")
                    
args = parser.parse_args()

# Create datasets from nc file
dset_fullres = netCDF4.Dataset(args.files[0])
dset_avg     = netCDF4.Dataset(args.files[1])


timestep_full = dset_fullres.dimensions['time'].size
#timestep_full = 12
timestep_avg = dset_avg.dimensions['time'].size

# Calculate the sample size used to create the time averaged netcdf file
calls_per_out = int( np.ceil( timestep_full / timestep_avg) )

print(f"{timestep_full} timesteps in {args.files[0]}")
print(f"{timestep_avg} timesteps in {args.files[1]}")
print(f" - assuming sampling period calls_per_out = {calls_per_out}\n")

# Initialize error code to 0 (no error)
exit_code = 0

# Loop over each variable, the datasets should contain the same variables
for var in dset_fullres.variables:

    # Only consider variables with 3 dimensions, these are the clubb vars we care about
    if dset_fullres[var].ndim == 3:

        print(f"comparing {var}")

        master_sum = 0.0
        samples = 0
        max_diff = 0.0

        # Calculate the scaling factor for this variable, this is to help  
        # normalize differences for large fields
        scale_factor = 1.0 / max( np.ceil( np.max(dset_avg[var]) ), 1.0 )

        for timestep in range( 0, timestep_full ):
                
            master_sum = master_sum + dset_fullres[var][timestep,:,:]
            samples = samples + 1

            if samples == calls_per_out or timestep == timestep_full-1:

                # Time index that we're comparing to in the averaged file
                avg_idx    = timestep // calls_per_out

                # Raw average of the master file over the sampling period
                master_avg = master_sum / samples

                # Average absolute difference of the master file averaged over the sampling period
                # and the average from the sample averaged netcdf file
                avg_abs_diff = np.average(abs(master_avg - dset_avg[var][avg_idx,:,:]))

                # Scale the differences 
                scaled_avg_abs_diff = avg_abs_diff * scale_factor

                # Report error if the scaled difference is above threshold
                if scaled_avg_abs_diff > abs_error_threshold:
                    print(f"{var}: master @ [{timestep+2-samples}:{timestep+1}] vs avg @ {avg_idx+1} : diff = {scaled_avg_abs_diff} ")
                    exit_code = 1

                # Record the maximum difference to report
                max_diff = max( scaled_avg_abs_diff, max_diff )
                
                # Reset sum and sample count
                master_sum = 0.0
                samples = 0

        print(f" - max diff {max_diff}")

sys.exit(exit_code)

            