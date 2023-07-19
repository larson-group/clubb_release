#!/usr/bin/python3

import netCDF4
import argparse
import sys
import os
import subprocess
import numpy as np
import tabulate

scriptPath = os.path.realpath(__file__)[0:-18]
#print(scriptPath)

# Threshold used to ignore field values for calculating % diff
field_threshold = 1.0e-7

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", action="store_true", help="add this option if you want to get output for each file compared put into a file made in the cwd.")
    parser.add_argument("-d", "--diffs", action="store_true", help="print the differences found in variables")
    parser.add_argument("-t", "--threshold", dest="threshold", type=float, action="store", help="print the differences found in variables")
    parser.add_argument("dirs", nargs=2, help="need 2 clubb output directories to diff. Usage: python run_bindiff_all.py dir_path1 dir_path2")

    args = parser.parse_args()

    paths_exist = os.path.exists(args.dirs[0]) and os.path.exists(args.dirs[1])
    
    if( paths_exist ):

        print("Directory 1 is", args.dirs[0])
        print("Directory 2 is", args.dirs[1])

        dir1_files = os.listdir(args.dirs[0])
        dir2_files = os.listdir(args.dirs[1])

        #print(dir1_files, dir2_files)

        cases = get_cases(dir1_files, dir2_files)

        if ( args.threshold is not None ):
          tot_abs_diff_thresh = args.threshold
        else:
          tot_abs_diff_thresh = 0.0

        print( "Using reporting threshold: ", tot_abs_diff_thresh, "\n" )

        diff_in_files = run_diff(cases, args.dirs[0], args.dirs[1], args.output, args.diffs, tot_abs_diff_thresh)

        if( args.output and diff_in_files ):
            print("\nOutput files made and placed in " + os.getcwd())

        if( diff_in_files ):
            print("\nThere were differences detected in netCDF (*.nc) files.")
            sys.exit(1)
        else:
            if( args.output ):
                print("\nNo differences detected, so no files were made!")
            else:
                print("\nNo differences detected!")
    else:
        print("Chosen directories do not exist. Please input valid directories")
        sys.exit(1)


def get_cases(dir1_files, dir2_files):

    cases = []

    with open(scriptPath + "RUN_CASES") as file:
        for line in file:
            #print(line.strip())
            case = line.strip()

            #If statement goes through every case in the RUN_CASES file and checks if the case is able to be run, 
            # and if that case has the files in the directories to diff it. If it passes all these conditions, 
            # the case is added to the runnable cases. 
            if(line[0] != "!" and len(line) > 1 and case + "_zt.nc" in dir1_files and case + "_zt.nc" in dir2_files
                    and case + "_zm.nc" in dir1_files and case + "_zm.nc" in dir2_files
                    and case + "_sfc.nc" in dir1_files and case + "_sfc.nc" in dir2_files):
                cases.append(case)
    #print(cases)
    return cases


def run_diff(cases, dir1_path, dir2_path, l_output, l_diffs, thresh):

    diff_in_files = False

    nc_data_formats = ["_zm.nc", "_zt.nc", "_sfc.nc"]

    # This for loop runs through all the cases you have the files to check, 
    # and each netcdf format. It runs the linux diff on the binary netcdf files 
    # to see which files should be looked at more closely
    for case in cases:

        print("DIFFING " + case + " netCDF (*.nc) files")

        for out_form in nc_data_formats:

            # Call linux diff command, if it shows a diff, check the netcdf values to confirm
            if( len(subprocess.getoutput("diff " + dir1_path + "/" + case + out_form + " " + dir2_path + "/" + case + out_form)) > 0 ):
                diff_in_files = run_py_diff(case+out_form, dir1_path, dir2_path, l_output, l_diffs, thresh) or diff_in_files

    return diff_in_files


def run_py_diff( test_file, dir1_path, dir2_path, l_save_output, l_print_diffs, tot_abs_diff_thresh ):

    diff_in_files = False
    any_diff_above_thresh = False

    #print("DIFFING " + test_file + "netCDF (*.nc) files")
    #print(dir1_path + test_file)
    #print(dir2_path + test_file)

    # Create datasets from nc files
    dset1 = netCDF4.Dataset(dir1_path+"/"+test_file)
    dset2 = netCDF4.Dataset(dir2_path+"/"+test_file)

    # Define table header, these are the values we're going to output
    table = [[  'Var', 
                'Max Abs Diff', 
                'Max % Diff', 
                'Total Abs Diff', 
                'Avg Abs Diff' ]]

    for var in dset1.variables:

        # The CLUBB netcdf variables we are interested in are all 
        # 4 dimensional (time, altitude, latitude, longitude),
        # but currently we don't actually have latitude or longitude in clubb, so those
        # dimensions are hardcoded to be 1. If in the future we remove those useless
        # dimensions (unlikely), then the variables of interested would be 2D. So 
        # for futureproofing, we will just check all variables with more than 1 dimension.
        if( dset1[var].ndim > 1 ):

            diff_in_files = True

            abs_diff = abs( dset1[var][:,:,:] - dset2[var][:,:,:] )

            # If the sum of all absolute differences is less than the threshold, then ignore this var
            if ( np.sum(abs_diff) < tot_abs_diff_thresh ):
              continue
            else:
              any_diff_above_thresh = True

            # Clip fields to ignore tiny values for the % diff
            field_1_clipped = np.clip( dset1[var][:,:,:], a_min = field_threshold, a_max = 9999999.0  )
            field_2_clipped = np.clip( dset2[var][:,:,:], a_min = field_threshold, a_max = 9999999.0 )

            # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
            percent_diff = 200.0 * ( field_1_clipped-field_2_clipped ) \
                                   / ( field_1_clipped+field_2_clipped )
                         
            # Append the table array with the values we want to print   
            table.append( [ var, 
                            np.max(abs_diff), 
                            np.max(percent_diff), 
                            np.sum(abs_diff), 
                            np.average(abs_diff) ] )

    # Print results 
    header = "************** Differences detected in " + test_file + "! **************"
    print( header )

    if( any_diff_above_thresh ):
      output = tabulate.tabulate(table, headers='firstrow')
    else:
      output = "Total absolute value of all differences are below threshold: " + str(tot_abs_diff_thresh)

    if( l_print_diffs ):
        # Print a very pretty table of the values
        print(output+"\n")

    if( l_save_output ):
        # Save all the output to a file
        f = open(test_file + "-diff_out", "w")
        f.write(header)
        f.write(output)
        f.close()

    return diff_in_files


if __name__ == "__main__":
    main()
