#!/usr/bin/python3

import netCDF4
import argparse
import sys
import os
import subprocess
import numpy as np
import tabulate

DEBUG = False

scriptPath = os.path.dirname(os.path.realpath(__file__))+'/'
if DEBUG:
    print(scriptPath)
outFilePath = os.path.normpath(scriptPath+"../output/bindiffs/")
filePostFix = "_diff{}.log"

## CONSTANTS ##
# Threshold used to ignore field values for calculating % diff
field_threshold = 1.0e-7

# NetCDF file postfixes we want to diff
nc_data_formats = ["_zm.nc", "_zt.nc", "_sfc.nc"]

# File containing all the case names
case_file = "RUN_CASES"

# Define output messages:
# 1. Sum of absolute differences below threshold for all variables in all files for a single case.
no_value_diff_out = ">>>Differences in linux diff but no differences in common variables detected for case {}.<<<"

def main():

    # Set up and parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--fileout", action="store", choices=["skip", "replace", "enumerate"], default=None,
                        help="Output a table of diffs into a text file, with one output file per each pair of diffed netcdf files.\n"+
                        "The text files will be named '<input_file>-diff_out' and will be stored in a subfolder of CLUBB's output folder.\n"+
                        "The choices for this option are:"+
                        "\n1. Flag omitted: No files will be created."+
                        "\n2. 'enumerate': If a diff log file for a case already exists\n\t\tthe new file will be created with an increasing numerical postfix."+
                        "\n3. 'skip': If a diff log file for a case already exists the creation of the log file will be skipped."+
                        "\n4. 'replace': If a diff log file for a case already exists the old file will be replaced with a new log file.")
    parser.add_argument("-v", "--verbose", action="store", type=int, default=1, help="Choose level of verbosity for outputs, i.e. what is printed to console.\n"+
                        "0: Output just a summary.\n"+
                        "1: Default. Add summarized results for each file.\n"+
                        "2: Add tables with detailed numerical differences in common variables for each file.")
    parser.add_argument("-t", "--threshold", dest="threshold", type=float, action="store", help="(float) Define the maximum absolute difference for an individual variable to be treated as different.")
    parser.add_argument("-g", "--ghostbuster", action="store_true", help="Perform a comparison that omits the 'ghost' level in '_zt.nc' output files.")
    parser.add_argument("dirs", nargs=2, help="Need 2 clubb output directories containing netCDF files with the same name to diff. Usage: python run_bindiff_all.py dir_path1 dir_path2")
    args = parser.parse_args()

    if DEBUG:
        print(args)

    # Check if folders exist
    paths_exist = os.path.exists(args.dirs[0]) and os.path.exists(args.dirs[1])

    if( paths_exist ):

        if args.verbose>=1:
            print("Directory 1 is", args.dirs[0])
            print("Directory 2 is", args.dirs[1])

        # Define difference detection threshold
        if ( args.threshold is not None ):
          tot_abs_diff_thresh = args.threshold
        else:
          tot_abs_diff_thresh = 0.0

        if args.verbose>=1:
            print( "Using reporting threshold: ", tot_abs_diff_thresh, "\n" )

        # Create folder in CLUBB output folder
        if args.fileout:
            if DEBUG:
                print("Creating output folder")
            if not os.path.exists(outFilePath):
                os.makedirs(outFilePath)

        linux_diff, diff_in_files, file_skipped = find_diffs_in_all_files(args.dirs[0], args.dirs[1], args.fileout, args.verbose, tot_abs_diff_thresh, args.ghostbuster)

        if DEBUG:
            print(linux_diff, diff_in_files, file_skipped)
        print("SUMMARY:")
        if linux_diff:
            if not diff_in_files:
                print("Differences detected by linux diff but no differences in the common variables of the netCDF files compared.")
            else:
                print("There were differences detected in the common variables in netCDF (*.nc) files.")
        else:
            print("Linux diff did not detect any differences in the compared files.")
        if args.fileout:
            if diff_in_files:
                if args.verbose>=1:
                    print("Output files made and placed in " + outFilePath)
                if file_skipped:
                    if args.verbose>=1:
                        print("The creation of some or all of the output files was skipped because files with the same name already existed. "+
                              "Use output options 'replace' or 'enumerate' to avoid this issue.")
            else:
                if args.verbose>=1:
                    print("No output files were created in " + outFilePath)
        if diff_in_files:
            sys.exit(1)

    else:
        print("Chosen directories do not exist. Please input valid directories")
        sys.exit(2)


def get_cases(dir1, dir2, verbose):
    # Generate and return a dict of all cases that have files present in the diff folders
    # This dict will have the following structure:
    # cases = {<case> : <folder_dict>, ...}
    # where folder_dict is another dict with two entries and the following structure:
    # folder_dict = {<dir>: [list, of, nc_files], ...}
    # where dir is either of the passed arguments and
    # [list, to, nc_files] contains the names of the existing netCDF files for <case> in the same order as nc_data_formats.
    # The entry will be None if a specific file does not exist.
    # Example: cases = {'arm': { <dir1>: ['arm_zm.nc', None, None], <dir2>: ['arm_zm.nc', None, None]},
    #                   'bomex': { <dir1>: ['bomex_zm.nc', 'bomex_zt.nc', 'bomex_sfc.nc']}, <dir2>: ['bomex_zm.nc', 'bomex_zt.nc', 'bomex_sfc.nc']},
    #                   'wangara': { <dir1>: [None, 'wangara_zt.nc',None]}, <dir2>: ['wangara_zm.nc', None, None]}
    #
    # Warnings will be printed for each missing file both on the console and in the output files.
    #####################################################################################

    # Get lists of files in each folder
    dir1_files = os.listdir(dir1)
    dir2_files = os.listdir(dir2)

    if DEBUG:
        print(dir1_files, dir2_files)

    # Create dict structure with case names as keys and a list for each folder of file names (or None if nonexistent) as values
    # This dict is returned from this function and will later be used as a basis for the comparison loop in 'find_diffs_in_all_files'
    cases = {}

    # Loop through all the case names listed in the RUN_CASES file and check if comparison is possible (i.e. at least one pair of case files is found in the folders)
    # First, open the RUN_CASES file
    with open(os.path.join(scriptPath, case_file)) as file:
        for line in file:
            if DEBUG:
                print(line)
            # Remove white space and the leading "!" in commented out lines
            # since we even want to compare commented out cases!
            #print(line.strip())
            case = line.strip()
            if DEBUG:
                print(case)
            # Skip line if empty
            if case == "" or case == "!":
                if DEBUG:
                    print("continue")
                continue
            # "Uncomment case name"
            if case[0] == "!":
                if DEBUG:
                    print("Commented")
                case = case[1:].strip()
                if DEBUG:
                    print(case)
            # Skip line if it contains multiple words
            if " " in case or "$" in case:
                if DEBUG:
                    print("continue")
                continue
            if DEBUG:
                print("Create entry")
            cases[case] = {dir1: [], dir2: []}

            # We want to compare all cases that have at least one common file in the folders, i.e. add that case to dict cases with a list of comparable files.
            # Note ANY missing files in a warning (also in log file)
            if DEBUG:
                print("iterating")
            for postfix in nc_data_formats:
                if DEBUG:
                    print(postfix)
                ncfname = case+postfix
                if DEBUG:
                    print(ncfname)
                # Check if files exist in either folder
                # Check dir1
                if ncfname in dir1_files:
                    if DEBUG:
                        print("Exists in " + dir1)
                    cases[case][dir1].append(ncfname)
                else:
                    if DEBUG:
                        print("Not in " + dir1)
                    cases[case][dir1].append(None)
                # Check dir2
                if ncfname in dir2_files:
                    if DEBUG:
                        print("Exists in " + dir2)
                    cases[case][dir2].append(ncfname)
                else:
                    if DEBUG:
                        print("Not in " + dir2)
                    cases[case][dir2].append(None)

            # Remove case from list if no FILES exist
            if not any(cases[case][dir1]) or not any(cases[case][dir2]):
                if DEBUG:
                    print("Removing case because no files exist.")
                    print(cases[case])
                del cases[case]
                if verbose>=2:
                    print("Warning! No files found for case {}. It will be skipped.".format(case))
            # Remove case from list if no file PAIRS exist
            elif not find_comparable_files(cases[case][dir1], cases[case][dir2]):
                if DEBUG:
                    print("Removing case because there are no pairs.")
                    print(cases[case])
                del cases[case]
                if verbose>=2:
                    print("Warning! There are files for case {} but no pairs for comparison were found. It will be skipped.".format(case))
    if verbose>=1:
        print("\nThe following cases will be compared: {}\n".format(list(cases.keys())))
    return cases

def find_comparable_files(l1, l2):
    return any([p1 and p2 for (p1,p2) in zip(l1,l2)])

def find_diffs_in_all_files(dir1, dir2, save_to_file, verbose, thresh, l_ghostbuster):
    # For each case with existing netCDF files in the diff folders:
    # 1. Create an output file if those are requested
    # 2. Loop through the netCDF files and call `find_diffs_in_common_vars` on each pair
    # 3. The cumulative boolean `diff_in_all_files` stores whether or not any differences have been found
    #
    # Parameters:
    #   dir1, dir2: Paths to folders containing netCDF files we want to compare
    #   save_to_file: Value of the --outfile command line argument
    #   verbose: Value of the --verbose command line argument, level of verbosity
    #   thresh: Detection threshold for differences in a variable
    #####################################################################################
    file_missing_msg = "Warning! The file {} does not exist in {}. No comparison possible.\n"

    # Find all cases that need to be compared (all files present in both folders)
    cases = get_cases(dir1, dir2, verbose)
    diff_in_all_files = False
    linux_diff = False
    file_skipped = False

    # This for loop runs through all the cases you have the files to check,
    # and each netcdf format. It runs the linux diff on the binary netcdf files
    # to see which files should be looked at more closely
    for case in cases:
        if verbose>=1:
            print("###DIFFING " + case + " netCDF (*.nc) files###")
        linux_diff_in_case = False
        diff_in_case = False

        # Write any output that should be written to file to string 'content'
        if save_to_file:
            content = "Output for 'run_bindiff_all.py' of the folders:\n"
            content += " - " + dir1
            content += "\n - " + dir2
            content += "\nDiffed case: " + case + "\n"
            content += "********************************************************************************************\n"

        # Looping over all files associated to <case>
        for i, postfix in enumerate(nc_data_formats):
            ncfname = case+postfix
            if not cases[case][dir1][i]:
                if verbose>=2:
                    print(file_missing_msg.format(ncfname, dir1))
                if save_to_file:
                    content += file_missing_msg.format(ncfname, dir1)
            if not cases[case][dir2][i]:
                if verbose>=2:
                    print(file_missing_msg.format(ncfname, dir2))
                if save_to_file:
                    content += file_missing_msg.format(ncfname, dir2)
            if not cases[case][dir1][i] or not cases[case][dir2][i]:
                if verbose>=1:
                    print(">File {} can not be compared since it does not exist both folders.<".format(ncfname))
                continue

            # Call linux diff command, if it shows a diff, check the netcdf values to confirm
            if( len(subprocess.getoutput("diff " + dir1 + "/" + ncfname + " " + dir2 + "/" + ncfname)) > 0 ):
                linux_diff_in_case = True
                if verbose>=1:
                    print(">The linux diff detected differences in " + ncfname + "<")
                if save_to_file:
                    content += ">The linux diff detected differences in " + ncfname + "<\n"
                # Update diff_in_all_files: if either it or the output of find_diffs_in_common_vars is true, we conclude that there is a difference. The printed messages will elaborate on the specific differences.
                diff_in_case, new_content = find_diffs_in_common_vars(ncfname, dir1, dir2, save_to_file, verbose, thresh, l_ghostbuster, postfix) or diff_in_case
                if verbose>=2:
                    print('')
                if save_to_file:
                    content += new_content
                    content += '\n'
            else:
                if verbose>=1:
                    print(">No differences detected by the linux diff in " + ncfname + "<")

        if diff_in_case:
            if verbose>=1:
                print(">>>Differences in common variables detected for case {}<<<".format(case))
            if save_to_file:
                # Create file to save diff log for <case> into
                diff_file_name = os.path.join(outFilePath, case + filePostFix)
                if DEBUG:
                    print(diff_file_name)
                # Check if we need to create a file
                if not os.path.exists(diff_file_name.format('')) or save_to_file=="replace":
                    # If file does not exist or should be overwritten, open file with default name and write content
                    with open(diff_file_name.format(''), "w") as caseLogFile:
                        caseLogFile.write(content)
                elif save_to_file == "enumerate":
                    # (File exists already) First find which file index is next and create the file with that index in its name
                    i = 1
                    diff_file_name = diff_file_name.format('{:02d}')
                    while os.path.exists(diff_file_name.format(i)):
                        i = i+1
                    with open(diff_file_name.format(i), "w") as caseLogFile:
                        caseLogFile.write(content)
                elif save_to_file == "skip":
                    file_skipped = True
                    if  verbose>=2:
                        print("Warning! Despite detected differences no log file was created for case {} since the file {} already exists and the 'skip' output option was used ".format(case, diff_file_name.format('')))
        elif linux_diff_in_case:
            if verbose>=1:
                print(no_value_diff_out.format(case))
        else:
            if verbose>=1:
                print(">>>The linux diff could not detect any differences in the file pairs for case {}.<<<".format(case))
        if verbose>=1:
            print('**********************************************************************************************************')

        # Update diff_in_all_files: Set to True if diff_in_case==True (differences in case found)
        linux_diff = linux_diff_in_case or linux_diff
        diff_in_all_files = diff_in_case or diff_in_all_files

    return (linux_diff, diff_in_all_files, file_skipped)

def find_diffs_in_common_vars( test_file, dir1, dir2, save_to_file, verbose, tot_abs_diff_thresh, l_ghostbuster, postfix ):
    # This is the integral function of this script!
    # Compare content of one specific pair of files with the same name in each folder:
    # 1. Find the variables that are present in only one of the files
    # 2. Find common variables and iterate through them
    # 3. Compare values between the files and detect differences
    # 4. Print findings
    # 5. Return TRUE if significant differences were detected and FALSE otherwise
    #####################################################################################

    # Assume no differences until proven otherwise
    diff_in_common_vars = False
    # Declare string containing all the output for <test_file>
    new_content = ""

    if DEBUG:
        print("DIFFING " + test_file + "netCDF (*.nc) files")
        print(os.path.join(dir1, test_file))
        print(os.path.join(dir2, test_file))

    # Create datasets from nc files
    dset1 = netCDF4.Dataset(os.path.join(dir1, test_file))
    dset2 = netCDF4.Dataset(os.path.join(dir2, test_file))

    # Find variables that are only present in ONE of the files
    diff1 = set(dset1.variables.keys()).difference(dset2.variables.keys())
    diff2 = set(dset2.variables.keys()).difference(dset1.variables.keys())
    # Print those variables
    if diff1:
        var_set_diff_out = "{} contains the following extra variables:".format(os.path.join(dir1, test_file))
        var_set_diff_out = var_set_diff_out + "\n" + "\n".join(diff1)
        if verbose>=1:
           print(var_set_diff_out)
        if save_to_file:
            new_content += var_set_diff_out + "\n"

    if diff2:
        var_set_diff_out = "{} contains the following extra variables:".format(os.path.join(dir2, test_file))
        var_set_diff_out = var_set_diff_out + "\n" + "\n".join(diff2)
        if verbose>=1:
            print(var_set_diff_out)
        if save_to_file:
            new_content += var_set_diff_out + "\n"
    if diff1 or diff2:
        if verbose>=1:
            print('')
        if save_to_file:
            new_content += "\n"

    # Define table header, these are the values we're going to output
    table = [[  'Var',
                'Max Abs Diff',
                'Max % Diff',
                'Total Abs Diff',
                'Avg Abs Diff' ]]

    # Create set of variables that are common to both files
    # These are the variables we can and want to compare the other variables are printed above
    vars_in_common = set(dset1.variables.keys()).intersection(dset2.variables.keys())
    for var in sorted(vars_in_common):
        # The CLUBB netcdf variables we are interested in are all
        # 4 dimensional (time, altitude, latitude, longitude),
        # but currently we don't actually have latitude or longitude in clubb, so those
        # dimensions are hardcoded to be 1. If in the future we remove those useless
        # dimensions (unlikely), then the variables of interested would be 2D. So
        # for futureproofing, we will just check all variables with more than 1 dimension.
        if( dset1[var].ndim > 1 and dset2[var].ndim > 1 ):

            if l_ghostbuster and postfix=="_zt.nc":
              z_levs_1=dset1[var].shape[1]
              z_levs_2=dset2[var].shape[1]
              if z_levs_1 == z_levs_2:
                abs_diff = abs( dset1[var][:,1:z_levs_1-1,:,:] - dset2[var][:,1:z_levs_2-1,:,:] )
              else:
                if z_levs_1 > z_levs_2:
                  abs_diff = abs( dset1[var][:,1:z_levs_1-1,:,:] - dset2[var][:,0:z_levs_2-1,:,:] )
                else:
                  abs_diff = abs( dset1[var][:,0:z_levs_1-1,:,:] - dset2[var][:,1:z_levs_2-1,:,:] )
            else:
              abs_diff = abs( dset1[var][:,:,:,:] - dset2[var][:,:,:,:] )

            # If the sum of all absolute differences is less than the threshold, then ignore this var
            if ( np.sum(abs_diff) <= tot_abs_diff_thresh ):
              continue
            else:
              diff_in_common_vars = True

            # Clip fields to ignore tiny values for the % diff
            if l_ghostbuster and postfix=="_zt.nc":
              if z_levs_1 == z_levs_2:
                field_1_clipped = np.clip( dset1[var][:,1:z_levs_1-1,:,:], a_min = field_threshold, a_max = 9999999.0  )
#yippee confetti field_2_clipped = np.clip( dset2[var][:,1:z_levs_2-1,:], a_min = field_threshold, a_max = 9999999.0 )
                field_2_clipped = np.clip( dset2[var][:,1:z_levs_2-1,:,:], a_min = field_threshold, a_max = 9999999.0 )
              else:
                if z_levs_1 > z_levs_2:
                  field_1_clipped = np.clip( dset1[var][:,1:z_levs_1-1,:,:], a_min = field_threshold, a_max = 9999999.0  )
                  field_2_clipped = np.clip( dset2[var][:,0:z_levs_2-1,:,:], a_min = field_threshold, a_max = 9999999.0 )
                else:
                  field_1_clipped = np.clip( dset1[var][:,0:z_levs_1-1,:,:], a_min = field_threshold, a_max = 9999999.0  )
                  field_2_clipped = np.clip( dset2[var][:,1:z_levs_2-1,:,:], a_min = field_threshold, a_max = 9999999.0 )
            else:
              field_1_clipped = np.clip( dset1[var][:,:,:,:], a_min = field_threshold, a_max = 9999999.0  )
              field_2_clipped = np.clip( dset2[var][:,:,:,:], a_min = field_threshold, a_max = 9999999.0 )

            # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
            percent_diff = 200.0 * ( field_1_clipped-field_2_clipped ) \
                                   / ( field_1_clipped+field_2_clipped )

            # Append the table array with the values we want to print
            table.append( [ var,
                            np.max(abs_diff),
                            np.max(percent_diff),
                            np.sum(abs_diff),
                            np.average(abs_diff) ] )
        else:
            if verbose>=2:
                print("The variable {} is not comparable because its dimension is not fit for comparison: {}, {}".format(var, dset1[var].ndim, dset2[var].ndim))


    if diff_in_common_vars:
      output = tabulate.tabulate(table, headers='firstrow')
    else:
      output = "Total absolute value of all differences are below threshold: " + str(tot_abs_diff_thresh)

    if verbose>=2:
        # Print a very pretty table of the values
        print(output+"\n")

    if save_to_file:
        # Save all the output to a file
        new_content += output + "\n"

    if not diff_in_common_vars:
        if verbose>=1:
            print(">>No differences detected in the common fields in file " + test_file + "<<")
    else:
        if verbose>=1:
            print(">>Differences above threshold were detected in the common fields in file " + test_file + "<<")
    return (diff_in_common_vars, new_content)


if __name__ == "__main__":
    main()
