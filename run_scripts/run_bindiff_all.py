#!/usr/bin/python3

import netCDF4
import argparse
import sys
import os
import subprocess
import numpy as np
import tabulate
import multiprocessing as mp
import io
import contextlib

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
nc_data_formats = [
    "_stats.nc",
    "_zm.nc",
    "_zt.nc",
    "_sfc.nc",
    "_lh_zt.nc",
    "_lh_sfc.nc",
    "_nl_lh_sample_points_2D.nc",
    "_u_lh_sample_points_2D.nc",
    "_rad_zt.nc",
    "_rad_zm.nc",
    "_multi_col_zm.nc",
    "_multi_col_zt.nc",
]

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
    parser.add_argument("-s", "--scale", action="store_true", help="Scale absolute differences by the average field value, using the same approach as check_multi_col_error.py.")
    parser.add_argument("dirs", nargs=2, help="Need 2 clubb output directories containing netCDF files with the same name to diff. Usage: python run_bindiff_all.py dir_path1 dir_path2")
    args = parser.parse_args()

    if DEBUG:
        print(args)

    # Check if folders exist
    paths_exist = os.path.exists(args.dirs[0]) and os.path.exists(args.dirs[1])

    if not paths_exist:
        print("Chosen directories do not exist. Please input valid directories")
        sys.exit(2)

    if os.path.samefile(args.dirs[0], args.dirs[1]):
        print("Input paths resolve to the same directory.")
        sys.exit(2)

    if args.verbose>=1:
        print("Directory 1 is", args.dirs[0])
        print("Directory 2 is", args.dirs[1])

    # Define difference detection threshold
    if ( args.threshold is not None ):
        abs_error_threshold = args.threshold
    else:
        abs_error_threshold = 0.0

    if args.verbose>=1:
        print( "Using reporting threshold: ", abs_error_threshold, "\n" )
        if args.scale:
            print("Absolute differences will be scaled by average field value.\n")

    # Create folder in CLUBB output folder
    if args.fileout:
        if DEBUG:
            print("Creating output folder")
        if not os.path.exists(outFilePath):
            os.makedirs(outFilePath)

    linux_diff, diff_in_files, file_skipped = find_diffs_in_all_files(
        args.dirs[0], args.dirs[1], args.fileout, args.verbose, abs_error_threshold, args.scale
    )

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



def get_cases(dir1, dir2, verbose):
    # Generate and return a dict of all cases that have files present in the diff folders
    # This dict will have the following structure:
    # cases = {<case> : [file1, file2, ...]}
    # Example: cases = {'arm': ['arm_stats.nc'],
    #                   'bomex': ['bomex_stats.nc'],
    #                   'wangara': ['wangara_stats.nc']}
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

            # Remove leading or trailing white space and the leading "!" in commented out lines
            # since we even want to compare commented out cases!
            case = line.strip().lstrip("!")

            # Skip line if empty
            if case == "" or case == "!":
                if DEBUG:
                    print("continue")
                continue

            # Skip line if it contains multiple words
            if " " in case or "$" in case:
                if DEBUG:
                    print("continue")
                continue

            if DEBUG:
                print(f"Adding {case} to case list")

            cases[case] = []

            common_files_found = False

            # We want to compare all cases that have at least one common file in the folders, 
            # i.e. add that case to dict cases with a list of comparable files.
            # Note ANY missing files in a warning (also in log file)
            for postfix in nc_data_formats:

                ncfname = case+postfix

                # Check if ncfname exists in both folders
                if ncfname in dir1_files and ncfname in dir2_files:
                    if DEBUG:
                        print(ncfname + " exists in " + dir1 + " and " + dir2)
                    cases[case].append(ncfname)
                    common_files_found = True
                else:
                    if DEBUG:
                        if ncfname not in dir1_files:
                            print(f"{ncfname} not in {dir1}")
                        if ncfname not in dir2_files:
                            print(f"{ncfname} not in {dir2}")

            # Remove the case if no common files were found 
            if not common_files_found:
                if DEBUG:
                    print("Removing case because there are no pairs to compare.")
                del cases[case]

    if DEBUG:
        print("cases dict :")
        print("\t", cases)
    if verbose>=1:
        print("\nThe following cases will be compared: {}\n".format(list(cases.keys())))

    return cases

def find_diffs_in_all_files(dir1, dir2, save_to_file, verbose, thresh, l_scale):
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

    case_order = list(cases.keys())

    # Determine worker count (nproc - 2, but at least 1 and no more than number of cases)
    cpu_total = os.cpu_count() or 1
    nproc = max(1, cpu_total - 2)
    nproc = min(nproc, max(1, len(case_order)))

    args_list = [
        (case, cases[case], dir1, dir2, save_to_file, verbose, thresh, l_scale)
        for case in case_order
    ]

    if nproc > 1:
        try:
            with mp.Pool(processes=nproc) as pool:
                results = pool.map(_diff_case, args_list)
        except (PermissionError, OSError):
            # Some CI/runtime environments disallow multiprocessing semaphores.
            # Fall back to serial comparison in that case.
            results = [_diff_case(args) for args in args_list]
    else:
        results = [_diff_case(args) for args in args_list]

    # Emit stdout in a stable order and aggregate results
    for result in results:
        if result["stdout"]:
            print(result["stdout"], end="")
        linux_diff = result["linux_diff_in_case"] or linux_diff
        diff_in_all_files = result["diff_in_case"] or diff_in_all_files
        file_skipped = result["file_skipped"] or file_skipped

    return (linux_diff, diff_in_all_files, file_skipped)


def _diff_case(args):
    case, files, dir1, dir2, save_to_file, verbose, thresh, l_scale = args
    linux_diff_in_case = False
    diff_in_case = False
    file_skipped = False

    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        if verbose >= 1:
            print("###DIFFING " + case + " netCDF (*.nc) files###")

        # Write any output that should be written to file to string 'content'
        if save_to_file:
            content = "Output for 'run_bindiff_all.py' of the folders:\n"
            content += " - " + dir1
            content += "\n - " + dir2
            content += "\nDiffed case: " + case + "\n"
            content += "********************************************************************************************\n"
        else:
            content = ""

        # Looping over all files associated to <case>
        for ncfname in files:

            # Call linux diff command, if it shows a diff, check the netcdf values to confirm
            if len(subprocess.getoutput("diff " + dir1 + "/" + ncfname + " " + dir2 + "/" + ncfname)) > 0:
                linux_diff_in_case = True
                if verbose >= 1:
                    print(">The linux diff detected differences in " + ncfname + "<")
                if save_to_file:
                    content += ">The linux diff detected differences in " + ncfname + "<\n"

                case_diff, new_content = find_diffs_in_common_vars(
                    ncfname, dir1, dir2, save_to_file, verbose, thresh, l_scale
                )
                diff_in_case = case_diff or diff_in_case

                if verbose >= 2:
                    print("")
                if save_to_file:
                    content += new_content
                    content += "\n"
            else:
                if verbose >= 1:
                    print(">No differences detected by the linux diff in " + ncfname + "<")

        if diff_in_case:
            if verbose >= 1:
                print(">>>Differences in common variables detected for case {}<<<".format(case))
            if save_to_file:
                # Create file to save diff log for <case> into
                diff_file_name = os.path.join(outFilePath, case + filePostFix)
                if DEBUG:
                    print(diff_file_name)
                # Check if we need to create a file
                if not os.path.exists(diff_file_name.format('')) or save_to_file == "replace":
                    # If file does not exist or should be overwritten, open file with default name and write content
                    with open(diff_file_name.format(''), "w") as caseLogFile:
                        caseLogFile.write(content)
                elif save_to_file == "enumerate":
                    # (File exists already) First find which file index is next and create the file with that index in its name
                    i = 1
                    diff_file_name = diff_file_name.format('{:02d}')
                    while os.path.exists(diff_file_name.format(i)):
                        i = i + 1
                    with open(diff_file_name.format(i), "w") as caseLogFile:
                        caseLogFile.write(content)
                elif save_to_file == "skip":
                    file_skipped = True
                    if verbose >= 2:
                        print("Warning! Despite detected differences no log file was created for case {} since the file {} already exists and the 'skip' output option was used ".format(case, diff_file_name.format('')))
        elif linux_diff_in_case:
            if verbose >= 1:
                print(no_value_diff_out.format(case))
        else:
            if verbose >= 1:
                print(">>>The linux diff could not detect any differences in the file pairs for case {}.<<<".format(case))
        if verbose >= 1:
            print('**********************************************************************************************************')

    return {
        "case": case,
        "linux_diff_in_case": linux_diff_in_case,
        "diff_in_case": diff_in_case,
        "file_skipped": file_skipped,
        "stdout": stdout.getvalue(),
    }

def find_diffs_in_common_vars( test_file, dir1, dir2, save_to_file, verbose, abs_error_threshold, l_scale ):
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

    # Create datasets from nc files.
    # Disable netCDF4 auto-masking so _FillValue=0.0 fields are compared as raw values.
    dset1 = netCDF4.Dataset(os.path.join(dir1, test_file))
    dset2 = netCDF4.Dataset(os.path.join(dir2, test_file))
    dset1.set_auto_mask(False)
    dset2.set_auto_mask(False)

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
        # Compare multidimensional fields only. Scalar/1D values are ignored.
        if( dset1[var].ndim > 1 and dset2[var].ndim > 1 ):
            data_1 = np.asarray(dset1[var][...])
            data_2 = np.asarray(dset2[var][...])

            try:
                abs_diff = abs(data_1 - data_2)
            except ValueError:
                if verbose >= 1:
                    print("Skipping variable {} because shapes are incompatible: {} vs {}".format(
                        var, data_1.shape, data_2.shape
                    ))
                continue

            # Match check_multi_col_error.py behavior:
            # scale absolute differences by average field magnitude.
            if l_scale:
              field_avg = ( np.average( data_1 ) + np.average( data_2 ) ) / 2.0
              abs_diff = abs_diff / np.ceil( np.abs( field_avg ) )

            # If the average absolute differences is less than the threshold, then ignore this var
            if ( np.average(abs_diff) <= abs_error_threshold ):
              continue
            else:
              diff_in_common_vars = True

            # Clip fields to ignore tiny values for the % diff
            field_1_clipped = np.clip( data_1, a_min = field_threshold, a_max = 9999999.0  )
            field_2_clipped = np.clip( data_2, a_min = field_threshold, a_max = 9999999.0 )

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
      output = "Total absolute value of all differences are below threshold: " + str(abs_error_threshold)

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
    dset1.close()
    dset2.close()
    return (diff_in_common_vars, new_content)


if __name__ == "__main__":
    main()
