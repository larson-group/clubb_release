#!/usr/bin/python3
"""
run_bindiff_all.py — Binary diff tool for CLUBB netCDF output.

Compares two directories of CLUBB simulation output and reports numerical
differences in the netCDF variables they share. Designed for regression
testing: run the same case(s) with two different builds or code versions,
then point this script at the two output directories.

How it works:
  1. Discovers cases by finding netCDF files common to both directories,
     grouping them by case name (the filename prefix before the stats suffix).
  2. For each case, does a fast byte-level comparison (filecmp) first.
     If files are byte-identical, no further work is needed.
  3. When bytes differ, opens both files with netCDF4, finds the variables
     common to both, and computes per-variable difference metrics (max
     absolute diff, percent diff, earliest diverging timestep, etc.).
  4. Variables whose average absolute difference falls below the active
     thresholds are silently skipped; the rest are reported in a table.
  5. Cases are compared in parallel (one worker per case) for speed;
     output is buffered per-case so results print in a stable order.

File format compatibility:
  The script recognizes both the legacy per-grid stats files (_zm.nc,
  _zt.nc, _sfc.nc, etc.) and the unified _stats.nc format introduced in
  Feb 2026. However, the two formats are not cross-comparable — a legacy
  directory should be diffed against another legacy directory, and a
  unified-stats directory against another unified-stats directory. The
  script matches files by exact filename, so mismatched formats will
  simply have no common files for that case and be skipped.

Usage:
  python run_bindiff_all.py <dir1> <dir2> [options]

  Options:
    -v {0,1,2}    Verbosity (0=summary, 1=per-file, 2=full diff tables)
    -t THRESHOLD   Min avg absolute diff to report a variable
    -pt THRESHOLD  Min avg absolute percent diff to report a variable
    -s             Scale diffs by average field magnitude
    -case CASE     Compare only the named case (e.g. 'bomex')
    -f {skip,replace,enumerate}  Write per-case diff logs to output/bindiffs/
"""

import netCDF4
import argparse
import sys
import os
import filecmp
import numpy as np
import tabulate
import multiprocessing as mp
import io
import contextlib

scriptPath = os.path.dirname(os.path.realpath(__file__))+'/'
outFilePath = os.path.normpath(scriptPath+"../output/bindiffs/")
filePostFix = "_diff{}.log"

## CONSTANTS ##
# Threshold used to ignore field values for calculating % diff
field_threshold = 1.0e-7

# NetCDF file postfixes we want to diff.
# CLUBB's stats system was overhauled (Feb 2026) to write a single unified
# _stats.nc file per case. Before that, each grid (zt, zm, sfc, rad, etc.)
# was written to its own file, so legacy output directories may contain any
# combination of these suffixes. We list them all so that old-vs-old and
# old-vs-new comparisons still work.
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

# Define output messages:
# 1. Sum of absolute differences below threshold for all variables in all files for a single case.
no_value_diff_out = ">>>Differences in linux diff but no differences in common variables detected for case {}.<<<"

# Candidate dimension names for each axis type. Legacy stats files used
# grid-specific dimension names (e.g. "zt", "zm", "lh_zt") while the
# unified _stats.nc uses generic names like "altitude" or "z". We accept
# all known variants so the location formatter works regardless of which
# generation of output is being compared.
TIME_DIM_NAMES = {"time", "t"}
VERTICAL_DIM_NAMES = {"zt", "zm", "z", "lh_zt", "altitude", "height", "lev"}
COLUMN_DIM_NAMES = {"ngrdcol", "col", "column", "columns"}

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
    parser.add_argument("-pt", "--percent_thresh", dest="percent_thresh", type=float, action="store", help="(float) Define the maximum average absolute percent difference for an individual variable to be treated as different.")
    parser.add_argument("-s", "--scale", action="store_true", help="Scale absolute differences by the average field value, using the same approach as check_multi_col_error.py.")
    parser.add_argument("-case", "--case", action="store", default=None, help="Compare only the specified case name (e.g. 'bomex'). When omitted, all cases found in both directories are compared.")
    parser.add_argument("dirs", nargs=2, help="Need 2 clubb output directories containing netCDF files with the same name to diff. Usage: python run_bindiff_all.py dir_path1 dir_path2")
    args = parser.parse_args()

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
    if ( args.percent_thresh is not None ):
        percent_error_threshold = args.percent_thresh
    else:
        percent_error_threshold = None

    if args.verbose>=1:
        print( "Using reporting threshold: ", abs_error_threshold, "\n" )
        if percent_error_threshold is not None:
            print( "Using percent reporting threshold: ", percent_error_threshold, "\n" )
        if args.scale:
            print("Absolute differences will be scaled by average field value.\n")

    # Create folder in CLUBB output folder
    if args.fileout:
        if not os.path.exists(outFilePath):
            os.makedirs(outFilePath)

    linux_diff, diff_in_files, file_skipped, passed_cases, failed_cases = find_diffs_in_all_files(
        args.dirs[0], args.dirs[1], args.fileout, args.verbose, abs_error_threshold, percent_error_threshold, args.scale,
        case_filter=args.case,
    )

    print("\nSUMMARY:")
    if linux_diff:
        if not diff_in_files:
            print("Differences detected by linux diff but no differences in the common variables of the netCDF files compared.")
        else:
            print("There were differences detected in the common variables in netCDF (*.nc) files.")
    else:
        print("Linux diff did not detect any differences in the compared files.")
    if args.verbose>=1:
        print("Pass: {} cases {}".format(len(passed_cases), passed_cases))
        print("Fail: {} cases {}".format(len(failed_cases), failed_cases))
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


def _find_axis_index(dimensions, candidate_names):
    """Map a variable's dimension tuple to a known axis type (time, vertical,
    column). Returns the positional index of the first dimension whose name
    matches any candidate, or None. This indirection is needed because legacy
    and modern stats files use different dimension names for the same axis."""
    for axis, dim_name in enumerate(dimensions):
        if dim_name.lower() in candidate_names:
            return axis
    return None


def _coord_value(dataset, dim_name, index):
    """Read the physical coordinate value at `index` along a named dimension.
    Returns None when the dimension has no corresponding coordinate variable
    (e.g. an unnamed index dimension). Clamps out-of-bounds indices so callers
    don't need to guard against edge cases in max-diff lookups."""
    if dim_name not in dataset.variables:
        return None
    values = np.asarray(dataset.variables[dim_name][...])
    if values.size == 0:
        return None
    safe_index = max(0, min(int(index), values.shape[0] - 1))
    return float(values[safe_index])


def _time_value_seconds(dataset, dim_name, index):
    """Convert a time coordinate value to seconds, normalizing across the
    different unit conventions that various stats files may use (days, hours,
    minutes, or already seconds)."""
    value = _coord_value(dataset, dim_name, index)
    if value is None:
        return None
    units = ""
    if dim_name in dataset.variables:
        units = getattr(dataset.variables[dim_name], "units", "").lower()
    if "day" in units:
        return value * 24.0 * 3600.0
    if "hour" in units:
        return value * 3600.0
    if "min" in units:
        return value * 60.0
    return value


def _height_value_meters(dataset, dim_name, index):
    """Convert a vertical coordinate value to meters, handling the unit
    variations that appear across legacy and modern output files."""
    value = _coord_value(dataset, dim_name, index)
    if value is None:
        return None
    units = ""
    if dim_name in dataset.variables:
        units = getattr(dataset.variables[dim_name], "units", "").lower()
    if units in {"km", "kilometer", "kilometers"}:
        return value * 1000.0
    if units in {"cm", "centimeter", "centimeters"}:
        return value / 100.0
    if units in {"mm", "millimeter", "millimeters"}:
        return value / 1000.0
    return value


def _format_max_diff_location(abs_diff, dimensions, dataset):
    """Build a human-readable string describing where the largest difference
    occurs (e.g. "t=5 (1800s), z=12 (500m), col=0"). This turns a raw
    multi-dimensional index into physical coordinates so the user can quickly
    locate the divergence point in the simulation output."""
    if abs_diff.size == 0:
        return ""

    max_index = np.unravel_index(np.argmax(abs_diff), abs_diff.shape)
    time_axis = _find_axis_index(dimensions, TIME_DIM_NAMES)
    vertical_axis = _find_axis_index(dimensions, VERTICAL_DIM_NAMES)
    column_axis = _find_axis_index(dimensions, COLUMN_DIM_NAMES)

    location_parts = []
    if time_axis is not None:
        time_index = int(max_index[time_axis])
        time_dim = dimensions[time_axis]
        time_sec = _time_value_seconds(dataset, time_dim, time_index)
        if time_sec is None:
            location_parts.append("t={}".format(time_index))
        else:
            location_parts.append("t={} ({:g}s)".format(time_index, time_sec))
    if vertical_axis is not None:
        z_index = int(max_index[vertical_axis])
        z_dim = dimensions[vertical_axis]
        z_m = _height_value_meters(dataset, z_dim, z_index)
        if z_m is None:
            location_parts.append("z={}".format(z_index))
        else:
            location_parts.append("z={} ({:g}m)".format(z_index, z_m))
    if column_axis is not None:
        location_parts.append("col={}".format(int(max_index[column_axis])))
    return ", ".join(location_parts)



def get_cases(dir1, dir2, verbose, case_filter=None):
    # Generate and return a dict of all comparable case files present in both diff folders.
    # The dict has the structure:
    # cases = {<case> : [file1, file2, ...]}
    #
    # Files are discovered directly from the output directories rather than RUN_CASES.
    # If case_filter is set, only return files belonging to that case.
    #####################################################################################

    dir1_files = set(os.listdir(dir1))
    dir2_files = set(os.listdir(dir2))
    common_files = sorted(dir1_files.intersection(dir2_files))

    cases = {}
    sorted_suffixes = sorted(nc_data_formats, key=len, reverse=True)
    for ncfname in common_files:
        matched_suffix = None
        for postfix in sorted_suffixes:
            if ncfname.endswith(postfix):
                matched_suffix = postfix
                break
        if matched_suffix is None:
            continue
        case = ncfname[: -len(matched_suffix)]
        if not case:
            continue
        if case_filter is not None and case != case_filter:
            continue
        cases.setdefault(case, []).append(ncfname)

    if verbose>=1:
        print("\nThe following cases will be compared: {}\n".format(list(cases.keys())))

    return cases

def find_diffs_in_all_files(dir1, dir2, save_to_file, verbose, thresh, percent_thresh, l_scale, case_filter=None):
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
    #   case_filter: If set, only compare this specific case name
    #####################################################################################
    # Find all cases that need to be compared (all files present in both folders)
    cases = get_cases(dir1, dir2, verbose, case_filter=case_filter)
    diff_in_all_files = False
    linux_diff = False
    file_skipped = False
    passed_cases = []
    failed_cases = []

    case_order = list(cases.keys())
    if not case_order:
        print("No comparable netCDF files were found in the provided directories.")
        sys.exit(2)

    # Determine worker count (nproc - 2, but at least 1 and no more than number of cases)
    cpu_total = os.cpu_count() or 1
    nproc = max(1, cpu_total)
    nproc = min(nproc, max(1, len(case_order)))

    args_list = [
        (case, cases[case], dir1, dir2, save_to_file, verbose, thresh, percent_thresh, l_scale)
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
        if result["diff_in_case"]:
            failed_cases.append(result["case"])
        else:
            passed_cases.append(result["case"])

    return (linux_diff, diff_in_all_files, file_skipped, passed_cases, failed_cases)


def _diff_case(args):
    """Process a single case (all its file pairs) in a worker-friendly way.
    Captures stdout to a StringIO buffer so that parallel workers don't
    interleave their output — the caller replays buffers in case order."""
    case, files, dir1, dir2, save_to_file, verbose, thresh, percent_thresh, l_scale = args
    linux_diff_in_case = False
    diff_in_case = False
    file_skipped = False

    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):

        if verbose >= 1:
            print(f"\n================================================ {case} ================================================")

        # Write any output that should be written to file to string 'content'
        if save_to_file:
            content = "Output for 'run_bindiff_all.py' of the folders:\n"
            content += " - " + dir1
            content += "\n - " + dir2
            content += "\nDiffed case: " + case + "\n"
            content += "==================================================================================\n"
        else:
            content = ""

        # Looping over all files associated to <case>
        for ncfname in files:
            file1 = os.path.join(dir1, ncfname)
            file2 = os.path.join(dir2, ncfname)

            # Call linux diff command, if it shows a diff, check the netcdf values to confirm
            if not filecmp.cmp(file1, file2, shallow=False):
                linux_diff_in_case = True
                if verbose >= 1:
                    print(">The linux diff detected differences in " + ncfname + "<")
                if save_to_file:
                    content += ">The linux diff detected differences in " + ncfname + "<\n"

                case_diff, new_content = find_diffs_in_common_vars(
                    ncfname, dir1, dir2, save_to_file, verbose, thresh, percent_thresh, l_scale
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
            print(f"\n=====================================================================================================")

    return {
        "case": case,
        "linux_diff_in_case": linux_diff_in_case,
        "diff_in_case": diff_in_case,
        "file_skipped": file_skipped,
        "stdout": stdout.getvalue(),
    }

def find_diffs_in_common_vars( test_file, dir1, dir2, save_to_file, verbose, abs_error_threshold, percent_error_threshold, l_scale ):
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

    # Create datasets from nc files.
    # Disable netCDF4 auto-masking so _FillValue=0.0 fields are compared as raw values.
    dset1 = netCDF4.Dataset(os.path.join(dir1, test_file))
    dset2 = netCDF4.Dataset(os.path.join(dir2, test_file))
    dset1.set_auto_mask(False)
    dset2.set_auto_mask(False)

    # Stop early if the file-pair does not have the same number of timesteps.
    if "time" in dset1.dimensions and "time" in dset2.dimensions:
        ntime_1 = len(dset1.dimensions["time"])
        ntime_2 = len(dset2.dimensions["time"])
        if ntime_1 != ntime_2:
            timestep_warning = (
                "\nWarning! {} and {} have different numbers of timesteps ({} vs {})."
            ).format(
                os.path.join(dir1, test_file),
                os.path.join(dir2, test_file),
                ntime_1,
                ntime_2,
            )
            if verbose >= 1:
                print(timestep_warning)
            if save_to_file:
                new_content += timestep_warning + "\n"
            dset1.close()
            dset2.close()
            return (True, new_content)

    # Find variables that are only present in ONE of the files
    diff1 = set(dset1.variables.keys()).difference(dset2.variables.keys())
    diff2 = set(dset2.variables.keys()).difference(dset1.variables.keys())
    # Print those variables
    if diff1:
        var_set_diff_out = "\n{} contains the following extra variables:".format(os.path.join(dir1, test_file))
        var_set_diff_out = var_set_diff_out + "\n" + "\n".join(diff1)
        if verbose>=1:
           print(var_set_diff_out)
        if save_to_file:
            new_content += var_set_diff_out + "\n"

    if diff2:
        var_set_diff_out = "\n{} contains the following extra variables:".format(os.path.join(dir2, test_file))
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
                'Max Diff Location',
                'Max % Diff',
                'Avg Abs Diff',
                'Avg % Diff',
                'Earliest Timestep' ]]
    n_vars_compared = 0
    n_vars_with_nonzero_diff = 0
    n_vars_exceeding_threshold = 0

    # Create set of variables that are common to both files
    # These are the variables we can and want to compare the other variables are printed above
    vars_in_common = set(dset1.variables.keys()).intersection(dset2.variables.keys())
    for var in sorted(vars_in_common):

        n_vars_compared += 1

        data_1 = np.asarray(dset1[var][...])
        data_2 = np.asarray(dset2[var][...])
        dimensions = dset1[var].dimensions

        if ( not np.issubdtype(data_1.dtype, np.number)
             or not np.issubdtype(data_2.dtype, np.number) ):
            if verbose >= 2:
                print("Skipping variable {} because it is not numeric: {} vs {}".format(
                    var, data_1.dtype, data_2.dtype
                ))
            continue

        try:
            abs_diff = abs(data_1 - data_2)
        except ValueError:
            print("The variable {} is not comparable because the shapes do not match: {}, {}".format(var, data_1.shape, data_2.shape))
            continue

        if not np.all(abs_diff == 0):
            n_vars_with_nonzero_diff += 1

        sum_abs_diff = np.sum(abs_diff)
        avg_abs_diff = sum_abs_diff / abs_diff.size

        # Match check_multi_col_error.py behavior:
        # scale absolute differences by average field magnitude.
        if l_scale:
            field_avg = ( np.average( data_1 ) + np.average( data_2 ) ) / 2.0
            abs_diff = abs_diff / np.ceil( np.abs( field_avg ) )
            sum_abs_diff = np.sum(abs_diff)
            avg_abs_diff = sum_abs_diff / abs_diff.size

        # Fast path: if the absolute threshold already rejects this variable,
        # we do not need clipping or percent-difference calculations.
        if avg_abs_diff <= abs_error_threshold:
            continue

        # Clip fields to ignore tiny values for the % diff
        field_1_clipped = np.clip( data_1, a_min = field_threshold, a_max = 9999999.0  )
        field_2_clipped = np.clip( data_2, a_min = field_threshold, a_max = 9999999.0 )

        # Calculate the percent difference, 100 * (a-b) / ((a+b)/2)
        percent_diff = 200.0 * ( field_1_clipped-field_2_clipped ) \
                                / ( field_1_clipped+field_2_clipped )
        avg_abs_percent_diff = np.average(np.abs(percent_diff))

        # Ignore the variable if it falls below either active reporting threshold.
        if percent_error_threshold is not None and avg_abs_percent_diff <= percent_error_threshold:
            continue

        diff_in_common_vars = True
        n_vars_exceeding_threshold += 1

        earliest_timestep = ""
        # Combine both thresholds: a grid point must exceed the absolute
        # threshold AND the percent threshold at the same timestep.
        exceed_mask = abs_diff > abs_error_threshold
        if percent_error_threshold is not None:
            exceed_mask = exceed_mask & (np.abs(percent_diff) > percent_error_threshold)
        if abs_diff.ndim >= 1 and np.any(exceed_mask):
            timestep_mask = np.any(exceed_mask, axis=tuple(range(1, abs_diff.ndim)))
            timestep_indices = np.where(timestep_mask)[0]
            if timestep_indices.size > 0:
                ts_idx = int(timestep_indices[0])
                time_axis = _find_axis_index(dimensions, TIME_DIM_NAMES)
                if time_axis is not None:
                    time_dim = dimensions[time_axis]
                    time_sec = _time_value_seconds(dset1, time_dim, ts_idx)
                    if time_sec is not None:
                        earliest_timestep = "{} ({:g}s)".format(ts_idx, time_sec)
                    else:
                        earliest_timestep = str(ts_idx)
                else:
                    earliest_timestep = str(ts_idx)
        max_diff_location = _format_max_diff_location(abs_diff, dimensions, dset1)

        # Append the table array with the values we want to print
        table.append( [ var,
                        np.max(abs_diff),
                        max_diff_location,
                        np.max(percent_diff),
                        avg_abs_diff,
                        avg_abs_percent_diff,
                        earliest_timestep ] )

    if diff_in_common_vars:
      output = tabulate.tabulate(table, headers='firstrow')
    else:
      output = "All common-variable differences were below the active thresholds."

    summary_output = (
        "Comparison summary:\n"
        + "Variables compared: " + str(n_vars_compared) + "\n"
        + "Variables with non-zero differences: " + str(n_vars_with_nonzero_diff) + "\n"
        + "Variables exceeding threshold: " + str(n_vars_exceeding_threshold)
    )

    if verbose>=2:
        # Print a very pretty table of the values
        print(output+"\n")
        print(summary_output + "\n")

    if save_to_file:
        # Save all the output to a file
        new_content += output + "\n"
        new_content += summary_output + "\n"

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
