"""
Run CLUBB with various flags are toggled.

In order to run this script you need to install the dependencies with 
    pip install -r run_bindiff_w_flags_requirements.txt

The simplest way to use this script would look like:
  python3 ./run_clubb_w_varying_flags.py -f <flag_config_file>

Inputs:
  1. JSON configuration file, describing what flags you want to toggle

All different input flag configuration files get generated and stored in the
input/tunable_parameters directory, depending on the input file passed with
--flag-config-file.  This input file defines what flags should be adjusted,
taking the configurable_model_flags.in file as the default.

An example for this file is the run_bindiff_w_flags_config_example.json found
in the run_scripts directory.

The code then gets run with all the different flag files.
"""

import shutil
import sys
import os
import re
import subprocess
import argparse
import json
import git


# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(",")


def get_cli_args():
    # Set up and parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--short-cases",
        action="store_true",
        default=False,
        help="Run short cases. This will omit\n\t\tthe gabls2, cloud_feedback_s6, "
           + "cloud_feedback_s11,\n\t\tcloud_feedback_s12, and twp_ice cases.",
    )
    parser.add_argument(
        "--priority-cases",
        action="store_true",
        default=False,
        help="Run priority cases. This will include\n\t\tonly the following cases if they are "
           + "in RUN_CASES: arm, atex,\n\t\tbomex, dycoms2_rf01, dycoms2_rf01_fixed_sst, "
           + "dycoms2_rf02_ds,\n\t\tdycoms2_rf02_nd, mpace_b, rico, wangara, arm_97,\n"
           + "\t\tcloud_feedback_s6, cloud_feedback_s11, cloud_feedback_s12,\n"
           + "\t\tgabls3_night, lba, and twp_ice.",
    )
    parser.add_argument(
        "--min-cases",
        action="store_true",
        default=False,
        help="Run a minimal set of cases (e.g. to\n\t\teconomize on output). "
           + "This will include only the following\n\t\tcases if they are in RUN_CASES: "
           + "arm, atex, bomex, dycoms2_rf01,\n\t\tdycoms2_rf02_ds, rico, wangara, arm_97, "
           + "gabls3_night, lba,\n\t\tand twp_ice.",
    )
    parser.add_argument(
        "--skip-default-flags",
        action="store_true",
        default=False,
        help="Skip the default flag configurations.",
    )
    parser.add_argument(
        "--central-run-script",
        action="store_true",
        default=False,
        help="Use the run_scm_all.bash script from the repository you running this script from. "
           + "Otherwise the script from the corresponding clones is taken.",
    )
    parser.add_argument(
        "-g",
        "--gg-test",
        action="store_true",
        default=False,
        help="A grid generalization test is being performed. Keep stats output files "
           + "to a minimal size",
    )
    parser.add_argument(
        "-f",
        "--flag-config-file",
        action="store",
        type=str,
        default="flag_config.json",
        help="Choose which file will be read to configure the different flag settings.\n",
    )
    args = parser.parse_args()
    if sum([args.short_cases, args.priority_cases, args.min_cases]) > 1:
        print(
            "Error: You can only enter one of --short-cases, --priority-cases or --min-cases."
        )
        exit(1)
    return args


def read_in_flag_settings(flag_config_file_path):
    print(f"Reading settings from {flag_config_file_path}...")
    if not flag_config_file_path.endswith(".json"):
        print(
            "Config file for different flag settings (--flag-config-file) "
            + "needs to be in the JSON format."
        )
        exit(1)
    with open(flag_config_file_path) as json_file:
        data = json.load(json_file)
    return data


def get_flag_file_names(flag_config_file, skip_default_flag, flags_to_change):
    """
    This function takes the flag config JSON file and returns a list of the names of
    all the flag .in files with the toggled flags from configurable_model_flags.in.
    """
    flag_file_names = {}
    default_name = "default"
    for run_name in flags_to_change:
        if run_name == default_name:
            print(
                f"Error: You are not allowed to use the name {default_name} "
                + "for one of the cases in your {flag_config_file}."
            )
            exit(1)
        flag_file_names[run_name] = run_name + "_tmp.in"
    if not skip_default_flag:
        flag_file_names[default_name] = "configurable_model_flags.in"
    return flag_file_names


def write_flag_change(flags_file, text_line, flags):
    for flag_name, flag_value in flags.items():
        pattern = r"\b" + re.escape(flag_name) + r"\b"
        flag_in_line = re.search(pattern, text_line)
        if flag_in_line:
            if isinstance(flag_value, bool):
                if flag_value:
                    text_line = re.sub(r"\..*\.", ".true.", text_line)
                else:
                    text_line = re.sub(r"\..*\.", ".false.", text_line)
            elif isinstance(flag_value, int):
                text_line = re.sub(r"= .*", f"= {flag_value}", text_line)
    flags_file.write(text_line)


def create_flag_files(abs_path_to_dirs, flags_to_change, flag_file_names):
    """
    This function creates the different input flag files from the flag configuration
    file content, given by flags_to_change.
    """
    print(
        f"Creating all input flag files in {abs_path_to_dirs}"
        + "input/tunable_parameters ..."
    )
    with open(
        abs_path_to_dirs
        + "input/tunable_parameters/configurable_model_flags.in",
        "r",
    ) as flags_file_default:

        for case_name, flag_file_name in flag_file_names.items():
            if flag_file_name != "configurable_model_flags.in":
                new_flags_file_path = (
                    abs_path_to_dirs
                    + "/input/tunable_parameters/"
                    + flag_file_name
                )

                # create and open new flag file
                with open(new_flags_file_path, "w") as flags_file:
                    for line in flags_file_default:
                        write_flag_change(flags_file, line, flags_to_change[case_name])

                # start again at beginning of file (set file pointer back to start)
                flags_file_default.seek(0)


def run_clubb_model_for_all_flag_settings(args, abs_path_to_dirs, flag_files):
    """
    This function takes all the names of the newly created flag .in files and
    runs CLUBB for all of the different flag files.
    """
    flags_to_add = []
    if args.priority_cases:
        flags_to_add.append("--priority-cases")
    elif args.min_cases:
        flags_to_add.append("--min-cases")
    elif args.short_cases:
        flags_to_add.append("--short-cases")

    # Initialize the flag_exit_code array.
    flag_exit_codes = []
    flag_exit_codes = [-999 for indx in range(len(flag_files.values()))]

    # Run every case in CLUBB for each flag file
    print(f"\nRunning CLUBB for {len(flag_files.values())} total flag files")
    indx = -1
    for _, flag_file in flag_files.items():
        print(f"\nRunning cases for {flag_file} ...")
        abs_clubb_path = f"{abs_path_to_dirs}"
        indx = indx + 1
        if args.central_run_script:
            result = subprocess.run(
                [
                    f"./run_scm_all.bash",
                    *flags_to_add,
                    "--clubb_exec_file",
                    f"{abs_clubb_path}/bin/clubb_standalone",
                    "-o",
                    f"{abs_clubb_path}/output/{flag_file.split('.')[0]}",
                    "--flags_file",
                    f"{abs_clubb_path}/input/tunable_parameters/{flag_file}",
                ],
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                universal_newlines = True
            )
            print(result.stdout)
            flag_exit_codes[indx] = result.returncode
        else:
            result = subprocess.run(
                [
                    f"{abs_clubb_path}/run_scripts/run_scm_all.bash",
                    *flags_to_add,
                    "--clubb_exec_file",
                    f"{abs_clubb_path}/bin/clubb_standalone",
                    "-o",
                    f"{abs_clubb_path}/output/{flag_file.split('.')[0]}",
                    "--flags_file",
                    f"{abs_clubb_path}/input/tunable_parameters/{flag_file}",
                ],
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                universal_newlines = True
            )
            print(result.stdout)
            flag_exit_codes[indx] = result.returncode

    exit_code = max(flag_exit_codes)
    return exit_code


def main():
    args = get_cli_args()

    path_to_runscript = os.path.dirname(os.path.realpath(__file__))
    dest_dir_path_expanded = path_to_runscript + "/.."
    abs_destination_dir_path = os.path.abspath(dest_dir_path_expanded) + "/"

    if not os.path.exists(abs_destination_dir_path):
        print(f"Directory {abs_destination_dir_path} does not exist.")
        print(f"Creating directory {abs_destination_dir_path}...")
        os.makedirs(abs_destination_dir_path)

    flags_to_change = read_in_flag_settings(args.flag_config_file)

    flag_file_names = get_flag_file_names(
        args.flag_config_file, args.skip_default_flags, flags_to_change
    )

    create_flag_files(
        abs_destination_dir_path, flags_to_change, flag_file_names
    )

    exit_code = run_clubb_model_for_all_flag_settings(
        args, abs_destination_dir_path, flag_file_names
    )

    if exit_code == 0:
        print(f"\nCLUBB run was successful for all cases and all files")
    else:
        print(f"\nThere was a failure in the CLUBB run")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
