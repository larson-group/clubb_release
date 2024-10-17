"""
Test whether two or more branches of CLUBB produce the same answer when various flags are toggled.

The simplest way to use this script would look like:
  python3 ./run_bindiff_w_flags.py -b <branch1>,<branch2>
                                   -f <flag_config_file> -d <destination directory>
This would then clone both branches into the <destination directory>.
Inputs:
  1. Two or more branches of CLUBB you want to compare
  2. JSON configuration file, describing what flags you want to toggle
  3. Path to directory where the clones should be stored

Output:
  1. Text summary of differences found in the output of the branches

Takes two or more branches or commit hashes that are on github
and clones each git repo into it's own directory.
For each of those clones all different input flag configuration files get generated,
depending on the input file passed with --flag-config-file.
This input file defines what flags should be adjusted, taking the configurable_model_flags.in
file in the corresponding clone as the default.
An example for this file is the run_bindiff_w_flags_config_example.json found
in the run_scripts directory.
Then for each of those clones, the code gets compiled and executed, with all the different
flag files.
Finally all the different outputs for the different flag files are compared for all possible
combinations of the passed branches and commit hashes.
If a directory with a clone is already existent in the directory passed with --destination-dir,
then the script will prompt the user, how it should continue.
The user can either overwrite that clone with a new clone, meaning that also the compilation and
cases will be run again, or the user can decide to skip the compilation and running of cases
for that specific clone.
For the latter case it is assumed that the output is already generated.

To change a flag setting afterwards without running everything again, you can do the following:
Suppose you ran the script like:
  python3 ./run_bindiff_w_flags -b master,ticket_42 -f flag_config.json -d ~/clubb_bindiff
1. Add the new flag configuration into flag_config.json
   in the clone you are running the script from.
2. Go into ~/clubb_bindiff/master and ~/clubb_bindiff/ticket_42
   and add the new flag configuration file with the name '<flag_case_name>_tmp.json'
   where the flag_case_name is the name of the case you entered in the flag_config.json.
3. Run CLUBB again for each branch with the new flag file and define the output directory
   in the run_scm script to be in the output directory in <flag_case_name>_tmp.
4. Run the run_bindiff_w_flags.py script again like before and enter 'no' if the script is
   prompting you with the question if this directory should be overwritten.

So if you wanted to manually change the code, you would need to recompile and run the model
again for all desired flag files.
"""

import shutil
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
        "-v",
        "--verbose",
        action="store",
        type=int,
        default=0,
        help="Choose level of verbosity for outputs, i.e. what is printed to console.\n"
           + "0: Output just a summary.\n"
           + "1: Default. Add summarized results for each file.\n"
           + "2: Add tables with detailed numerical differences in common variables for each file.",
    )
    parser.add_argument(
        "-f",
        "--flag-config-file",
        action="store",
        type=str,
        default="flag_config.json",
        help="Choose which file will be read to configure the different flag settings.\n",
    )
    parser.add_argument(
        "-d",
        "--destination-dir",
        action="store",
        type=str,
        default="~/clubb_bindiff",
        help="Choose in which directory the .\n",
    )
    parser.add_argument(
        "-b",
        "--branches",
        action="store",
        type=list_of_strings,
        help="Enter a comma seperated list of two or more branches, commit hashes or any sort "
        + "of git reference from github, that you want to compare.",
    )
    args = parser.parse_args()
    if sum([args.short_cases, args.priority_cases, args.min_cases]) > 1:
        print(
            "Error: You can only enter one of --short-cases, --priority-cases or --min-cases."
        )
        exit(1)
    if len(args.branches) < 2:
        print(
            "Error: You need to enter at least two branches or commits with -b or --branches."
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


def create_dir_for_branch_and_compile(branch, abs_path_to_dirs):
    """
    This function takes all the branches, clones them into the passed abs_path_to_dirs directory
    and compiles the model for each of the clones.
    """
    skip_compilation_and_run = False
    if not os.path.exists(abs_path_to_dirs + branch):
        os.makedirs(abs_path_to_dirs + branch)
    else:
        print(f"{abs_path_to_dirs + branch} does already exist.")
        ans = ""
        while ans != "yes" and ans != "no":
            ans = input(
                "Should the directory be overwritten (if not the already existent directory "
                + "will be taken, assuming it has the executable and all flag files) [yes/no]: "
            )
        if ans == "yes":
            shutil.rmtree(abs_path_to_dirs + branch)
            os.makedirs(abs_path_to_dirs + branch)
        elif ans == "no":
            print(f"Taking the already existent directory {branch}.")
            skip_compilation_and_run = True

    if not skip_compilation_and_run:
        print(f"Cloning CLUBB into {abs_path_to_dirs + branch}...")
        git.Git(abs_path_to_dirs + branch).clone(
            "git@github.com:larson-group/clubb.git"
        )
        git_obj = git.Git(abs_path_to_dirs + branch + "/clubb")

        git_obj.checkout(branch)

        print(f"Compiling CLUBB in {abs_path_to_dirs + branch}...")
        git_obj.execute([f"{abs_path_to_dirs + branch}/clubb/compile/compile.bash"])
    else:
        print(f"Skipping cloning and compilation for {branch}.")

    return skip_compilation_and_run


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


def create_flag_files(branch, abs_path_to_dirs, flags_to_change, flag_file_names):
    """
    This function creates the different input flag files from the flag configuration
    file content, given by flags_to_change.
    """
    print(
        f"Creating all input flag files in {abs_path_to_dirs + branch}/clubb/"
        + "input/tunable_parameters ..."
    )
    with open(
        abs_path_to_dirs
        + branch
        + "/clubb/input/tunable_parameters/configurable_model_flags.in",
        "r",
    ) as flags_file_default:

        for case_name, flag_file_name in flag_file_names.items():
            if flag_file_name != "configurable_model_flags.in":
                new_flags_file_path = (
                    abs_path_to_dirs
                    + branch
                    + "/clubb/input/tunable_parameters/"
                    + flag_file_name
                )

                # create and open new flag file
                with open(new_flags_file_path, "w") as flags_file:
                    for line in flags_file_default:
                        write_flag_change(flags_file, line, flags_to_change[case_name])

                # start again at beginning of file (set file pointer back to start)
                flags_file_default.seek(0)


def run_clubb_model_for_all_flag_settings(args, branch, abs_path_to_dirs, flag_files):
    """
    This function takes all the branches and names of the newly created flag .in files and
    runs CLUBB for each of the branch clones for all different flag files.
    """
    flags_to_add = []
    if args.priority_cases:
        flags_to_add.append("--priority-cases")
    elif args.min_cases:
        flags_to_add.append("--min-cases")
    elif args.short_cases:
        flags_to_add.append("--short-cases")

    for _, flag_file in flag_files.items():
        if args.verbose == 0:
            opt_kwargs = {"stdout": subprocess.DEVNULL}
        else:
            opt_kwargs = {}

        print(f"\nRunning cases for {flag_file} for {branch}...")
        abs_clubb_path = f"{abs_path_to_dirs}{branch}/clubb"
        if args.central_run_script:
            subprocess.call(
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
                **opt_kwargs,
            )
        else:
            subprocess.call(
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
                **opt_kwargs,
            )


def compare_outputs(branches, abs_path_to_dirs, flag_files, verbose_level):
    """
    This function takes all the branches and names of the newly created flag configuration .in files
    and compares the output for all those different flag files for all possible combinations
    of branches.
    """
    branch_combinations = {
        frozenset({branch1, branch2})
        for branch1 in branches
        for branch2 in branches
        if branch1 != branch2
    }
    for branch1, branch2 in branch_combinations:
        for _, flag_file in flag_files.items():
            print(f"\nComparing {branch1} and {branch2} for {flag_file}...")
            subprocess.call(
                [
                    "python3",
                    "./run_bindiff_all.py",
                    f"{abs_path_to_dirs}{branch1}/clubb/output/{flag_file.split('.')[0]}",
                    f"{abs_path_to_dirs}{branch2}/clubb/output/{flag_file.split('.')[0]}",
                    "--verbose",
                    f"{verbose_level}",
                ]
            )


def main():
    args = get_cli_args()

    dest_dir_path_expanded = os.path.expanduser(args.destination_dir)
    abs_destination_dir_path = os.path.abspath(dest_dir_path_expanded) + "/"

    if not os.path.exists(abs_destination_dir_path):
        print(f"Directory {abs_destination_dir_path} does not exist.")
        print(f"Creating directory {abs_destination_dir_path}...")
        os.makedirs(abs_destination_dir_path)

    flags_to_change = read_in_flag_settings(args.flag_config_file)

    flag_file_names = get_flag_file_names(
        args.flag_config_file, args.skip_default_flags, flags_to_change
    )
    for branch in args.branches:
        skip_run = create_dir_for_branch_and_compile(branch, abs_destination_dir_path)
        if not skip_run:
            create_flag_files(
                branch, abs_destination_dir_path, flags_to_change, flag_file_names
            )
            run_clubb_model_for_all_flag_settings(
                args, branch, abs_destination_dir_path, flag_file_names
            )
        else:
            print(f"Skipping run for {branch}.")

    compare_outputs(
        args.branches, abs_destination_dir_path, flag_file_names, args.verbose
    )


if __name__ == "__main__":
    main()
