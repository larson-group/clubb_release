#!/usr/bin/env python3
"""
Run CLUBB with various flag configurations.

This script:
  1. Reads a JSON configuration file describing flag sets to test.
  2. Creates modified versions of configurable_model_flags.in.
  3. Runs CLUBB for each flag configuration.
  4. Reports failures immediately and exits with an aggregated code.

Modes:
  - Single-case mode:
        ./run_clubb_flags.py rico
    → runs run_scm.py rico for each flag configuration.

  - Multi-case mode (no case name):
        ./run_clubb_flags.py --short-cases
    → runs run_scm_all.py -short_cases for each flag configuration.
"""

import sys
import os
import re
import subprocess
import argparse
import json


# ------------------------------------------------------------------------------
# CLI ARGUMENT PROCESSING
# ------------------------------------------------------------------------------

def get_cli_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )

    # case_name is now OPTIONAL
    parser.add_argument(
        "case_name",
        type=str,
        nargs="?",
        help=(
            "Name of the CLUBB case to run (optional).\n"
            "If omitted, you must specify exactly one of:\n"
            "  --all, --short-cases, --priority-cases, --min-cases."
        )
    )

    parser.add_argument(
        "--all",
        action="store_true",
        default=False,
        help="Run all cases (multi-case mode, via run_scm_all.py).",
    )
    parser.add_argument(
        "--short-cases", action="store_true", default=False,
        help="Run short cases only (multi-case mode)."
    )
    parser.add_argument(
        "--priority-cases", action="store_true", default=False,
        help="Run priority cases only (multi-case mode)."
    )
    parser.add_argument(
        "--min-cases", action="store_true", default=False,
        help="Run minimal set of cases (multi-case mode)."
    )
    parser.add_argument(
        "--skip-default-flags", action="store_true", default=False,
        help="Do not run the default flag configuration."
    )
    parser.add_argument(
        "-f", "--flag-config-file", type=str, default="flag_config.json",
        help="JSON file describing alternate flag settings."
    )

    args = parser.parse_args()

    # Count subset flags
    subset_flags_count = sum([
        args.all,
        args.short_cases,
        args.priority_cases,
        args.min_cases,
    ])

    if args.case_name is None:
        # Multi-case mode: must have only one subset flag
        if subset_flags_count >= 1:
            print("\nError: only one of the following may be specified:\n"
                  "  --all, --short-cases, --priority-cases, --min-cases\n")
            sys.exit(1)
    else:
        # Single-case mode: subset flags are not allowed
        if subset_flags_count > 0:
            print("\nError: When providing a case_name, you may not also use "
                  "--all, --short-cases, --priority-cases, or --min-cases.\n")
            sys.exit(1)

    return args


# ------------------------------------------------------------------------------
# READ JSON FLAG SETTINGS
# ------------------------------------------------------------------------------

def read_flag_settings(path):
    if not path.endswith(".json"):
        print("Error: Flag config file must be a JSON file.")
        sys.exit(1)

    print(f"Reading settings from {path}...")
    with open(path) as f:
        return json.load(f)


# ------------------------------------------------------------------------------
# DETERMINE FLAG FILE NAMES
# ------------------------------------------------------------------------------

def get_flag_file_names(skip_default, flag_dict):
    names = {}
    default_name = "default"

    for case_name in flag_dict:
        if case_name == default_name:
            print(f"Error: '{default_name}' may not be used as a flag set name.")
            sys.exit(1)
        names[case_name] = f"{case_name}_tmp.in"

    if not skip_default:
        names[default_name] = "configurable_model_flags.in"

    return names


# ------------------------------------------------------------------------------
# FLAG MODIFYING HELPER
# ------------------------------------------------------------------------------

def write_flag_change(output_file, line, case_flags):
    """
    Replace the value of a flag on a given line.
    Boolean replacement is safe; integer replacement replaces RHS of '='.
    """
    for flag_name, flag_value in case_flags.items():
        if re.search(rf"\b{re.escape(flag_name)}\b", line):
            if isinstance(flag_value, bool):
                value = ".true." if flag_value else ".false."
                line = re.sub(r"=\s*\.\w+\.", f"= {value}", line)
            elif isinstance(flag_value, int):
                line = re.sub(r"=\s*.*", f"= {flag_value}", line)

    output_file.write(line)


# ------------------------------------------------------------------------------
# CREATE FLAG FILES
# ------------------------------------------------------------------------------

def create_flag_files(root, flag_dict, flag_file_names):
    print("Creating input flag files...")

    base_path = os.path.join(root, "input/tunable_parameters")
    default_file = os.path.join(base_path, "configurable_model_flags.in")

    with open(default_file, "r") as original:
        for case_name, filename in flag_file_names.items():
            if filename == "configurable_model_flags.in":
                continue

            out_path = os.path.join(base_path, filename)
            with open(out_path, "w") as newf:
                for line in original:
                    write_flag_change(newf, line, flag_dict[case_name])
                original.seek(0)   # rewind for next flag file


# ------------------------------------------------------------------------------
# RUN CLUBB FOR EACH FLAG CONFIGURATION
# ------------------------------------------------------------------------------

def run_clubb(root, flag_files, args):
    """
    If args.case_name is set:
        Runs run_scm.py <flags> <case_name> for each flag configuration.

    If args.case_name is None:
        Runs run_scm_all.py <subset-flag> <flags> for each flag configuration,
        where <subset-flag> is one of: -all, -short_cases, -priority_cases, -min_cases.
    """
    is_single_case_mode = args.case_name is not None

    if is_single_case_mode:
        script_path = os.path.join(root, "run_scripts", "run_scm.py")
        case_name = args.case_name
        print(f"\nRunning CLUBB case '{case_name}' for {len(flag_files)} flag configurations...")
        subset_flag = None  # not used in this mode
    else:
        script_path = os.path.join(root, "run_scripts", "run_scm_all.py")
        print(f"\nRunning CLUBB multi-case mode for {len(flag_files)} flag configurations...")

        # Determine which subset flag to pass to run_scm_all.py
        if args.all:
            subset_flag = "-all"
        elif args.short_cases:
            subset_flag = "-short_cases"
        elif args.priority_cases:
            subset_flag = "-priority_cases"
        elif args.min_cases:
            subset_flag = "-min_cases"
        else:
            # use whatever default run_scm_all.py uses
            subset_flag = None

    exit_codes = []

    for case_label, flag_file in flag_files.items():
        print(f"\n--- Running flag set: {case_label} ({flag_file}): ", end="")

        out_dir = os.path.join(root, "output", case_label)
        flag_path = os.path.join(root, "input/tunable_parameters", flag_file)

        cmd = [script_path]

        # In multi-case mode, pass the subset flag (-all, -short_cases, etc.)
        if not is_single_case_mode and subset_flag is not None:
            cmd.append(subset_flag)

        cmd += [
            "-out_dir", out_dir,
            "-tout", "0",
            "-debug", "0",
            "-flags", flag_path,
        ]

        # In single-case mode, append the case name at the end
        if is_single_case_mode:
            cmd.append(case_name)

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )

        exit_codes.append(result.returncode)

        if result.returncode != 0:
            print(result.stdout)   # Print ONLY if failed
            print(f"\n!===== FAILURE for flag set '{case_label}' (file: {flag_file})")
        else:
            print("PASS")

    return max(exit_codes)


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

def main():
    args = get_cli_args()

    # Paths relative to script location
    script_dir = os.path.dirname(os.path.realpath(__file__))
    root = os.path.abspath(os.path.join(script_dir, ".."))

    flags = read_flag_settings(args.flag_config_file)
    flag_files = get_flag_file_names(args.skip_default_flags, flags)

    create_flag_files(root, flags, flag_files)

    exit_code = run_clubb(root, flag_files, args)

    if exit_code == 0:
        print("\n🎉 All CLUBB runs succeeded.")
    else:
        print("\n⚠️ Some CLUBB runs failed. See messages above.")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
