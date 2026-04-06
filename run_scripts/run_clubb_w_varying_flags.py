#!/usr/bin/env python3
"""
Run CLUBB with various flag configurations.

This script:
  1. Reads a JSON configuration file describing flag sets to test.
  2. Expands the requested cases and flag sets into a list of run_scm.py tasks.
  3. Runs those tasks with bounded parallelism.
  4. Reports failures immediately and exits with an aggregated code.

Modes:
  - Single-case mode:
        ./run_clubb_flags.py rico
    → runs run_scm.py rico for each flag configuration.

  - Multi-case mode (no case name):
        ./run_clubb_flags.py --short-cases
    → runs run_scm.py once per (case, flag set) pair for the selected case list.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from run_scm_all import (
    ALL_CASES,
    MIN_CASES,
    PRIORITY_CASES,
    SHORT_CASES,
    STANDARD_CASES,
    positive_int,
)

DEFAULT_MAX_WORKERS = 8


# ------------------------------------------------------------------------------
# CLI ARGUMENT PROCESSING
# ------------------------------------------------------------------------------


def get_cli_args():
    """Parse wrapper arguments and validate single-case vs multi-case mode."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "case_name",
        type=str,
        nargs="?",
        help=(
            "Name of the CLUBB case to run (optional).\n"
            "If omitted, you may specify at most one of:\n"
            "  --all, --short-cases, --priority-cases, --min-cases."
        )
    )

    parser.add_argument(
        "--all",
        action="store_true",
        default=False,
        help="Run all cases (multi-case mode).",
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
    parser.add_argument(
        "--max-iters", "-max_iters",
        dest="max_iters",
        type=int,
        default=None,
        help="Maximum number of iterations to pass to run_scm.py.",
    )
    parser.add_argument(
        "--tout", "-tout",
        dest="tout",
        type=int,
        default=None,
        help="Stats output interval in seconds to pass to run_scm.py.",
    )
    parser.add_argument(
        "--nproc", "-nproc",
        dest="nproc",
        type=positive_int,
        default=DEFAULT_MAX_WORKERS,
        help=(
            "Maximum number of concurrent run_scm.py workers "
            f"(default: {DEFAULT_MAX_WORKERS})."
        ),
    )

    args = parser.parse_args()

    subset_flags_count = sum([
        args.all,
        args.short_cases,
        args.priority_cases,
        args.min_cases,
    ])

    if args.case_name is None:
        if subset_flags_count > 1:
            print("\nError: only one of the following may be specified:\n"
                  "  --all, --short-cases, --priority-cases, --min-cases\n")
            sys.exit(1)
    else:
        if subset_flags_count > 0:
            print("\nError: When providing a case_name, you may not also use "
                  "--all, --short-cases, --priority-cases, or --min-cases.\n")
            sys.exit(1)

    return args


# ------------------------------------------------------------------------------
# CONFIGURATION / TASK BUILDING
# ------------------------------------------------------------------------------


def read_flag_settings(path):
    """Load the JSON mapping from flag-set name to override values."""
    if not path.endswith(".json"):
        print("Error: Flag config file must be a JSON file.")
        sys.exit(1)

    print(f"Reading settings from {path}...")
    with open(path) as f:
        return json.load(f)


def get_flag_sets(skip_default, flag_dict):
    """Normalize JSON flag sets and optionally inject the unmodified default run."""
    flag_sets = {}

    if not skip_default:
        flag_sets["default"] = None

    for flag_set_name, overrides in flag_dict.items():
        if flag_set_name == "default":
            print("Error: 'default' may not be used as a flag set name.")
            sys.exit(1)
        flag_sets[flag_set_name] = overrides

    return flag_sets


def determine_run_cases(args):
    """Resolve the requested case selection into a concrete list of case names."""
    if args.case_name is not None:
        return [args.case_name]
    if args.short_cases:
        return SHORT_CASES
    if args.all:
        return ALL_CASES
    if args.priority_cases:
        return PRIORITY_CASES
    if args.min_cases:
        return MIN_CASES
    return STANDARD_CASES


def format_override_value(value):
    """Render Python values into Fortran-friendly override strings."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, (int, float)):
        return str(value)
    return str(value)


def build_override_arg(overrides):
    """Serialize one flag set into the comma-delimited -override format."""
    if not overrides:
        return None
    return ",".join(
        f"{key}={format_override_value(value)}"
        for key, value in overrides.items()
    )


def build_tasks(root, flag_sets, run_cases, args):
    """Expand flag sets and cases into independent run_scm.py task records."""
    run_scm_script = os.path.join(root, "run_scripts", "run_scm.py")
    tasks = []

    # Each task is one case under one flag set. Keeping tasks independent lets
    # the worker pool schedule the full case list (all cases and flag sets)
    for flag_set_name, overrides in flag_sets.items():
        override_arg = build_override_arg(overrides)
        out_dir = os.path.join(root, "output", flag_set_name)

        for case_name in run_cases:
            cmd = [
                run_scm_script,
                "-out_dir", out_dir,
                "-debug", "0",
            ]

            if args.tout is not None:
                cmd += ["-tout", str(args.tout)]

            if args.max_iters is not None:
                cmd += ["-max_iters", str(args.max_iters)]

            if override_arg is not None:
                cmd += ["-override", override_arg]

            cmd.append(case_name)

            tasks.append({
                "flag_set": flag_set_name,
                "case": case_name,
                "cmd": cmd,
            })

    return tasks


# ------------------------------------------------------------------------------
# TASK EXECUTION
# ------------------------------------------------------------------------------


def run_task(task):
    """Execute one run_scm.py task and capture failure output for reporting."""
    fd, tmp_out = tempfile.mkstemp(
        prefix=f"run_scm_{task['flag_set']}_{task['case']}_",
        suffix=".log",
    )
    os.close(fd)

    try:
        with open(tmp_out, "w") as log:
            result = subprocess.run(
                task["cmd"],
                stdout=log,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )

        error_output = ""
        if result.returncode != 0:
            with open(tmp_out) as f:
                error_output = f.read()

        return task["flag_set"], task["case"], result.returncode, error_output
    finally:
        Path(tmp_out).unlink(missing_ok=True)


def run_clubb_tasks(tasks, args):
    """Run the tasks (all cases and flag sets) with bounded concurrency and summarize failures."""
    if not tasks:
        print("No tasks to run.")
        return 0

    max_workers = min(args.nproc, len(tasks))
    print(f"Running {len(tasks)} run_scm.py task(s) with up to {max_workers} worker(s)...")

    exit_status = 0
    failures_by_flag_set = {}

    # Use one global pool across all (flag set, case) tasks so workers stay
    # busy even when some cases finish faster than others.
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(run_task, task): (task["flag_set"], task["case"])
            for task in tasks
        }

        # Consume completed tasks as they finish so errors are surfaced early
        # and we can build a per-flag-set failure summary at the end.
        for future in as_completed(futures):
            flag_set_name, case_name = futures[future]
            try:
                _, _, code, error_output = future.result()
            except Exception as exc:
                code = 1
                error_output = f"Error running {flag_set_name}/{case_name}: {exc}"

            if code == 0:
                print(f"COMPLETE -- {flag_set_name} / {case_name}")
            else:
                print(f"ERROR -- {flag_set_name} / {case_name}")
                if error_output:
                    print(error_output)
                failures_by_flag_set.setdefault(flag_set_name, []).append(case_name)
                exit_status = 1

    print("\n=================== Runs Complete ===================")
    if failures_by_flag_set:
        for flag_set_name in sorted(failures_by_flag_set):
            failed_cases = ", ".join(sorted(failures_by_flag_set[flag_set_name]))
            print(f"{flag_set_name} failures: {failed_cases}")
    else:
        print("All flag sets ran to completion.")

    return exit_status


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------


def main():
    """Wire together argument parsing, task expansion, execution, and exit code."""
    args = get_cli_args()

    script_dir = os.path.dirname(os.path.realpath(__file__))
    root = os.path.abspath(os.path.join(script_dir, ".."))

    # Build the full task matrix up front so the execution phase is only
    # responsible for scheduling and reporting.
    flag_dict = read_flag_settings(args.flag_config_file)
    flag_sets = get_flag_sets(args.skip_default_flags, flag_dict)
    run_cases = determine_run_cases(args)
    tasks = build_tasks(root, flag_sets, run_cases, args)

    exit_code = run_clubb_tasks(tasks, args)

    if exit_code == 0:
        print("\nAll CLUBB runs succeeded.")
    else:
        print("\nSome CLUBB runs failed. See messages above.")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
