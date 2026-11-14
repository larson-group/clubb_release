#!/usr/bin/env python3
"""
Run CLUBB with various flag configurations.

This script:
  1. Resolves flag groups or lists from input/clubb_config_flags.jsonc.
  2. Expands the requested cases and flag sets into a list of run_scm.py tasks.
  3. Runs those tasks with bounded parallelism.
  4. Reports failures immediately and exits with an aggregated code.

Modes:
  - Single-case mode:
        ./run_clubb_w_varying_flags.py rico
    → runs run_scm.py rico for each flag configuration.

  - Multi-case mode (no case name):
        ./run_clubb_w_varying_flags.py -cases short_cases
    → runs run_scm.py once per (case, flag set) pair for the selected case list.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_ROOT))

from input.scripts.input_catalog import positive_int, resolve_selection

from input.scripts.flag_sets import DEFAULT_FLAG_GROUP, build_override_arg, get_flag_sets, read_flag_settings  # noqa: E402

from tuner.system_defaults import default_max_workers as default_workers
from input.scripts.create_case_namelist import parse_forwarded_args

DEFAULT_MAX_WORKERS = default_workers()


# ------------------------------------------------------------------------------
# CLI ARGUMENT PROCESSING
# ------------------------------------------------------------------------------



def get_cli_args():
    """Read sweep selection and forward model/runtime options unchanged."""
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help")
    parser.add_argument("case_name", nargs="?", help="Case, group or run list")
    parser.add_argument("-cases", nargs="+", help="Cases, groups or lists from the catalog")
    parser.add_argument("-flag_group", default=DEFAULT_FLAG_GROUP, help="Flag group, list or preset from clubb_config_flags.jsonc")
    parser.add_argument("-skip_default_flags", action="store_true", help="Omit the unmodified baseline")
    parser.add_argument("-workers", dest="nproc", type=positive_int, default=DEFAULT_MAX_WORKERS)
    parser.add_argument("-clubb_root", default=str(REPO_ROOT), help="Model checkout to run (default: this checkout)")
    parser.add_argument("-override", action="append", help="Final assignments or JSON edits")
    args, extra = parse_forwarded_args(parser)
    if args.cases and args.case_name:
        parser.error("Use either a positional selection or -cases")
    for value in args.override or []:
        extra.extend(["-override", value])
    args.run_scm_extra_args = extra
    return args


def determine_run_cases(args):
    """Expand the requested case selection before constructing the sweep matrix."""
    return resolve_selection("cases", args.cases or args.case_name or "standard_cases")


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
                "-output_dir", out_dir,
                "-debug", "0",
            ]

            if override_arg is not None:
                cmd += ["-override", override_arg]

            cmd += args.run_scm_extra_args
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

    root = str(Path(args.clubb_root).resolve())

    # Build the full task matrix up front so the execution phase is only
    # responsible for scheduling and reporting.
    print(f"Reading settings from {args.flag_group}...")
    try:
        flag_dict = read_flag_settings(args.flag_group)
        flag_sets = get_flag_sets(args.skip_default_flags, flag_dict)
    except ValueError as exc:
        print(f"Error: {exc}")
        sys.exit(1)
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
