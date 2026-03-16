#!/usr/bin/env python3

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# Directory where this script lives, assumes clubb/run_scripts, which is important
# since this is used to find CLUBB_ROOT
RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT  = os.path.join(RUN_SCRIPTS, "..") 

run_scm_script = os.path.join(RUN_SCRIPTS, "run_scm.py") 

ALL_CASES = [
    "arm", "arm_3year", "arm_97", "astex_a209", "atex", "bomex", "cgils_s6", "cgils_s6_p2k", 
    "cgils_s11", "cgils_s11_p2k", "cgils_s12", "cgils_s12_p2k", "clex9_nov02", "clex9_oct14", 
    "cloud_feedback_s6", "cloud_feedback_s6_p2k", "cloud_feedback_s11", "cloud_feedback_s11_p2k", 
    "cloud_feedback_s12", "cloud_feedback_s12_p2k", "cobra", "dycoms2_rf01", "dycoms2_rf01_fixed_sst", 
    "dycoms2_rf02_do", "dycoms2_rf02_ds", "dycoms2_rf02_morr", "dycoms2_rf02_nd", "dycoms2_rf02_so", 
    "fire", "gabls2", "gabls3", "gabls3_night", "jun25_altocu", "lba", "mc3e", "mpace_a", "mpace_b", 
    "mpace_b_silhs", "nov11_altocu", "rico", "rico_silhs", "twp_ice", "wangara"
]

STANDARD_CASES = [
    "arm", "arm_97", "astex_a209", "atex", "bomex", "cgils_s6", "cgils_s11", "cgils_s12",
    "clex9_nov02", "clex9_oct14", "dycoms2_rf01", "dycoms2_rf01_fixed_sst",
    "dycoms2_rf02_do", "dycoms2_rf02_ds", "dycoms2_rf02_nd", "dycoms2_rf02_so",
    "fire", "gabls2", "gabls3", "gabls3_night", "jun25_altocu", "lba", "mc3e", "mpace_a", 
    "mpace_b", "mpace_b_silhs", "nov11_altocu", "rico", "rico_silhs", "twp_ice", "wangara"
]

PRIORITY_CASES = [
    "arm", "atex", "bomex", "dycoms2_rf01", "dycoms2_rf01_fixed_sst", 
    "dycoms2_rf02_ds", "dycoms2_rf02_nd", "mpace_b", "rico", "wangara", "arm_97",
    "cgils_s6", "cgils_s11", "cgils_s12", "gabls3_night", "lba", "twp_ice"
]

MIN_CASES = [
    "arm", "atex", "bomex", "dycoms2_rf01", "dycoms2_rf02_ds", 
    "rico", "wangara", "arm_97", "gabls3_night", "lba", "twp_ice"
]

SHORT_CASES = [
    "gabls2", "cgils_s6", "cgils_s11", "cgils_s12", "cloud_feedback_s6", 
    "cloud_feedback_s11", "cloud_feedback_s12", "twp_ice"
]

def positive_int(value):
    """Argparse type checker for strictly positive integers."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def run_case(case, options, verbose=False):
    """Run one SCM case with run_scm.py and return (case, code, error_output)."""
    cmd = [run_scm_script] + options + [case]

    try:
        if verbose:
            print(f"\n[VERBOSE] Running {case}: {' '.join(cmd)}\n")
            result = subprocess.run(cmd)
            return case, result.returncode, ""

        fd, tmp_out = tempfile.mkstemp(prefix=f"run_scm_{case}_", suffix=".log")
        os.close(fd)
        try:
            with open(tmp_out, "w") as log:
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

            error_output = ""
            if result.returncode != 0:
                with open(tmp_out) as f:
                    error_output = f.read()

            return case, result.returncode, error_output
        finally:
            Path(tmp_out).unlink(missing_ok=True)

    except Exception as e:
        return case, 1, f"Error running {case}: {e}"


def main():
    #=================== Argument parsing ===================

    parser = argparse.ArgumentParser(
        description="Simplified CLUBB SCM runner (no nightly mode, predefined case lists)"
    )
    parser.add_argument("-all", action="store_true", help="Run ALL cases, even unmaintained ones")
    parser.add_argument("-short_cases", action="store_true", help="Run short cases only")
    parser.add_argument("-priority_cases", action="store_true", help="Run priority cases only")
    parser.add_argument("-min_cases", action="store_true", help="Run minimal case set")
    parser.add_argument(
        "-nproc",
        type=positive_int,
        default=8,
        metavar="N",
        help="Number of processes to use (default: 8)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Show output from each run_scm.py call")

    args, extra_opts = parser.parse_known_args()

    #=================== Determine case list ===================

    if args.short_cases:
        run_cases = SHORT_CASES
        print("Performing short-cases run")
    elif args.all:
        run_cases = ALL_CASES
        print("Running all cases run")
        print("WARNING: this runs unsupported cases")
    elif args.priority_cases:
        run_cases = PRIORITY_CASES
        print("Performing priority-cases run")
    elif args.min_cases:
        run_cases = MIN_CASES
        print("Performing min-cases run")
    else:
        run_cases = STANDARD_CASES
        print("Performing standard run")

    #=================== Run all cases ===================

    max_workers = min(args.nproc, len(run_cases))
    print(f"Using {max_workers} process(es)")
    print("Cases to run: " + ", ".join(run_cases))
    exit_codes = {}

    # Use a process pool so multiple independent SCM cases can run in parallel.
    # max_workers limits how many case processes are active at the same time.
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}

        # Submit every requested case to the worker pool and keep a mapping from
        # each Future back to its case name so completion messages stay readable.
        for case in run_cases:
            futures[executor.submit(run_case, case, extra_opts, args.verbose)] = case

        # Handle results as workers finish rather than in submission order. This
        # lets fast cases report completion immediately while slower cases continue.
        for future in as_completed(futures):
            case = futures[future]
            try:
                _, code, error_output = future.result()
            except Exception as exc:
                code = 1
                error_output = f"Error running {case}: {exc}"

            exit_codes[case] = code
            if code == 0:
                print(f"COMPLETE -- {case}")
            else:
                print(f"ERROR -- {case}")
                if error_output:
                    print(error_output)

    #=================== Print summary ===================

    exit_status = 0
    print("\n=================== Runs Complete ===================")
    for case in run_cases:
        code = exit_codes.get(case, 1)
        if code != 0:
            print(f"{case} failure")
            exit_status = 1

    if exit_status == 0:
        print("All cases ran to completion.")

    return exit_status


if __name__ == "__main__":
    sys.exit(main())
