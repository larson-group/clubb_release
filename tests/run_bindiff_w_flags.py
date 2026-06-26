#!/usr/bin/env python3
"""
Test whether two or more CLUBB refs produce the same answer when model flags
are toggled.

This script clones each requested branch, tag, commit hash, or other git ref
into a destination directory, compiles each clone, runs the selected SCM case
set for each JSON flag configuration, and compares the resulting output trees.
Comparison is run even if some cases fail. Failed runs are represented by the
files that each run left in the output tree, and the comparison reports any
missing files or differing NetCDF values.

The wrapper uses the current Python run scripts:

  1. run_clubb_w_varying_flags.py runs each selected case for each flag set.
  2. run_bindiff_all.py performs the per-case NetCDF comparisons for matching
     flag-set output directories.

Inputs:
  1. Two or more git refs to compare.
  2. A JSON configuration file describing flag sets to test.
  3. A destination directory where cloned repos and outputs are stored.

Output:
  1. Per-flag-set bindiff output from run_bindiff_all.py.
  2. A final nonzero exit status when compared output trees differ.

The JSON flag config has the same shape used by run_clubb_w_varying_flags.py:

  {
    "flag_set_name": {
      "l_some_flag": true,
      "some_integer_option": 2
    }
  }

The unmodified default flag set is included unless --skip-default-flags is
provided. Alternate flag sets are applied with run_scm.py -override; this
script no longer creates temporary configurable_model_flags.in files.

Only bindiff-specific options are consumed here. Case-selection options and
run_scm.py options are forwarded to run_clubb_w_varying_flags.py, and unknown
options are then forwarded to run_scm.py. For example:

  python3 tests/run_bindiff_w_flags.py \\
      -b master,my_branch \\
      -f input/flag_sets/run_bindiff_w_flags_config_example.json \\
      -d ~/clubb_bindiff \\
      --priority-cases -nproc 4 -max_iters 3

In that example, --priority-cases and -nproc are consumed by
run_clubb_w_varying_flags.py, while -max_iters is forwarded to run_scm.py.

If a clone directory already exists, the script prompts whether to overwrite
it. Answering "no" reuses the existing directory and assumes the expected
output has already been generated. Use --overwrite-existing for noninteractive
test jobs such as Jenkins.
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


DEFAULT_REPO_URL = "https://github.com/larson-group/clubb.git"
STATS_FILE_SUFFIX = "_stats.nc"
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FLAG_CONFIG_FILE = REPO_ROOT / "input" / "flag_sets" / "run_bindiff_w_flags_config_core_flags.json"
RUN_BINDIFF_ALL = REPO_ROOT / "run_scripts" / "run_bindiff_all.py"


def list_of_strings(arg):
    """Parse a comma-separated list and reject empty entries."""
    values = [item.strip() for item in arg.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("must contain at least one value")
    return values


def safe_dir_name(git_ref):
    """Make a filesystem-safe directory name for branch names like feature/foo."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", git_ref).strip("_") or "ref"


def get_cli_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Run CLUBB bindiffs across multiple git refs and multiple model-flag "
            "configurations."
        ),
        epilog=(
            "Any unrecognized options are forwarded to run_clubb_w_varying_flags.py. "
            "That script consumes case-selection options such as --priority-cases "
            "and forwards remaining options such as -max_iters to run_scm.py."
        ),
    )
    parser.add_argument(
        "-b",
        "--branches",
        required=True,
        type=list_of_strings,
        help=(
            "Comma-separated list of two or more branches, commit hashes, tags, "
            "or other git refs to compare."
        ),
    )
    parser.add_argument(
        "-f",
        "--flag-config-file",
        default=str(DEFAULT_FLAG_CONFIG_FILE),
        help="JSON flag-set config passed to run_clubb_w_varying_flags.py.",
    )
    parser.add_argument(
        "-d",
        "--destination-dir",
        default="~/clubb_bindiff",
        help="Directory where branch clones and outputs are stored.",
    )
    parser.add_argument(
        "--skip-default-flags",
        action="store_true",
        default=False,
        help="Do not run the unmodified default flag configuration.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=1,
        help="Verbosity level passed to run_bindiff_all.py.",
    )
    parser.add_argument(
        "--repo-url",
        default=DEFAULT_REPO_URL,
        help=f"Git repository URL to clone. Default: {DEFAULT_REPO_URL}",
    )
    parser.add_argument(
        "--compile-arg",
        action="append",
        default=[],
        help=(
            "Extra argument passed to compile.py. May be repeated. "
            "Use --compile-arg=-debug for values that begin with '-'."
        ),
    )
    parser.add_argument(
        "--no-compile",
        action="store_true",
        default=False,
        help="Skip compilation for newly cloned or overwritten refs.",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        default=False,
        help="Overwrite existing ref clone directories without prompting.",
    )

    args, run_args = parser.parse_known_args()

    if len(args.branches) < 2:
        parser.error("at least two branches or refs must be provided with -b/--branches")

    args.run_args = run_args
    args.flag_config_file = str(Path(args.flag_config_file).expanduser().resolve())
    args.destination_dir = Path(args.destination_dir).expanduser().resolve()
    return args


def run_checked(cmd, cwd=None, stdout=None, env=None):
    """Run a subprocess and raise on failure."""
    printable_cwd = f" (cwd={cwd})" if cwd else ""
    print(f"Running: {' '.join(map(str, cmd))}{printable_cwd}")
    subprocess.run(
        cmd,
        cwd=cwd,
        stdout=stdout,
        stderr=subprocess.STDOUT,
        env=env,
        check=True,
    )


def compile_environment():
    """Return an environment suitable for compile.py in a plain shell."""
    env = os.environ.copy()
    if "FC" not in env and "LMOD_FAMILY_COMPILER" not in env and shutil.which("gfortran"):
        env["FC"] = "gfortran"
    return env


def prepare_clone(git_ref, args):
    """Clone and compile one requested ref unless an existing checkout is reused."""
    ref_dir = args.destination_dir / safe_dir_name(git_ref)
    clone_dir = ref_dir / "clubb"
    skip_run = False

    if ref_dir.exists():
        print(f"{ref_dir} already exists.")
        if args.overwrite_existing:
            print(f"Overwriting existing directory for {git_ref}.")
            shutil.rmtree(ref_dir)
        elif not clone_dir.is_dir():
            print(f"{ref_dir} does not contain a clubb checkout. Overwriting it.")
            shutil.rmtree(ref_dir)
        else:
            answer = ""
            while answer not in {"yes", "no"}:
                answer = input(
                    "Overwrite it? If not, existing output is reused and runs are skipped "
                    "[yes/no]: "
                ).strip().lower()

            if answer == "yes":
                shutil.rmtree(ref_dir)
            else:
                print(f"Reusing existing directory for {git_ref}.")
                skip_run = True

    if skip_run:
        return clone_dir, True

    ref_dir.mkdir(parents=True, exist_ok=True)

    print(f"Cloning {git_ref} into {clone_dir}...")
    run_checked(["git", "clone", args.repo_url, str(clone_dir)])
    run_checked(["git", "checkout", git_ref], cwd=clone_dir)

    if args.no_compile:
        print(f"Skipping compilation for {git_ref}.")
        return clone_dir, False

    compile_script = clone_dir / "compile.py"
    if not compile_script.is_file():
        raise FileNotFoundError(
            f"{compile_script} was not found. This wrapper expects refs with "
            "the Python compile script."
        )

    print(f"Compiling CLUBB for {git_ref}...")
    run_checked(
        [sys.executable, str(compile_script), *args.compile_arg],
        cwd=clone_dir,
        env=compile_environment(),
    )
    return clone_dir, False


def build_varying_flags_command(clone_dir, args):
    """Build the run_clubb_w_varying_flags.py command for one clone."""
    command = [
        sys.executable,
        str(clone_dir / "run_scripts" / "run_clubb_w_varying_flags.py"),
        "-f",
        args.flag_config_file,
    ]

    if args.skip_default_flags:
        command.append("--skip-default-flags")

    command.extend(args.run_args)
    return command


def run_clubb_model_for_all_flag_settings(git_ref, clone_dir, args):
    """Run all requested cases and flag sets through the Python varying-flags runner."""
    print(f"\nRunning varying-flag cases for {git_ref}...")
    stdout = subprocess.DEVNULL if args.verbose == 0 else None
    result = subprocess.run(
        build_varying_flags_command(clone_dir, args),
        stdout=stdout,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        print(
            f"Run phase for {git_ref} exited with code {result.returncode}; "
            "continuing to output comparison."
        )
    return result.returncode


def get_flag_dirs(root):
    """Return child directories keyed by flag-set name."""
    return {path.name: path for path in root.iterdir() if path.is_dir()}


def get_nc_files(root):
    """Return relative NetCDF paths under a flag-set output directory."""
    return {str(path.relative_to(root)) for path in root.rglob("*.nc")}


def case_from_nc_file(path):
    """Infer case name from a *_stats.nc relative path."""
    filename = Path(path).name
    if filename.endswith(STATS_FILE_SUFFIX):
        return filename[: -len(STATS_FILE_SUFFIX)]
    return None


def get_nc_cases(root):
    """Return case names with stats NetCDF files under a flag-set output directory."""
    return {
        case
        for case in (case_from_nc_file(path) for path in get_nc_files(root))
        if case
    }


def print_name_diff(kind, names):
    """Print a sorted list of names for missing flag sets or files."""
    print(f"{kind}:")
    for name in sorted(names):
        print(f"  {name}")


def add_summary(summary, flag_set, case_name, reason):
    """Record why one flag-set/case pair differed."""
    summary.setdefault(flag_set, {}).setdefault(case_name, set()).add(reason)


def summarize_differences(summary):
    """Print the compact flag-set/case difference summary."""
    if not summary:
        print("\nDIFFERENCE SUMMARY: no flag sets had differing cases.")
        return

    print("\nDIFFERENCE SUMMARY:")
    for flag_set in sorted(summary):
        print(f"  {flag_set}:")
        for case_name in sorted(summary[flag_set]):
            reasons = ", ".join(sorted(summary[flag_set][case_name]))
            print(f"    {case_name}: {reasons}")


def compare_varying_flags_outputs(dir1, dir2, verbose):
    """Compare two run_clubb_w_varying_flags.py output directories."""
    dir1 = Path(dir1).resolve()
    dir2 = Path(dir2).resolve()

    if not dir1.is_dir() or not dir2.is_dir():
        print("Both inputs must be existing directories.")
        return 2

    if dir1 == dir2:
        print("Input paths resolve to the same directory.")
        return 2

    flag_dirs_1 = get_flag_dirs(dir1)
    flag_dirs_2 = get_flag_dirs(dir2)

    if not flag_dirs_1 or not flag_dirs_2:
        print("Expected both inputs to contain per-flag output subdirectories.")
        return 2

    only_in_1 = set(flag_dirs_1) - set(flag_dirs_2)
    only_in_2 = set(flag_dirs_2) - set(flag_dirs_1)
    common = sorted(set(flag_dirs_1) & set(flag_dirs_2))

    difference_summary = {}
    had_difference = False

    if only_in_1:
        print_name_diff(f"Flag directories only in {dir1}", only_in_1)
        had_difference = True

    if only_in_2:
        print_name_diff(f"Flag directories only in {dir2}", only_in_2)
        had_difference = True

    for flag_name in common:
        flag_dir_1 = flag_dirs_1[flag_name]
        flag_dir_2 = flag_dirs_2[flag_name]
        nc_files_1 = get_nc_files(flag_dir_1)
        nc_files_2 = get_nc_files(flag_dir_2)

        missing_in_2 = nc_files_1 - nc_files_2
        missing_in_1 = nc_files_2 - nc_files_1

        if missing_in_2:
            print_name_diff(
                f"NetCDF files in {flag_dir_1} but not {flag_dir_2}",
                missing_in_2,
            )
            had_difference = True
            for nc_file in missing_in_2:
                case_name = case_from_nc_file(nc_file)
                if case_name:
                    add_summary(difference_summary, flag_name, case_name, f"NetCDF missing from {flag_dir_2}")

        if missing_in_1:
            print_name_diff(
                f"NetCDF files in {flag_dir_2} but not {flag_dir_1}",
                missing_in_1,
            )
            had_difference = True
            for nc_file in missing_in_1:
                case_name = case_from_nc_file(nc_file)
                if case_name:
                    add_summary(difference_summary, flag_name, case_name, f"NetCDF missing from {flag_dir_1}")

        common_cases = sorted(get_nc_cases(flag_dir_1) & get_nc_cases(flag_dir_2))
        if not common_cases:
            print(f"\nNo common successful cases to compare for flag set: {flag_name}")
            continue

        print(f"\nComparing flag set: {flag_name}")
        result = subprocess.run(
            [
                sys.executable,
                str(RUN_BINDIFF_ALL),
                "-v",
                str(verbose),
                str(flag_dir_1),
                str(flag_dir_2),
            ],
            check=False,
        )

        if result.returncode != 0:
            had_difference = True
            for case_name in common_cases:
                case_result = subprocess.run(
                    [
                        sys.executable,
                        str(RUN_BINDIFF_ALL),
                        "-v",
                        "0",
                        "-case",
                        case_name,
                        str(flag_dir_1),
                        str(flag_dir_2),
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                if case_result.returncode != 0:
                    add_summary(difference_summary, flag_name, case_name, "NetCDF values differ")

    summarize_differences(difference_summary)

    if had_difference:
        print("\nSUMMARY: differences were detected across varying-flags output.")
        return 1

    print("\nSUMMARY: no differences detected across varying-flags output.")
    return 0


def compare_outputs(ref_to_clone, args):
    """Compare each pair of cloned refs across matching flag-set output directories."""
    exit_status = 0

    for ref1, ref2 in itertools.combinations(ref_to_clone, 2):
        output1 = ref_to_clone[ref1] / "output"
        output2 = ref_to_clone[ref2] / "output"
        print(f"\nComparing {ref1} and {ref2}...")
        if compare_varying_flags_outputs(output1, output2, args.verbose) != 0:
            exit_status = 1

    return exit_status


def main():
    args = get_cli_args()

    if not Path(args.flag_config_file).is_file():
        print(f"Flag config file does not exist: {args.flag_config_file}")
        return 2

    if not RUN_BINDIFF_ALL.is_file():
        print(f"Bindiff script does not exist: {RUN_BINDIFF_ALL}")
        return 2

    args.destination_dir.mkdir(parents=True, exist_ok=True)

    ref_to_clone = {}
    for git_ref in args.branches:
        try:
            clone_dir, skip_run = prepare_clone(git_ref, args)
            ref_to_clone[git_ref] = clone_dir
            if skip_run:
                print(f"Skipping run for {git_ref}.")
            else:
                run_clubb_model_for_all_flag_settings(git_ref, clone_dir, args)
        except subprocess.CalledProcessError as exc:
            print(f"Command failed for {git_ref} with exit code {exc.returncode}.")
            return exc.returncode
        except Exception as exc:
            print(f"Error preparing {git_ref}: {exc}")
            return 1

    return compare_outputs(ref_to_clone, args)


if __name__ == "__main__":
    sys.exit(main())
