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

The unmodified default flag set is included unless -skip_default_flags is
provided. Alternate flag sets are applied with run_scm.py -override; this
script no longer creates temporary configurable_model_flags.in files.

Only bindiff-specific options are consumed here. Case-selection options and
run_scm.py options are forwarded to run_clubb_w_varying_flags.py, and unknown
options are then forwarded to run_scm.py. For example:

  python3 tests/run_bindiff_w_flags.py \\
      -branches master,my_branch \\
      -flag_config_file input/flag_sets/run_bindiff_w_flags_config_example.json \\
      -output_root ~/clubb_bindiff \\
      -priority_cases -workers 4 -max_iters 3

In that example, -priority_cases and -workers are consumed by
run_clubb_w_varying_flags.py, while -max_iters is forwarded to run_scm.py.

If a clone directory already exists, the script prompts whether to overwrite
it. Answering "no" reuses the existing directory and assumes the expected
output has already been generated. Use -overwrite_existing for noninteractive
test jobs such as Jenkins.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


DEFAULT_REPO_URL = "https://github.com/larson-group/clubb.git"
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
            "That script consumes case-selection options such as -priority_cases "
            "and forwards remaining options such as -max_iters to run_scm.py."
        ),
        add_help=False, allow_abbrev=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        '-branches', dest='branches',
        required=True,
        type=list_of_strings,
        help=(
            "Comma-separated list of two or more branches, commit hashes, tags, "
            "or other git refs to compare."
        ),
    )
    parser.add_argument(
        '-flag_config_file', dest='flag_config_file',
        default=str(DEFAULT_FLAG_CONFIG_FILE),
        help="JSON flag-set config passed to run_clubb_w_varying_flags.py.",
    )
    parser.add_argument(
        '-output_root', dest='destination_dir',
        default="~/clubb_bindiff",
        help="Directory where branch clones and outputs are stored.",
        metavar='DIR',
    )
    parser.add_argument(
        '-skip_default_flags', dest='skip_default_flags',
        action="store_true",
        default=False,
        help="Do not run the unmodified default flag configuration.",
    )
    parser.add_argument(
        '-verbose', dest='verbose',
        type=int,
        default=1,
        help="Verbosity level passed to run_bindiff_all.py.",
    )
    parser.add_argument(
        '-repo_url', dest='repo_url',
        default=DEFAULT_REPO_URL,
        help=f"Git repository URL to clone. Default: {DEFAULT_REPO_URL}",
    )
    parser.add_argument(
        '-compile_arg', dest='compile_arg',
        action="append",
        default=[],
        help=(
            "Extra argument passed to compile.py. May be repeated. "
            "Use -compile_arg=-debug for values that begin with '-'."
        ),
    )
    parser.add_argument(
        '-no_compile', dest='no_compile',
        action="store_true",
        default=False,
        help="Skip compilation for newly cloned or overwritten refs.",
    )
    parser.add_argument(
        '-overwrite_existing', dest='overwrite_existing',
        action="store_true",
        default=False,
        help="Overwrite existing ref clone directories without prompting.",
    )

    args, run_args = parser.parse_known_args()

    if len(args.branches) < 2:
        parser.error("at least two branches or refs must be provided with -branches")

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
        [sys.executable, str(compile_script),
         *clone_arguments(clone_dir, args.compile_arg, ["compile.py"])],
        cwd=clone_dir,
        env=compile_environment(),
    )
    return clone_dir, False



def clone_arguments(clone_dir, arguments, scripts):
    """Translate canonical options when comparing refs with older script interfaces."""
    supported = set()
    for script in scripts:
        path = clone_dir / script
        if not path.is_file():
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument":
                supported.update(arg.value for arg in node.args
                                 if isinstance(arg, ast.Constant) and isinstance(arg.value, str))
    legacy = {
        "-workers": ("-nproc", "--nproc"),
        "-flag_config_file": ("--flag-config-file", "-f"),
        "-skip_default_flags": ("--skip-default-flags",),
        "-all": ("--all",),
        "-short_cases": ("--short-cases",),
        "-priority_cases": ("--priority-cases",),
        "-min_cases": ("--min-cases",),
        "-output_dir": ("-out_dir",),
        "-params_file": ("-params",),
        "-silhs_params_file": ("-silhs_params",),
        "-install_dir": ("-install",),
    }
    translated = []
    for token in arguments:
        option, separator, value = token.partition("=")
        if option not in supported:
            option = next((old for old in legacy.get(option, ()) if old in supported), option)
        translated.append(option + separator + value)
    return translated


def build_varying_flags_command(clone_dir, args):
    """Build the run_clubb_w_varying_flags.py command for one clone."""
    command = [
        sys.executable,
        str(clone_dir / "run_scripts" / "run_clubb_w_varying_flags.py"),
        "-flag_config_file",
        args.flag_config_file,
    ]

    if args.skip_default_flags:
        command.append("-skip_default_flags")

    command.extend(args.run_args)
    scripts = ["run_scripts/run_clubb_w_varying_flags.py", "run_scripts/run_scm.py",
               "utilities/create_case_namelist.py"]
    return command[:2] + clone_arguments(clone_dir, command[2:], scripts)


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


def compare_outputs(ref_to_clone, args):
    """Compare each pair of cloned refs across matching flag-set output directories."""
    exit_status = 0

    for ref1, ref2 in itertools.combinations(ref_to_clone, 2):
        output1 = ref_to_clone[ref1] / "output"
        output2 = ref_to_clone[ref2] / "output"
        print(f"\nComparing {ref1} and {ref2}...")
        result = subprocess.run(
            [
                sys.executable,
                str(RUN_BINDIFF_ALL),
                "-flag_sets",
                "-verbose",
                str(args.verbose),
                str(output1),
                str(output2),
            ],
            check=False,
        )
        if result.returncode != 0 and exit_status == 0:
            exit_status = result.returncode

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
