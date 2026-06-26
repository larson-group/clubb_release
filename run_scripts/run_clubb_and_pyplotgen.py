#!/usr/bin/env python3
#####################################################################
#
# Script that runs one or more CLUBB cases (BOMEX, RICO, etc.), 
# outputs the results to directory output/$NAME,
# plots the results using pyplotgen, and stores the plots in  
# directory output/pyplots_$NAME.
#
#####################################################################

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


# These are all the cloud cases that will be run.
CASES = [
    "arm_97",
    "arm",
    "bomex",
    "rico",
    "cgils_s6",
    "cgils_s12",
    "dycoms2_rf01",
    "dycoms2_rf02_ds",
    "lba",
    "gabls3_night",
]

# This is the default run name. It writes to output/NAME and output/pyplots_NAME.
DEFAULT_RUN_NAME = "new"

# Plotting options for pyplotgen (see postprocessing/pyplotgen/README.md)
PYPLOTGEN_BASE_OPTIONS = ["--plot-budgets", "-l"]

# Code begins ------------------------------------------------------


def get_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run selected CLUBB cases and plot the output with pyplotgen."
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=CASES,
        help="Case names to run and plot. Defaults to the case list in this script.",
    )
    parser.add_argument(
        "--name",
        default=DEFAULT_RUN_NAME,
        help="Run name. Writes CLUBB output to output/NAME and plots to output/pyplots_NAME.",
    )
    parser.add_argument(
        "--output-root",
        help="Directory containing CLUBB output subdirectories. Default: repo output/.",
    )
    parser.add_argument(
        "--compare",
        action="append",
        default=[],
        help=(
            "Existing output subdirectory or absolute path to include in the plot. "
            "May be repeated."
        ),
    )
    parser.add_argument(
        "--skip-compile",
        action="store_true",
        help="Skip the compile step and use the existing CLUBB executable.",
    )
    parser.add_argument(
        "--max-iters",
        type=int,
        help="Forward -max_iters to run_scm.py for each case.",
    )
    parser.add_argument(
        "--enable-plot-multithreading",
        action="store_true",
        help="Allow pyplotgen to use multiprocessing. Disabled by default for predictable wrapper behavior.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would run without creating output.",
    )
    return parser.parse_args()


def format_command(command: list[str | Path]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run_checked(command: list[str | Path], *, stdout=None, stderr=None, cwd=None, env=None, dry_run=False) -> None:
    if dry_run:
        cwd_text = f" (cwd={cwd})" if cwd is not None else ""
        print(f"Dry run: {format_command(command)}{cwd_text}", flush=True)
        return
    subprocess.run(
        [str(part) for part in command],
        check=True,
        stdout=stdout,
        stderr=stderr,
        cwd=cwd,
        env=env,
    )


def resolve_compare_output_dirs(output_root: Path, compare_output_dirs: list[str]) -> list[Path]:
    resolved = []
    for compare_output_dir in compare_output_dirs:
        if not compare_output_dir:
            continue
        compare_path = Path(compare_output_dir).expanduser()
        if not compare_path.is_absolute():
            compare_path = output_root / compare_path
        compare_path = compare_path.resolve()
        if compare_path.is_dir():
            resolved.append(compare_path)
        else:
            print(f"Warning: comparison output directory does not exist and will be skipped: {compare_path}", flush=True)
    return resolved


def main() -> int:
    args = get_cli_args()

    # Figure out the directory where the script is located
    scriptPath = Path(__file__).resolve().parent
    clubb_root = scriptPath.parent

    # Store the current directory location so it can be restored
    restoreDir = Path.cwd()

    try:
        # Change directories to the one that this script is located in,
        #   assumed to be directory run_scripts
        os.chdir(scriptPath)

        cases = list(args.cases)
        output_root = (
            Path(args.output_root).expanduser().resolve()
            if args.output_root
            else (clubb_root / "output").resolve()
        )
        run_output_path = output_root / args.name
        pyplots_output_path = output_root / f"pyplots_{args.name}"

        # Clear existing output before running.
        if not args.dry_run:
            if run_output_path.is_dir():
                print(f"Removing existing directory {run_output_path}", flush=True)
                shutil.rmtree(run_output_path)
            if pyplots_output_path.is_dir():
                print(f"Removing existing directory {pyplots_output_path}", flush=True)
                shutil.rmtree(pyplots_output_path)
            run_output_path.mkdir(parents=True)

        # Compile CLUBB
        if not args.skip_compile:
            run_checked([sys.executable, clubb_root / "compile.py"], cwd=clubb_root, dry_run=args.dry_run)

        # Run all the cases specified in variable CASES.
        for case in cases:
            print(f"Running case {case} . . . ", flush=True)
            run_scm_command: list[str | Path] = [
                sys.executable,
                scriptPath / "run_scm.py",
                "-out_dir",
                run_output_path,
            ]
            if args.max_iters is not None:
                run_scm_command.extend(["-max_iters", str(args.max_iters)])
            run_scm_command.append(case)
            run_checked(
                run_scm_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=scriptPath,
                dry_run=args.dry_run,
            )

        print("\nAll cases have finished running.", flush=True)

        # Create plots using pyplotgen.
        print("\nCreating plots using pyplotgen . . .\n", flush=True)
        pyplotgen = clubb_root / "postprocessing" / "pyplotgen" / "pyplotgen.py"
        plot_inputs = [*resolve_compare_output_dirs(output_root, args.compare), run_output_path]
        pyplotgen_options = [*PYPLOTGEN_BASE_OPTIONS, "--cases", *cases, "-c"]
        if not args.enable_plot_multithreading:
            pyplotgen_options.insert(0, "--disable-multithreading")

        pyplotgen_env = os.environ.copy()
        mpl_config_dir = output_root / ".matplotlib"
        if not args.dry_run:
            mpl_config_dir.mkdir(parents=True, exist_ok=True)
        pyplotgen_env.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

        run_checked(
            [sys.executable, pyplotgen, *pyplotgen_options, *plot_inputs, "-o", pyplots_output_path],
            env=pyplotgen_env,
            dry_run=args.dry_run,
        )
    finally:
        os.chdir(restoreDir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
