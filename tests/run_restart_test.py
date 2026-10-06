#!/usr/bin/env python3
"""Check restart consistency with the native or JAX standalone runner.

Run a complete case, preserve its output in restart/, then initialize another
run from the saved interior record nearest the effective midpoint. Compare
the selected final statistic bit for bit. Forwarded timestep/duration and
physics settings apply to both runs.

This existing native workflow exercises src/clubb_driver.F90 (init_clubb_case
and restart_clubb) and src/Input_fields/input_fields.F90 (stat_fields_reader).
-jax uses the corresponding JAX ports and compares every saved column; the
native path retains its original first-column comparison. Native runs require
a compiled standalone executable; JAX uses the managed launcher environment.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

import numpy as np
from netCDF4 import Dataset
from run_scripts.run_scm import extract_jax_options


TESTS_DIR = Path(__file__).resolve().parent
CLUBB_ROOT = TESTS_DIR.parent
RUN_SCRIPTS = CLUBB_ROOT / "run_scripts"
RUN_SCM = RUN_SCRIPTS / "run_scm.py"
OUTPUT_DIR = CLUBB_ROOT / "output"
RESTART_DIR = CLUBB_ROOT / "restart"


def _read_model_times(model_file: Path) -> tuple[float, float]:
    values: dict[str, float] = {}
    line_re = re.compile(r"^\s*([a-zA-Z_]\w*)\s*=\s*([-+0-9.eEdD]+)")
    with model_file.open(encoding="utf-8") as f:
        for line in f:
            m = line_re.match(line.split("!")[0])
            if m:
                key, val = m.groups()
                try:
                    values[key.lower()] = float(val.replace("D", "E").replace("d", "e"))
                except ValueError:
                    pass
    if "time_initial" not in values or "time_final" not in values:
        raise ValueError(f"Could not parse time_initial/time_final from {model_file}")
    return values["time_initial"], values["time_final"]


def _saved_restart_time(case_name: str, time_initial: float, time_final: float) -> float:
    """Choose the saved interior record nearest the effective midpoint."""
    split_files = [OUTPUT_DIR / f"{case_name}_{grid}.nc" for grid in ("zt", "zm", "sfc")]
    files = split_files if split_files[0].is_file() else [OUTPUT_DIR / f"{case_name}_stats.nc"]
    saved_times = None
    for path in files:
        with Dataset(path) as ds:
            times = np.asarray(ds["time"][:].compressed(), dtype=float)
        saved_times = times if saved_times is None else np.intersect1d(saved_times, times)
    saved_times = saved_times[(saved_times > time_initial) & (saved_times < time_final)]
    if not saved_times.size:
        raise ValueError("Restart test requires a saved record strictly inside the run")
    halfway_time = 0.5 * (time_initial + time_final)
    return float(saved_times[np.argmin(np.abs(saved_times - halfway_time))])


def _run_scm(
    case_name: str, extra_opts: list[str], override: str | None = None
) -> tuple[int, str]:
    cmd = [str(RUN_SCM), *extra_opts]
    if override is not None:
        cmd.extend(["-override", override])
    cmd.append(case_name)
    result = subprocess.run(
        cmd,
        cwd=CLUBB_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    return result.returncode, result.stdout


def _cleanup_case_outputs(case_name: str) -> None:
    for path in glob.glob(str(OUTPUT_DIR / f"{case_name}*")):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def _compare_final_timestep(case_name: str, var_name: str, all_columns: bool = False) -> bool:
    restart_file = RESTART_DIR / f"{case_name}_stats.nc"
    output_file = OUTPUT_DIR / f"{case_name}_stats.nc"
    if not restart_file.is_file() or not output_file.is_file():
        return False

    with Dataset(restart_file) as ds_restart, Dataset(output_file) as ds_output:
        if var_name not in ds_restart.variables or var_name not in ds_output.variables:
            return False

        var_restart = np.asarray(ds_restart.variables[var_name][...])
        var_output = np.asarray(ds_output.variables[var_name][...])

        if var_restart.ndim < 1 or var_output.ndim < 1:
            return False

        # Native input_netcdf:get_var fixes the stored column index at 1.
        # The legacy test checks only that column even for multicolumn runs.
        # JAX reads each saved column separately, so compare all of them.
        if not all_columns and var_restart.ndim == 3 and var_output.ndim == 3:
            var_restart = var_restart[-1, :, 0]
            var_output = var_output[-1, :, 0]
        else:
            var_restart = var_restart[-1, ...]
            var_output = var_output[-1, ...]

        if var_restart.shape != var_output.shape:
            return False

        return np.array_equal(var_restart, var_output)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run CLUBB restart bit-for-bit test.",
        #RawDescriptionHelpFormatter so the line breaks are kept.
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Unrecognized options are forwarded to run_scm.py for both runs,\n"
            "e.g. run_restart_test.py bomex -multicol 4 (case_name must come first).\n"
            "Caveats:\n"
            "  -override    is applied before the test's restart settings.\n"
            "  -output_dir, -tout, -stats, etc. change the output location or contents,\n"
            "               so the test may not find or compare the results."
        ),
        add_help=False, allow_abbrev=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument("case_name", help="Case name (e.g., bomex, rico_silhs)")
    parser.add_argument(
        "-jax", action="store_true",
        help="Use the JAX runner for both runs; accepts attached -jax=cpu or -jax=gpu.",
    )
    parser.add_argument(
        '-var', dest='var', default="thlm", help="Variable name to compare (default: thlm)"
    )
    parser.add_argument(
        '-keep_artifacts', dest='keep_artifacts',
        action="store_true",
        help="Keep output/<case>* and restart/ files after the test finishes.",
    )
    # Unrecognized options are forwarded to run_scm.py (e.g. -multicol 4)
    argv, jax_options, jax_occurrences = extract_jax_options(sys.argv[1:])
    if jax_occurrences > 1:
        parser.error("-jax may be specified only once.")
    args, extra_opts = parser.parse_known_args(argv)
    if args.jax:
        extra_opts.insert(0, "-jax" if jax_options is None else f"-jax={jax_options}")

    model_file = CLUBB_ROOT / "input" / "case_setups" / f"{args.case_name}_model.in"
    if not model_file.is_file():
        print(f"{model_file} does not exist")
        return 1

    try:
        time_initial, time_final = _read_model_times(model_file)
    except ValueError as exc:
        print(str(exc))
        return 1
    halfway_time = 0.5 * (time_initial + time_final)

    try:
        shutil.rmtree(RESTART_DIR, ignore_errors=True)
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        _cleanup_case_outputs(args.case_name)

        print(f"Running standard {args.case_name} case... ", end="", flush=True)
        ret, out = _run_scm(args.case_name, extra_opts, override="l_restart=.false.")
        if ret != 0:
            print("FAILED")
            print(out, end="")
            return 1
        print("Done!")

        # run_scm writes the effective namelist, including forwarded timestep,
        # duration and override settings. Use its times for shortened tests too.
        generated_input = OUTPUT_DIR / f"{args.case_name}.in"
        if generated_input.is_file():
            time_initial, time_final = _read_model_times(generated_input)
        try:
            halfway_time = _saved_restart_time(args.case_name, time_initial, time_final)
        except (OSError, KeyError, ValueError) as exc:
            print(f"Could not select restart time: {exc}")
            return 1

        RESTART_DIR.mkdir(parents=True, exist_ok=True)
        run_files = list(OUTPUT_DIR.glob(f"{args.case_name}*"))
        if not run_files:
            print(f"No standard output files found for {args.case_name}")
            return 1
        for path in run_files:
            shutil.move(str(path), str(RESTART_DIR / path.name))

        override = (
            f"l_restart = .true.,"
            f"restart_path_case = 'restart/{args.case_name}',"
            f"time_restart = {halfway_time}"
        )
        print(
            f"Running restart {args.case_name} case from halfway point... ",
            end="",
            flush=True,
        )
        ret, out = _run_scm(args.case_name, extra_opts, override=override)
        if ret != 0:
            print("FAILED")
            print(out, end="")
            return 1
        print("Done!")

        if _compare_final_timestep(args.case_name, args.var, all_columns=args.jax):
            print("Results bit-for-bit.")
            return 0

        print("Results not bit-for-bit!")
        return 1
    finally:
        if not args.keep_artifacts:
            shutil.rmtree(RESTART_DIR, ignore_errors=True)
            _cleanup_case_outputs(args.case_name)


if __name__ == "__main__":
    sys.exit(main())
