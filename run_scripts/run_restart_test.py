#!/usr/bin/env python3
"""
Restart consistency test:
1) Run a full case.
2) Move that run's outputs into restart/.
3) Rerun from halfway time using restart input.
4) Compare final-timestep values of a target variable bit-for-bit.
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

import numpy as np
from netCDF4 import Dataset


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
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


def _run_scm(case_name: str, override: str | None = None) -> tuple[int, str]:
    cmd = [str(RUN_SCM)]
    if not (CLUBB_ROOT / "install" / "latest" / "clubb_standalone").is_file() and (
        CLUBB_ROOT / "bin" / "clubb_standalone"
    ).is_file():
        cmd.append("-legacy")
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


def _compare_final_timestep(case_name: str, var_name: str) -> bool:
    restart_file = RESTART_DIR / f"{case_name}_zt.nc"
    output_file = OUTPUT_DIR / f"{case_name}_zt.nc"
    if not restart_file.is_file() or not output_file.is_file():
        return False

    with Dataset(restart_file) as ds_restart, Dataset(output_file) as ds_output:
        if var_name not in ds_restart.variables or var_name not in ds_output.variables:
            return False

        var_restart = np.asarray(ds_restart.variables[var_name][...])
        var_output = np.asarray(ds_output.variables[var_name][...])

        if var_restart.ndim < 1 or var_output.ndim < 1:
            return False

        var_restart = var_restart[-1, ...]
        var_output = var_output[-1, ...]

        if var_restart.shape != var_output.shape:
            return False

        return np.array_equal(var_restart, var_output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run CLUBB restart bit-for-bit test.")
    parser.add_argument("case_name", help="Case name (e.g., bomex, rico_silhs)")
    parser.add_argument(
        "-v", "--var", default="thlm", help="Variable name to compare (default: thlm)"
    )
    parser.add_argument(
        "--keep_artifacts",
        action="store_true",
        help="Keep output/<case>* and restart/ files after the test finishes.",
    )
    args = parser.parse_args()

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
        ret, out = _run_scm(args.case_name)
        if ret != 0:
            print("FAILED")
            print(out, end="")
            return 1
        print("Done!")

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
        ret, out = _run_scm(args.case_name, override=override)
        if ret != 0:
            print("FAILED")
            print(out, end="")
            return 1
        print("Done!")

        if _compare_final_timestep(args.case_name, args.var):
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
