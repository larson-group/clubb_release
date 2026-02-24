#!/usr/bin/env python3
"""
Run a CLUBB timestep convergence test.

This is a Python port of run_clubb_conv_test.bash that prefers run_scm.py
runtime options/overrides over editing case setup files.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
RUN_SCM_PY = SCRIPT_DIR / "run_scm.py"
OUTPUT_DIR = REPO_ROOT / "output"
RESTART_DIR = REPO_ROOT / "restart"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CLUBB timestep convergence test.")
    parser.add_argument(
        "-p",
        "--plot-result",
        action="store_true",
        help="Generate and save convergence plot in output/<case>_convergence.png.",
    )
    parser.add_argument(
        "--case",
        default="bomex",
        help="Case name to run (default: bomex).",
    )
    parser.add_argument(
        "--var",
        default="rcm",
        help="Variable to test from <case>_stats.nc (default: rcm).",
    )
    return parser.parse_args()


def remove_case_outputs(case_name: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUTPUT_DIR.glob(f"{case_name}*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def move_case_outputs(case_name: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    moved_any = False
    for src in OUTPUT_DIR.glob(f"{case_name}*"):
        moved_any = True
        dst = destination / src.name
        if dst.exists():
            if dst.is_dir():
                shutil.rmtree(dst)
            else:
                dst.unlink()
        shutil.move(str(src), str(dst))

    if not moved_any:
        raise RuntimeError(f"No output files found for case '{case_name}' in {OUTPUT_DIR}.")


def run_scm(case_name: str, dt: int, time_initial: float, time_final: float, l_restart: bool) -> None:
    override_parts = [
        f"time_initial={time_initial}",
        f"time_final={time_final}",
        f"stats_tsamp={dt}",
    ]

    if l_restart:
        restart_path = f'restart/{case_name}'
        override_parts.extend(
            [
                "l_restart=.true.",
                f'restart_path_case="{restart_path}"',
                "time_restart=5400.0",
            ]
        )
    else:
        override_parts.append("l_restart=.false.")

    cmd = [
        sys.executable,
        str(RUN_SCM_PY),
        "-dt_main",
        str(dt),
        "-dt_rad",
        str(dt),
        "-tout",
        str(dt),
        "-override",
        ",".join(override_parts),
        case_name,
    ]

    result = subprocess.run(
        cmd,
        cwd=SCRIPT_DIR,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"run_scm.py failed for case '{case_name}' (exit {result.returncode}).")


def extract_profile(stats_file: Path, var_name: str) -> np.ndarray:
    if not stats_file.exists():
        raise RuntimeError(f"Missing expected file: {stats_file}")

    with netCDF4.Dataset(stats_file) as ds:
        if var_name not in ds.variables:
            raise RuntimeError(f"Variable '{var_name}' not found in {stats_file}")
        data = np.asarray(ds.variables[var_name][:])

    if data.ndim < 2:
        raise RuntimeError(
            f"Variable '{var_name}' has invalid shape {data.shape} in {stats_file}"
        )

    # Assume stats dimensions are time-first and repeatedly pick first index on
    # remaining non-profile axes (e.g. first column).
    last_timestep = data[-1]
    while last_timestep.ndim > 1:
        last_timestep = last_timestep[..., 0]

    if last_timestep.ndim != 1:
        raise RuntimeError(
            f"Failed to reduce '{var_name}' from {stats_file} to a 1D profile; shape={last_timestep.shape}"
        )

    return np.asarray(last_timestep)


def compute_slope(test_vars: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    # Preserve legacy behavior from the bash script:
    # these are log(2,4,8,16,32,64,128,256) for the 8 non-reference runs.
    log_timesteps = np.log(np.array([2, 4, 8, 16, 32, 64, 128, 256], dtype=float))

    rmse_values = []
    for i in range(1, test_vars.shape[1]):
        rmse = np.sqrt(np.var(test_vars[:, 0] - test_vars[:, i]))
        rmse_values.append(np.log(rmse))

    rmse_values = np.array(rmse_values, dtype=float)
    slope, intercept = np.polyfit(log_timesteps[:4], rmse_values[:4], 1)
    _ = intercept
    return float(slope), log_timesteps, rmse_values


def maybe_plot(case_name: str, var_name: str, slope: float, log_timesteps: np.ndarray, rmse_values: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    fitted_points = [rmse_values[0]]
    slope_1_line = [rmse_values[0]]
    for i in range(1, 4):
        fitted_points.append(slope * (log_timesteps[i] - log_timesteps[0]) + rmse_values[0])
        slope_1_line.append((log_timesteps[i] - log_timesteps[0]) + rmse_values[0])

    plt.plot(log_timesteps[0:4], slope_1_line[0:4], "k--", label="slope-1 line")
    plt.plot(
        log_timesteps[0:4],
        fitted_points[0:4],
        "k:",
        label=f"{case_name} fit (m = {slope:.2f}, 4 pts.)",
    )
    plt.scatter(log_timesteps, rmse_values, label=f"{case_name} data points")
    plt.xlabel("log(dt)")
    plt.ylabel("log(RMSE)")
    plt.title(f"Convergence test for {case_name.upper()} ({var_name})")
    plt.legend()
    out_png = OUTPUT_DIR / f"{case_name}_convergence.png"
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()


def main() -> int:
    args = parse_args()
    case_name = args.case
    var_name = args.var
    timesteps = [1, 2, 4, 8, 16, 30, 60, 120, 240]

    timestep_dirs = [OUTPUT_DIR / f"timestep_{dt}" for dt in timesteps]

    if not RUN_SCM_PY.exists():
        print(f"Missing script: {RUN_SCM_PY}")
        return 1

    try:
        # Clean prior artifacts used by this test.
        shutil.rmtree(RESTART_DIR, ignore_errors=True)
        for d in timestep_dirs:
            shutil.rmtree(d, ignore_errors=True)

        print(f"Creating initial {case_name} 90-min run... ", end="", flush=True)
        remove_case_outputs(case_name)
        run_scm(case_name=case_name, dt=1, time_initial=0.0, time_final=5400.0, l_restart=False)
        move_case_outputs(case_name, RESTART_DIR)
        print("Done!")

        profiles: list[np.ndarray] = []
        for dt in timesteps:
            print(f"Restarting {case_name} with time step {dt} seconds... ", end="", flush=True)
            remove_case_outputs(case_name)
            run_scm(case_name=case_name, dt=dt, time_initial=5400.0, time_final=9000.0, l_restart=True)
            curr_out = OUTPUT_DIR / f"timestep_{dt}"
            move_case_outputs(case_name, curr_out)

            stats_file = curr_out / f"{case_name}_stats.nc"
            profile = extract_profile(stats_file, var_name)
            profiles.append(profile)
            print("Done!")

        test_vars = np.column_stack(profiles)
        slope, log_timesteps, rmse_values = compute_slope(test_vars)
        print("Generating plot...")
        print(f"slope = {slope}")

        if args.plot_result:
            maybe_plot(case_name, var_name, slope, log_timesteps, rmse_values)

        if slope > 0.5:
            print("CLUBB converged!")
            return 0

        print("CLUBB failed to converge!")
        return 1

    except Exception as err:
        print(err)
        print("CLUBB failed to converge!")
        return 1
    finally:
        # Match legacy cleanup behavior.
        shutil.rmtree(RESTART_DIR, ignore_errors=True)
        for d in timestep_dirs:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
