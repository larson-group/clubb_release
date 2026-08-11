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


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
RUN_SCRIPTS = REPO_ROOT / "run_scripts"
RUN_SCM_PY = RUN_SCRIPTS / "run_scm.py"
OUTPUT_DIR = REPO_ROOT / "output"
RESTART_DIR = REPO_ROOT / "restart"


def parse_args() -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(description="Run CLUBB timestep convergence test.")
    parser.add_argument(
        "-p",
        "--plot-result",
        action="store_true",
        help="Generate convergence and final-profile comparison plots in output/.",
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
    args, run_scm_args = parser.parse_known_args()
    return args, run_scm_args


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


def run_scm(
    case_name: str,
    dt: int,
    time_initial: float,
    time_final: float,
    l_restart: bool,
    run_scm_args: list[str],
) -> None:
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
        *run_scm_args,
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
        cwd=RUN_SCRIPTS,
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


def fit_convergence_line(x_values: np.ndarray, y_values: np.ndarray) -> tuple[float, float, float]:
    """Fit y = slope*x + intercept and return its coefficient of determination."""
    slope, intercept = np.polyfit(x_values, y_values, 1)
    fitted_values = slope * x_values + intercept
    residual_sum_squares = np.sum((y_values - fitted_values) ** 2)
    total_sum_squares = np.sum((y_values - np.mean(y_values)) ** 2)
    r_squared = 1.0 if total_sum_squares == 0.0 else 1.0 - residual_sum_squares / total_sum_squares
    return float(slope), float(intercept), float(r_squared)


def compute_fits(
    test_vars: np.ndarray,
    timesteps: list[int],
) -> tuple[float, float, float, float, float, float, np.ndarray, np.ndarray]:
    rmse_values = []
    for i in range(1, test_vars.shape[1]):
        rmse = np.sqrt(np.var(test_vars[:, 0] - test_vars[:, i]))
        rmse_values.append(np.log(rmse))

    rmse_values = np.array(rmse_values, dtype=float)
    non_reference_timesteps = np.asarray(timesteps[1:], dtype=float)
    if non_reference_timesteps.size != rmse_values.size:
        raise RuntimeError(
            "The number of non-reference timesteps must match the number of convergence errors"
        )
    log_timesteps = np.log(non_reference_timesteps)
    first_four_slope, first_four_intercept, first_four_r_squared = fit_convergence_line(
        log_timesteps[:4], rmse_values[:4]
    )
    all_points_slope, all_points_intercept, all_points_r_squared = fit_convergence_line(
        log_timesteps, rmse_values
    )
    return (
        first_four_slope,
        first_four_intercept,
        first_four_r_squared,
        all_points_slope,
        all_points_intercept,
        all_points_r_squared,
        log_timesteps,
        rmse_values,
    )


def maybe_plot(
    case_name: str,
    var_name: str,
    first_four_slope: float,
    first_four_intercept: float,
    first_four_r_squared: float,
    all_points_slope: float,
    all_points_intercept: float,
    all_points_r_squared: float,
    log_timesteps: np.ndarray,
    rmse_values: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt

    fit_x = log_timesteps[:4]
    first_four_fit = first_four_slope * fit_x + first_four_intercept
    all_points_fit = all_points_slope * log_timesteps + all_points_intercept
    slope_1_line = (fit_x - fit_x[0]) + rmse_values[0]

    plt.plot(fit_x, slope_1_line, "k--", label="slope-1 line")
    plt.plot(
        fit_x,
        first_four_fit,
        "k:",
        label=f"first 4: m = {first_four_slope:.2f}, R² = {first_four_r_squared:.3f}",
    )
    plt.plot(
        log_timesteps,
        all_points_fit,
        "C1-.",
        label=f"all 8: m = {all_points_slope:.2f}, R² = {all_points_r_squared:.3f}",
    )
    plt.scatter(log_timesteps, rmse_values, label=f"{case_name} data points")
    plt.xlabel("log(dt)")
    plt.ylabel("log(RMSE)")
    plt.title(f"Convergence test for {case_name.upper()} ({var_name})")
    plt.legend()
    out_png = OUTPUT_DIR / f"{case_name}_convergence.png"
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()
    print(f"Generated {out_png}")


def extract_vertical_coordinate(stats_file: Path, profile_size: int) -> np.ndarray:
    """Return the vertical coordinate that matches a final profile's length."""
    with netCDF4.Dataset(stats_file) as ds:
        for coordinate_name in ("zt", "zm", "lh_zt", "rad_zt", "rad_zm"):
            if coordinate_name not in ds.variables:
                continue
            coordinate = np.asarray(ds.variables[coordinate_name][:])
            if coordinate.ndim == 1 and coordinate.size == profile_size:
                return coordinate

    raise RuntimeError(
        f"No vertical coordinate in {stats_file} matches profile length {profile_size}"
    )


def plot_profile_comparisons(
    case_name: str,
    var_name: str,
    timesteps: list[int],
    profiles: list[np.ndarray],
    heights: np.ndarray,
) -> None:
    """Save final vertical profiles for every timestep in the convergence run."""
    import matplotlib.pyplot as plt

    if heights.ndim != 1:
        raise RuntimeError(f"Expected a 1D vertical coordinate; got shape {heights.shape}")

    fig, ax = plt.subplots()
    for dt, profile in zip(timesteps, profiles, strict=True):
        if profile.shape != heights.shape:
            raise RuntimeError(
                f"Profile shape {profile.shape} does not match vertical coordinate shape {heights.shape}"
            )
        ax.plot(profile, heights, label=f"dt = {dt} s")

    ax.set_xlabel(var_name)
    ax.set_ylabel("height (m)")
    ax.set_title(f"Final {var_name} profiles for {case_name.upper()}")
    ax.legend()
    out_png = OUTPUT_DIR / f"{case_name}_{var_name}_timestep_profiles.png"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Generated {out_png}")


def main() -> int:
    args, run_scm_args = parse_args()
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
        run_scm(
            case_name=case_name,
            dt=1,
            time_initial=0.0,
            time_final=5400.0,
            l_restart=False,
            run_scm_args=run_scm_args,
        )
        move_case_outputs(case_name, RESTART_DIR)
        print("Done!")

        profiles: list[np.ndarray] = []
        heights: np.ndarray | None = None
        for dt in timesteps:
            print(f"Restarting {case_name} with time step {dt} seconds... ", end="", flush=True)
            remove_case_outputs(case_name)
            run_scm(
                case_name=case_name,
                dt=dt,
                time_initial=5400.0,
                time_final=9000.0,
                l_restart=True,
                run_scm_args=run_scm_args,
            )
            curr_out = OUTPUT_DIR / f"timestep_{dt}"
            move_case_outputs(case_name, curr_out)

            stats_file = curr_out / f"{case_name}_stats.nc"
            profile = extract_profile(stats_file, var_name)
            profiles.append(profile)
            if heights is None:
                heights = extract_vertical_coordinate(stats_file, profile.size)
            print("Done!")

        test_vars = np.column_stack(profiles)
        (
            first_four_slope,
            first_four_intercept,
            first_four_r_squared,
            all_points_slope,
            all_points_intercept,
            all_points_r_squared,
            log_timesteps,
            rmse_values,
        ) = compute_fits(test_vars, timesteps)
        print("Generating plot...")
        print(f"first 4 points: slope = {first_four_slope}, R^2 = {first_four_r_squared}")
        print(f"all 8 points: slope = {all_points_slope}, R^2 = {all_points_r_squared}")

        if args.plot_result:
            maybe_plot(
                case_name,
                var_name,
                first_four_slope,
                first_four_intercept,
                first_four_r_squared,
                all_points_slope,
                all_points_intercept,
                all_points_r_squared,
                log_timesteps,
                rmse_values,
            )
            if heights is None:
                raise RuntimeError("No profiles were generated for comparison plotting")
            plot_profile_comparisons(case_name, var_name, timesteps, profiles, heights)

        if first_four_slope > 0.5:
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
