#!/usr/bin/env python3
"""
Run a SILHS convergence check by comparing AKm and lh_AKm for a small and large
number of SILHS samples.
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
# all_stats.in includes LH microphysics fields needed by this test, including
# variables referenced by silhs_noninteractive_stats.
DEFAULT_STATS_FILE = REPO_ROOT / "input" / "stats" / "all_stats.in"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run SILHS convergence test (Python port of run_silhs_test.bash)."
    )
    parser.add_argument(
        "--case",
        default="rico",
        help="Case name to run (default: rico).",
    )
    parser.add_argument(
        "--n-small",
        type=int,
        default=8,
        help="Small number of SILHS samples (default: 8).",
    )
    parser.add_argument(
        "--n-large",
        type=int,
        default=1000,
        help="Large number of SILHS samples (default: 1000).",
    )
    parser.add_argument(
        "--keep-outputs",
        action="store_true",
        help="Keep output/small and output/large on exit.",
    )
    parser.add_argument(
        "--stats",
        default=str(DEFAULT_STATS_FILE),
        help="Stats file passed to run_scm.py (default: input/stats/all_stats.in).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show run_scm.py output while running.",
    )
    return parser.parse_args()


def remove_case_outputs(case_name: str) -> None:
    for path in OUTPUT_DIR.glob(f"{case_name}*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def copy_case_outputs(case_name: str, destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)

    copied = False
    for src in OUTPUT_DIR.glob(f"{case_name}*"):
        copied = True
        dst = destination / src.name
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)

    if not copied:
        raise RuntimeError(f"No outputs found for case '{case_name}' in {OUTPUT_DIR}.")


def run_case(case_name: str, stats_file: Path, n_samples: int, verbose: bool) -> None:
    override_parts = [
        'lh_microphys_type="non-interactive"',
        f"lh_num_samples={n_samples}",
        "lh_sequence_length=1",
        "l_silhs_KK_convergence_adj_mean=.true.",
        "l_local_kk=.false.",
    ]
    cmd = [
        sys.executable,
        str(RUN_SCM_PY),
        "-stats",
        str(stats_file),
        "-override",
        ",".join(override_parts),
        case_name,
    ]
    if verbose:
        result = subprocess.run(cmd, cwd=SCRIPT_DIR, check=False)
    else:
        result = subprocess.run(
            cmd,
            cwd=SCRIPT_DIR,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"run_scm.py failed for case '{case_name}' (exit {result.returncode}).")


def load_2d_var(dataset: netCDF4.Dataset, var_name: str) -> np.ndarray:
    if var_name not in dataset.variables:
        raise RuntimeError(f"Variable '{var_name}' not found in {dataset.filepath()}.")

    data = np.asarray(dataset.variables[var_name][:])
    data = np.squeeze(data)

    if data.ndim != 2 or data.shape[0] == 0 or data.shape[1] == 0:
        raise RuntimeError(
            f"Variable '{var_name}' must be 2D and non-empty after squeeze; got shape {data.shape}."
        )
    return data


def compute_rmse(lh_zt_file: Path) -> tuple[float, tuple[int, int]]:
    if not lh_zt_file.exists():
        raise RuntimeError(f"Required file not found: {lh_zt_file}")

    with netCDF4.Dataset(lh_zt_file, mode="r") as ds:
        akm = load_2d_var(ds, "AKm")
        lh_akm = load_2d_var(ds, "lh_AKm")

    if akm.shape != lh_akm.shape:
        raise RuntimeError(f"AKm and lh_AKm shape mismatch: {akm.shape} vs {lh_akm.shape}")

    nlev = akm.shape[1]
    mse = sum((np.mean(lh_akm[:, k]) - np.mean(akm[:, k])) ** 2 for k in range(nlev)) / nlev
    rmse = float(np.sqrt(mse))
    return rmse, akm.shape


def evaluate_convergence(rmse_small: float, rmse_large: float, n_small: int, n_large: int) -> bool:
    if rmse_small < 0 or rmse_large < 0 or rmse_small > 1e-5 or rmse_large > 1e-5:
        return False

    if rmse_small == 0.0:
        return rmse_large == 0.0

    threshold = 4.0 * np.sqrt(n_small) / np.sqrt(n_large)
    return (rmse_large / rmse_small) < threshold


def main() -> int:
    args = parse_args()
    if args.n_small <= 0 or args.n_large <= 0:
        print("n-small and n-large must be positive integers.")
        return 1

    case_model = REPO_ROOT / "input" / "case_setups" / f"{args.case}_model.in"
    if not case_model.exists():
        print(f"Case model file not found: {case_model}")
        return 1
    if not RUN_SCM_PY.exists():
        print(f"Missing script: {RUN_SCM_PY}")
        return 1
    stats_file = Path(args.stats).resolve()
    if not stats_file.exists():
        print(f"Stats file not found: {stats_file}")
        return 1

    output_small = OUTPUT_DIR / "small"
    output_large = OUTPUT_DIR / "large"

    try:
        print(f"Running {args.case} case with small number of sample points...")
        remove_case_outputs(args.case)
        run_case(args.case, stats_file, args.n_small, args.verbose)
        copy_case_outputs(args.case, output_small)
        print("Done!")

        print(f"Running {args.case} case with large number of sample points...")
        remove_case_outputs(args.case)
        run_case(args.case, stats_file, args.n_large, args.verbose)
        copy_case_outputs(args.case, output_large)
        print("Done!")

        small_file = output_small / f"{args.case}_lh_zt.nc"
        large_file = output_large / f"{args.case}_lh_zt.nc"

        rmse_small, shape_small = compute_rmse(small_file)
        rmse_large, shape_large = compute_rmse(large_file)

        if shape_small != shape_large:
            print(f"Shape mismatch: small={shape_small}, large={shape_large}")
            return 1

        print(f"N_small = {args.n_small}")
        print(f"N_large = {args.n_large}")
        print(f"RMSE_small = {rmse_small}")
        print(f"RMSE_large = {rmse_large}")

        if evaluate_convergence(rmse_small, rmse_large, args.n_small, args.n_large):
            print("SILHS converged!")
            return 0

        threshold = 4.0 * np.sqrt(args.n_small) / np.sqrt(args.n_large)
        ratio = np.inf if rmse_small == 0.0 else rmse_large / rmse_small
        print(
            "FAIL: RMSE_large/RMSE_small = "
            f"{ratio} is greater than or equal to 4*sqrt(N_SMALL)/sqrt(N_LARGE) = {threshold}"
        )
        print("SILHS failed to converge!")
        return 1

    except Exception as err:
        print(err)
        print("SILHS failed to converge!")
        return 1
    finally:
        if not args.keep_outputs:
            shutil.rmtree(output_small, ignore_errors=True)
            shutil.rmtree(output_large, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
