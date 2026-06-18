#!/usr/bin/env python3
"""
Run consistency checks for CLUBB unified stats output.

High-level picture
------------------
This script runs the same standalone CLUBB case several ways, then compares the
generated unified stats NetCDF files. The goal is to make sure stats output is
only affected by the physical/model configuration, not by bookkeeping choices
inside the stats system.

The workflow has four checks:

1. Batch-size consistency:
   Run the same multicol case with different runtime batch sizes. The NetCDF
   payloads should be exactly identical because batching is only how the driver
   chunks columns internally.

2. Output-interval consistency:
   Run one fine-output reference, then one or more coarser stats_tout runs. For
   each coarse output record, manually average the matching fine records using
   time_bnds and compare that reconstruction to the coarse NetCDF value.

3. Stats-window consistency:
   Run with stats_tstart/stats_tend set to a smaller time window. The windowed
   fine-output file should exactly equal the subset of records from the full
   fine-output reference over that same time_bnds interval.

4. Stats-window plus output-interval consistency:
   Run the same stats_tstart/stats_tend window with a coarse stats_tout, then
   manually reconstruct that windowed/coarse output from the full fine-output
   reference.

The default configuration uses BOMEX, standard_stats.in, and four deterministic
multicol columns. It does not shorten the case with max_iters; by default it
uses the case's own time_initial/time_final from the model file. The windowing
checks derive a middle-of-run window from that default case time range unless
-window_start/-window_end are provided explicitly.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass

import netCDF4
import numpy as np


RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.abspath(os.path.join(RUN_SCRIPTS, ".."))
RUN_SCM = os.path.join(RUN_SCRIPTS, "run_scm.py")

# Default case setup:
#
# - bomex is the default because it is small, common, deterministic, and has a
#   60 s dt_main in the checked-in model file.
# - standard_stats.in keeps the default stats registry compact while still
#   exercising both zt and zm profile output. Use -stats to broaden coverage
#   with standard_stats.in or another registry.
# - C8/0.2:0.8/4 creates four deterministic parameter columns, which is enough
#   to exercise stats batching without making the test large.
# - DEFAULT_BATCH_SIZES is ordered: the first value is the reference batch size,
#   and later values are exact-compared against that reference.
# - DEFAULT_FINE_TOUT should normally match the case's stats_tsamp/dt_main so
#   the reference file has one output record per sample.
# - DEFAULT_COARSE_TOUTS are manually reconstructed from DEFAULT_FINE_TOUT
#   records. The first coarse value is also used for the window+average test and
#   for deriving the automatic stats window alignment.
DEFAULT_CASE = "bomex"
DEFAULT_OUT_ROOT = os.path.join(CLUBB_ROOT, "output", "stats_output_consistency")
DEFAULT_STATS_FILE = os.path.join(CLUBB_ROOT, "input", "stats", "standard_stats.in")
DEFAULT_MULTICOL = "C8/0.2:0.8/4"
DEFAULT_BATCH_SIZES = "4,2,1"
DEFAULT_FINE_TOUT = 60
DEFAULT_COARSE_TOUTS = "300,600"

# Tolerances:
#
# Batch-size comparisons and fine-output window subset comparisons use exact
# array equality. Those paths should be bit-for-bit because they compare values
# that CLUBB wrote through the same stats_tout cadence.
#
# Manual averaging comparisons use np.allclose(actual, reconstructed,
# atol=ABS_TOL, rtol=REL_TOL). A point fails only when:
#
#   abs(actual - reconstructed) > ABS_TOL + REL_TOL * abs(reconstructed)
#
# This is intentionally tight. It allows tiny roundoff differences from
# reconstructing a Fortran time average in Python while still catching real
# output-window or averaging regressions. The printed max_abs/max_rel values are
# diagnostics over all checked points; they are not separate failure thresholds.
# In particular, max_rel can look large for near-zero reconstructed values even
# when the absolute difference is harmless.
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-12


@dataclass(frozen=True)
class RunSpec:
    """Describe one standalone run used by the consistency check."""

    label: str
    batch_size: int
    tout: int
    output_dir: str
    stats_window: tuple[int, int] | None = None


def print_banner(title: str) -> None:
    """Print a simple section banner."""
    print(f"\n{'=' * 24} {title} {'=' * 24}", flush=True)


def print_status(label: str, passed: bool, detail: str) -> None:
    """Print a uniform PASS/FAIL line."""
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] {label}: {detail}", flush=True)


def parse_csv_ints(value: str, option_name: str) -> list[int]:
    """Parse a comma-separated positive-integer option."""
    values: list[int] = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            parsed = int(item)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"{option_name} must contain only integers"
            ) from exc
        if parsed <= 0:
            raise argparse.ArgumentTypeError(
                f"{option_name} values must be positive integers"
            )
        values.append(parsed)

    if not values:
        raise argparse.ArgumentTypeError(f"{option_name} must include at least one value")
    return values


def abs_path(path: str | None) -> str | None:
    """Return an absolute path while preserving missing optional arguments."""
    return os.path.abspath(path) if path else None


def model_file_path(case_name: str) -> str:
    """Return the default model file path for one SCM case."""
    return os.path.join(CLUBB_ROOT, "input", "case_setups", f"{case_name}_model.in")


def read_model_timing(case_name: str) -> dict[str, int]:
    """Read integer-valued time settings from a case model namelist."""
    path = model_file_path(case_name)
    if not os.path.isfile(path):
        raise RuntimeError(f"Case model file not found: {path}")

    content = open(path, encoding="utf-8").read()
    values: dict[str, int] = {}
    for key in ("time_initial", "time_final", "dt_main"):
        match = re.search(rf"(?im)^\s*{key}\s*=\s*([-+]?\d+(?:\.\d*)?(?:[eEdD][-+]?\d+)?)", content)
        if not match:
            raise RuntimeError(f"Could not find {key} in {path}")

        raw_value = float(match.group(1).replace("D", "E").replace("d", "e"))
        int_value = int(round(raw_value))
        if not math.isclose(raw_value, int_value, abs_tol=1.0e-8, rel_tol=0.0):
            raise RuntimeError(
                f"{key}={raw_value} in {path} is not integer-valued. "
                "Pass an explicit aligned stats window or use integer case timing."
            )
        values[key] = int_value

    return values


def lcm_many(values: list[int]) -> int:
    """Return the least common multiple of one or more positive integers."""
    lcm_value = 1
    for value in values:
        if value <= 0:
            raise RuntimeError(f"Cannot compute LCM with non-positive value {value}")
        lcm_value = math.lcm(lcm_value, value)
    return lcm_value


def align_up(value: int, origin: int, interval: int) -> int:
    """Round value up to the next interval boundary anchored at origin."""
    offset = value - origin
    return origin + int(math.ceil(offset / interval)) * interval


def align_down(value: int, origin: int, interval: int) -> int:
    """Round value down to the previous interval boundary anchored at origin."""
    offset = value - origin
    return origin + int(math.floor(offset / interval)) * interval


def derive_stats_window(args: argparse.Namespace) -> tuple[int, int]:
    """Pick or validate the stats window used by the windowing tests."""
    timing = read_model_timing(args.case_name)
    time_initial = timing["time_initial"]
    time_final = timing["time_final"]
    dt_main = args.dt_main if args.dt_main is not None else timing["dt_main"]
    window_average_tout = args.coarse_touts[0]
    alignment = lcm_many([dt_main, args.fine_tout, window_average_tout])

    if args.window_start is not None or args.window_end is not None:
        if args.window_start is None or args.window_end is None:
            raise RuntimeError("-window_start and -window_end must be provided together")
        window_start = args.window_start
        window_end = args.window_end
    else:
        duration = time_final - time_initial
        if duration <= 0:
            raise RuntimeError(
                f"Invalid default case duration: time_initial={time_initial}, "
                f"time_final={time_final}"
            )
        window_start = align_up(time_initial + duration // 4, time_initial, alignment)
        window_end = align_down(time_initial + (3 * duration) // 4, time_initial, alignment)

    if window_start < time_initial or window_end > time_final or window_end <= window_start:
        raise RuntimeError(
            "Invalid stats window "
            f"[{window_start}, {window_end}] for case range [{time_initial}, {time_final}]"
        )
    if (window_start - time_initial) % dt_main != 0:
        raise RuntimeError(
            f"stats window start {window_start} does not align to dt_main={dt_main}"
        )
    if (window_end - time_initial) % dt_main != 0:
        raise RuntimeError(f"stats window end {window_end} does not align to dt_main={dt_main}")
    if (window_start - time_initial) % args.fine_tout != 0:
        raise RuntimeError(
            f"stats window start {window_start} does not align to "
            f"fine_tout={args.fine_tout}"
        )
    if (window_end - window_start) % args.fine_tout != 0:
        raise RuntimeError(
            f"stats window length {window_end - window_start} is not divisible by "
            f"fine_tout={args.fine_tout}"
        )
    if (window_end - window_start) % window_average_tout != 0:
        raise RuntimeError(
            f"stats window length {window_end - window_start} is not divisible by "
            f"window average tout={window_average_tout}"
        )

    return window_start, window_end


def window_label(stats_window: tuple[int, int]) -> str:
    """Return a compact path-safe label for a stats window."""
    return f"window{stats_window[0]}_{stats_window[1]}"


def run_subprocess(cmd: list[str], *, cwd: str | None = None) -> None:
    """Run one child process, streaming output directly to the terminal."""
    print(f"+ {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit {result.returncode}: {' '.join(cmd)}"
        )


def clear_dir(path: str) -> None:
    """Remove a previous run directory so stale files cannot affect the check."""
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def stats_path(run_dir: str, case_name: str) -> str:
    """Return the expected unified stats file path for one run."""
    return os.path.join(run_dir, f"{case_name}_stats.nc")


def build_run_command(args: argparse.Namespace, spec: RunSpec) -> list[str]:
    """Construct one run_scm.py invocation."""
    cmd = [
        sys.executable,
        RUN_SCM,
        "-stats",
        os.path.abspath(args.stats),
        "-out_dir",
        spec.output_dir,
        "-multicol",
        args.multicol,
        "-batch_size",
        str(spec.batch_size),
        "-tout",
        str(spec.tout),
    ]

    forwarded_paths = (
        ("-config", args.config),
        ("-params", args.params),
        ("-flags", args.flags),
        ("-silhs_params", args.silhs_params),
        ("-exe", args.exe),
    )
    for opt, value in forwarded_paths:
        if value is not None:
            cmd.extend([opt, os.path.abspath(value)])

    if args.debug is not None:
        cmd.extend(["-debug", str(args.debug)])
    if args.dt_main is not None:
        cmd.extend(["-dt_main", str(args.dt_main)])
    if args.dt_rad is not None:
        cmd.extend(["-dt_rad", str(args.dt_rad)])
    if spec.stats_window is not None:
        window_start, window_end = spec.stats_window
        cmd.extend([
            "-stats_tstart",
            str(window_start),
            "-stats_tend",
            str(window_end),
        ])

    cmd.append(args.case_name)
    return cmd


def time_axis_index(var: netCDF4.Variable) -> int:
    """Find the time axis in a stats variable using dimension names."""
    for axis, dim_name in enumerate(var.dimensions):
        if dim_name.lower() == "time":
            return axis
    raise RuntimeError(f"Variable {var.name} has no time dimension: {var.dimensions}")


def read_dataset_time_bounds(dataset: netCDF4.Dataset) -> np.ndarray:
    """Load time bounds as an array shaped (ntime, 2)."""
    if "time_bnds" not in dataset.variables:
        raise RuntimeError(f"{dataset.filepath()} is missing time_bnds")

    time_bnds_var = dataset.variables["time_bnds"]
    semantics = getattr(time_bnds_var, "interval_semantics", "").strip()
    if semantics != "(start, end]":
        raise RuntimeError(
            f"{dataset.filepath()} has unsupported time_bnds interval semantics "
            f"'{semantics}'. Expected '(start, end]'."
        )

    time_bnds = np.asarray(time_bnds_var[:], dtype=np.float64)
    dims = tuple(dim.lower() for dim in time_bnds_var.dimensions)
    if len(dims) == 2 and "time" in dims:
        if dims.index("time") == 1:
            time_bnds = np.transpose(time_bnds)

    if time_bnds.ndim != 2 or time_bnds.shape[1] != 2:
        raise RuntimeError(
            f"{dataset.filepath()} has unexpected time_bnds shape {time_bnds.shape}"
        )
    return time_bnds


def select_exact_tiling_records(
    time_bnds: np.ndarray,
    window_start: float,
    window_end: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Select fine records whose time_bnds exactly tile one coarse window."""
    selected = np.where(
        (time_bnds[:, 0] >= window_start - ABS_TOL)
        & (time_bnds[:, 1] <= window_end + ABS_TOL)
    )[0]
    if selected.size == 0:
        raise RuntimeError(
            f"No fine stats records fall inside window [{window_start}, {window_end}]"
        )

    selected_bnds = time_bnds[selected]
    if not np.isclose(selected_bnds[0, 0], window_start, atol=ABS_TOL, rtol=0.0):
        raise RuntimeError(
            f"Fine stats records start at {selected_bnds[0, 0]}, "
            f"expected {window_start}"
        )
    if not np.isclose(selected_bnds[-1, 1], window_end, atol=ABS_TOL, rtol=0.0):
        raise RuntimeError(
            f"Fine stats records end at {selected_bnds[-1, 1]}, expected {window_end}"
        )

    for left, right in zip(selected_bnds[:-1], selected_bnds[1:]):
        if not np.isclose(left[1], right[0], atol=ABS_TOL, rtol=0.0):
            raise RuntimeError(
                "Fine stats records do not tile contiguously: "
                f"{left.tolist()} then {right.tolist()}"
            )

    weights = selected_bnds[:, 1] - selected_bnds[:, 0]
    if np.any(weights <= 0.0):
        raise RuntimeError(f"Encountered non-positive time_bnds widths: {weights}")
    return selected, weights


def variable_payload(var: netCDF4.Variable) -> tuple[np.ndarray, np.ndarray]:
    """Read one NetCDF variable while preserving any mask separately."""
    data = var[:]
    if np.ma.isMaskedArray(data):
        mask = np.ma.getmaskarray(data)
        try:
            values = np.asarray(data.filled(0))
        except TypeError:
            values = np.asarray(data.filled(""))
    else:
        values = np.asarray(data)
        mask = np.zeros(values.shape, dtype=bool)
    return values, mask


def arrays_exact_equal(var_a: netCDF4.Variable, var_b: netCDF4.Variable) -> bool:
    """Return True when two variable payloads are exactly identical."""
    values_a, mask_a = variable_payload(var_a)
    values_b, mask_b = variable_payload(var_b)
    return np.array_equal(mask_a, mask_b) and np.array_equal(values_a, values_b)


def numeric_max_abs(var_a: netCDF4.Variable, var_b: netCDF4.Variable) -> float | None:
    """Return a helpful max-abs difference for numeric variables."""
    if not np.issubdtype(var_a.dtype, np.number) or not np.issubdtype(var_b.dtype, np.number):
        return None
    values_a, _ = variable_payload(var_a)
    values_b, _ = variable_payload(var_b)
    if values_a.shape != values_b.shape or values_a.size == 0:
        return None
    return float(np.max(np.abs(values_a.astype(np.float64) - values_b.astype(np.float64))))


def compare_dimensions_exact(
    reference: netCDF4.Dataset,
    test: netCDF4.Dataset,
) -> list[str]:
    """Compare dimension names, sizes, and unlimited flags."""
    mismatches: list[str] = []
    reference_dims = set(reference.dimensions)
    test_dims = set(test.dimensions)
    if reference_dims != test_dims:
        mismatches.append(
            "dimension set mismatch: "
            f"only reference={sorted(reference_dims - test_dims)}, "
            f"only test={sorted(test_dims - reference_dims)}"
        )
        return mismatches

    for dim_name in sorted(reference_dims):
        ref_dim = reference.dimensions[dim_name]
        test_dim = test.dimensions[dim_name]
        if len(ref_dim) != len(test_dim):
            mismatches.append(
                f"dimension {dim_name} length mismatch: {len(ref_dim)} vs {len(test_dim)}"
            )
        if ref_dim.isunlimited() != test_dim.isunlimited():
            mismatches.append(
                f"dimension {dim_name} unlimited flag mismatch: "
                f"{ref_dim.isunlimited()} vs {test_dim.isunlimited()}"
            )
    return mismatches


def compare_exact_stats(reference_path: str, test_path: str) -> tuple[list[str], str]:
    """Compare two stats files exactly at the variable-payload level."""
    mismatches: list[str] = []
    variables_checked = 0

    with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(test_path) as test:
        mismatches.extend(compare_dimensions_exact(reference, test))

        reference_vars = set(reference.variables)
        test_vars = set(test.variables)
        if reference_vars != test_vars:
            mismatches.append(
                "variable set mismatch: "
                f"only reference={sorted(reference_vars - test_vars)}, "
                f"only test={sorted(test_vars - reference_vars)}"
            )
            return mismatches, "variable sets differ"

        for var_name in sorted(reference_vars):
            ref_var = reference.variables[var_name]
            test_var = test.variables[var_name]
            variables_checked += 1

            if ref_var.dimensions != test_var.dimensions:
                mismatches.append(
                    f"{var_name}: dimension mismatch "
                    f"{ref_var.dimensions} vs {test_var.dimensions}"
                )
                continue
            if ref_var.shape != test_var.shape:
                mismatches.append(
                    f"{var_name}: shape mismatch {ref_var.shape} vs {test_var.shape}"
                )
                continue
            if not arrays_exact_equal(ref_var, test_var):
                max_abs = numeric_max_abs(ref_var, test_var)
                if max_abs is None:
                    mismatches.append(f"{var_name}: exact payload mismatch")
                else:
                    mismatches.append(
                        f"{var_name}: exact payload mismatch, max_abs={max_abs:.6e}"
                    )

    return mismatches, f"checked {variables_checked} variables exactly"


def compare_non_time_dimensions(
    fine: netCDF4.Dataset,
    coarse: netCDF4.Dataset,
) -> list[str]:
    """Compare dimensions except time, which is expected to differ."""
    mismatches: list[str] = []
    fine_dims = set(fine.dimensions) - {"time"}
    coarse_dims = set(coarse.dimensions) - {"time"}
    if fine_dims != coarse_dims:
        mismatches.append(
            "non-time dimension set mismatch: "
            f"only fine={sorted(fine_dims - coarse_dims)}, "
            f"only coarse={sorted(coarse_dims - fine_dims)}"
        )
        return mismatches

    for dim_name in sorted(fine_dims):
        if len(fine.dimensions[dim_name]) != len(coarse.dimensions[dim_name]):
            mismatches.append(
                f"dimension {dim_name} length mismatch: "
                f"{len(fine.dimensions[dim_name])} vs {len(coarse.dimensions[dim_name])}"
            )
    return mismatches


def weighted_time_mean(
    dataset: netCDF4.Dataset,
    var_name: str,
    record_indices: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Rebuild one coarse-window mean from fine stats records."""
    var = dataset.variables[var_name]
    axis = time_axis_index(var)
    data = np.asarray(var[:], dtype=np.float64)
    data = np.moveaxis(data, axis, 0)
    selected = data[record_indices, ...]
    return np.tensordot(weights, selected, axes=(0, 0)) / np.sum(weights)


def read_time_record(
    dataset: netCDF4.Dataset,
    var_name: str,
    record_index: int,
) -> np.ndarray:
    """Read one record from a time-dependent variable."""
    var = dataset.variables[var_name]
    axis = time_axis_index(var)
    data = np.asarray(var[:], dtype=np.float64)
    data = np.moveaxis(data, axis, 0)
    return data[record_index, ...]


def validate_time_bounds_window(
    time_bnds: np.ndarray,
    expected_window: tuple[int, int],
    label: str,
) -> list[str]:
    """Confirm time_bnds exactly covers the requested stats window."""
    mismatches: list[str] = []
    window_start, window_end = expected_window
    if time_bnds.size == 0:
        return [f"{label}: no time_bnds records were written"]
    if not np.isclose(time_bnds[0, 0], window_start, atol=ABS_TOL, rtol=0.0):
        mismatches.append(
            f"{label}: first time_bnds start is {time_bnds[0, 0]}, "
            f"expected {window_start}"
        )
    if not np.isclose(time_bnds[-1, 1], window_end, atol=ABS_TOL, rtol=0.0):
        mismatches.append(
            f"{label}: final time_bnds end is {time_bnds[-1, 1]}, "
            f"expected {window_end}"
        )
    if np.any(time_bnds[:, 0] < window_start - ABS_TOL) or np.any(
        time_bnds[:, 1] > window_end + ABS_TOL
    ):
        mismatches.append(
            f"{label}: at least one time_bnds record falls outside "
            f"[{window_start}, {window_end}]"
        )
    for left, right in zip(time_bnds[:-1], time_bnds[1:]):
        if not np.isclose(left[1], right[0], atol=ABS_TOL, rtol=0.0):
            mismatches.append(
                f"{label}: time_bnds do not tile contiguously: "
                f"{left.tolist()} then {right.tolist()}"
            )
            break
    return mismatches


def time_dependent_payload(var: netCDF4.Variable) -> tuple[np.ndarray, np.ndarray]:
    """Read a time-dependent variable with time moved to axis 0."""
    values, mask = variable_payload(var)
    axis = time_axis_index(var)
    return np.moveaxis(values, axis, 0), np.moveaxis(mask, axis, 0)


def array_max_abs(values_a: np.ndarray, values_b: np.ndarray) -> float | None:
    """Return max-abs difference for numeric arrays with matching shapes."""
    if values_a.shape != values_b.shape or values_a.size == 0:
        return None
    if not np.issubdtype(values_a.dtype, np.number) or not np.issubdtype(
        values_b.dtype, np.number
    ):
        return None
    return float(np.max(np.abs(values_a.astype(np.float64) - values_b.astype(np.float64))))


def compare_window_subset(
    full_path: str,
    window_path: str,
    expected_window: tuple[int, int],
) -> tuple[list[str], str]:
    """Compare a windowed fine-tout file to the matching full-run records."""
    mismatches: list[str] = []
    variables_checked = 0

    with netCDF4.Dataset(full_path) as full, netCDF4.Dataset(window_path) as windowed:
        mismatches.extend(compare_non_time_dimensions(full, windowed))

        full_vars = set(full.variables)
        window_vars = set(windowed.variables)
        if full_vars != window_vars:
            mismatches.append(
                "variable set mismatch: "
                f"only full={sorted(full_vars - window_vars)}, "
                f"only windowed={sorted(window_vars - full_vars)}"
            )
            return mismatches, "variable sets differ"

        full_bnds = read_dataset_time_bounds(full)
        window_bnds = read_dataset_time_bounds(windowed)
        mismatches.extend(
            validate_time_bounds_window(window_bnds, expected_window, window_path)
        )
        selected, _ = select_exact_tiling_records(
            full_bnds,
            float(expected_window[0]),
            float(expected_window[1]),
        )
        selected_bnds = full_bnds[selected]
        if window_bnds.shape != selected_bnds.shape or not np.allclose(
            window_bnds,
            selected_bnds,
            atol=ABS_TOL,
            rtol=0.0,
        ):
            mismatches.append(
                f"time_bnds mismatch for selected window: "
                f"full selected shape={selected_bnds.shape}, "
                f"windowed shape={window_bnds.shape}"
            )

        for var_name in sorted(full_vars):
            full_var = full.variables[var_name]
            window_var = windowed.variables[var_name]
            variables_checked += 1

            if full_var.dimensions != window_var.dimensions:
                mismatches.append(
                    f"{var_name}: dimension mismatch "
                    f"{full_var.dimensions} vs {window_var.dimensions}"
                )
                continue

            if "time" not in full_var.dimensions:
                if not arrays_exact_equal(full_var, window_var):
                    mismatches.append(f"{var_name}: non-time metadata mismatch")
                continue

            if var_name == "time_bnds":
                continue

            full_values, full_mask = time_dependent_payload(full_var)
            window_values, window_mask = time_dependent_payload(window_var)
            selected_values = full_values[selected, ...]
            selected_mask = full_mask[selected, ...]

            if selected_values.shape != window_values.shape:
                mismatches.append(
                    f"{var_name}: selected window shape mismatch "
                    f"{selected_values.shape} vs {window_values.shape}"
                )
                continue

            if not np.array_equal(selected_mask, window_mask) or not np.array_equal(
                selected_values,
                window_values,
            ):
                max_abs = array_max_abs(selected_values, window_values)
                if max_abs is None:
                    mismatches.append(f"{var_name}: selected window payload mismatch")
                else:
                    mismatches.append(
                        f"{var_name}: selected window payload mismatch, "
                        f"max_abs={max_abs:.6e}"
                    )

    detail = f"checked {variables_checked} variables against window {expected_window}"
    return mismatches, detail


def compare_tout_average(
    fine_path: str,
    coarse_path: str,
    expected_window: tuple[int, int] | None = None,
) -> tuple[list[str], str]:
    """Compare a coarse tout file to manual averages from a fine tout file."""
    mismatches: list[str] = []
    compared_vars = 0
    max_abs_seen = 0.0
    max_rel_seen = 0.0

    with netCDF4.Dataset(fine_path) as fine, netCDF4.Dataset(coarse_path) as coarse:
        mismatches.extend(compare_non_time_dimensions(fine, coarse))

        fine_vars = set(fine.variables)
        coarse_vars = set(coarse.variables)
        if fine_vars != coarse_vars:
            mismatches.append(
                "variable set mismatch: "
                f"only fine={sorted(fine_vars - coarse_vars)}, "
                f"only coarse={sorted(coarse_vars - fine_vars)}"
            )
            return mismatches, "variable sets differ"

        fine_bnds = read_dataset_time_bounds(fine)
        coarse_bnds = read_dataset_time_bounds(coarse)
        if expected_window is not None:
            mismatches.extend(
                validate_time_bounds_window(coarse_bnds, expected_window, coarse_path)
            )

        if "time" not in coarse.variables:
            mismatches.append(f"{coarse_path} is missing time")
        else:
            coarse_time = np.asarray(coarse.variables["time"][:], dtype=np.float64)
            if coarse_time.shape != (coarse_bnds.shape[0],):
                mismatches.append(
                    f"coarse time shape {coarse_time.shape} does not match "
                    f"time_bnds records {coarse_bnds.shape[0]}"
                )
            elif not np.allclose(coarse_time, coarse_bnds[:, 1], atol=ABS_TOL, rtol=0.0):
                mismatches.append("coarse time values do not equal time_bnds window ends")

        for var_name in sorted(fine_vars):
            fine_var = fine.variables[var_name]
            coarse_var = coarse.variables[var_name]

            if "time" not in fine_var.dimensions:
                if not arrays_exact_equal(fine_var, coarse_var):
                    mismatches.append(f"{var_name}: non-time metadata mismatch")
                continue

            if var_name in {"time", "time_bnds"}:
                continue
            if "time" not in coarse_var.dimensions:
                mismatches.append(f"{var_name}: fine has time dimension but coarse does not")
                continue
            if not np.issubdtype(fine_var.dtype, np.number):
                mismatches.append(f"{var_name}: time-dependent nonnumeric variable is unsupported")
                continue

            compared_vars += 1
            for coarse_record, bounds in enumerate(coarse_bnds):
                selected, weights = select_exact_tiling_records(
                    fine_bnds,
                    float(bounds[0]),
                    float(bounds[1]),
                )
                reconstructed = weighted_time_mean(fine, var_name, selected, weights)
                actual = read_time_record(coarse, var_name, coarse_record)

                if reconstructed.shape != actual.shape:
                    mismatches.append(
                        f"{var_name}[record={coarse_record}]: shape mismatch "
                        f"{reconstructed.shape} vs {actual.shape}"
                    )
                    continue

                diff = actual - reconstructed
                max_abs = float(np.max(np.abs(diff))) if diff.size else 0.0
                denom = np.maximum(np.abs(reconstructed), 1.0e-30)
                max_rel = float(np.max(np.abs(diff) / denom)) if diff.size else 0.0
                max_abs_seen = max(max_abs_seen, max_abs)
                max_rel_seen = max(max_rel_seen, max_rel)

                if not np.allclose(actual, reconstructed, atol=ABS_TOL, rtol=REL_TOL):
                    mismatches.append(
                        f"{var_name}[record={coarse_record}, "
                        f"window=({bounds[0]}, {bounds[1]}]] mismatch "
                        f"max_abs={max_abs:.6e}, max_rel={max_rel:.6e}"
                    )

    detail = (
        f"checked {compared_vars} time-dependent variables, "
        f"max_abs={max_abs_seen:.6e}, max_rel={max_rel_seen:.6e}"
    )
    return mismatches, detail


def run_specs(
    args: argparse.Namespace,
    stats_window: tuple[int, int],
) -> tuple[RunSpec, list[RunSpec], list[RunSpec], RunSpec, RunSpec]:
    """Build the reference, batch-check, and coarse-tout run specifications."""
    case_root = os.path.join(os.path.abspath(args.out_root), args.case_name)
    window_name = window_label(stats_window)
    reference = RunSpec(
        label=f"bs{args.batch_sizes[0]}_tout{args.fine_tout}",
        batch_size=args.batch_sizes[0],
        tout=args.fine_tout,
        output_dir=os.path.join(case_root, f"bs{args.batch_sizes[0]}_tout{args.fine_tout}"),
    )
    batch_runs = [
        RunSpec(
            label=f"bs{batch_size}_tout{args.fine_tout}",
            batch_size=batch_size,
            tout=args.fine_tout,
            output_dir=os.path.join(case_root, f"bs{batch_size}_tout{args.fine_tout}"),
        )
        for batch_size in args.batch_sizes[1:]
    ]
    coarse_runs = [
        RunSpec(
            label=f"bs{args.batch_sizes[0]}_tout{tout}",
            batch_size=args.batch_sizes[0],
            tout=tout,
            output_dir=os.path.join(case_root, f"bs{args.batch_sizes[0]}_tout{tout}"),
        )
        for tout in args.coarse_touts
    ]
    window_fine_run = RunSpec(
        label=f"bs{args.batch_sizes[0]}_tout{args.fine_tout}_{window_name}",
        batch_size=args.batch_sizes[0],
        tout=args.fine_tout,
        output_dir=os.path.join(
            case_root,
            f"bs{args.batch_sizes[0]}_tout{args.fine_tout}_{window_name}",
        ),
        stats_window=stats_window,
    )
    window_average_run = RunSpec(
        label=f"bs{args.batch_sizes[0]}_tout{args.coarse_touts[0]}_{window_name}",
        batch_size=args.batch_sizes[0],
        tout=args.coarse_touts[0],
        output_dir=os.path.join(
            case_root,
            f"bs{args.batch_sizes[0]}_tout{args.coarse_touts[0]}_{window_name}",
        ),
        stats_window=stats_window,
    )
    return reference, batch_runs, coarse_runs, window_fine_run, window_average_run


def parse_args() -> argparse.Namespace:
    """Define the CLI for the stats consistency test."""
    parser = argparse.ArgumentParser(
        description=(
            "Run one CLUBB case several ways and verify unified stats output is "
            "consistent across batch sizes and stats_tout averaging intervals."
        )
    )
    parser.add_argument(
        "case_name",
        nargs="?",
        default=DEFAULT_CASE,
        help=f"Case name to run. Default: {DEFAULT_CASE}",
    )
    parser.add_argument(
        "-out_root",
        default=DEFAULT_OUT_ROOT,
        help="Top-level output directory. Default: output/stats_output_consistency",
    )
    parser.add_argument(
        "-stats",
        default=DEFAULT_STATS_FILE,
        help=(
            "Stats registry file to use. Default: "
            f"{os.path.relpath(DEFAULT_STATS_FILE, CLUBB_ROOT)}"
        ),
    )
    parser.add_argument(
        "-multicol",
        default=DEFAULT_MULTICOL,
        help=f"run_scm.py multicol spec. Default: {DEFAULT_MULTICOL}",
    )
    parser.add_argument(
        "-batch_sizes",
        default=DEFAULT_BATCH_SIZES,
        help=(
            "Comma-separated batch sizes. The first value is the reference run. "
            f"Default: {DEFAULT_BATCH_SIZES}"
        ),
    )
    parser.add_argument(
        "-fine_tout",
        type=int,
        default=DEFAULT_FINE_TOUT,
        help=f"Reference stats_tout in seconds. Default: {DEFAULT_FINE_TOUT}",
    )
    parser.add_argument(
        "-coarse_touts",
        default=DEFAULT_COARSE_TOUTS,
        help=f"Comma-separated coarse stats_tout values. Default: {DEFAULT_COARSE_TOUTS}",
    )
    parser.add_argument(
        "-window_start",
        type=int,
        help=(
            "Optional stats_tstart for the windowing checks. "
            "Must be provided with -window_end."
        ),
    )
    parser.add_argument(
        "-window_end",
        type=int,
        help=(
            "Optional stats_tend for the windowing checks. "
            "Must be provided with -window_start."
        ),
    )
    parser.add_argument("-config", help="Optional config directory forwarded to run_scm.py")
    parser.add_argument("-params", help="Optional params file forwarded to run_scm.py")
    parser.add_argument("-flags", help="Optional model flags file forwarded to run_scm.py")
    parser.add_argument(
        "-silhs_params",
        help="Optional SILHS params file forwarded to run_scm.py",
    )
    parser.add_argument("-exe", help="Optional CLUBB executable forwarded to run_scm.py")
    parser.add_argument("-debug", type=int, help="Optional debug level forwarded to run_scm.py")
    parser.add_argument("-dt_main", type=int, help="Optional dt_main forwarded to run_scm.py")
    parser.add_argument("-dt_rad", type=int, help="Optional dt_rad forwarded to run_scm.py")
    args = parser.parse_args()

    args.batch_sizes = parse_csv_ints(args.batch_sizes, "-batch_sizes")
    args.coarse_touts = parse_csv_ints(args.coarse_touts, "-coarse_touts")
    if args.fine_tout <= 0:
        parser.error("-fine_tout must be positive")
    args.stats = abs_path(args.stats)
    args.config = abs_path(args.config)
    args.params = abs_path(args.params)
    args.flags = abs_path(args.flags)
    args.silhs_params = abs_path(args.silhs_params)
    args.exe = abs_path(args.exe)

    if not os.path.isfile(args.stats):
        parser.error(f"stats file not found: {args.stats}")
    return args


def main() -> None:
    """Run the configured stats consistency workflow."""
    args = parse_args()
    stats_window = derive_stats_window(args)
    reference, batch_runs, coarse_runs, window_fine_run, window_average_run = run_specs(
        args,
        stats_window,
    )
    all_runs = [reference] + batch_runs + coarse_runs + [window_fine_run, window_average_run]

    print_banner("Test Setup")
    print(f"Case: {args.case_name}")
    print(f"Stats file: {args.stats}")
    print(f"Multicol spec: {args.multicol}")
    print(f"Reference batch size: {reference.batch_size}")
    print(f"Batch sizes checked: {', '.join(str(run.batch_size) for run in batch_runs)}")
    print(f"Fine tout: {args.fine_tout}")
    print(f"Coarse touts: {', '.join(str(tout) for tout in args.coarse_touts)}")
    print(f"Window checks stats_tstart/stats_tend: {stats_window[0]}, {stats_window[1]}")
    print(f"Output root: {os.path.join(os.path.abspath(args.out_root), args.case_name)}")

    print("\nPlanned runs:")
    for spec in all_runs:
        window_detail = ""
        if spec.stats_window is not None:
            window_detail = (
                f", stats_tstart={spec.stats_window[0]}, stats_tend={spec.stats_window[1]}"
            )
        print(
            f" - {spec.label}: batch_size={spec.batch_size}, "
            f"tout={spec.tout}{window_detail}, output_dir={spec.output_dir}"
        )

    print_banner("Run Phase")
    for spec in all_runs:
        clear_dir(spec.output_dir)
        cmd = build_run_command(args, spec)
        print(f"Running {spec.label}")
        try:
            run_subprocess(cmd, cwd=RUN_SCRIPTS)
            path = stats_path(spec.output_dir, args.case_name)
            if not os.path.isfile(path):
                raise RuntimeError(f"Expected stats file not found: {path}")
            print_status(spec.label, True, f"stats file: {path}")
        except Exception as exc:
            print_status(spec.label, False, str(exc))
            raise

    mismatches: list[str] = []
    reference_path = stats_path(reference.output_dir, args.case_name)

    print_banner("Batch-Size Checks")
    for spec in batch_runs:
        test_path = stats_path(spec.output_dir, args.case_name)
        run_mismatches, detail = compare_exact_stats(reference_path, test_path)
        label = f"{spec.label}_vs_{reference.label}"
        if run_mismatches:
            mismatches.extend(f"{label}: {message}" for message in run_mismatches)
            print_status(label, False, "; ".join(run_mismatches[:5]))
        else:
            print_status(label, True, detail)

    print_banner("Tout Averaging Checks")
    for spec in coarse_runs:
        coarse_path = stats_path(spec.output_dir, args.case_name)
        run_mismatches, detail = compare_tout_average(reference_path, coarse_path)
        label = f"{spec.label}_manual_average_vs_{reference.label}"
        if run_mismatches:
            mismatches.extend(f"{label}: {message}" for message in run_mismatches)
            print_status(label, False, "; ".join(run_mismatches[:5]))
        else:
            print_status(label, True, detail)

    print_banner("Windowing Checks")
    window_fine_path = stats_path(window_fine_run.output_dir, args.case_name)
    run_mismatches, detail = compare_window_subset(
        reference_path,
        window_fine_path,
        stats_window,
    )
    label = f"{window_fine_run.label}_selected_window_vs_{reference.label}"
    if run_mismatches:
        mismatches.extend(f"{label}: {message}" for message in run_mismatches)
        print_status(label, False, "; ".join(run_mismatches[:5]))
    else:
        print_status(label, True, detail)

    window_average_path = stats_path(window_average_run.output_dir, args.case_name)
    run_mismatches, detail = compare_tout_average(
        reference_path,
        window_average_path,
        expected_window=stats_window,
    )
    label = f"{window_average_run.label}_manual_average_vs_{reference.label}"
    if run_mismatches:
        mismatches.extend(f"{label}: {message}" for message in run_mismatches)
        print_status(label, False, "; ".join(run_mismatches[:5]))
    else:
        print_status(label, True, detail)

    print_banner("Summary")
    if mismatches:
        print(f"Found {len(mismatches)} mismatch(es):")
        for message in mismatches[:20]:
            print(f" - {message}")
        if len(mismatches) > 20:
            print(f" - ... {len(mismatches) - 20} additional mismatch(es)")
        raise SystemExit(1)

    print("All stats consistency checks passed.")


if __name__ == "__main__":
    main()
