#!/usr/bin/env python3
"""
Check that the loss driver saves and reports internally consistent output.

This script checks two loss-driver contracts:

1. Saved stats consistency:
   - Build one shared 4-column parameter file.
   - Read the configured loss window for the requested case.
   - Run normal `run_scm.py` with stats settings that match `run_scm_loss.py`.
   - Run `run_scm_loss.py` for the same case, fields, and parameters.
   - Compare the requested fields in the two single-record stats files.

   The normal run is configured with the same loss-window stats cadence that
   `run_scm_loss.py` writes: samples are accumulated every dt_main, then one
   output record is written at the end of the requested loss window.

2. Reported metric consistency:
   - Parse the loss metric table printed by `run_scm_loss.py`.
   - Recompute those metrics from the saved loss stats file and normalized
     benchmark file.
   - Compare the reported values, including `best_col`, against the recomputed
     values.

   This matters because the loss metrics are printed from in-memory Fortran
   stats, not by reading the saved NetCDF file. Recomputing them from the saved
   artifacts catches cases where the reported loss and saved profiles drift
   apart.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

import netCDF4
import numpy as np


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.abspath(os.path.join(TESTS_DIR, ".."))
RUN_SCRIPTS = os.path.join(CLUBB_ROOT, "run_scripts")
if CLUBB_ROOT not in sys.path:
    sys.path.insert(0, CLUBB_ROOT)

from tuner.case_defaults import DEFAULT_LOSS_FIELDS, read_case_defaults  # noqa: E402
from utilities.loss_metrics import calculate_column_loss_metrics  # noqa: E402
from utilities.create_case_namelist import read_model_times  # noqa: E402


DEFAULT_OUT_ROOT = os.path.join(CLUBB_ROOT, "output", "loss_output_consistency")
DEFAULT_CONFIG_DIR = os.path.join(CLUBB_ROOT, "input", "tunable_parameters")
DEFAULT_STATS_FILE = os.path.join(CLUBB_ROOT, "input", "stats", "tuning_stats.in")
CREATE_MULTI_COL_PARAMS = os.path.join(CLUBB_ROOT, "utilities", "create_multi_col_params.py")
RUN_SCM = os.path.join(RUN_SCRIPTS, "run_scm.py")
RUN_SCM_LOSS = os.path.join(RUN_SCRIPTS, "run_scm_loss.py")

NUM_COLS = 4
HR_SPEC = "C8/0.2:0.8/4"

# The two runs should write the same values. Keep this tight, but allow a tiny
# margin for different executable paths or NetCDF/Python readback details.
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-12

# The loss log prints metrics with 8 digits after the decimal in scientific
# notation, so the log-vs-file metric check cannot use the tighter BFB
# tolerance used for the NetCDF-vs-NetCDF comparison above.
METRIC_ABS_TOL = 5.0e-8
METRIC_REL_TOL = 5.0e-8
TIME_MATCH_TOL = 1.0e-6

LOSS_METRIC_PREFIXES = {
    "scaled_rmse": "scaled_rmse_col",
    "corr": "corr_col",
    "std_ratio": "std_ratio_col",
    "crmse_norm": "crmse_norm_col",
    "bias_norm": "bias_norm_col",
}


def print_banner(title: str) -> None:
    """Print a simple section banner."""
    print(f"\n{'=' * 24} {title} {'=' * 24}", flush=True)


def print_status(label: str, passed: bool, detail: str) -> None:
    """Print a uniform PASS/FAIL line."""
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] {label}: {detail}\n", flush=True)


def parse_csv_values(values: list[str] | None, label: str) -> list[str]:
    """Parse one or more comma-separated CLI values."""
    parsed: list[str] = []
    for value in values or []:
        for item in str(value).split(","):
            item = item.strip()
            if item:
                parsed.append(item)
    if not parsed:
        raise RuntimeError(f"{label} must include at least one value")
    return parsed


def abs_path(path: str | None) -> str | None:
    """Return an absolute path when the user provided one."""
    return os.path.abspath(path) if path else None


def run_subprocess(cmd: list[str], *, cwd: str | None = None) -> None:
    """Run a child process and fail if it exits nonzero."""
    print(f"+ {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit {result.returncode}: {' '.join(cmd)}"
        )


def clear_dir(path: str) -> None:
    """Remove stale run output before writing a fresh result."""
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def resolve_base_params(args: argparse.Namespace) -> str:
    """Choose the scalar tunable-parameter file used to seed the 4 columns."""
    if args.params:
        return args.params
    if args.config:
        return os.path.join(args.config, "tunable_parameters.in")
    return os.path.join(DEFAULT_CONFIG_DIR, "tunable_parameters.in")


def build_shared_multicol_params(args: argparse.Namespace, out_root: str) -> str:
    """Create the shared multicol parameter file used by both runners."""
    params_file = resolve_base_params(args)
    if not os.path.isfile(params_file):
        raise RuntimeError(f"Base params file not found: {params_file}")

    shared_params = os.path.join(out_root, f"{args.case_name}_shared_multicol_params.in")
    cmd = [
        sys.executable,
        CREATE_MULTI_COL_PARAMS,
        "-param_file",
        params_file,
        "-out_file",
        shared_params,
        "-hr",
        HR_SPEC,
        "-batch_size",
        str(NUM_COLS),
    ]
    run_subprocess(cmd, cwd=RUN_SCRIPTS)
    return shared_params


def loss_window(case_name: str) -> tuple[float, float]:
    """Read the configured absolute loss window for a case."""
    case_defaults = read_case_defaults(case_name)
    window_start, window_end = [
        float(value) for value in case_defaults["time_average_range"]
    ]
    if window_end <= window_start:
        raise RuntimeError(
            f"Invalid time_average_range for {case_name}: {window_start}, {window_end}"
        )
    return window_start, window_end


def case_dt_main(case_name: str) -> float:
    """Read the model timestep used by the loss-driver namelist builder."""
    model_file = os.path.join(CLUBB_ROOT, "input", "case_setups", f"{case_name}_model.in")
    values = read_model_times(model_file)
    dt_main = values.get("dt_main")
    if dt_main is None:
        raise RuntimeError(f"Could not determine dt_main from {model_file}")
    return float(dt_main)


def common_forwarded_options(args: argparse.Namespace) -> list[str]:
    """Return CLI options both run_scm.py and run_scm_loss.py understand."""
    options: list[str] = []
    if args.config:
        options.extend(["-config", args.config])
    if args.flags:
        options.extend(["-flags", args.flags])
    return options


def normal_run_command(
    args: argparse.Namespace,
    shared_params: str,
    output_dir: str,
    window_start: float,
    window_end: float,
) -> list[str]:
    """Build a normal run_scm.py command that matches loss-run stats settings."""
    window_length = window_end - window_start
    dt_main = case_dt_main(args.case_name)
    override = (
        f"stats_tsamp={dt_main:.1f},"
        f"time_final={window_end:.1f},"
        "debug_level=-1"
    )
    return [
        sys.executable,
        RUN_SCM,
        *common_forwarded_options(args),
        "-params",
        shared_params,
        "-stats",
        DEFAULT_STATS_FILE,
        "-out_dir",
        output_dir,
        "-tout",
        str(int(window_length)),
        "-stats_tstart",
        f"{window_start:.1f}",
        "-stats_tend",
        f"{window_end:.1f}",
        "-override",
        override,
        args.case_name,
    ]


def loss_run_command(
    args: argparse.Namespace,
    shared_params: str,
    output_dir: str,
    fields: list[str],
) -> list[str]:
    """Build the loss-driver run command."""
    return [
        sys.executable,
        RUN_SCM_LOSS,
        *common_forwarded_options(args),
        "-params",
        shared_params,
        "-out_dir",
        output_dir,
        "-fields",
        ",".join(fields),
        args.case_name,
    ]


def stats_path(output_dir: str, case_name: str) -> str:
    """Return the expected unified stats file path."""
    return os.path.join(output_dir, f"{case_name}_stats.nc")


def loss_log_path(output_dir: str, case_name: str) -> str:
    """Return the run_scm_loss.py log path."""
    return os.path.join(output_dir, f"{case_name}_log")


def normalized_benchmark_path(output_dir: str, case_name: str) -> str:
    """Return the normalized benchmark file generated by run_scm_loss.py."""
    return os.path.join(output_dir, f"{case_name}_normalized_benchmark.nc")


def time_axis_index(var: netCDF4.Variable) -> int:
    """Find the time axis in a stats variable by dimension name."""
    for axis, dim_name in enumerate(var.dimensions):
        if dim_name.lower() == "time":
            return axis
    raise RuntimeError(f"Variable {var.name} has no time dimension: {var.dimensions}")


def read_dataset_time_bounds(dataset: netCDF4.Dataset) -> np.ndarray:
    """Load `time_bnds` and require the stats interval convention we compare."""
    if "time_bnds" not in dataset.variables:
        raise RuntimeError(f"{dataset.filepath()} is missing time_bnds")
    time_bnds_var = dataset.variables["time_bnds"]
    semantics = getattr(time_bnds_var, "interval_semantics", "").strip()
    if semantics != "(start, end]":
        raise RuntimeError(
            f"{dataset.filepath()} has unsupported time_bnds interval semantics "
            f"'{semantics}'. Expected '(start, end]'."
        )
    return np.asarray(time_bnds_var[:], dtype=np.float64)


def single_record(dataset: netCDF4.Dataset, var_name: str) -> np.ndarray:
    """Extract the single stats record expected in this focused test."""
    if var_name not in dataset.variables:
        raise RuntimeError(f"Variable {var_name} not found in {dataset.filepath()}")

    var = dataset.variables[var_name]
    axis = time_axis_index(var)
    data = np.asarray(var[:], dtype=np.float64)
    data = np.moveaxis(data, axis, 0)
    if data.shape[0] != 1:
        raise RuntimeError(
            f"Expected one record for {var_name} in {dataset.filepath()}, "
            f"found {data.shape[0]}"
        )
    return data[0, ...]


def compare_loss_against_normal(
    normal_path: str,
    loss_path: str,
    fields: list[str],
    window_start: float,
    window_end: float,
) -> tuple[list[str], list[str]]:
    """Compare loss stats fields against an equivalent normal stats file."""
    mismatches: list[str] = []
    diagnostics: list[str] = []
    expected_bnds = np.array([[window_start, window_end]], dtype=np.float64)

    with netCDF4.Dataset(normal_path) as normal_ds, netCDF4.Dataset(loss_path) as loss_ds:
        normal_bnds = read_dataset_time_bounds(normal_ds)
        loss_bnds = read_dataset_time_bounds(loss_ds)

        for label, bnds in (("normal", normal_bnds), ("loss", loss_bnds)):
            if bnds.shape != expected_bnds.shape:
                mismatches.append(
                    f"{label} time_bnds shape is {bnds.shape}, expected {expected_bnds.shape}"
                )
            elif not np.allclose(bnds, expected_bnds, atol=ABS_TOL, rtol=0.0):
                mismatches.append(
                    f"{label} time_bnds {bnds.tolist()} do not match "
                    f"{expected_bnds.tolist()}"
                )

        if "clubb_params" in normal_ds.variables and "clubb_params" in loss_ds.variables:
            normal_params = np.asarray(normal_ds.variables["clubb_params"][:], dtype=np.float64)
            loss_params = np.asarray(loss_ds.variables["clubb_params"][:], dtype=np.float64)
            if normal_params.shape != loss_params.shape or not np.array_equal(
                normal_params,
                loss_params,
            ):
                mismatches.append("clubb_params metadata differ between normal and loss runs")

        for var_name in fields:
            normal_values = single_record(normal_ds, var_name)
            loss_values = single_record(loss_ds, var_name)
            if normal_values.shape != loss_values.shape:
                mismatches.append(
                    f"{var_name}: shape mismatch, normal {normal_values.shape} "
                    f"vs loss {loss_values.shape}"
                )
                continue

            diff = loss_values - normal_values
            max_abs = float(np.max(np.abs(diff)))
            denom = np.maximum(np.abs(normal_values), 1.0e-30)
            max_rel = float(np.max(np.abs(diff) / denom))
            diagnostics.append(
                f"{var_name}(max_abs={max_abs:.6e}, max_rel={max_rel:.6e})"
            )

            if not np.allclose(loss_values, normal_values, atol=ABS_TOL, rtol=REL_TOL):
                mismatches.append(
                    f"{var_name}: mismatch "
                    f"(max_abs={max_abs:.6e}, max_rel={max_rel:.6e})"
                )

    return mismatches, diagnostics


def loss_metric_columns(header: list[str], prefix: str) -> list[int]:
    """Return header indices for a per-column loss metric, sorted by column."""
    columns: list[tuple[int, int]] = []
    for idx, name in enumerate(header):
        if not name.startswith(prefix):
            continue
        col_text = name[len(prefix):]
        try:
            col_num = int(col_text)
        except ValueError as exc:
            raise RuntimeError(f"Could not parse loss metric column name: {name}") from exc
        columns.append((col_num, idx))
    if not columns:
        raise RuntimeError(f"Loss log header is missing metric columns for {prefix}")
    return [idx for _, idx in sorted(columns)]


def parse_loss_log_metrics(output_dir: str, case_name: str) -> dict[str, dict[str, object]]:
    """Parse the row-oriented loss metric table printed by run_scm_loss.py."""
    path = loss_log_path(output_dir, case_name)
    if not os.path.isfile(path):
        raise RuntimeError(f"Expected loss log not found: {path}")

    parsed: dict[str, dict[str, object]] = {}
    header: list[str] | None = None
    metric_indices: dict[str, list[int]] = {}
    variable_idx = -1
    window_idx = -1
    best_col_idx = -1

    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            tokens = raw_line.split()
            if not tokens:
                continue

            if tokens[0] == "variable" and "best_col" in tokens:
                header = tokens
                variable_idx = header.index("variable")
                window_idx = header.index("window") if "window" in header else -1
                best_col_idx = header.index("best_col")
                metric_indices = {
                    metric_name: loss_metric_columns(header, prefix)
                    for metric_name, prefix in LOSS_METRIC_PREFIXES.items()
                }
                continue

            if header is None or len(tokens) != len(header):
                continue

            best_col = tokens[best_col_idx]
            if not best_col.startswith("col"):
                continue
            if window_idx >= 0 and tokens[window_idx] != "window1":
                raise RuntimeError(
                    "run_loss_output_consistency.py currently expects one loss window; "
                    f"found {tokens[window_idx]}"
                )

            var_name = tokens[variable_idx]
            parsed[var_name] = {
                metric_name: np.asarray(
                    [float(tokens[idx]) for idx in indices],
                    dtype=np.float64,
                )
                for metric_name, indices in metric_indices.items()
            }
            parsed[var_name]["best_col"] = best_col

    if not parsed:
        raise RuntimeError(f"Could not find a loss metric table in {path}")
    return parsed


def time_values_seconds(dataset: netCDF4.Dataset) -> np.ndarray:
    """Read a NetCDF time coordinate and convert simple units to seconds."""
    if "time" not in dataset.variables:
        raise RuntimeError(f"{dataset.filepath()} is missing time")

    time_var = dataset.variables["time"]
    units = getattr(time_var, "units", "seconds").strip().lower()
    first_word = units.split()[0] if units else "seconds"
    if first_word.startswith("hour"):
        multiplier = 3600.0
    elif first_word.startswith("minute"):
        multiplier = 60.0
    elif first_word.startswith("second"):
        multiplier = 1.0
    else:
        raise RuntimeError(
            f"{dataset.filepath()} has unsupported time units '{units}'"
        )
    return np.asarray(time_var[:], dtype=np.float64) * multiplier


def stats_profile_matrix(
    dataset: netCDF4.Dataset,
    var_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return one loss stats field as (level, column) plus its model grid."""
    if var_name not in dataset.variables:
        raise RuntimeError(f"Variable {var_name} not found in {dataset.filepath()}")

    var = dataset.variables[var_name]
    record = single_record(dataset, var_name)
    dims = list(var.dimensions)
    time_axis = time_axis_index(var)
    dims.pop(time_axis)

    col_axes = [idx for idx, dim_name in enumerate(dims) if dim_name.lower() == "col"]
    if len(col_axes) != 1:
        raise RuntimeError(
            f"Loss metric check requires one col dimension for {var_name}; "
            f"found {var.dimensions}"
        )
    col_axis = col_axes[0]

    vertical_axes = [
        idx
        for idx, dim_name in enumerate(dims)
        if idx != col_axis and dim_name.lower() in ("zt", "zm")
    ]
    if len(vertical_axes) != 1:
        raise RuntimeError(
            f"Loss metric check requires a zt or zm dimension for {var_name}; "
            f"found {var.dimensions}"
        )
    vertical_axis = vertical_axes[0]
    vertical_dim = dims[vertical_axis]

    if vertical_dim not in dataset.variables:
        raise RuntimeError(f"{dataset.filepath()} is missing grid variable {vertical_dim}")

    matrix = np.moveaxis(record, [vertical_axis, col_axis], [0, 1])
    if matrix.ndim != 2:
        raise RuntimeError(
            f"Expected {var_name} to reduce to (level, col); found {matrix.shape}"
        )

    model_z = np.asarray(dataset.variables[vertical_dim][:], dtype=np.float64)
    if model_z.ndim != 1 or model_z.shape[0] != matrix.shape[0]:
        raise RuntimeError(
            f"Grid {vertical_dim} shape {model_z.shape} does not match {var_name} "
            f"shape {matrix.shape}"
        )
    return matrix, model_z


def benchmark_profile_on_model_grid(
    benchmark_ds: netCDF4.Dataset,
    var_name: str,
    model_z: np.ndarray,
    window_start: float,
    window_end: float,
) -> np.ndarray:
    """Rebuild the benchmark profile used by the loss driver."""
    if var_name not in benchmark_ds.variables:
        raise RuntimeError(f"Variable {var_name} not found in {benchmark_ds.filepath()}")

    var = benchmark_ds.variables[var_name]
    dims = list(var.dimensions)
    try:
        time_axis = next(idx for idx, dim_name in enumerate(dims) if dim_name.lower() == "time")
    except StopIteration as exc:
        raise RuntimeError(f"Benchmark variable {var_name} has no time dimension") from exc

    vertical_candidates = [
        idx
        for idx, dim_name in enumerate(dims)
        if idx != time_axis
        and dim_name in benchmark_ds.variables
        and np.asarray(benchmark_ds.variables[dim_name][:]).ndim == 1
        and len(benchmark_ds.dimensions[dim_name]) > 1
    ]
    if len(vertical_candidates) != 1:
        raise RuntimeError(
            f"Could not identify one vertical coordinate for benchmark variable "
            f"{var_name}: {var.dimensions}"
        )
    vertical_axis = vertical_candidates[0]
    vertical_dim = dims[vertical_axis]

    source_z = np.asarray(benchmark_ds.variables[vertical_dim][:], dtype=np.float64)
    values = np.asarray(var[:], dtype=np.float64)
    values = np.moveaxis(values, [time_axis, vertical_axis], [0, 1])
    if values.ndim > 2:
        trailing_axes = tuple(range(2, values.ndim))
        if any(values.shape[axis] != 1 for axis in trailing_axes):
            raise RuntimeError(
                f"Benchmark variable {var_name} has unsupported trailing shape "
                f"{values.shape[2:]}"
            )
        values = np.squeeze(values, axis=trailing_axes)

    if values.shape[1] != source_z.shape[0]:
        raise RuntimeError(
            f"Benchmark grid {vertical_dim} shape {source_z.shape} does not match "
            f"{var_name} shape {values.shape}"
        )

    benchmark_times = time_values_seconds(benchmark_ds)
    time_mask = (
        (benchmark_times > window_start + TIME_MATCH_TOL)
        & (benchmark_times <= window_end + TIME_MATCH_TOL)
    )
    if not np.any(time_mask):
        raise RuntimeError(
            f"No benchmark records for {var_name} in ({window_start}, {window_end}]"
        )

    if model_z[0] < source_z[0] - TIME_MATCH_TOL or model_z[-1] > source_z[-1] + TIME_MATCH_TOL:
        raise RuntimeError(
            f"Model grid for {var_name} extends outside benchmark grid "
            f"({source_z[0]}, {source_z[-1]})"
        )

    if np.any(np.diff(source_z) < 0.0):
        sort_idx = np.argsort(source_z)
        source_z = source_z[sort_idx]
        values = values[:, sort_idx]

    interpolated = np.empty((int(np.count_nonzero(time_mask)), model_z.shape[0]), dtype=np.float64)
    for row_idx, profile in enumerate(values[time_mask, :]):
        interpolated[row_idx, :] = np.interp(model_z, source_z, profile)
    return np.mean(interpolated, axis=0)


def calculate_loss_metrics(
    model_matrix: np.ndarray,
    benchmark_profile: np.ndarray,
) -> dict[str, object]:
    """Mirror the loss-driver per-field metrics for all parameter columns."""
    raw_metrics = calculate_column_loss_metrics(model_matrix, benchmark_profile)
    metrics = {
        "scaled_rmse": raw_metrics["scaled_rmse"],
        "corr": raw_metrics["correlation"],
        "std_ratio": raw_metrics["std_ratio"],
        "crmse_norm": raw_metrics["centered_rmse_norm"],
        "bias_norm": raw_metrics["bias_norm"],
    }
    best_col_idx = int(np.argmin(metrics["scaled_rmse"])) + 1
    metrics["best_col"] = f"col{best_col_idx}"
    return metrics


def recompute_loss_metrics_from_netcdf(
    stats_file: str,
    benchmark_file: str,
    case_name: str,
    fields: list[str],
    window_start: float,
    window_end: float,
) -> dict[str, dict[str, object]]:
    """Recompute reported loss metrics from saved stats and benchmark files."""
    case_defaults = read_case_defaults(case_name)
    z_min, z_max = [
        float(value) for value in case_defaults["altitude_comparison_range"]
    ]
    recomputed: dict[str, dict[str, object]] = {}

    with netCDF4.Dataset(stats_file) as stats_ds, netCDF4.Dataset(benchmark_file) as benchmark_ds:
        for var_name in fields:
            model_matrix_all, model_z_all = stats_profile_matrix(stats_ds, var_name)
            selected = np.nonzero((model_z_all >= z_min) & (model_z_all <= z_max))[0]
            if selected.size == 0:
                raise RuntimeError(
                    f"Altitude range [{z_min}, {z_max}] selected no levels for {var_name}"
                )

            model_matrix = model_matrix_all[selected, :]
            model_z = model_z_all[selected]
            benchmark_profile = benchmark_profile_on_model_grid(
                benchmark_ds,
                var_name,
                model_z,
                window_start,
                window_end,
            )
            recomputed[var_name] = calculate_loss_metrics(model_matrix, benchmark_profile)

    return recomputed


def compare_reported_metrics_against_netcdf(
    loss_dir: str,
    case_name: str,
    fields: list[str],
    window_start: float,
    window_end: float,
) -> tuple[list[str], list[str]]:
    """Compare printed loss metrics with metrics recomputed from NetCDF output."""
    benchmark_file = normalized_benchmark_path(loss_dir, case_name)
    loss_stats_file = stats_path(loss_dir, case_name)
    if not os.path.isfile(benchmark_file):
        raise RuntimeError(f"Expected normalized benchmark file not found: {benchmark_file}")

    reported = parse_loss_log_metrics(loss_dir, case_name)
    recomputed = recompute_loss_metrics_from_netcdf(
        loss_stats_file,
        benchmark_file,
        case_name,
        fields,
        window_start,
        window_end,
    )

    mismatches: list[str] = []
    diagnostics: list[str] = []

    for var_name in fields:
        if var_name not in reported:
            mismatches.append(f"{var_name}: missing from printed loss table")
            continue
        if var_name not in recomputed:
            mismatches.append(f"{var_name}: missing from recomputed metrics")
            continue

        var_diagnostics: list[str] = []
        for metric_name in LOSS_METRIC_PREFIXES:
            reported_values = reported[var_name][metric_name]
            recomputed_values = recomputed[var_name][metric_name]
            assert isinstance(reported_values, np.ndarray)
            assert isinstance(recomputed_values, np.ndarray)
            if reported_values.shape != recomputed_values.shape:
                mismatches.append(
                    f"{var_name}:{metric_name}: shape mismatch, reported "
                    f"{reported_values.shape} vs recomputed {recomputed_values.shape}"
                )
                continue

            diff = reported_values - recomputed_values
            max_abs = float(np.max(np.abs(diff)))
            denom = np.maximum(np.abs(recomputed_values), 1.0e-30)
            max_rel = float(np.max(np.abs(diff) / denom))
            var_diagnostics.append(
                f"{metric_name}(max_abs={max_abs:.6e}, max_rel={max_rel:.6e})"
            )

            if not np.allclose(
                reported_values,
                recomputed_values,
                atol=METRIC_ABS_TOL,
                rtol=METRIC_REL_TOL,
            ):
                mismatches.append(
                    f"{var_name}:{metric_name}: mismatch "
                    f"(max_abs={max_abs:.6e}, max_rel={max_rel:.6e})"
                )

        reported_best = reported[var_name]["best_col"]
        recomputed_best = recomputed[var_name]["best_col"]
        if reported_best != recomputed_best:
            mismatches.append(
                f"{var_name}: best_col mismatch, reported {reported_best} "
                f"vs recomputed {recomputed_best}"
            )
        diagnostics.append(f"{var_name}({'; '.join(var_diagnostics)})")

    return mismatches, diagnostics


def parse_args() -> argparse.Namespace:
    """Define the CLI for the loss consistency test."""
    parser = argparse.ArgumentParser(
        description=(
            "Run normal CLUBB and the loss driver for one case, then verify the "
            "loss-driver stats output and printed metrics are consistent."
        )
    )
    parser.add_argument("case_name", help="Case name, e.g. arm or bomex")
    parser.add_argument(
        "-out_root",
        default=DEFAULT_OUT_ROOT,
        help="Top-level output directory. Default: output/loss_output_consistency",
    )
    parser.add_argument(
        "-config",
        help=(
            "Optional config directory forwarded to run_scm.py/run_scm_loss.py. "
            "Its tunable_parameters.in seeds the shared multicol params unless "
            "-params is also supplied."
        ),
    )
    parser.add_argument(
        "-params",
        help="Optional scalar tunable parameters file used to seed the shared params.",
    )
    parser.add_argument(
        "-flags",
        help="Optional configurable_model_flags.in file forwarded to both runners.",
    )
    parser.add_argument(
        "-fields",
        nargs="+",
        help=(
            "Comma-separated CLUBB fields to compare. Defaults to the loss-driver "
            "field default."
        ),
    )
    args = parser.parse_args()

    args.config = abs_path(args.config)
    args.params = abs_path(args.params)
    args.flags = abs_path(args.flags)
    if args.config and not os.path.isdir(args.config):
        parser.error(f"config directory not found: {args.config}")
    if args.params and not os.path.isfile(args.params):
        parser.error(f"params file not found: {args.params}")
    if args.flags and not os.path.isfile(args.flags):
        parser.error(f"flags file not found: {args.flags}")
    if not os.path.isfile(DEFAULT_STATS_FILE):
        parser.error(f"stats file not found: {DEFAULT_STATS_FILE}")
    return args


def main() -> None:
    """Run the focused loss consistency workflow."""
    args = parse_args()
    fields = parse_csv_values(args.fields, "-fields") if args.fields else list(DEFAULT_LOSS_FIELDS)
    window_start, window_end = loss_window(args.case_name)

    case_root = os.path.join(os.path.abspath(args.out_root), args.case_name)
    normal_dir = os.path.join(case_root, "normal")
    loss_dir = os.path.join(case_root, "loss")
    clear_dir(case_root)

    shared_params = build_shared_multicol_params(args, case_root)

    print_banner("Test Setup")
    print(f"Case: {args.case_name}")
    print(f"Fields: {', '.join(fields)}")
    print(f"Loss window: [{window_start}, {window_end}]")
    print(f"Shared params: {shared_params}")
    print(f"Output root: {case_root}")

    print_banner("Run Phase")
    clear_dir(normal_dir)
    normal_cmd = normal_run_command(
        args,
        shared_params,
        normal_dir,
        window_start,
        window_end,
    )
    print("Running normal CLUBB")
    run_subprocess(normal_cmd, cwd=RUN_SCRIPTS)
    normal_path = stats_path(normal_dir, args.case_name)
    if not os.path.isfile(normal_path):
        raise RuntimeError(f"Expected normal stats file not found: {normal_path}")
    print_status("normal_run", True, normal_path)

    clear_dir(loss_dir)
    loss_cmd = loss_run_command(args, shared_params, loss_dir, fields)
    print("Running loss driver")
    run_subprocess(loss_cmd, cwd=RUN_SCRIPTS)
    loss_path = stats_path(loss_dir, args.case_name)
    if not os.path.isfile(loss_path):
        raise RuntimeError(f"Expected loss stats file not found: {loss_path}")
    print_status("loss_run", True, loss_path)

    print_banner("Comparison")
    mismatches, diagnostics = compare_loss_against_normal(
        normal_path,
        loss_path,
        fields,
        window_start,
        window_end,
    )
    if mismatches:
        print_status(
            "loss_vs_equivalent_normal",
            False,
            "; ".join(mismatches + diagnostics),
        )
        print_banner("Summary")
        print(f"Found {len(mismatches)} mismatch(es):")
        for message in mismatches:
            print(f" - {message}")
        raise SystemExit(1)

    print_status(
        "loss_vs_equivalent_normal",
        True,
        "; ".join(diagnostics),
    )

    mismatches, diagnostics = compare_reported_metrics_against_netcdf(
        loss_dir,
        args.case_name,
        fields,
        window_start,
        window_end,
    )
    if mismatches:
        print_status(
            "reported_loss_metrics_vs_netcdf",
            False,
            "; ".join(mismatches + diagnostics),
        )
        print_banner("Summary")
        print(f"Found {len(mismatches)} mismatch(es):")
        for message in mismatches:
            print(f" - {message}")
        raise SystemExit(1)

    print_status(
        "reported_loss_metrics_vs_netcdf",
        True,
        "; ".join(diagnostics),
    )

    print_banner("Summary")
    print("Loss-driver stats output and reported metrics are internally consistent.")


if __name__ == "__main__":
    main()
