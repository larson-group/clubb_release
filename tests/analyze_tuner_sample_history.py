#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Tuner sample-history validation
#
# The Dash/CLI tuner writes sample-history files so we can inspect every
# parameter set that was evaluated, not just the final best result.  Each tuner
# job writes one or more compressed NumPy chunks named
# sample_history_000001.npz, sample_history_000002.npz, and so on.  The chunks
# contain:
#
#   - all_params: every sampled parameter vector
#   - loss_metrics: every metric for every case/field/time-window observation
#   - sample_id and batch_id: row identifiers for the sample axis
#   - metadata arrays that describe parameters, cases, fields, metrics, and
#     the observation axis
#
# This file is a Jenkins-friendly test for that contract.  It intentionally
# fails loudly if history is missing, chunks are malformed, metadata does not
# match across chunks, sample IDs are incomplete, numeric arrays contain NaNs or
# infinities, or the observation axis does not cover every expected
# case/field/window combination exactly once.
#
# When the input path is a single job directory or one .npz file, only that
# history is validated.  When the input path is a parent output directory, the
# test recursively discovers and validates every tuner job with sample-history
# chunks below it.  After validation, it prints a small SVD diagnostic that is
# useful for humans but is not the main pass/fail mechanism.
# -----------------------------------------------------------------------------
"""Validate and summarize tuner sample-history NumPy chunks."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tuner.sample_history import CHUNK_PREFIX, SAMPLE_HISTORY_DIR, sample_history_paths  # noqa: E402


METADATA_ARRAYS = (
    "schema_version",
    "param_names",
    "case_names",
    "field_names",
    "metric_names",
    "case_window_counts",
    "obs_case",
    "obs_field",
    "obs_window",
    "obs_window_start_seconds",
    "obs_window_end_seconds",
)
SAMPLE_ARRAYS = (
    "all_params",
    "loss_metrics",
    "sample_id",
    "batch_id",
)
REQUIRED_ARRAYS = METADATA_ARRAYS + SAMPLE_ARRAYS
EXPECTED_SCHEMA_VERSION = "sample_history_observation_v1"


def history_paths(path: Path) -> list[Path]:
    """Return sample-history chunks from a job directory or explicit file."""
    path = Path(path)
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise RuntimeError(f"Sample-history path does not exist: {path}")
    chunks = sample_history_paths(path)
    if not chunks:
        raise RuntimeError(f"No {CHUNK_PREFIX}*.npz chunks found in {path} or {path / SAMPLE_HISTORY_DIR}")
    return chunks


def _chunk_glob() -> str:
    return f"{CHUNK_PREFIX}[0-9][0-9][0-9][0-9][0-9][0-9].npz"


def history_targets(path: Path) -> list[Path]:
    """Return every sample-history target represented by a file, job dir, or parent dir."""
    path = Path(path)
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise RuntimeError(f"Sample-history path does not exist: {path}")

    if sample_history_paths(path):
        return [path]

    targets = []
    seen = set()
    for chunk_path in sorted(path.rglob(_chunk_glob())):
        target = chunk_path.parent.parent if chunk_path.parent.name == SAMPLE_HISTORY_DIR else chunk_path.parent
        resolved = target.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        targets.append(target)
    if not targets:
        raise RuntimeError(
            f"No {CHUNK_PREFIX}*.npz chunks found in {path} or any child sample-history directory"
        )
    return targets


def _chunk_index(path: Path) -> int | None:
    stem = path.stem
    if not stem.startswith(CHUNK_PREFIX):
        return None
    suffix = stem[len(CHUNK_PREFIX) :]
    return int(suffix) if suffix.isdigit() else None


def check_chunk_paths(paths: list[Path]) -> list[str]:
    """Return failures for chunk presence, naming, and sequence checks."""
    failures = []
    if not paths:
        return ["no sample-history chunks were found"]
    indices = [_chunk_index(path) for path in paths]
    if any(index is None for index in indices):
        failures.append(f"chunk names must match {CHUNK_PREFIX}NNNNNN.npz")
        return failures
    expected = list(range(1, len(paths) + 1))
    if indices != expected:
        failures.append(f"chunk indices are not sequential from 1: found {indices}, expected {expected}")
    return failures


def _same_array(left: np.ndarray, right: np.ndarray) -> bool:
    if left.shape != right.shape:
        return False
    if np.issubdtype(left.dtype, np.floating) or np.issubdtype(right.dtype, np.floating):
        return np.array_equal(left, right, equal_nan=True)
    return np.array_equal(left, right)


def load_history(path: Path) -> dict:
    """Load chunks and concatenate arrays along the sample axis."""
    paths = history_paths(path)
    chunk_path_failures = check_chunk_paths(paths)
    if chunk_path_failures:
        raise RuntimeError("; ".join(chunk_path_failures))
    chunks = [np.load(chunk_path) for chunk_path in paths]
    first = chunks[0]
    for path, chunk in zip(paths, chunks):
        missing = [name for name in REQUIRED_ARRAYS if name not in chunk.files]
        if missing:
            raise RuntimeError(f"{path} is missing required arrays: {', '.join(missing)}")
        if chunk["all_params"].shape[0] == 0:
            raise RuntimeError(f"{path} contains no samples")

    for chunk in chunks[1:]:
        for name in METADATA_ARRAYS:
            if not _same_array(first[name], chunk[name]):
                raise RuntimeError(f"Metadata mismatch across chunks for {name}")

    history = {
        name: first[name].copy()
        for name in METADATA_ARRAYS
    }
    history["all_params"] = np.concatenate([chunk["all_params"] for chunk in chunks], axis=0)
    history["loss_metrics"] = np.concatenate([chunk["loss_metrics"] for chunk in chunks], axis=0)
    history["sample_id"] = np.concatenate([chunk["sample_id"] for chunk in chunks], axis=0)
    history["batch_id"] = np.concatenate([chunk["batch_id"] for chunk in chunks], axis=0)
    history["chunk_paths"] = paths
    for chunk in chunks:
        chunk.close()
    return history


def check_history(history: dict, *, min_samples: int = 1, require_param_variation: bool = False) -> list[str]:
    """Return a list of consistency-check failure messages."""
    failures = []
    all_params = history["all_params"]
    loss_metrics = history["loss_metrics"]
    sample_id = history["sample_id"]
    batch_id = history["batch_id"]

    expected_param_count = len(history["param_names"])
    expected_obs_count = len(history["obs_case"])
    expected_metric_count = len(history["metric_names"])
    sample_count = all_params.shape[0]

    schema_version = history["schema_version"]
    if schema_version.shape != (1,) or str(schema_version[0]) != EXPECTED_SCHEMA_VERSION:
        failures.append(
            f"schema_version is not {EXPECTED_SCHEMA_VERSION}: "
            f"{schema_version.tolist() if hasattr(schema_version, 'tolist') else schema_version}"
        )
    if sample_count < int(min_samples):
        failures.append(f"sample count {sample_count} is less than required minimum {int(min_samples)}")
    if expected_param_count < 1:
        failures.append("param_names is empty")
    if expected_metric_count < 1:
        failures.append("metric_names is empty")
    if len(history["case_names"]) < 1:
        failures.append("case_names is empty")
    if len(history["field_names"]) < 1:
        failures.append("field_names is empty")
    if expected_obs_count < 1:
        failures.append("observation axis is empty")
    if len(set(map(str, history["param_names"]))) != expected_param_count:
        failures.append("param_names contains duplicates")
    if len(set(map(str, history["metric_names"]))) != expected_metric_count:
        failures.append("metric_names contains duplicates")
    if all_params.ndim != 2 or all_params.shape[1] != expected_param_count:
        failures.append("all_params shape does not match param_names")
    if loss_metrics.shape != (sample_count, expected_obs_count, expected_metric_count):
        failures.append("loss_metrics shape does not match sample/observation/metric axes")
    if sample_id.shape != (sample_count,):
        failures.append("sample_id length does not match sample count")
    if batch_id.shape != (sample_count,):
        failures.append("batch_id length does not match sample count")
    if len(np.unique(sample_id)) != sample_count:
        failures.append("sample_id contains duplicates")
    if not np.all(np.isfinite(all_params)):
        failures.append("all_params contains non-finite values")
    if not np.all(np.isfinite(loss_metrics)):
        failures.append("loss_metrics contains non-finite values")
    if not np.array_equal(np.sort(sample_id), np.arange(sample_count)):
        failures.append("sample_id values are not exactly 0..sample_count-1")
    if not np.all(batch_id >= 0):
        failures.append("batch_id contains negative values")
    if (
        require_param_variation
        and all_params.ndim == 2
        and all_params.shape[0] >= 2
        and not np.any(np.ptp(all_params, axis=0) > 0.0)
    ):
        failures.append("all parameter columns are constant despite --require-param-variation")
    failures.extend(check_observation_axis(history))
    return failures


def check_observation_axis(history: dict) -> list[str]:
    """Return failures for observation metadata and case/field/window coverage."""
    failures = []
    obs_arrays = [
        "obs_case",
        "obs_field",
        "obs_window",
        "obs_window_start_seconds",
        "obs_window_end_seconds",
    ]
    lengths = {name: len(history[name]) for name in obs_arrays}
    if len(set(lengths.values())) != 1:
        failures.append(f"observation metadata lengths do not match: {lengths}")
        return failures

    case_count = len(history["case_names"])
    field_count = len(history["field_names"])
    case_window_counts = np.asarray(history["case_window_counts"], dtype=int)
    if case_window_counts.shape != (case_count,):
        failures.append("case_window_counts length does not match case_names")
        return failures
    if np.any(case_window_counts < 1):
        failures.append("case_window_counts contains values < 1")

    obs_case = np.asarray(history["obs_case"], dtype=int)
    obs_field = np.asarray(history["obs_field"], dtype=int)
    obs_window = np.asarray(history["obs_window"], dtype=int)
    starts = np.asarray(history["obs_window_start_seconds"], dtype=float)
    ends = np.asarray(history["obs_window_end_seconds"], dtype=float)

    if np.any((obs_case < 0) | (obs_case >= case_count)):
        failures.append("obs_case contains out-of-range case indices")
    if np.any((obs_field < 0) | (obs_field >= field_count)):
        failures.append("obs_field contains out-of-range field indices")
    for obs_idx, (case_idx, window_idx) in enumerate(zip(obs_case, obs_window)):
        if 0 <= case_idx < case_count and not (0 <= window_idx < case_window_counts[case_idx]):
            failures.append(f"obs_window[{obs_idx}] is out of range for its case")
            break

    finite_start = np.isfinite(starts)
    finite_end = np.isfinite(ends)
    if not np.array_equal(finite_start, finite_end):
        failures.append("observation window starts/ends must be both finite or both NaN")
    if np.any(finite_start & (ends <= starts)):
        failures.append("finite observation windows must have end > start")

    expected = {
        (case_idx, field_idx, window_idx)
        for case_idx, window_count in enumerate(case_window_counts)
        for field_idx in range(field_count)
        for window_idx in range(int(window_count))
    }
    actual = set(zip(obs_case.tolist(), obs_field.tolist(), obs_window.tolist()))
    if actual != expected:
        failures.append(
            f"observation axis does not exactly cover case/field/window grid "
            f"(actual={len(actual)}, expected={len(expected)})"
        )
    if len(actual) != len(obs_case):
        failures.append("observation axis contains duplicate case/field/window entries")
    return failures


def observation_labels(history: dict, *, metric_name: str | None = None) -> list[str]:
    """Return readable labels for observation columns."""
    case_names = list(history["case_names"])
    field_names = list(history["field_names"])
    starts = history["obs_window_start_seconds"]
    ends = history["obs_window_end_seconds"]
    labels = []
    for obs_idx, (case_idx, field_idx, window_idx) in enumerate(
        zip(history["obs_case"], history["obs_field"], history["obs_window"])
    ):
        start = starts[obs_idx]
        end = ends[obs_idx]
        if np.isfinite(start) and np.isfinite(end):
            window_text = f"w{int(window_idx) + 1}:{int(start)}-{int(end)}s"
        else:
            window_text = f"w{int(window_idx) + 1}"
        base = f"{case_names[int(case_idx)]}/{field_names[int(field_idx)]}/{window_text}"
        labels.append(base if metric_name is None else f"{base}/{metric_name}")
    return labels


def svd_matrix(history: dict, metric: str) -> tuple[np.ndarray, list[str]]:
    """Return the loss matrix and labels used for SVD."""
    metric_names = list(history["metric_names"])
    loss_metrics = history["loss_metrics"]
    if metric == "all":
        labels = []
        pieces = []
        for metric_idx, metric_name in enumerate(metric_names):
            pieces.append(loss_metrics[:, :, metric_idx])
            labels.extend(observation_labels(history, metric_name=str(metric_name)))
        return np.concatenate(pieces, axis=1), labels

    if metric not in metric_names:
        raise RuntimeError(f"Unknown metric '{metric}'. Use one of {metric_names} or 'all'.")
    metric_idx = metric_names.index(metric)
    return loss_metrics[:, :, metric_idx], observation_labels(history, metric_name=metric)


def standardized_svd(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return SVD of centered/scaled finite columns and the retained column mask."""
    finite_columns = np.all(np.isfinite(matrix), axis=0)
    work = matrix[:, finite_columns]
    if work.size == 0:
        raise RuntimeError("No finite SVD columns remain")
    work = work - np.mean(work, axis=0, keepdims=True)
    scale = np.std(work, axis=0, keepdims=True)
    nonconstant_mask = scale[0] > 0.0
    if not np.any(nonconstant_mask):
        raise RuntimeError("No non-constant SVD columns remain")
    work = work[:, nonconstant_mask]
    scale = scale[:, nonconstant_mask]
    work = work / scale
    u, singular_values, vt = np.linalg.svd(work, full_matrices=False)
    retained_columns = finite_columns.copy()
    finite_indices = np.flatnonzero(finite_columns)
    retained_columns[finite_indices[~nonconstant_mask]] = False
    return u, singular_values, vt, retained_columns


def top_param_correlations(params: np.ndarray, param_names: list[str], scores: np.ndarray, limit: int) -> list[str]:
    """Return parameter names with largest absolute correlation to one score vector."""
    centered_params = params - np.mean(params, axis=0, keepdims=True)
    centered_scores = scores - np.mean(scores)
    param_std = np.std(centered_params, axis=0)
    score_std = float(np.std(centered_scores))
    if score_std == 0.0:
        return []
    correlations = np.zeros(params.shape[1], dtype=float)
    active = param_std > 0.0
    correlations[active] = (
        np.mean(centered_params[:, active] * centered_scores[:, None], axis=0)
        / (param_std[active] * score_std)
    )
    active_indices = np.flatnonzero(active)
    top_indices = active_indices[np.argsort(np.abs(correlations[active_indices]))[::-1][:limit]]
    return [
        f"{param_names[idx]}={correlations[idx]:+.3f}"
        for idx in top_indices
        if np.isfinite(correlations[idx])
    ]


def print_check_descriptions(*, min_samples: int, require_param_variation: bool) -> None:
    """Print the validation contract enforced before diagnostics run."""
    print("Checks performed")
    print("  - sample-history chunks exist and are sequentially numbered")
    print("  - every chunk contains required metadata, sample, and loss arrays")
    print("  - metadata arrays match exactly across chunks")
    print("  - schema version matches the current observation-axis format")
    print("  - sample IDs are unique and cover 0..sample_count-1")
    print("  - parameter and loss arrays have the expected shapes and finite values")
    print("  - observation metadata covers every case/field/window exactly once")
    print(f"  - at least {int(min_samples)} sample(s) are present")
    if require_param_variation:
        print("  - at least one tuned parameter varies across samples")


def print_summary(
    history: dict,
    *,
    metric: str,
    modes: int,
    top: int,
    min_samples: int,
    require_param_variation: bool,
) -> None:
    """Print consistency and SVD summaries."""
    chunk_paths = history["chunk_paths"]
    sample_count = history["all_params"].shape[0]
    metric_names = list(history["metric_names"])
    print("Tuner Sample History")
    print(f"  job: {chunk_paths[0].parent.parent if chunk_paths[0].parent.name == SAMPLE_HISTORY_DIR else chunk_paths[0].parent}")
    print(f"  chunks: {len(chunk_paths)}")
    print(f"  first_chunk: {chunk_paths[0]}")
    print(f"  samples: {sample_count}")
    print(f"  params: {len(history['param_names'])}")
    print(f"  cases: {len(history['case_names'])}")
    print(f"  fields: {len(history['field_names'])}")
    print(f"  observations: {len(history['obs_case'])}")
    print(f"  metrics: {', '.join(str(name) for name in metric_names)}")
    print_check_descriptions(
        min_samples=min_samples,
        require_param_variation=require_param_variation,
    )

    failures = check_history(
        history,
        min_samples=min_samples,
        require_param_variation=require_param_variation,
    )
    if failures:
        print("  checks: FAIL")
        for failure in failures:
            print(f"    - {failure}")
        raise RuntimeError("Sample-history consistency checks failed")
    print("  checks: OK")

    if sample_count < 2:
        print("SVD skipped: at least two samples are required.")
        return

    matrix, labels = svd_matrix(history, metric)
    try:
        u, singular_values, vt, retained_columns = standardized_svd(matrix)
    except RuntimeError as exc:
        print(f"SVD skipped: {exc}")
        return
    labels = [label for label, keep in zip(labels, retained_columns) if keep]
    variance = singular_values**2
    explained = variance / np.sum(variance) if np.sum(variance) > 0.0 else np.zeros_like(variance)

    print(f"SVD ({metric})")
    print(f"  matrix_shape: {matrix.shape[0]} samples x {len(labels)} retained columns")
    for mode_idx in range(min(modes, len(singular_values))):
        print(
            f"  mode {mode_idx + 1}: singular={singular_values[mode_idx]:.6g}, "
            f"explained={explained[mode_idx] * 100.0:.2f}%"
        )
        loading = vt[mode_idx]
        top_obs = np.argsort(np.abs(loading))[::-1][:top]
        print("    strongest observations:")
        for obs_idx in top_obs:
            print(f"      {labels[obs_idx]} loading={loading[obs_idx]:+.3f}")
        scores = u[:, mode_idx] * singular_values[mode_idx]
        param_lines = top_param_correlations(
            history["all_params"],
            list(history["param_names"]),
            scores,
            top,
        )
        if param_lines:
            print("    strongest parameter correlations:")
            for line in param_lines:
                print(f"      {line}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate and analyze tuner sample-history .npz chunks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python tests/analyze_tuner_sample_history.py output/tuner/job\n"
            "  python tests/analyze_tuner_sample_history.py output/tuner/new_tuner_test "
            "-metric all -modes 5 -top 10 -min_samples 8 --require_param_variation\n"
        ),
    )
    parser.add_argument(
        "path",
        help="Tuner job directory, parent output directory, or one sample_history_*.npz file.",
    )
    parser.add_argument(
        "-metric",
        default="scaled_rmse",
        help="Metric to analyze, or 'all'. Default: scaled_rmse.",
    )
    parser.add_argument("-modes", type=int, default=3, help="Number of SVD modes to print. Default: 3.")
    parser.add_argument("-top", type=int, default=5, help="Top observations/params per mode. Default: 5.")
    parser.add_argument("-min_samples", type=int, default=1, help="Minimum sample rows required. Default: 1.")
    parser.add_argument(
        "--require_param_variation",
        action="store_true",
        help="Fail if every tuned parameter column is constant across the sample history.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        targets = history_targets(Path(args.path))
        print(f"Discovered {len(targets)} sample-history target(s).")
        for target_idx, target in enumerate(targets, start=1):
            if len(targets) > 1:
                print("")
                print(f"=== Sample history {target_idx}/{len(targets)}: {target} ===")
            history = load_history(target)
            print_summary(
                history,
                metric=args.metric,
                modes=max(1, args.modes),
                top=max(1, args.top),
                min_samples=max(1, args.min_samples),
                require_param_variation=bool(args.require_param_variation),
            )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
