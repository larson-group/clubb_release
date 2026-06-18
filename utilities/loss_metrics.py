"""Profile metric helpers shared by loss-driver tests and Dash plots.

These functions mirror the raw profile math in `src/clubb_loss_driver.F90`.
They intentionally do not read NetCDF files, choose time windows, interpolate
grids, or rank tuner results; callers own that context and pass in profiles.
"""

from __future__ import annotations

import numpy as np


LOSS_METRIC_NAMES = (
    "scaled_rmse",
    "correlation",
    "std_ratio",
    "centered_rmse_norm",
    "bias_norm",
)


def benchmark_range_norm(benchmark_profile: np.ndarray, fallback: float = 1.0) -> float:
    """Return the loss-driver range normalization for one benchmark profile."""
    benchmark_arr = np.asarray(benchmark_profile, dtype=np.float64)
    norm = float(np.max(benchmark_arr) - np.min(benchmark_arr))
    if not np.isfinite(norm) or norm <= 0.0:
        return float(fallback)
    return norm


def calculate_taylor_metrics(model_profile: np.ndarray, benchmark_profile: np.ndarray) -> dict[str, float]:
    """Mirror `clubb_loss_driver.F90` Taylor metrics for one profile pair."""
    model_arr = np.asarray(model_profile, dtype=np.float64)
    benchmark_arr = np.asarray(benchmark_profile, dtype=np.float64)
    num_levels = model_arr.size
    if num_levels <= 0 or num_levels != benchmark_arr.size:
        raise RuntimeError("Taylor metric profiles must be non-empty and matching")

    model_mean = float(np.sum(model_arr) / num_levels)
    benchmark_mean = float(np.sum(benchmark_arr) / num_levels)
    model_centered = model_arr - model_mean
    benchmark_centered = benchmark_arr - benchmark_mean

    model_centered_sumsq = float(np.sum(model_centered**2))
    benchmark_centered_sumsq = float(np.sum(benchmark_centered**2))
    covariance_sum = float(np.sum(model_centered * benchmark_centered))
    centered_diff_sumsq = float(np.sum((model_centered - benchmark_centered) ** 2))

    model_stddev = float(np.sqrt(model_centered_sumsq / num_levels))
    benchmark_stddev = float(np.sqrt(benchmark_centered_sumsq / num_levels))
    centered_rmse = float(np.sqrt(centered_diff_sumsq / num_levels))
    bias = model_mean - benchmark_mean

    if model_stddev > 0.0 and benchmark_stddev > 0.0:
        correlation = covariance_sum / np.sqrt(model_centered_sumsq * benchmark_centered_sumsq)
        correlation = float(max(-1.0, min(1.0, correlation)))
    else:
        correlation = 0.0

    if benchmark_stddev > 0.0:
        std_ratio = model_stddev / benchmark_stddev
        centered_rmse_norm = centered_rmse / benchmark_stddev
        bias_norm = bias / benchmark_stddev
    else:
        std_ratio = 0.0
        centered_rmse_norm = centered_rmse
        bias_norm = bias

    return {
        "correlation": float(correlation),
        "std_ratio": float(std_ratio),
        "centered_rmse_norm": float(centered_rmse_norm),
        "bias_norm": float(bias_norm),
    }


def calculate_scaled_rmse(
    model_profile: np.ndarray,
    benchmark_profile: np.ndarray,
    *,
    norm: float | None = None,
) -> float:
    """Return the loss-driver scaled squared profile error."""
    model_arr = np.asarray(model_profile, dtype=np.float64)
    benchmark_arr = np.asarray(benchmark_profile, dtype=np.float64)
    if model_arr.size <= 0 or model_arr.size != benchmark_arr.size:
        raise RuntimeError("Loss metric profiles must be non-empty and matching")
    norm_value = benchmark_range_norm(benchmark_arr) if norm is None else float(norm)
    diff = model_arr - benchmark_arr
    return float(np.sum((diff / norm_value) ** 2))


def calculate_profile_loss_metrics(model_profile: np.ndarray, benchmark_profile: np.ndarray) -> dict[str, float]:
    """Return all raw loss-driver metrics for one profile pair."""
    metrics = {
        "scaled_rmse": calculate_scaled_rmse(model_profile, benchmark_profile),
    }
    metrics.update(calculate_taylor_metrics(model_profile, benchmark_profile))
    return metrics


def calculate_column_loss_metrics(model_matrix: np.ndarray, benchmark_profile: np.ndarray) -> dict[str, np.ndarray]:
    """Return all raw loss-driver metrics for each parameter column."""
    model_arr = np.asarray(model_matrix, dtype=np.float64)
    benchmark_arr = np.asarray(benchmark_profile, dtype=np.float64)
    if model_arr.ndim != 2:
        raise RuntimeError("model_matrix must be two-dimensional")
    if model_arr.shape[0] <= 0 or model_arr.shape[0] != benchmark_arr.size:
        raise RuntimeError("Model matrix levels must be non-empty and match the benchmark profile")

    num_cols = model_arr.shape[1]
    metrics = {
        name: np.zeros(num_cols, dtype=np.float64)
        for name in LOSS_METRIC_NAMES
    }

    for col_idx in range(num_cols):
        col_metrics = calculate_profile_loss_metrics(model_arr[:, col_idx], benchmark_arr)
        for metric_name in LOSS_METRIC_NAMES:
            metrics[metric_name][col_idx] = col_metrics[metric_name]

    return metrics
