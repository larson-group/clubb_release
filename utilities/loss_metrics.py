"""Profile metrics and loss-table formatting shared by runners and Dash plots.

Public helpers: normalize_loss_weights validates positive mean/std/corr weights;
normalize_named_weights validates field inclusion weights; profile_comparison
returns four profile metrics and the squared Taylor contributions; root_mean_square
and case_rmse aggregate errors. calculate_profile_loss_metrics, calculate_rmse
and calculate_column_loss_metrics provide NumPy convenience results;
print_loss_matrix formats four metrics and the separate case loss.
loss_metric_values_valid checks permitted undefined diagnostics;
loss_result_for_json converts unavailable numeric values to JSON nulls.

The array backend is NumPy by default; JAX callers supply jnp and sqrt_clipped
for zero norm derivatives. Contribution weights are fixed before evaluation.
These functions mirror the raw profile math in `src/clubb_loss_driver.F90`.
They intentionally do not read NetCDF files, choose time windows, interpolate
grids, or rank tuner results; callers own that context and pass in profiles.
"""

from __future__ import annotations

import numpy as np


LOSS_WEIGHT_NAMES = ("mean", "std", "corr")
DEFAULT_LOSS_WEIGHTS = dict.fromkeys(LOSS_WEIGHT_NAMES, 1.0)


def normalize_loss_weights(weights=None) -> dict[str, float]:
    """Validate fixed contribution weights without rescaling their magnitude.

    Missing entries use one, so unit weights reproduce normalized RMSE. Require
    finite strictly positive values for all three contributions.
    """
    if weights is not None and not isinstance(weights, dict):
        raise ValueError("loss_weights must map mean, std and corr to numbers")
    weights = dict(weights or {})
    unknown = set(weights) - set(LOSS_WEIGHT_NAMES)
    if unknown:
        raise ValueError("Unknown loss weight(s): " + ", ".join(sorted(unknown)))
    result = {}
    for name in LOSS_WEIGHT_NAMES:
        try:
            value = float(weights.get(name, 1.0))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} loss weight must be a finite positive number") from exc
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} loss weight must be a finite positive number")
        result[name] = value
    return result


def normalize_named_weights(weights, names, label) -> dict[str, float]:
    """Validate fixed case or field weights, preserving the supplied name order.

    Missing names use one. Zero excludes an item; at least one item must remain.
    The returned complete mapping is suitable for request serialization.
    """
    if weights is not None and not isinstance(weights, dict):
        raise ValueError(f"{label} weights must be a mapping")
    supplied = dict(weights or {})
    if set(supplied) - set(names):
        raise ValueError(f"Unknown {label} weight key(s): " + ", ".join(sorted(set(supplied) - set(names))))
    result = {}
    for name in names:
        try:
            value = float(supplied.get(name, 1.0))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} weight for {name} must be a finite nonnegative number") from exc
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{label} weight for {name} must be a finite nonnegative number")
        result[name] = value
    if not any(result.values()):
        raise ValueError(f"At least one {label} weight must be positive")
    return result


def profile_comparison(model_profile, benchmark_profile, loss_weights=(1.0, 1.0, 1.0),
                       *, array_module=np, sqrt=None):
    """Return four profile metrics and their three squared-error contributions.

    The final array axis contains paired levels; leading axes broadcast. Fixed,
    validated loss_weights are ordered mean/std/corr. NumPy and JAX share this
    arithmetic through array_module. JAX callers supply the existing sqrt_clipped
    kernel to choose zero norm derivatives at exactly flat profiles.

    Bias uses the benchmark standard deviation as its scale, with the existing
    scale-one convention for flat benchmarks. Amplitude and pattern errors use
    actual moments even when the exported correlation/std ratio are undefined
    and reported as NaN. The pattern covariance term retains a
    valid total-error derivative at flat model profiles.
    """
    xp = array_module
    root = xp.sqrt if sqrt is None else sqrt
    model = xp.asarray(model_profile)
    benchmark = xp.asarray(benchmark_profile)
    if model.shape[-1] < 1 or model.shape[-1] != benchmark.shape[-1]:
        raise ValueError("Profile comparisons require nonempty matching level axes")
    model, benchmark = xp.broadcast_arrays(model, benchmark)
    model_mean = xp.mean(model, axis=-1)
    benchmark_mean = xp.mean(benchmark, axis=-1)
    model_centered = model - model_mean[..., None]
    benchmark_centered = benchmark - benchmark_mean[..., None]
    model_variance = xp.mean(model_centered**2, axis=-1)
    benchmark_variance = xp.mean(benchmark_centered**2, axis=-1)
    # A rounded mean can leave spurious variance in exactly constant profiles.
    # Equality detects only that case; small but real variation stays defined.
    model_variance = xp.where(xp.all(model == model[..., :1], axis=-1), 0.0, model_variance)
    benchmark_variance = xp.where(xp.all(benchmark == benchmark[..., :1], axis=-1), 0.0, benchmark_variance)
    model_std = root(model_variance)
    benchmark_std = root(benchmark_variance)
    norm = xp.where(xp.isfinite(benchmark_std) & (benchmark_std > 0.0), benchmark_std, 1.0)
    bias = (model_mean - benchmark_mean) / norm
    centered_error = xp.mean(((model_centered - benchmark_centered) / norm[..., None])**2, axis=-1)
    amplitude_error = ((model_std - benchmark_std) / norm)**2
    covariance = xp.mean((model_centered / norm[..., None]) *
                         (benchmark_centered / norm[..., None]), axis=-1)
    pattern_from_covariance = 2.0 * ((model_std / norm) * (benchmark_std / norm) - covariance)
    # Unit-pattern differences retain small shape errors near a perfect match.
    # At zero model variance, use the covariance expression so its linear
    # derivative survives even though the norm derivative is chosen as zero.
    unit_model = model_centered / xp.where(model_std > 0.0, model_std, 1.0)[..., None]
    unit_benchmark = benchmark_centered / xp.where(benchmark_std > 0.0, benchmark_std, 1.0)[..., None]
    pattern_from_differences = ((model_std / norm) * (benchmark_std / norm)
                                * xp.mean((unit_model - unit_benchmark)**2, axis=-1))
    pattern_error = xp.where((model_std > 0.0) & (benchmark_std > 0.0),
                             pattern_from_differences, pattern_from_covariance)
    # The true pattern contribution is nonnegative. Preserve its covariance
    # derivative at zero; maximum(x, 0) would halve that derivative at a tie.
    pattern_error = xp.where(pattern_error >= 0.0, pattern_error, 0.0)

    mean_weight, std_weight, corr_weight = loss_weights
    if mean_weight == std_weight == corr_weight:
        # This expression avoids cancellation and undefined shape derivatives
        # at flat model profiles while recovering ordinary RMSE at unit weights.
        squared_loss = mean_weight * xp.mean(((model - benchmark) / norm[..., None])**2, axis=-1)
    elif std_weight >= corr_weight:
        squared_loss = (mean_weight * bias**2 + corr_weight * centered_error
                        + (std_weight - corr_weight) * amplitude_error)
    else:
        squared_loss = (mean_weight * bias**2 + std_weight * centered_error
                        + (corr_weight - std_weight) * pattern_error)

    unit_correlation = xp.mean(unit_model * unit_benchmark, axis=-1)
    correlation = xp.where((benchmark_std > 0.0) & (model_std > 0.0),
                            xp.clip(unit_correlation, -1.0, 1.0), xp.nan)
    std_ratio = xp.where(benchmark_std > 0.0, model_std / norm, xp.nan)
    return {
        "metrics": {"rmse": root(squared_loss), "correlation": correlation,
                    "std_ratio": std_ratio, "bias_norm": bias},
        "components": {"mean": bias**2, "std": amplitude_error, "corr": pattern_error},
    }


def root_mean_square(values, weights=None, axis=None, *, array_module=np, sqrt=None):
    """Combine errors by a weighted RMS over the requested axes.

    Weights must be validated nonnegative constants that broadcast to values.
    Scaling before squaring keeps finite invalid-profile penalties finite even
    in float32. Zero-weight items are excluded before arithmetic. JAX callers
    supply sqrt_clipped for the zero-loss derivative convention.
    """
    xp = array_module
    root = xp.sqrt if sqrt is None else sqrt
    values = xp.asarray(values)
    weights = xp.ones_like(values) if weights is None else xp.broadcast_to(xp.asarray(weights), values.shape)
    weight_scale = xp.max(weights, axis=axis, keepdims=True)
    weights = weights / weight_scale
    values = xp.where(weights > 0.0, values, 0.0)
    scale = xp.max(xp.abs(values), axis=axis, keepdims=True)
    safe_scale = xp.where(scale > 0.0, scale, 1.0)
    denominator = xp.sum(weights, axis=axis, keepdims=True)
    squared = xp.sum(weights * (values / safe_scale)**2, axis=axis, keepdims=True) / denominator
    return xp.squeeze(safe_scale * root(squared), axis=axis)


def case_rmse(metric_rmse, field_weights=None, *, array_module=np, sqrt=None):
    """Reduce (window, field, candidate) RMSE values to one loss per candidate."""
    xp = array_module
    values = xp.asarray(metric_rmse)
    weights = None if field_weights is None else xp.asarray(field_weights)[None, :, None]
    return root_mean_square(values, weights, axis=(0, 1), array_module=xp, sqrt=sqrt)


LOSS_METRIC_NAMES = (
    "rmse",
    "correlation",
    "std_ratio",
    "bias_norm",
)


def loss_metric_values_valid(name, values, *, array_module=np):
    """Accept finite values and undefined correlation/std ratios, rejecting infinities.

    Numeric loss outputs use NaN only for undefined diagnostics; RMSE and bias
    must remain finite. This predicate does not change candidate penalties.
    """
    xp = array_module
    values = xp.asarray(values)
    return xp.isfinite(values) | (xp.isnan(values) if name in ('correlation', 'std_ratio') else False)


def loss_result_for_json(value):
    """Convert numeric loss results to strict JSON values; unavailable values become null.

    Keep NaNs in numeric arrays for analysis, and use null at JSON boundaries
    for both undefined diagnostics and unavailable gradients. Input is an array,
    scalar or nested dict/list/tuple; the original object is not mutated.
    """
    if isinstance(value, dict):
        return {key: loss_result_for_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [loss_result_for_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return loss_result_for_json(value.tolist())
    if isinstance(value, np.generic):
        return loss_result_for_json(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def calculate_profile_loss_metrics(model_profile: np.ndarray, benchmark_profile: np.ndarray,
                                   loss_weights=None) -> dict[str, float]:
    """Return the four weighted profile metrics using fixed positive weights."""
    weights = normalize_loss_weights(loss_weights)
    result = profile_comparison(np.asarray(model_profile, dtype=np.float64),
                                np.asarray(benchmark_profile, dtype=np.float64),
                                tuple(weights[name] for name in LOSS_WEIGHT_NAMES))
    return {name: float(value) for name, value in result['metrics'].items()}


def calculate_rmse(model_profile: np.ndarray, benchmark_profile: np.ndarray,
                   loss_weights=None) -> float:
    """Return the weighted Taylor RMSE normalized by benchmark standard deviation."""
    return calculate_profile_loss_metrics(model_profile, benchmark_profile, loss_weights)['rmse']


def calculate_column_loss_metrics(model_matrix: np.ndarray, benchmark_profile: np.ndarray, loss_weights=None) -> dict[str, np.ndarray]:
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
        col_metrics = calculate_profile_loss_metrics(model_arr[:, col_idx], benchmark_arr, loss_weights)
        for metric_name in LOSS_METRIC_NAMES:
            metrics[metric_name][col_idx] = col_metrics[metric_name]

    return metrics


def print_loss_matrix(
    clubb_var_names,
    rmse: np.ndarray,
    correlation: np.ndarray,
    std_ratio: np.ndarray,
    bias_norm: np.ndarray,
    case_loss: np.ndarray,
):
    """Print the loss table in the same row-oriented format as the Fortran tools."""
    header = ["variable"]
    if rmse.shape[0] > 1:
        header.append("window")
    header.extend(f"rmse_col{i + 1}" for i in range(rmse.shape[2]))
    header.extend(f"corr_col{i + 1}" for i in range(rmse.shape[2]))
    header.extend(f"std_ratio_col{i + 1}" for i in range(rmse.shape[2]))
    header.extend(f"bias_norm_col{i + 1}" for i in range(rmse.shape[2]))
    header.append("best_col")
    print(" ".join(header))

    for window_idx in range(rmse.shape[0]):
        for row_idx, var_name in enumerate(clubb_var_names):
            best_col_idx = int(np.argmin(rmse[window_idx, row_idx, :])) + 1
            row = [str(var_name).strip()]
            if rmse.shape[0] > 1:
                row.append(f"window{window_idx + 1}")
            row.extend(f"{value:16.8E}" for value in rmse[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in correlation[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in std_ratio[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in bias_norm[window_idx, row_idx, :])
            row.append(f"col{best_col_idx}")
            print(" ".join(row))
    print("case_loss " + " ".join(f"{value:16.8E}" for value in case_loss))
