"""Shared Taylor diagnostic helpers returned by the loss driver.

The Fortran loss driver supplies profile diagnostics, but the tuner ranking is a
Python policy decision.  This module keeps that policy versioned and explicit so
old runs remain reproducible when formulas or finite fallback penalties change.
"""

from __future__ import annotations

import math

TAYLOR_METRIC_NAMES = (
    "correlation",
    "std_ratio",
    "centered_rmse_norm",
    "bias_norm",
)

LOSS_METRIC_NAMES = (
    "scaled_rmse",
    *TAYLOR_METRIC_NAMES,
)

LOSS_POLICY_VERSION = "taylor_metrics_v1"

# Prefer the shape-sensitive objective for fresh Tune experiments.  Existing
# saved requests retain their explicit loss mode, so this only changes new
# workspaces and requests that omit the field.
DEFAULT_LOSS_MODE = "shape_first"
LOSS_MODE_NAMES = (
    "scaled_rmse",
    "centered_rmse_bias",
    "taylor_components",
    "taylor_components_squared",
    "shape_first",
    "bias_light_taylor",
    "decomposed_taylor",
)

DEFAULT_AGGREGATION_MODE = "quantile_weighted"
AGGREGATION_MODE_NAMES = (
    "quantile_weighted",
    "mean_max",
    "mean_worst_quantile",
)

# Four equally populated bins are formed after sorting active losses from best
# to worst.  The values control the importance of those bins, not the number of
# samples admitted to them.  They are deliberately normalized at use time so a
# user can enter intuitive relative weights in either the CLI or Dash.
DEFAULT_AGGREGATION_WEIGHTS = (0.1, 0.4, 0.4, 0.1)
AGGREGATION_WEIGHT_COUNT = len(DEFAULT_AGGREGATION_WEIGHTS)
TIME_WINDOW_AGGREGATION_SCOPES = ("overall", "by_case")
DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE = "overall"

DEFAULT_NUM_TIME_WINDOWS = 1

MEAN_MAX_MEAN_WEIGHT = 0.6
MEAN_MAX_MAX_WEIGHT = 0.4
WORST_QUANTILE_MEAN_WEIGHT = 0.7
WORST_QUANTILE_WEIGHT = 0.3
WORST_QUANTILE_FRACTION = 0.25

INVALID_CORRELATION_PENALTY = 2.0
INVALID_LOG_STD_RATIO_PENALTY = 5.0
INVALID_SCALED_RMSE_PENALTY = 1.0e30
CORRELATION_ROUNDOFF_TOLERANCE = 1.0e-6

LOSS_POLICY_CONSTANTS = {
    "invalid_correlation_penalty": INVALID_CORRELATION_PENALTY,
    "invalid_log_std_ratio_penalty": INVALID_LOG_STD_RATIO_PENALTY,
    "invalid_scaled_rmse_penalty": INVALID_SCALED_RMSE_PENALTY,
    "correlation_roundoff_tolerance": CORRELATION_ROUNDOFF_TOLERANCE,
    "mean_max_mean_weight": MEAN_MAX_MEAN_WEIGHT,
    "mean_max_max_weight": MEAN_MAX_MAX_WEIGHT,
    "worst_quantile_mean_weight": WORST_QUANTILE_MEAN_WEIGHT,
    "worst_quantile_weight": WORST_QUANTILE_WEIGHT,
    "worst_quantile_fraction": WORST_QUANTILE_FRACTION,
    "default_quantile_weights": list(DEFAULT_AGGREGATION_WEIGHTS),
}


def _finite_float(value: float | int | None) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result):
        return None
    return result


def _finite_or_zero(value: float | int | None) -> float:
    result = _finite_float(value)
    if result is None:
        return 0.0
    return result


def _sanitize_taylor_components(metrics: dict) -> dict:
    """Return finite Taylor components and flags used by all loss modes.

    `centered_rmse_norm` measures profile-shape and amplitude mismatch after the
    mean is removed.  `bias_norm` measures the removed mean offset, so the loss
    can separately distinguish shape errors from bulk vertical shifts.

    Correlation should capture vertical-structure agreement, but degenerate
    profiles can make it undefined.  Small finite roundoff excursions outside
    [-1, 1] are clamped; non-finite or clearly out-of-range values receive a
    named finite penalty and are flagged instead of being silently converted.

    `abs(log(std_ratio))` measures vertical-variability amplitude error.  Log
    space is symmetric: a profile with twice the benchmark variability and one
    with half the benchmark variability have equal-magnitude penalties.  Invalid
    or nonpositive ratios receive a named finite penalty so degenerate Taylor
    metrics stay rankable and visible in diagnostics.
    """
    centered_rmse_raw = _finite_float(metrics.get("centered_rmse_norm"))
    bias_raw = _finite_float(metrics.get("bias_norm"))
    correlation_raw = _finite_float(metrics.get("correlation"))
    std_ratio_raw = _finite_float(metrics.get("std_ratio"))

    centered_rmse_norm = max(0.0, 0.0 if centered_rmse_raw is None else centered_rmse_raw)
    abs_bias_norm = abs(0.0 if bias_raw is None else bias_raw)

    correlation_clamped = False
    correlation_invalid = False
    correlation_out_of_range = False
    if correlation_raw is None:
        sanitized_correlation = 0.0
        correlation_penalty = INVALID_CORRELATION_PENALTY
        correlation_invalid = True
    elif -1.0 - CORRELATION_ROUNDOFF_TOLERANCE <= correlation_raw <= 1.0 + CORRELATION_ROUNDOFF_TOLERANCE:
        sanitized_correlation = max(-1.0, min(1.0, correlation_raw))
        correlation_clamped = sanitized_correlation != correlation_raw
        correlation_penalty = max(0.0, 1.0 - sanitized_correlation)
    else:
        sanitized_correlation = max(-1.0, min(1.0, correlation_raw))
        correlation_penalty = INVALID_CORRELATION_PENALTY
        correlation_invalid = True
        correlation_out_of_range = True

    std_ratio_invalid = False
    std_ratio_nonpositive = False
    if std_ratio_raw is None:
        sanitized_std_ratio = 0.0
        abs_log_std_ratio = INVALID_LOG_STD_RATIO_PENALTY
        std_ratio_invalid = True
    elif std_ratio_raw <= 0.0:
        sanitized_std_ratio = std_ratio_raw
        abs_log_std_ratio = INVALID_LOG_STD_RATIO_PENALTY
        std_ratio_invalid = True
        std_ratio_nonpositive = True
    else:
        sanitized_std_ratio = std_ratio_raw
        abs_log_std_ratio = abs(math.log(std_ratio_raw))

    return {
        "centered_rmse_norm": float(centered_rmse_norm),
        "abs_bias_norm": float(abs_bias_norm),
        "correlation_penalty": float(correlation_penalty),
        "abs_log_std_ratio": float(abs_log_std_ratio),
        "sanitized_correlation": float(sanitized_correlation),
        "sanitized_std_ratio": float(sanitized_std_ratio),
        "centered_rmse_invalid": centered_rmse_raw is None,
        "bias_norm_invalid": bias_raw is None,
        "correlation_invalid": correlation_invalid,
        "correlation_clamped": correlation_clamped,
        "correlation_out_of_range": correlation_out_of_range,
        "std_ratio_invalid": std_ratio_invalid,
        "std_ratio_nonpositive": std_ratio_nonpositive,
        "possible_degenerate_taylor_metrics": (
            centered_rmse_raw is None
            or bias_raw is None
            or correlation_invalid
            or std_ratio_invalid
        ),
    }


def _linear_loss_components(components: dict, weights: dict[str, float]) -> dict[str, float]:
    return {
        component_name: float(weight * components[component_name])
        for component_name, weight in weights.items()
    }


def _squared_loss_components(components: dict, weights: dict[str, float]) -> dict[str, float]:
    return {
        component_name: float(weight * components[component_name] ** 2)
        for component_name, weight in weights.items()
    }


def compute_field_loss_diagnostics(metrics: dict) -> dict:
    """Compute every supported per-field Taylor loss and its contributions.

    The named modes are fixed policy choices:

    `centered_rmse_bias` is the historical Python smart loss.  It rewards low
    centered profile mismatch plus low mean bias.

    `taylor_components` keeps centered RMSE dominant, but adds mild explicit
    penalties for vertical-shape correlation and vertical-variability amplitude.

    `taylor_components_squared` uses the same physical components, squared, so
    one very bad component is rejected more aggressively.

    `shape_first` intentionally gives profile shape and variability more direct
    influence when RMS-like objectives improve magnitude but leave physically
    wrong vertical structure.

    `bias_light_taylor` reduces the mean-bias term and prioritizes centered
    vertical structure while still keeping bias constrained.

    `decomposed_taylor` uses only the decomposed Taylor components plus bias:
    abs(bias_norm) + (1 - R) + abs(log(std_ratio)).  It intentionally removes
    centered RMSE because centered RMSE is determined by correlation and
    std_ratio.
    """
    components = _sanitize_taylor_components(metrics)
    scaled_rmse_raw = _finite_float(metrics.get("scaled_rmse"))
    scaled_rmse = INVALID_SCALED_RMSE_PENALTY if scaled_rmse_raw is None else max(0.0, scaled_rmse_raw)
    components["scaled_rmse"] = float(scaled_rmse)
    components["scaled_rmse_invalid"] = scaled_rmse_raw is None
    component_weights_by_mode = {
        "scaled_rmse": {
            "scaled_rmse": 1.0,
        },
        "centered_rmse_bias": {
            "centered_rmse_norm": 1.0,
            "abs_bias_norm": 1.0,
        },
        "taylor_components": {
            "centered_rmse_norm": 1.0,
            "abs_bias_norm": 0.5,
            "correlation_penalty": 0.25,
            "abs_log_std_ratio": 0.25,
        },
        "taylor_components_squared": {
            "centered_rmse_norm": 1.0,
            "abs_bias_norm": 0.5,
            "correlation_penalty": 0.25,
            "abs_log_std_ratio": 0.25,
        },
        "shape_first": {
            "centered_rmse_norm": 0.5,
            "abs_bias_norm": 0.5,
            "correlation_penalty": 1.0,
            "abs_log_std_ratio": 0.5,
        },
        "bias_light_taylor": {
            "centered_rmse_norm": 1.0,
            "abs_bias_norm": 0.25,
            "correlation_penalty": 0.5,
            "abs_log_std_ratio": 0.5,
        },
        "decomposed_taylor": {
            "abs_bias_norm": 1.0,
            "correlation_penalty": 1.0,
            "abs_log_std_ratio": 1.0,
        },
    }

    loss_components_by_mode = {}
    per_field_losses = {}
    for mode_name in LOSS_MODE_NAMES:
        weights = component_weights_by_mode[mode_name]
        if mode_name == "taylor_components_squared":
            contributions = _squared_loss_components(components, weights)
        else:
            contributions = _linear_loss_components(components, weights)
        loss_components_by_mode[mode_name] = contributions
        per_field_losses[mode_name] = float(sum(contributions.values()))

    return {
        **components,
        "per_field_losses": per_field_losses,
        "loss_components_by_mode": loss_components_by_mode,
    }


def field_smart_loss(metrics: dict, loss_mode: str = DEFAULT_LOSS_MODE) -> float:
    """Return the selected per-field smart loss.

    This wrapper preserves the original call pattern while allowing callers to
    select one of the versioned Python loss modes.
    """
    return float(compute_field_loss_diagnostics(metrics)["per_field_losses"][loss_mode])


def _active_loss_pairs(losses: list[float], weights: list[float] | None) -> tuple[list[tuple[float, float]], int]:
    if weights is None:
        weights = [1.0] * len(losses)

    active_pairs = []
    excluded_count = 0
    for loss, weight in zip(losses, weights):
        loss_value = _finite_or_zero(loss)
        weight_value = _finite_or_zero(weight)
        if weight_value <= 0.0:
            excluded_count += 1
            continue
        active_pairs.append((loss_value, weight_value))
    excluded_count += max(0, len(losses) - len(weights))
    return active_pairs, excluded_count


def normalize_aggregation_weights(values=None) -> list[float]:
    """Return four finite, nonnegative quantile weights with unit sum.

    The public request keeps the user's relative values for provenance, but the
    numerical aggregation always uses this normalized form.  Requiring exactly
    four entries avoids a silent change in the meaning of the best-to-worst
    quartile controls.
    """
    raw_values = DEFAULT_AGGREGATION_WEIGHTS if values is None else values
    if isinstance(raw_values, str):
        raw_values = [item.strip() for item in raw_values.split(",") if item.strip()]
    try:
        weights = [float(value) for value in raw_values]
    except (TypeError, ValueError) as exc:
        raise ValueError("aggregation_weights must contain four finite nonnegative numbers") from exc
    if len(weights) != AGGREGATION_WEIGHT_COUNT:
        raise ValueError(f"aggregation_weights must contain exactly {AGGREGATION_WEIGHT_COUNT} numbers")
    if any(not math.isfinite(value) or value < 0.0 for value in weights):
        raise ValueError("aggregation_weights must contain four finite nonnegative numbers")
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("aggregation_weights must contain at least one positive number")
    return [float(value / total) for value in weights]


def aggregate_losses(
    losses: list[float],
    weights: list[float] | None = None,
    aggregation_mode: str = DEFAULT_AGGREGATION_MODE,
    aggregation_weights=None,
) -> dict:
    """Aggregate active case-field losses and return diagnostics.

    `quantile_weighted` sorts active losses best-to-worst, partitions them into
    four equally populated bins, and combines their weighted means using the
    supplied best-to-worst quantile weights.  Empty bins are omitted and the
    remaining weights are renormalized, so one-window requests retain their
    actual loss rather than being artificially scaled.

    `mean_max` is the historical blend: weighted mean rewards broad improvement,
    while the hard maximum rejects a parameter set that is very bad for one
    active case-field pair.

    `mean_worst_quantile` softens the hard maximum by selecting the worst 25% of
    active items by raw loss value, always including at least one item, and then
    computing the weighted mean over that subset.
    """
    active_pairs, excluded_count = _active_loss_pairs(losses, weights)
    if aggregation_mode not in AGGREGATION_MODE_NAMES:
        raise ValueError(f"Unknown aggregation mode: {aggregation_mode}")

    base = {
        "aggregation_mode": aggregation_mode,
        "active_item_count": len(active_pairs),
        "excluded_item_count": excluded_count,
        "weighted_mean": 0.0,
    }
    if not active_pairs:
        if aggregation_mode == "quantile_weighted":
            return {
                **base,
                "quantile_weights": normalize_aggregation_weights(aggregation_weights),
                "quantile_bin_means": [],
                "loss": 0.0,
            }
        if aggregation_mode == "mean_worst_quantile":
            return {**base, "worst_quantile_mean": 0.0, "loss": 0.0}
        return {**base, "max_loss": 0.0, "loss": 0.0}

    weighted_sum = sum(loss * weight for loss, weight in active_pairs)
    weight_sum = sum(weight for _loss, weight in active_pairs)
    weighted_mean = weighted_sum / weight_sum

    if aggregation_mode == "quantile_weighted":
        quantile_weights = normalize_aggregation_weights(aggregation_weights)
        sorted_pairs = sorted(active_pairs, key=lambda pair: pair[0])
        bins = [[] for _ in range(AGGREGATION_WEIGHT_COUNT)]
        item_count = len(sorted_pairs)
        for item_idx, pair in enumerate(sorted_pairs):
            # Integer boundaries keep the bins contiguous and differ by at
            # most one item, even when the number of time windows is not four.
            bin_idx = min(AGGREGATION_WEIGHT_COUNT - 1, item_idx * AGGREGATION_WEIGHT_COUNT // item_count)
            bins[bin_idx].append(pair)
        bin_diagnostics = []
        for bin_idx, bin_pairs in enumerate(bins):
            if not bin_pairs:
                continue
            bin_weight_sum = sum(weight for _loss, weight in bin_pairs)
            bin_mean = sum(loss * weight for loss, weight in bin_pairs) / bin_weight_sum
            bin_diagnostics.append(
                {
                    "bin_index": bin_idx,
                    "label": ("best", "lower_middle", "upper_middle", "worst")[bin_idx],
                    "item_count": len(bin_pairs),
                    "mean": float(bin_mean),
                    "requested_weight": float(quantile_weights[bin_idx]),
                }
            )
        active_quantile_weight = sum(item["requested_weight"] for item in bin_diagnostics)
        loss = sum(item["mean"] * item["requested_weight"] for item in bin_diagnostics) / active_quantile_weight
        return {
            **base,
            "weighted_mean": float(weighted_mean),
            "quantile_weights": quantile_weights,
            "quantile_bin_means": bin_diagnostics,
            "active_quantile_weight": float(active_quantile_weight),
            "loss": float(loss),
        }

    if aggregation_mode == "mean_worst_quantile":
        worst_count = max(1, math.ceil(len(active_pairs) * WORST_QUANTILE_FRACTION))
        worst_pairs = sorted(active_pairs, key=lambda pair: pair[0], reverse=True)[:worst_count]
        worst_weighted_sum = sum(loss * weight for loss, weight in worst_pairs)
        worst_weight_sum = sum(weight for _loss, weight in worst_pairs)
        worst_quantile_mean = worst_weighted_sum / worst_weight_sum
        loss = WORST_QUANTILE_MEAN_WEIGHT * weighted_mean + WORST_QUANTILE_WEIGHT * worst_quantile_mean
        return {
            **base,
            "weighted_mean": float(weighted_mean),
            "worst_quantile_mean": float(worst_quantile_mean),
            "worst_quantile_fraction": WORST_QUANTILE_FRACTION,
            "worst_quantile_item_count": worst_count,
            "loss": float(loss),
        }

    max_loss = max(loss for loss, _weight in active_pairs)
    loss = MEAN_MAX_MEAN_WEIGHT * weighted_mean + MEAN_MAX_MAX_WEIGHT * max_loss
    return {
        **base,
        "weighted_mean": float(weighted_mean),
        "max_loss": float(max_loss),
        "loss": float(loss),
    }


def blended_mean_max_loss(losses: list[float], weights: list[float] | None = None) -> float:
    """Backward-compatible wrapper for the historical mean/max aggregation."""
    return float(aggregate_losses(losses, weights, "mean_max")["loss"])
