"""Combine prepared case losses using the tuner's cross-case policy.

Per-profile error construction and field/window RMS live below the tuner in the
loss driver. aggregate_losses() retains the existing case-ranking choices;
normalize_aggregation_weights() validates their four best-to-worst bin weights.
blended_mean_max_loss() forwards to the mean/max cross-case policy.
This module does not reconstruct a case objective from profile diagnostics.
"""

from __future__ import annotations

import math

from utilities.loss_metrics import LOSS_METRIC_NAMES, DEFAULT_LOSS_WEIGHTS

LOSS_POLICY_VERSION = "weighted_taylor_rmse_v1"

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

DEFAULT_NUM_TIME_WINDOWS = 1

MEAN_MAX_MEAN_WEIGHT = 0.6
MEAN_MAX_MAX_WEIGHT = 0.4
WORST_QUANTILE_MEAN_WEIGHT = 0.7
WORST_QUANTILE_WEIGHT = 0.3
WORST_QUANTILE_FRACTION = 0.25

INVALID_RMSE_PENALTY = 1.0e30
LOSS_POLICY_CONSTANTS = {
    "invalid_rmse_penalty": INVALID_RMSE_PENALTY,
    "normalization": "benchmark_standard_deviation",
    "flat_reference_scale": 1.0,
    "case_aggregation": "rms_fields_and_windows",
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
    """Aggregate active prepared case losses and return diagnostics.

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
