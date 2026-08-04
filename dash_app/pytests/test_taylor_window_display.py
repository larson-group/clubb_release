"""Regression coverage for the Taylor diagnostic subwindow display selector."""

import math

from dash_app.tune_tab.callbacks_display import (
    _collect_taylor_points,
    _parameter_group_specs,
    _window_display_metrics,
)


def _metrics():
    return {
        "loss": 2.4,
        "scaled_rmse": 1.5,
        "correlation": -0.4,
        "std_ratio": 0.4,
        "centered_rmse_norm": 1.3,
        "bias_norm": -0.2,
        "subwindows": [
            {
                "window_index": 1,
                "loss": 3.0,
                "scaled_rmse": 2.0,
                "correlation": -0.4,
                "std_ratio": 0.4,
                "centered_rmse_norm": 1.3,
                "bias_norm": -0.2,
            },
            {
                "window_index": 2,
                "loss": 1.0,
                "scaled_rmse": 0.5,
                "correlation": 0.8,
                "std_ratio": 0.9,
                "centered_rmse_norm": 0.4,
                "bias_norm": 0.1,
            },
            {
                "window_index": 3,
                "loss": 2.0,
                "scaled_rmse": 1.0,
                "correlation": 0.4,
                "std_ratio": 0.7,
                "centered_rmse_norm": 0.8,
                "bias_norm": 0.0,
            },
        ],
    }


def test_taylor_display_selects_best_and_worst_by_subwindow_smart_loss():
    worst, worst_info = _window_display_metrics(_metrics(), "worst")
    best, best_info = _window_display_metrics(_metrics(), "best")

    assert worst_info["window_index"] == 1
    assert worst["correlation"] == -0.4
    assert worst["aggregate_field_loss"] == 2.4
    assert best_info["window_index"] == 2
    assert best["correlation"] == 0.8
    assert best["aggregate_field_loss"] == 2.4


def test_taylor_display_average_means_raw_metrics_but_preserves_aggregate_loss_in_points():
    average, info = _window_display_metrics(_metrics(), "average")
    assert info["window_index"] is None
    assert info["window_count"] == 3
    assert math.isclose(average["correlation"], (0.8 + 0.4 - 0.4) / 3.0)

    results = [{"rank": 1, "total_loss": 3.2, "field_metrics": {"arm": {"cloud_frac": _metrics()}}}]
    point = _collect_taylor_points(results, "average")[0]
    assert math.isclose(point["correlation"], average["correlation"])
    assert math.isclose(point["field_loss"], 2.4)
    assert point["window_index"] is None


def test_parameter_box_selector_is_limited_to_aggregate_or_all_cases():
    aggregate = _parameter_group_specs([{"rank": 1}], {"arm": [{"rank": 1}], "bomex": [{"rank": 1}]}, "aggregate")
    all_groups = _parameter_group_specs([{"rank": 1}], {"arm": [{"rank": 1}], "bomex": [{"rank": 1}]}, "all")

    assert [item[0] for item in aggregate] == ["aggregate"]
    assert [item[0] for item in all_groups] == ["aggregate", "case:arm", "case:bomex"]
