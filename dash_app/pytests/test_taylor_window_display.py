"""Regression coverage for the Taylor diagnostic subwindow display selector."""

import math

from dash_app.tune_tab.callbacks_display import (
    _parameter_box_context,
    _collect_taylor_points,
    _parameter_group_specs,
    _window_display_metrics,
    build_parameter_box_figure,
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


def test_parameter_box_range_normalizes_each_parameter_and_marks_config_defaults():
    request = {
        "config": "experiment",
        "parameter_ranges": [
            {"name": "C1", "targets": ["C1"], "min": 1.0, "max": 3.0},
            {"name": "C8", "targets": ["C8"], "min": 0.0, "max": 10.0},
        ],
        "settings_resolution": {
            "tunable_default_ranges": {
                "C1": {"default": 2.0},
                "C8": {"default": 2.5},
            }
        },
    }
    results = [
        {"rank": 1, "total_loss": 1.0, "params": {"C1": 1.5, "C8": 5.0}},
        {"rank": 2, "total_loss": 2.0, "params": {"C1": 2.5, "C8": 10.0}},
    ]

    figure = build_parameter_box_figure(
        results,
        selected_groups="aggregate",
        scale_mode="normalized",
        parameter_context=_parameter_box_context(request),
    )

    boxes = {trace.x[0]: list(trace.y) for trace in figure.data if trace.type == "box"}
    assert boxes == {"C1": [0.25, 0.75], "C8": [0.5, 1.0]}
    defaults = next(trace for trace in figure.data if trace.name == "Config default")
    assert list(defaults.x) == ["C1", "C8"]
    assert list(defaults.y) == [0.5, 0.25]
    assert defaults.marker.symbol == "diamond"
    assert figure.layout.yaxis.tickformat == ".0%"


def test_parameter_box_unscaled_mode_preserves_raw_values_without_default_markers():
    results = [{"rank": 1, "total_loss": 1.0, "params": {"C1": 2.5}}]
    figure = build_parameter_box_figure(
        results,
        selected_groups="aggregate",
        scale_mode="unscaled",
        parameter_context={"ranges": {"C1": (1.0, 3.0)}, "defaults": {"C1": 2.0}},
    )

    box = next(trace for trace in figure.data if trace.type == "box")
    assert list(box.y) == [2.5]
    assert all(trace.name != "Config default" for trace in figure.data)
    assert figure.layout.yaxis.title.text == "parameter value"
