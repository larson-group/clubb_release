"""Refresh-hydration tests for Plot-tab state helpers."""

from types import SimpleNamespace

from dash import Patch

from dash_app.plot_tab.callbacks_case import _is_positive_click, _is_same_case
from dash_app.plot_tab.plot_types import shared


def test_dynamic_remove_button_mount_is_not_a_click():
    assert not _is_positive_click(None)
    assert not _is_positive_click(0)
    assert _is_positive_click(1)


def test_same_case_compares_requested_output_directories_before_reload():
    current = {"name": "arm", "output_dirs": ["/tmp/output"]}

    assert _is_same_case(current, "arm", ["/tmp/output"])
    assert not _is_same_case(current, "arm", ["/tmp/other-output"])
    assert not _is_same_case(current, "bomex", ["/tmp/output"])


def test_nominal_height_spacing_uses_native_vertical_levels(monkeypatch):
    dataset = SimpleNamespace(
        var_info={"w_1": {"z_dim": "zt"}, "wp2": {"z_dim": "zt"}}
    )
    collection = SimpleNamespace(datasets=[dataset])
    monkeypatch.setattr(
        shared,
        "get_z_values",
        lambda _dataset, _z_dim: [20.0, 60.0, 100.0, 140.0],
    )

    assert shared.collection_nominal_height_spacing(collection) == 40.0


def test_time_patch_preserves_manual_x_zoom_but_autoscales_unzoomed_view():
    manual_patch = Patch()
    shared.apply_patch_x_range(
        manual_patch,
        [2.0, 8.0],
        {"xaxis.range[0]": 3.0, "xaxis.range[1]": 4.0},
    )
    assert "layout" not in manual_patch.to_plotly_json()

    automatic_patch = Patch()
    shared.apply_patch_x_range(automatic_patch, [2.0, 8.0], None)
    assert automatic_patch.to_plotly_json()["operations"] == [
        {
            "operation": "Assign",
            "location": ["layout", "xaxis", "range"],
            "params": {"value": [2.0, 8.0]},
        }
    ]
