"""Refresh-hydration tests for Plot-tab state helpers."""

from types import SimpleNamespace

from dash import Dash, Patch, no_update

from dash_app.plot_tab.callbacks_case import (
    _catalog_tracking_paths,
    _delete_output_directory,
    _is_positive_click,
    _is_same_case,
    _normalize_output_dirs,
    _update_output_dirs,
    register_case_callbacks,
)
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


def test_selection_normalization_preserves_temporarily_missing_directories():
    assert _normalize_output_dirs(["/tmp/missing", "/tmp/other", "/tmp/missing"]) == [
        "/tmp/missing",
        "/tmp/other",
    ]


def test_adding_or_removing_one_output_preserves_every_other_selection():
    current = ["/tmp/one", "/tmp/two"]

    assert _update_output_dirs(current, "add", "/tmp/three") == [
        "/tmp/one",
        "/tmp/two",
        "/tmp/three",
    ]
    assert _update_output_dirs(current, "remove", "/tmp/one") == ["/tmp/two"]


def test_output_changes_commit_only_after_picker_closes():
    app = Dash(__name__)
    register_case_callbacks(app)
    entry = next(
        value
        for key, value in app.callback_map.items()
        if "plots-loaded-output-dirs.data" in key
    )
    commit = entry["callback"].__wrapped__

    loaded, warning = commit(True, ["/tmp/one", "/tmp/two"], ["/tmp/one"])
    assert loaded is no_update
    assert warning == "Close dropdown to load changes"

    loaded, warning = commit(False, ["/tmp/one", "/tmp/two"], ["/tmp/one"])
    assert loaded == ["/tmp/one", "/tmp/two"]
    assert warning == ""


def test_case_loading_consumes_committed_not_draft_outputs():
    app = Dash(__name__)
    register_case_callbacks(app)
    consumers = [
        entry
        for entry in app.callback_map.values()
        if "plots-loaded-output-dirs" in {item["id"] for item in entry["inputs"]}
    ]

    assert len(consumers) == 3
    assert all(
        "plots-output-dirs" not in {item["id"] for item in entry["inputs"]}
        for entry in consumers
    )


def test_output_catalog_refresh_is_tab_scoped_and_ten_second_driven():
    app = Dash(__name__)
    register_case_callbacks(app)

    refresh = app.callback_map["plots-output-catalog.data"]
    assert [(item["id"], item["property"]) for item in refresh["inputs"]] == [
        ("dashboard-tabs", "value"),
        ("plots-output-refresh-interval", "n_intervals"),
    ]
    assert [(item["id"], item["property"]) for item in refresh["state"]] == [
        ("plots-output-dirs", "data"),
        ("plots-output-catalog", "data"),
    ]


def test_catalog_is_only_state_for_the_callback_that_mutates_active_outputs():
    app = Dash(__name__)
    register_case_callbacks(app)
    selection = next(
        entry
        for entry in app.callback_map.values()
        if '{"path":["ALL"],"type":"plots-add-output-dir"}' in {
            item["id"] for item in entry["inputs"]
        }
    )

    assert "plots-output-catalog" not in {item["id"] for item in selection["inputs"]}
    assert ("plots-output-catalog", "data") in {
        (item["id"], item["property"]) for item in selection["state"]
    }
    click_properties = {item["property"] for item in selection["inputs"]}
    assert click_properties == {"n_clicks_timestamp"}


def test_output_delete_is_limited_to_subdirectories_under_output(monkeypatch, tmp_path):
    from dash_app.plot_tab import callbacks_case

    output_root = tmp_path / "output"
    target = output_root / "run"
    target.mkdir(parents=True)
    (target / "arm_stats.nc").write_bytes(b"CDF")
    monkeypatch.setattr(callbacks_case, "OUTPUT_ROOT", output_root)

    _delete_output_directory(target)

    assert not target.exists()
    outside = tmp_path / "outside"
    outside.mkdir()
    try:
        _delete_output_directory(outside)
    except ValueError as exc:
        assert "inside output/" in str(exc)
    else:
        raise AssertionError("outside directory was accepted for deletion")
    assert outside.exists()


def test_catalog_refresh_only_drives_output_selection_ui():
    app = Dash(__name__)
    register_case_callbacks(app)
    consumers = [
        key
        for key, entry in app.callback_map.items()
        if "plots-output-catalog" in {item["id"] for item in entry["inputs"]}
    ]

    assert consumers == [
        "..plots-available-output-list.children...plots-active-output-list.children...plots-output-menu.className.."
    ]


def test_removed_external_output_remains_tracked_for_the_available_menu():
    catalog = [
        {"path": "/repo/output/run", "catalog_origin": "output"},
        {"path": "/tmp/external-run", "catalog_origin": "external"},
    ]

    assert _catalog_tracking_paths(catalog, ["/repo/output/run"]) == [
        "/repo/output/run",
        "/tmp/external-run",
    ]


def test_case_scan_uses_union_and_keeps_matching_output_order(tmp_path, monkeypatch):
    first = str((tmp_path / "first").resolve())
    second = str((tmp_path / "second").resolve())
    case_maps = {
        first: {"arm": f"{first}/arm_stats.nc"},
        second: {
            "arm": f"{second}/arm_stats.nc",
            "rico": f"{second}/rico_stats.nc",
        },
    }
    monkeypatch.setattr(shared, "_scan_cases_in_directory", case_maps.__getitem__)

    cases = shared.scan_output_cases([first, second])

    assert cases == {
        "arm": [f"{first}/arm_stats.nc", f"{second}/arm_stats.nc"],
        "rico": [f"{second}/rico_stats.nc"],
    }


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
