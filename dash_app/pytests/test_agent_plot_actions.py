"""Focused checks for agent-controlled Plot-tab windows and PNG export."""

from pathlib import Path

import plotly.graph_objects as go
import pytest

from dash_app.shared import actions, activity
from dash_app.services.jobs import JobStore
from dash_app.plot_tab.benchmark_overlay import resolve_enabled_sources, toggle_enabled_source


def _case_data():
    return {
        "profile_vars": [{"label": "wp2", "value": "wp2"}],
        "default_height_range": [0.0, 1000.0],
        "time_slider_start_min_seconds": 0.0,
        "default_time_start_seconds": 1200.0,
        "time_slider_duration_min_minutes": 10.0,
        "time_slider_duration_max_minutes": 120.0,
        "benchmarks": {"available_sources": ["sam", "coamps"]},
    }


def test_plot_request_carries_custom_physical_controls(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "_profile_selection", lambda *_args, **_kwargs: (_case_data(), {"time_start_seconds": 1800.0, "average_minutes": 30.0}))

    result = actions.plot_profiles("arm", ["wp2"], time_start_seconds=1800, average_minutes=30)

    assert result["time_start_seconds"] == 1800.0
    request = activity.read_activity()["plot_request"]
    assert request["time_start_seconds"] == 1800.0
    assert request["average_minutes"] == 30.0


def test_plot_request_can_select_one_immutable_scm_run(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    store = JobStore(Path(tmp_path) / "jobs.json")
    record, _created = store.submit(
        "scm",
        "plot-run-selection-123",
        {"case": "bomex"},
    )
    output_directory = str(Path(tmp_path) / "isolated-output")
    store.update(record["job_id"], output_directory=output_directory)
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    captured = {}

    def fake_selection(*_args, **kwargs):
        captured.update(kwargs)
        return _case_data(), {
            "time_start_seconds": 0.0,
            "average_minutes": 1.0,
        }

    monkeypatch.setattr(actions, "_profile_selection", fake_selection)

    result = actions.plot_profiles(
        "bomex",
        ["wp2"],
        run_id=record["run_id"],
    )

    assert captured["output_dirs"] == [output_directory]
    assert result["run_id"] == record["run_id"]
    assert activity.read_activity()["plot_request"]["output_dir"] == output_directory


def test_plot_request_can_select_sam_overlay_for_immutable_run(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(
        actions,
        "_profile_selection",
        lambda *_args, **_kwargs: (_case_data(), {"time_start_seconds": 0.0, "average_minutes": 1.0}),
    )

    result = actions.plot_profiles("arm", ["wp2"], benchmark_sources=["sam"])

    assert result["benchmark_sources"] == ["sam"]
    assert activity.read_activity()["plot_request"]["benchmark_sources"] == ["sam"]


def test_plot_request_rejects_unavailable_benchmark_source(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(
        actions,
        "_profile_selection",
        lambda *_args, **_kwargs: (_case_data(), {}),
    )

    with pytest.raises(ValueError, match="unavailable"):
        actions.plot_profiles("arm", ["wp2"], benchmark_sources=["wrf"])


def test_ui_and_mcp_benchmark_selection_share_the_same_source_transition():
    case_data = {"benchmarks": {"available_sources": ["sam", "coamps"]}}

    assert resolve_enabled_sources(case_data, ["sam"]) == ["sam"]
    assert toggle_enabled_source(case_data, [], "sam") == ["sam"]
    assert toggle_enabled_source(case_data, ["sam"], "sam") == []


def test_png_export_uses_the_profile_builder_and_reports_saved_files(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "_profile_selection", lambda *_args, **_kwargs: (_case_data(), {"window_preset": "loss", "time_start_seconds": 1200.0, "average_minutes": 30.0}))
    monkeypatch.setattr(actions, "_PROFILE_EXPORT_DIR", Path(tmp_path) / "exports")
    captured = {}
    monkeypatch.setattr(actions.profile_plot, "build_figure", lambda state, context: captured.update({"state": state, "context": context}) or object())

    def fake_write(_figure, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"png")

    monkeypatch.setattr(actions, "_write_profile_figure_png", fake_write)

    result = actions.save_profile_png("arm", ["wp2"], window_preset="loss")

    assert captured["state"] == {"var": "wp2", "size": "normal"}
    assert captured["context"]["time_range"] == 30.0
    assert result["status"] == "saved"
    assert Path(result["paths"][0]).is_file()


def test_profile_png_renderer_writes_a_real_png(tmp_path):
    figure = go.Figure(go.Scatter(x=[0.0, 1.0], y=[0.0, 1000.0], name="CLUBB"))
    figure.update_layout(title="wp2", xaxis_title="wp2", yaxis_title="Height")
    path = Path(tmp_path) / "profile.png"

    actions._write_profile_figure_png(figure, path)

    assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
