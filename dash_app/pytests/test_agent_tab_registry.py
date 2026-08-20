"""Tests for the generic tab-operation interface used by every agent adapter."""

from pathlib import Path

from dash_app.shared import actions, activity


def _isolated_activity(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()


def _live_dashboard(monkeypatch):
    monkeypatch.setattr(
        actions.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "available", "dashboard_id": "test-dashboard"},
    )


def test_manifest_covers_every_top_level_dashboard_tab():
    manifest = actions.inspect_dashboard()
    tabs = {item["tab"]: item for item in manifest["tabs"]}

    assert set(tabs) == {
        "tutorial",
        "compile",
        "run",
        "profile",
        "tune",
        "plots",
        "reports",
        "misc",
    }
    assert all(any(operation["name"] == "navigate" for operation in item["operations"]) for item in tabs.values())
    assert {operation["name"] for operation in tabs["plots"]["operations"]} >= {"set_view", "list", "remove"}
    assert {operation["name"] for operation in tabs["tune"]["operations"]} == {"navigate"}
    assert {operation["name"] for operation in tabs["reports"]["operations"]} >= {"open_report"}


def test_generic_invoke_routes_plot_state_to_existing_tab_handler(monkeypatch):
    _live_dashboard(monkeypatch)
    captured = {}

    def fake_plot(case, variables, **kwargs):
        captured.update({"case": case, "variables": variables, **kwargs})
        return {"status": "requested"}

    monkeypatch.setattr(actions, "plot_profiles", fake_plot)

    result = actions.invoke_dashboard(
        "plots",
        "set_view",
        {"case": "arm", "variables": ["wp2"], "output_dir": "dash_default", "window_preset": "loss", "benchmark_sources": ["sam"]},
    )

    assert result == {"tab": "plots", "operation": "set_view", "status": "requested"}
    assert captured == {
        "case": "arm",
        "variables": ["wp2"],
        "run_id": None,
        "output_dir": "dash_default",
        "output_dirs": None,
        "time_start_seconds": None,
        "average_minutes": None,
        "window_preset": "loss",
        "benchmark_sources": ["sam"],
    }


def test_generic_invoke_rejects_scientific_run_actions(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    captured = {}

    def fake_start(case, stats, overrides, cli_options, config):
        captured.update(
            {
                "case": case,
                "stats": stats,
                "overrides": overrides,
                "cli_options": cli_options,
                "config": config,
            }
        )
        return {"pid": 123, "log": "/tmp/arm.log", "start_time": 1.0}

    monkeypatch.setattr(actions, "start_case_process", fake_start)
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    import pytest

    with pytest.raises(ValueError, match="typed MCP domain operation"):
        actions.invoke_dashboard(
            "run",
            "start",
            {
                "case": "arm",
                "overrides": "iiPDF_type=1 l_predict_upwp_vpwp=false",
            },
        )
    assert captured == {}


def test_browser_handoff_returns_structured_unavailable_without_live_dashboard(monkeypatch):
    monkeypatch.setattr(
        actions.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "unavailable", "reason": "test"},
    )

    result = actions.invoke_dashboard("run", "navigate")

    assert result["error"]["code"] == "DASHBOARD_UNAVAILABLE"
    assert result["dashboard"]["status"] == "unavailable"


def test_tutorial_operation_emits_one_unified_browser_request(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    _live_dashboard(monkeypatch)

    tutorial = actions.invoke_dashboard("tutorial", "open_lesson", {"lesson": "tutorial-adg1-explorer"})
    tutorial_request = activity.read_activity()["ui_request"]
    assert tutorial["lesson"] == "tutorial-adg1-explorer"
    assert tutorial_request["tab"] == "tutorial"
    assert tutorial_request["operation"] == "open_lesson"

def test_static_report_operation_opens_the_published_bundle(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    _live_dashboard(monkeypatch)

    result = actions.invoke_dashboard(
        "reports",
        "open_report",
        {"report_id": "sam-heatmap-coloring-examples-20h30"},
    )
    request = activity.read_activity()["ui_request"]

    assert result["status"] == "requested"
    assert result["report_id"] == "sam-heatmap-coloring-examples-20h30"
    assert request["tab"] == "reports"
    assert request["operation"] == "open_report"
    assert request["report_id"] == "sam-heatmap-coloring-examples-20h30"
