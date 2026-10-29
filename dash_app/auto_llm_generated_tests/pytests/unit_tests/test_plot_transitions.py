"""Focused coverage for the shared Plot state-transition service."""

from pathlib import Path

import pytest

from dash_app.services import plots
from dash_app.shared import actions, activity


def _case_data():
    return {
        "name": "arm",
        "profile_vars": [{"label": "wp2", "value": "wp2"}],
        "budget_groups": [
            {"label": "wp2 budget", "value": "wp2"},
            {"label": "rtm budget", "value": "rtm"},
        ],
        "time_slider_start_min_seconds": 0.0,
        "time_slider_start_max_seconds": 7200.0,
        "default_time_start_seconds": 1200.0,
        "time_slider_duration_min_minutes": 10.0,
        "time_slider_duration_max_minutes": 120.0,
        "loss_time_start_seconds": 600.0,
        "loss_time_duration_minutes": 20.0,
        "benchmarks": {"available_sources": ["sam"]},
    }


def test_budget_request_defaults_to_wp2_and_applies_one_grid_instance():
    case_data = _case_data()
    request = plots.validate_plot_request(
        case_data,
        {
            "operation": "add_budget",
            "case": "arm",
            "time_start_seconds": 1800.0,
            "average_minutes": 30.0,
            "benchmark_sources": ["sam"],
        },
    )

    transition = plots.apply_plot_request(case_data, request, [], {}, 0)

    assert request["budget_group"] == "wp2"
    assert request["time_start_seconds"] == 1800.0
    assert request["average_minutes"] == 30.0
    assert request["benchmark_sources"] == ["sam"]
    assert transition.plot_order == [0]
    assert transition.plot_state["0"] == {"plot_type": "budget", "var": "wp2", "size": "normal"}
    assert transition.next_id == 1


def test_budget_group_validation_is_case_and_output_specific():
    with pytest.raises(ValueError, match="unavailable"):
        plots.validate_plot_request(
            _case_data(),
            {"operation": "add_budget", "case": "arm", "budget_group": "wp3"},
        )


def test_native_append_transition_uses_the_same_budget_validator():
    transition = plots.append_plot_instance(
        _case_data(),
        [4],
        {"4": {"plot_type": "profile", "var": "wp2", "size": "normal"}},
        5,
        plot_type="budget",
        variable="wp2",
    )

    assert transition.plot_order == [4, 5]
    assert transition.plot_state["5"]["plot_type"] == "budget"
    assert transition.plot_state["5"]["var"] == "wp2"


def test_plot_instances_have_stable_ids_and_remove_preserves_next_id():
    order = [2, 5]
    state = {
        "2": {"plot_type": "profile", "var": "wp2", "size": "normal"},
        "5": {"plot_type": "budget", "var": "rtm", "size": "large"},
    }

    assert plots.list_plot_instances(order, state) == [
        {"id": 2, "plot_type": "profile", "selection": "wp2", "size": "normal"},
        {"id": 5, "plot_type": "budget", "selection": "rtm", "size": "large"},
    ]
    transition = plots.remove_plot_instance(2, order, state, 6)

    assert transition.plot_order == [5]
    assert transition.plot_state == {"5": state["5"]}
    assert transition.next_id == 6


def test_plot_remove_rejects_an_id_not_owned_by_current_state():
    with pytest.raises(ValueError, match="not an active Plot card"):
        plots.remove_plot_instance(9, [2], {"2": {"plot_type": "profile", "var": "wp2"}}, 3)


def test_budget_handoff_accepts_exact_physical_values_with_a_named_preset():
    request = plots.validate_plot_request(
        _case_data(),
        {
            "operation": "add_budget",
            "case": "arm",
            "budget_group": "wp2",
            "time_start_seconds": 600.0,
            "average_minutes": 20.0,
            "window_preset": "loss",
        },
    )

    assert request["window_preset"] == "loss"
    assert request["time_start_seconds"] == 600.0
    assert request["average_minutes"] == 20.0


def test_typed_budget_action_publishes_one_validated_handoff(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "_profile_case_data", lambda *_args, **_kwargs: _case_data())

    result = actions.add_budget_plot("arm")
    request = activity.read_activity()["plot_request"]

    assert result["status"] == "requested"
    assert result["plot_type"] == "budget"
    assert result["budget_group"] == "wp2"
    assert request["operation"] == "add_budget"
    assert request["budget_group"] == "wp2"
    assert request["tab"] == "plots"


def test_typed_budget_action_selects_nested_output_instead_of_legacy_root(tmp_path, monkeypatch):
    nested = tmp_path / "output" / "dash_default"
    nested.mkdir(parents=True)
    monkeypatch.setattr(actions, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(actions, "_validate_case", lambda value: str(value))
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()
    captured = {}

    def fake_case_data(*_args, **kwargs):
        captured["output_dirs"] = kwargs["output_dirs"]
        return _case_data()

    monkeypatch.setattr(actions, "_profile_case_data", fake_case_data)

    result = actions.add_budget_plot("arm", output_dir="dash_default")

    assert captured["output_dirs"] == [str(nested.resolve())]
    assert result["output_directory"] == str(nested.resolve())
    assert activity.read_activity()["plot_request"]["output_dir"] == str(nested.resolve())


def test_plot_output_selector_rejects_outside_and_conflicting_paths(tmp_path, monkeypatch):
    from dash_app.services.jobs import JobStore

    output = tmp_path / "output" / "dash_default"
    other = tmp_path / "output" / "other"
    output.mkdir(parents=True)
    other.mkdir(parents=True)
    monkeypatch.setattr(actions, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(actions, "_validate_case", lambda value: str(value))
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "_profile_case_data", lambda *_args, **_kwargs: _case_data())

    with pytest.raises(ValueError, match="inside the repository output"):
        actions.add_budget_plot("arm", output_dir=str(tmp_path / "outside"))

    store = JobStore(tmp_path / "jobs.json")
    record, _created = store.submit("scm", "plot-selector-request", {"case": "arm"})
    store.update(record["job_id"], output_directory=str(output))
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    with pytest.raises(ValueError, match="select different output directories"):
        actions.add_budget_plot("arm", run_id=record["run_id"], output_dir="other")

    accepted = actions.add_budget_plot("arm", run_id=record["run_id"], output_dir=str(output))
    assert accepted["output_directory"] == str(output.resolve())


def test_plot_manifest_and_invoker_expose_budget_as_a_separate_typed_operation(monkeypatch):
    monkeypatch.setattr(
        actions.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "available", "dashboard_id": "test-dashboard"},
    )
    captured = {}

    def fake_add(case, budget_group="wp2", **kwargs):
        captured.update({"case": case, "budget_group": budget_group, **kwargs})
        return {"status": "requested", "plot_type": "budget", "budget_group": budget_group}

    monkeypatch.setattr(actions, "add_budget_plot", fake_add)
    manifest = actions.inspect_dashboard("plots")
    result = actions.invoke_dashboard(
        "plots",
        "add_budget",
        {"case": "arm", "budget_group": "wp2"},
    )

    assert "add_budget" in {item["name"] for item in manifest["tabs"][0]["operations"]}
    assert result["status"] == "requested"
    assert captured["case"] == "arm"
    assert captured["budget_group"] == "wp2"


def test_typed_plot_list_and_remove_validate_current_snapshot(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    activity.set_plot_instances(
        {
            "case": "arm",
            "output_dirs": ["output"],
            "plots": [{"id": 4, "plot_type": "budget", "selection": "wp2", "size": "normal"}],
            "next_id": 5,
        }
    )
    monkeypatch.setattr(
        actions.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "available", "dashboard_id": "test-dashboard"},
    )

    listed = actions.invoke_dashboard("plots", "list")
    removed = actions.invoke_dashboard("plots", "remove", {"plot_id": 4})
    request = activity.read_activity()["ui_request"]

    assert listed["plots"] == [{"id": 4, "plot_type": "budget", "selection": "wp2", "size": "normal"}]
    assert removed["status"] == "requested"
    assert removed["plot_id"] == 4
    assert request["operation"] == "remove"
    assert request["plot_id"] == 4

    with pytest.raises(ValueError, match="not an active Plot card"):
        actions.invoke_dashboard("plots", "remove", {"plot_id": 9})


def test_typed_budget_action_rejects_unavailable_benchmark_source(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "_profile_case_data", lambda *_args, **_kwargs: _case_data())

    # The budget family has no benchmark control in its typed schema.  The
    # shared service still guarantees that any future common overlay field is
    # validated against the loaded case rather than accepted blindly.
    with pytest.raises(ValueError):
        plots.validate_plot_request(
            _case_data(),
            {"operation": "add_budget", "case": "arm", "benchmark_sources": ["coamps"]},
        )


def test_ui_and_mcp_use_the_same_benchmark_source_resolver():
    case_data = {"benchmarks": {"available_sources": ["sam"]}}
    assert plots.validate_plot_request(
        {**_case_data(), **case_data},
        {"operation": "add_budget", "case": "arm", "benchmark_sources": ["sam"]},
    )["benchmark_sources"] == plots.resolve_benchmark_sources(case_data, ["sam"], strict=True)
    assert plots.toggle_benchmark_source(case_data, [], "sam") == ["sam"]
