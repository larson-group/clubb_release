"""Unchanged/hidden views should not repeat expensive data work."""

import json
from types import SimpleNamespace

from dash import Dash, no_update

from dash_app.profile_tab import callbacks as profile
from dash_app.profile_tab.runtime import profile_results_signature
from utilities import timing_profiles


def callback(app, name):
    return next(
        entry["callback"].__wrapped__
        for entry in app.callback_map.values()
        if entry.get("callback") and entry["callback"].__name__ == name
    )


def test_profile_poll_only_publishes_changed_results(tmp_path, monkeypatch):
    directory = tmp_path / "example"
    directory.mkdir()
    timing_profiles.write_profile_manifest(directory, {"run_id": "example"})
    results = directory / "timings.csv"
    results.write_text("first")
    reads = []

    def read(job):
        reads.append(results.read_text())
        return job["run_id"], [{"timer": reads[-1]}]

    monkeypatch.setattr(profile, "read_profile_results", read)
    monkeypatch.setattr(profile, "callback_context", SimpleNamespace(triggered_id="profile-interval"))
    app = Dash(__name__)
    profile.register_profile_callbacks(app)
    refresh = callback(app, "refresh_profile")
    job = {"output": str(tmp_path), "run_id": "example", "state": "running", "log_tail": "started"}
    first = refresh(1, {}, {"profile": job}, "profile", job, {})
    assert first[0] is no_update
    assert first[1]["row_count"] == 1
    assert reads == ["first"]

    # Heartbeats and logs can change without invalidating figures or reading CSVs.
    latest = {**job, "log_tail": "still running", "updated_at": 5}
    same = refresh(2, {}, {"profile": latest}, "profile", job, first[1])
    assert same[0] == latest
    assert same[1] is no_update
    assert same[4] == "still running"
    assert len(reads) == 1

    # Atomic replacement with equal size/row count must still refresh the data.
    replacement = directory / "replacement"
    replacement.write_text("other")
    replacement.replace(results)
    changed = refresh(3, {}, {"profile": latest}, "profile", latest, first[1])
    assert changed[1] is not no_update
    assert changed[1]["row_count"] == 1
    assert reads == ["first", "other"]

    finished = {**latest, "state": "finished"}
    terminal = refresh(4, {}, {"profile": finished}, "profile", latest, changed[1])
    assert terminal[1]["state"] == "finished"
    assert terminal[2] == "Complete"
    assert len(reads) == 2
    hidden = refresh(5, {}, {"profile": finished}, "welcome", finished, terminal[1])
    assert all(value is no_update for value in hidden)


def test_profile_signature_tracks_creation_deletion_and_direct_profile_paths(tmp_path):
    directory = tmp_path / "example"
    directory.mkdir()
    timing_profiles.write_profile_manifest(directory, {"run_id": "example"})
    job = {"output": str(tmp_path), "run_id": "example"}
    before = profile_results_signature(job)
    assert json.loads(json.dumps(before), parse_int=float) == before
    assert profile_results_signature({"output": str(tmp_path)}) is None
    assert profile_results_signature({**job, "output": str(directory)}) == before
    (directory / "batches.csv").write_text("new data")
    created = profile_results_signature(job)
    assert created != before
    (directory / "batches.csv").unlink()
    assert profile_results_signature(job) == before


def test_profile_hidden_graphs_and_unchanged_choices_do_not_trigger_work(monkeypatch):
    def unexpected(*args):
        raise AssertionError("Hidden Profile must not load figure data")

    monkeypatch.setattr(profile, "load_profile_plot_data", unexpected)
    app = Dash(__name__)
    profile.register_profile_callbacks(app)
    refresh = callback(app, "refresh_profile_figures")
    assert all(value is no_update for value in refresh(*([None] * 12), "welcome", {}))
    timer = callback(app, "update_timer_choices")
    options, value = timer([{"timer_name": "example"}], "example")
    assert "example" in {item["value"] for item in options}
    assert value is no_update


def test_profile_load_derives_both_tables_once_per_selected_profile(tmp_path, monkeypatch):
    calls = []

    def derive(directory):
        calls.append(directory)
        return [{"summary": 1}], [{"process": 2}]

    monkeypatch.setattr(timing_profiles, "_derived_rows", derive)
    summaries, processes = timing_profiles.load_profiles(tmp_path, ["example", "example"])
    assert calls == [tmp_path / "example"]
    assert summaries == [{"summary": 1, "profile_id": "example", "profile_label": "example"}]
    assert processes == [{"process": 2, "profile_id": "example", "profile_label": "example"}]


def test_profile_library_refresh_reloads_unchanged_selection(tmp_path, monkeypatch):
    monkeypatch.setattr(profile, "discover_profile_library", lambda _: [{"run_id": "example", "label": "Example"}])
    context = SimpleNamespace(triggered_id="profile-active-results")
    monkeypatch.setattr(profile, "callback_context", context)
    app = Dash(__name__)
    profile.register_profile_callbacks(app)
    refresh = callback(app, "refresh_profile_library")
    args = (str(tmp_path), {"kind": "refreshed"}, {}, ["example"], "example", "example", str(tmp_path), {})
    unchanged = refresh(*args)
    assert unchanged[2] is no_update
    assert unchanged[4] is no_update
    assert unchanged[6] is no_update
    context.triggered_id = "profile-library-action"
    assert refresh(*args)[2] == ["example"]


def test_report_catalog_only_scans_while_visible(monkeypatch):
    from dash_app.reports_tab import tab

    app = Dash(__name__)
    tab.build_tab(app)
    scans = []
    monkeypatch.setattr(tab, "discover_reports", lambda: scans.append(True) or [])
    refresh = callback(app, "refresh_report_catalog")
    assert refresh(1, "welcome", "old", None) == (no_update,) * 3
    assert scans == []
    changed = refresh(1, "reports", "old", None)
    assert scans == [True]
    assert changed[1] == tab.catalog_token([])
    assert refresh(2, "reports", changed[1], None) == (no_update,) * 3


def test_tune_activity_only_polls_while_visible_and_suppresses_equal_results(monkeypatch):
    from dash_app.shared import broker_client
    from dash_app.tune_tab.callbacks_workspaces import register_workspace_callbacks

    calls = []
    records = [{"id": "example", "state": "running"}]

    def perform(*args, **kwargs):
        calls.append((args, kwargs))
        return {"activity": records}

    monkeypatch.setattr(broker_client, "perform_action", perform)
    app = Dash(__name__)
    register_workspace_callbacks(app)
    refresh = callback(app, "refresh_workspace_activity")
    assert refresh(1, "welcome", []) is no_update
    assert not calls
    assert refresh(1, "tune", []) == records
    assert calls[-1][1]["ensure_running"] is False
    assert refresh(2, "tune", records) is no_update
    records = [{"id": "example", "state": "finished"}]
    assert refresh(3, "tune", []) == records
