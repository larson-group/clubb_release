import json
from pathlib import Path

from dash_app.shared import activity


def test_activity_stream_keeps_events_and_latest_plot_request(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")

    activity.reset_activity()
    first = activity.publish_event("run", "Starting ARM", "iiPDF_type=1", status="running")
    second = activity.publish_plot_request("arm", ["rtp3", "thlp3"], time_start_seconds=42000, average_minutes=30)
    snapshot = activity.read_activity()

    assert [event["id"] for event in snapshot["events"]] == [first["id"], second["id"]]
    assert snapshot["plot_request"]["tab"] == "plots"
    assert snapshot["plot_request"]["time_start_seconds"] == 42000.0


def test_ui_request_queue_preserves_rapid_handoffs(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")

    activity.reset_activity()
    first = activity.publish_plot_request("arm", ["rtp3"])
    second = activity.publish_plot_request("arm", ["thlp3"])
    third = activity.publish_budget_request("arm", "wp2")

    queued = activity.read_activity()["ui_request_queue"]
    assert [item["id"] for item in queued] == [first["id"], second["id"], third["id"]]

    assert activity.claim_ui_request(0)["id"] == first["id"]
    assert activity.claim_ui_request(0) is None
    assert activity.claim_ui_request(first["id"]) is None
    assert activity.acknowledge_ui_request(first["id"])
    assert activity.claim_ui_request(first["id"])["id"] == second["id"]
    assert activity.acknowledge_ui_request(second["id"])
    assert activity.claim_ui_request(second["id"])["id"] == third["id"]
    assert activity.acknowledge_ui_request(third["id"])
    assert activity.claim_ui_request(third["id"]) is None
    assert activity.claim_ui_request(0) is None


def test_global_broker_snapshot_excludes_detailed_scm_history(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()

    activity.set_broker_job(
        "compile",
        {"state": "running", "job": {"pid": 17, "log": "/tmp/compile.log"}, "log_tail": "building\n"},
    )
    activity.set_broker_job(
        "profile",
        {"state": "running", "pid": 19, "log": "/tmp/profile.log", "log_tail": "timing\n"},
    )
    jobs = activity.broker_jobs()
    assert jobs["compile"]["state"] == "running"
    assert jobs["compile"]["log_tail"] == "building\n"
    assert jobs["profile"]["log_tail"] == "timing\n"
    assert "runs" not in jobs
    assert "scm_batches" not in jobs
    assert set(jobs["run_summary"]) == {
        "state", "revision", "active_count", "queued", "running", "stopping"
    }


def test_active_job_count_is_generic_across_broker_job_groups():
    jobs = {
        "compile": {"state": "finished"},
        "runs": {
            "arm": {"state": "running"},
            "bomex": {"state": "queued"},
            "rico": {"state": "error"},
        },
        "future_kind": [
            {"state": "submitting"},
            {"state": "stopping"},
            {"state": "success"},
        ],
    }

    assert activity.count_active_jobs(jobs) == 4


def test_active_job_count_treats_an_active_record_as_one_job():
    aggregate_job = {
        "state": "running",
        "metadata": {"state": "running"},
    }

    assert activity.count_active_jobs(aggregate_job) == 1


def test_legacy_agent_state_is_dropped_on_read(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.ACTIVITY_PATH.write_text(
        '{"version": 5, "agents": [{"id": "stale"}], "messages": [{"body": "old"}]}',
        encoding="utf-8",
    )

    snapshot = activity.read_activity()

    assert "agents" not in snapshot
    assert "messages" not in snapshot
    assert snapshot["version"] == 10


def test_legacy_run_handoff_is_dropped_without_stealing_tab_focus(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    run_request = {"id": 4, "type": "run", "operation": "start", "tab": "run"}
    plot_request = {"id": 5, "type": "profile", "operation": "set_view", "tab": "plots"}
    activity.ACTIVITY_PATH.write_text(
        json.dumps(
            {
                "ui_request": run_request,
                "ui_request_in_flight": run_request,
                "ui_request_queue": [run_request, plot_request],
            }
        ),
        encoding="utf-8",
    )

    snapshot = activity.read_activity()

    assert snapshot["ui_request"] is None
    assert snapshot["ui_request_in_flight"] is None
    assert snapshot["ui_request_queue"] == [plot_request]
