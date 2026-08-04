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


def test_broker_job_snapshot_survives_activity_reset_hook(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()

    activity.set_broker_job(
        "compile",
        {"state": "running", "job": {"pid": 17, "log": "/tmp/compile.log"}, "log_tail": "building\n"},
    )
    activity.set_broker_run("arm", {"state": "running", "proc_data": {"pid": 18, "log": "/tmp/arm.log"}})
    activity.update_broker_run("arm", log_tail="running arm\n")

    jobs = activity.broker_jobs()
    assert jobs["compile"]["state"] == "running"
    assert jobs["compile"]["log_tail"] == "building\n"
    assert jobs["runs"]["arm"]["log_tail"] == "running arm\n"


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
    assert snapshot["version"] == 6
