"""Broker result-run metadata must retain its ID while reaching terminal state.

The ordinary Dash result buttons use this same watcher as typed leaderboard
reruns. The real metadata includes run_id, which must not collide with the
status updater's selector argument.
"""

from dash_app.shared import actions, activity


def test_replay_watcher_records_completion_with_run_id_in_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()
    run = {"run_id": "completed-run", "pid": 123, "state": "running"}
    activity.set_broker_loss_run(run["run_id"], run)
    monkeypatch.setattr(actions, "poll_loss_runs", lambda runs: (
        {"completed-run": {**run, "state": "success", "returncode": 0}}, False,
    ))
    actions._watch_tuning_loss_run(run)
    assert activity.broker_jobs()["loss_runs"][run["run_id"]]["state"] == "success"
