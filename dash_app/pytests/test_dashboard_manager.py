from __future__ import annotations

import os


def test_manager_restart_policy_is_ten_seconds_for_five_minutes():
    from dash_app import manager

    assert manager.RESTART_INTERVAL_SECONDS == 10.0
    assert manager.RESTART_DEADLINE_SECONDS == 300.0


def test_active_job_labels_cover_every_broker_job_family():
    from dash_app import manager

    jobs = {
        "compile": {"state": "running"},
        "profile": {"state": "running"},
        "run_summary": {"active_count": 1},
        "tune": {"state": "stopping"},
        "loss_runs": {"loss_1": {"state": "queued"}},
    }

    assert manager._active_job_labels(jobs) == [
        "compile",
        "Profile",
        "Tune",
        "SCM:1 active",
        "Tune-result:loss_1",
    ]


def test_managed_broker_connection_is_required_in_manager_environment(monkeypatch):
    from dash_app.shared import broker
    from dash_app.shared.manager_lease import MANAGER_REQUIRED_ENV

    monkeypatch.setenv(MANAGER_REQUIRED_ENV, "1")
    assert broker._connection_is_stale({"runtime_fingerprint": "same"}, "same")
    assert not broker._connection_is_stale(
        {"runtime_fingerprint": "same", "manager_required": True}, "same"
    )


def test_manager_lease_validates_timestamp_and_process_identity(monkeypatch):
    from dash_app.shared import manager_lease

    class Process:
        def __init__(self, pid):
            assert pid == os.getpid()

        def create_time(self):
            return 123.0

    monkeypatch.setattr(manager_lease.psutil, "Process", Process)
    monkeypatch.setattr(manager_lease, "LEASE_TIMEOUT_SECONDS", 30.0)
    record = {"pid": os.getpid(), "started_at": 123.0, "updated_at": 100.0}

    assert manager_lease.manager_lease_is_live(record, now=129.9)
    assert not manager_lease.manager_lease_is_live(record, now=130.1)
    assert not manager_lease.manager_lease_is_live(
        {**record, "started_at": 321.0}, now=101.0
    )


def test_broker_lease_expiry_stops_work_before_server(monkeypatch):
    from dash_app.shared import actions, broker

    calls = []

    class StopEvent:
        stopped = False

        def wait(self, _timeout):
            return self.stopped

        def set(self):
            self.stopped = True

    class Server:
        def shutdown(self):
            calls.append("server")

    monkeypatch.setattr(broker, "manager_lease_is_live", lambda: False)
    monkeypatch.setattr(broker, "count_active_jobs", lambda _jobs: 0)
    monkeypatch.setattr(broker, "broker_jobs", lambda: {})
    monkeypatch.setattr(
        broker.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "unavailable"},
    )
    monkeypatch.setattr(broker, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        actions,
        "stop_all_broker_work",
        lambda **_kwargs: calls.append("work") or {"errors": []},
    )

    broker._stop_after_manager_lease_expires(Server(), StopEvent())

    assert calls == ["work", "server"]


def test_broker_lease_expiry_stops_orphaned_dash_group(monkeypatch):
    from dash_app.shared import actions, broker

    calls = []

    class StopEvent:
        stopped = False

        def wait(self, _timeout):
            return self.stopped

        def set(self):
            self.stopped = True

    class Server:
        def shutdown(self):
            calls.append("server")

    monkeypatch.setattr(broker, "manager_lease_is_live", lambda: False)
    monkeypatch.setattr(broker, "count_active_jobs", lambda _jobs: 0)
    monkeypatch.setattr(broker, "broker_jobs", lambda: {})
    monkeypatch.setattr(broker, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        broker.dashboard_registry,
        "dashboard_status",
        lambda: {"status": "available", "pid": 456, "started_at": 123.0},
    )
    monkeypatch.setattr(
        broker.dashboard_registry,
        "process_identity_is_live",
        lambda *_args: False,
    )
    monkeypatch.setattr(broker.os, "getpgid", lambda pid: pid)
    monkeypatch.setattr(
        broker.os,
        "killpg",
        lambda pid, signum: calls.append(("dash", pid, signum)),
    )
    monkeypatch.setattr(
        actions,
        "stop_all_broker_work",
        lambda **_kwargs: calls.append("work") or {"errors": []},
    )

    broker._stop_after_manager_lease_expires(Server(), StopEvent())

    assert calls == [("dash", 456, broker.signal.SIGTERM), "work", "server"]


def test_stop_all_work_uses_durable_tune_metadata(monkeypatch):
    from dash_app.shared import actions

    stopped = []

    class Store:
        def list_kind(self, _kind):
            return []

        def update(self, *_args, **_kwargs):
            return None

    monkeypatch.setattr(
        actions,
        "broker_jobs",
        lambda: {
            "compile": None,
            "runs": {},
            "tune": {"state": "running", "job": {"pid": 123, "job_dir": "/tmp/tune"}},
            "loss_runs": {},
        },
    )
    monkeypatch.setattr(actions, "stop_tuning_job", lambda job: stopped.append(job))
    monkeypatch.setattr(actions, "update_broker_job", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(actions, "_JOB_STORE", Store())
    actions._BROKER_SHUTTING_DOWN.clear()
    try:
        result = actions.stop_all_broker_work(reason="test shutdown")
    finally:
        actions._BROKER_SHUTTING_DOWN.clear()

    assert stopped == [{"pid": 123, "job_dir": "/tmp/tune"}]
    assert result["requested"] == ["Tune"]
    assert result["errors"] == []
