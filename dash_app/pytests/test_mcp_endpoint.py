import json
from pathlib import Path

from dash_app.agent_integration import endpoint
from dash_app.shared.runtime import atomic_write_json


def _record(instance_id: str, *, dashboard_pid=10, endpoint_pid=20, token="owned"):
    return {
        "version": endpoint.ENDPOINT_SCHEMA_VERSION,
        "instance_id": instance_id,
        "endpoint_url": "http://127.0.0.1:41234/mcp",
        "endpoint_token": "endpoint-secret",
        "endpoint_pid": endpoint_pid,
        "endpoint_started_at": 2.0,
        "dashboard": {"pid": dashboard_pid, "started_at": 1.0},
        "broker": {"pid": 30, "url": "http://127.0.0.1:41235/api", "token": token},
    }


def test_endpoint_record_requires_exact_process_and_broker_identity(monkeypatch):
    instance_id = "a" * 32
    live_processes = {(10, 1.0), (20, 2.0)}
    monkeypatch.setattr(endpoint, "process_identity_is_live", lambda pid, started: (pid, started) in live_processes)
    monkeypatch.setattr(endpoint, "broker_connection_is_live", lambda connection: connection["token"] == "owned")

    assert endpoint.endpoint_record_is_live(_record(instance_id))
    assert not endpoint.endpoint_record_is_live(_record(instance_id, dashboard_pid=11))
    assert not endpoint.endpoint_record_is_live(_record(instance_id, token="other"))


def test_reconcile_removes_stale_private_endpoint_records(tmp_path, monkeypatch):
    monkeypatch.setattr(endpoint, "ENDPOINT_ROOT", tmp_path)
    instance_id = "b" * 32
    directory = tmp_path / instance_id
    directory.mkdir()
    atomic_write_json(directory / "endpoint.json", _record(instance_id))
    atomic_write_json(directory / "endpoint-config.json", {"instance_id": instance_id})
    monkeypatch.setattr(endpoint, "endpoint_record_is_live", lambda _record: False)

    assert endpoint.reconcile_instance_records() == []
    assert not directory.exists()


def test_start_endpoint_publishes_owner_bound_manual_connection(monkeypatch, tmp_path):
    monkeypatch.setattr(endpoint, "ENDPOINT_ROOT", tmp_path)
    monkeypatch.setattr(endpoint, "reconcile_instance_records", lambda: [])
    monkeypatch.setattr(endpoint, "_local_port", lambda: 41236)
    monkeypatch.setattr(endpoint, "process_started_at", lambda _pid: 1.0)
    broker = {"url": "http://127.0.0.1:41235/api", "token": "broker-secret", "pid": 30, "version": 1}

    class FakeProcess:
        pid = 77

        def __init__(self):
            self.stopped = False

        def poll(self):
            return 0 if self.stopped else None

    process = FakeProcess()

    def fake_popen(command, **_kwargs):
        config = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        record = {
            "version": config["version"],
            "instance_id": config["instance_id"],
            "endpoint_url": "http://127.0.0.1:41236/mcp",
            "endpoint_token": config["endpoint_token"],
            "endpoint_pid": process.pid,
            "endpoint_started_at": 2.0,
            "dashboard": config["dashboard"],
            "broker": config["broker"],
        }
        atomic_write_json(Path(config["registry_path"]), record)
        return process

    monkeypatch.setattr(endpoint.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(endpoint.os, "kill", lambda _pid, _signal: setattr(process, "stopped", True))
    handle = endpoint.start_dashboard_endpoint(broker, dashboard_pid=10, dashboard_started_at=1.0)

    details = handle.public_details()
    assert details["instance_id"] == handle.instance_id
    assert details["endpoint_url"] == "http://127.0.0.1:41236/mcp"
    assert details["endpoint_token"] == handle.record["endpoint_token"]
    assert details["dashboard_pid"] == 10
    assert details["broker_pid"] == 30

    handle.stop()
    assert not handle.config_path.exists()
    assert not (tmp_path / handle.instance_id / "endpoint.json").exists()
