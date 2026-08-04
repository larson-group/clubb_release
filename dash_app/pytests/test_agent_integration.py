import json
import stat
from pathlib import Path
from urllib.error import URLError

from flask import Flask

from dash_app.agent_integration import client as agent_client
from dash_app.shared import activity, broker, gateway
from dash_app.shared.provenance import runtime_source_fingerprint


def test_runtime_fingerprint_tracks_runtime_content_and_untracked_modules(tmp_path):
    runtime = tmp_path / "dash_app" / "shared"
    runtime.mkdir(parents=True)
    source = runtime / "runtime_module.py"
    source.write_text("VALUE = 1\n")
    first = runtime_source_fingerprint(tmp_path)

    source.write_text("VALUE = 2\n")
    second = runtime_source_fingerprint(tmp_path)
    assert second != first

    test_file = tmp_path / "dash_app" / "pytests" / "test_only.py"
    test_file.parent.mkdir()
    test_file.write_text("VALUE = 3\n")
    assert runtime_source_fingerprint(tmp_path) == second


def test_local_gateway_uses_transient_authenticated_actions(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    monkeypatch.setattr(gateway, "CONNECTION_PATH", Path(tmp_path) / "connection.json")
    activity.reset_activity()
    server = Flask(__name__)
    gateway.configure_gateway(server, 24567)
    connection = json.loads(gateway.CONNECTION_PATH.read_text())
    assert stat.S_IMODE(gateway.CONNECTION_PATH.stat().st_mode) == 0o600
    headers = {"X-CLUBB-Agent-Token": connection["token"]}
    client = server.test_client()

    assert client.post("/api/agent-integration/connect", headers=headers).status_code == 404
    status = client.get("/api/agent-integration/status", headers=headers)
    assert status.status_code == 200
    assert "agents" not in status.json

    manifest = client.post(
        "/api/agent-integration/actions",
        json={"action": "inspect_dashboard", "payload": {"tab": "tutorial"}},
        headers=headers,
    )
    assert manifest.status_code == 202
    assert manifest.json["tabs"][0]["tab"] == "tutorial"
    assert any(operation["name"] == "open_lesson" for operation in manifest.json["tabs"][0]["operations"])

    blocked = client.post(
        "/api/agent-integration/actions",
        json={"action": "run_scm", "payload": {"case": "arm"}},
        headers=headers,
    )
    assert blocked.status_code == 400
    assert blocked.json["error"]["code"] == "DEPRECATED_DIRECT_ACTION"
    assert client.post("/api/agent-integration/heartbeat", headers=headers).status_code == 404
    assert client.get("/api/agent-integration/messages", headers=headers).status_code == 404


def test_gateway_status_exposes_broker_snapshot(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    monkeypatch.setattr(gateway, "CONNECTION_PATH", Path(tmp_path) / "connection.json")
    activity.reset_activity()
    activity.set_broker_metadata(pid=77, state="running")
    activity.set_broker_job("compile", {"state": "running", "job": {"pid": 78}, "log_tail": "still building"})
    server = Flask(__name__)
    gateway.configure_gateway(server, 24568)
    connection = json.loads(gateway.CONNECTION_PATH.read_text())
    response = server.test_client().get(
        "/api/agent-integration/status",
        headers={"X-CLUBB-Agent-Token": connection["token"]},
    )

    assert response.status_code == 200
    assert response.json["broker"]["pid"] == 77
    assert response.json["jobs"]["compile"]["log_tail"] == "still building"


def test_gateway_rejects_oversized_payload(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    monkeypatch.setattr(gateway, "CONNECTION_PATH", Path(tmp_path) / "connection.json")
    activity.reset_activity()
    server = Flask(__name__)
    gateway.configure_gateway(server, 24569)
    connection = json.loads(gateway.CONNECTION_PATH.read_text())
    response = server.test_client().post(
        "/api/agent-integration/actions",
        json={"action": "inspect_dashboard", "payload": {"value": "x" * (300 * 1024)}},
        headers={"X-CLUBB-Agent-Token": connection["token"]},
    )

    assert response.status_code == 413
    assert response.json["error"]["code"] == "REQUEST_ENTITY_TOO_LARGE"


def test_gateway_rejects_non_loopback_even_with_connection_token(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    monkeypatch.setattr(gateway, "CONNECTION_PATH", Path(tmp_path) / "connection.json")
    activity.reset_activity()
    server = Flask(__name__)
    gateway.configure_gateway(server, 24570)
    connection = json.loads(gateway.CONNECTION_PATH.read_text())
    response = server.test_client().get(
        "/api/agent-integration/status",
        headers={"X-CLUBB-Agent-Token": connection["token"]},
        environ_overrides={"REMOTE_ADDR": "192.0.2.10"},
    )

    assert response.status_code == 403


def test_broker_stop_targets_only_the_recorded_sidecar(monkeypatch):
    calls = []
    monkeypatch.setattr(broker, "read_connection", lambda: {"pid": 901})
    states = iter((True, False, False))
    monkeypatch.setattr(broker, "_pid_is_alive", lambda _pid: next(states))
    monkeypatch.setattr(broker.os, "kill", lambda pid, signal: calls.append((pid, signal)))

    broker.stop_broker()

    assert calls == [(901, broker.signal.SIGTERM)]


def test_broker_passes_process_creation_time_to_stable_endpoint(monkeypatch):
    observed = {}

    class FakeServer:
        server_port = 24571

        def serve_forever(self):
            return None

        def shutdown(self):
            return None

    monkeypatch.setattr(broker, "make_server", lambda *_args, **_kwargs: FakeServer())
    monkeypatch.setattr(broker, "broker_process_started_at", lambda pid: observed.__setitem__("pid", pid) or 123.5)
    monkeypatch.setattr(
        broker,
        "write_connection",
        lambda _port, pid=None, **_kwargs: {"url": "http://127.0.0.1:24571/api", "token": "token", "pid": pid},
    )
    monkeypatch.setattr(broker, "install_gateway_routes", lambda *_args: None)
    monkeypatch.setattr(broker, "set_broker_metadata", lambda **_kwargs: None)
    monkeypatch.setattr(broker, "_recover_compile_monitoring", lambda: None)
    monkeypatch.setattr(broker, "_recover_tune_keepalive", lambda: None)
    monkeypatch.setattr(broker, "_recover_scm_monitoring", lambda: None)
    monkeypatch.setattr(
        broker,
        "start_broker_endpoint",
        lambda _connection, **kwargs: (
            observed.__setitem__("started_at", kwargs["owner_started_at"]) or {}
        ),
    )
    monkeypatch.setattr(broker, "update_connection_endpoint", lambda _endpoint: None)
    monkeypatch.setattr(broker, "stop_broker_endpoint", lambda: None)

    broker.serve()

    assert observed["pid"] == broker.os.getpid()
    assert observed["started_at"] == 123.5


def test_broker_protocol_marks_old_sidecar_incompatible():
    assert broker._connection_is_compatible({"version": broker.BROKER_PROTOCOL_VERSION})
    assert not broker._connection_is_compatible({"version": broker.BROKER_PROTOCOL_VERSION - 1})
    assert not broker._connection_is_compatible({})


def test_broker_runtime_fingerprint_staleness_requires_replacement():
    assert not broker._connection_is_stale({"runtime_fingerprint": "same"}, "same")
    assert broker._connection_is_stale({"runtime_fingerprint": "old"}, "new")
    assert broker._connection_is_stale({}, "new")


def test_broker_active_work_guard_includes_queued_and_submitting_jobs(monkeypatch):
    class Response:
        def __init__(self, payload):
            self.payload = payload

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            return json.dumps(self.payload).encode("utf-8")

    monkeypatch.setattr(
        broker,
        "urlopen",
        lambda *_args, **_kwargs: Response({"jobs": {"runs": {"run-1": {"state": "queued"}}}}),
    )
    assert broker._broker_has_active_work({"url": "http://127.0.0.1:1", "token": "token"})


def test_require_fresh_preserves_active_stale_broker(monkeypatch, tmp_path):
    connection = {"version": broker.BROKER_PROTOCOL_VERSION, "runtime_fingerprint": "old", "pid": 901}
    monkeypatch.setattr(broker, "_existing_live_connection", lambda: connection)
    monkeypatch.setattr(broker, "_connection_is_compatible", lambda _connection: True)
    monkeypatch.setattr(broker, "_connection_is_stale", lambda _connection, _expected: True)
    monkeypatch.setattr(broker.dashboard_registry, "dashboard_is_live", lambda: False)
    monkeypatch.setattr(broker, "_broker_has_active_work", lambda _connection: True)
    monkeypatch.setattr(broker, "BROKER_LOCK_PATH", tmp_path / "broker.lock")

    import pytest

    with pytest.raises(RuntimeError, match="work is active"):
        broker.ensure_broker(expected_fingerprint="new", require_fresh=True)


def test_gateway_discovers_connection_from_another_private_runtime_root(tmp_path, monkeypatch):
    primary = tmp_path / "primary" / "connection.json"
    alternate = tmp_path / "alternate" / "connection.json"
    alternate.parent.mkdir(parents=True)
    alternate.write_text(json.dumps({"url": "http://127.0.0.1:24571/api", "token": "alternate"}))
    alternate.chmod(0o600)
    monkeypatch.setattr(gateway, "CONNECTION_PATH", primary)
    monkeypatch.setattr(gateway, "readable_private_paths", lambda _root, _name: [alternate])

    assert gateway.read_connection()["token"] == "alternate"


def test_transient_client_connect_does_not_persist_session(monkeypatch):
    monkeypatch.setattr(agent_client, "request_gateway", lambda *args, **kwargs: {"broker": {"state": "running"}})

    assert agent_client.connect() == {
        "broker": {"state": "running"},
        "status": "connected",
        "allowed_actions": ["inspect_dashboard", "invoke_dashboard"],
    }


def test_gateway_client_retries_when_broker_republishes_connection(monkeypatch):
    connections = iter(
        (
            {"url": "http://127.0.0.1:1111/api", "token": "old"},
            {"url": "http://127.0.0.1:2222/api", "token": "new"},
        )
    )
    requests = []

    class _Response:
        status = 200

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            return json.dumps({"status": "ok"}).encode("utf-8")

    def _urlopen(request, timeout):
        requests.append(request.full_url)
        if len(requests) == 1:
            raise URLError("broker replacement")
        return _Response()

    monkeypatch.setattr(agent_client, "read_connection", lambda: next(connections))
    monkeypatch.setattr(agent_client, "urlopen", _urlopen)

    assert agent_client.request_gateway("GET", "status") == {"status": "ok"}
    assert requests == ["http://127.0.0.1:1111/api/status", "http://127.0.0.1:2222/api/status"]
