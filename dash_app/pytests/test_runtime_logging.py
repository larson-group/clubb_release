import os
import stat
from pathlib import Path

from dash_app.shared import runtime_logging


def test_dashboard_log_path_uses_private_runtime(monkeypatch, tmp_path):
    monkeypatch.setattr(runtime_logging, "private_path", lambda _root, name: tmp_path / name)
    assert runtime_logging.dashboard_log_path(tmp_path) == tmp_path / "dash.log"


def test_relay_preserves_stream_and_rotates_private_log(tmp_path, monkeypatch):
    path = tmp_path / "dash.log"
    monkeypatch.setattr(runtime_logging, "APP_LOG_MAX_BYTES", 4)
    monkeypatch.setattr(runtime_logging, "APP_LOG_BACKUPS", 2)
    runtime_logging._write_chunk(path, b"old!")
    runtime_logging._write_chunk(path, b"new")
    assert path.read_bytes() == b"new"
    assert path.with_name("dash.log.1").read_bytes() == b"old!"
    assert stat.S_IMODE(path.stat().st_mode) == 0o600


def test_connection_log_metadata_contains_only_paths(tmp_path, monkeypatch):
    from dash_app.shared import gateway

    connection_path = tmp_path / "connection.json"
    monkeypatch.setattr(gateway, "CONNECTION_PATH", connection_path)
    gateway.write_connection(24567)
    updated = gateway.update_connection_logs({"app": str(tmp_path / "dash.log")})
    assert updated["log_paths"] == {"app": str(tmp_path / "dash.log")}
    assert "token" in updated
    assert "payload" not in updated
    assert stat.S_IMODE(connection_path.stat().st_mode) == 0o600


def test_connection_publishes_runtime_fingerprint_metadata(tmp_path, monkeypatch):
    from dash_app.shared import gateway

    connection_path = tmp_path / "connection.json"
    monkeypatch.setattr(gateway, "CONNECTION_PATH", connection_path)
    connection = gateway.write_connection(24567, runtime_fingerprint="abc123")

    assert connection["runtime_fingerprint"] == "abc123"
    assert connection["runtime_fingerprint_scope"]
