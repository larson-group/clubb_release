"""Legacy rotating endpoint helpers retained only for migration tests.

Production uses :mod:`broker_endpoint`, which owns one stable endpoint per
checkout. These records are reconciled and removed when the broker starts.
"""

from __future__ import annotations

import json
import os
import re
import secrets
import signal
import socket
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

import psutil

from dash_app.shared.activity import REPO_ROOT
from dash_app.shared.runtime import atomic_write_json, private_path


ENDPOINT_SCHEMA_VERSION = 1
ENDPOINT_ROOT = private_path(REPO_ROOT, "mcp-instances")
ENDPOINT_PATH = "/mcp"
ENDPOINT_START_TIMEOUT_SECONDS = 8.0
ENDPOINT_STOP_TIMEOUT_SECONDS = 3.0
_INSTANCE_ID_RE = re.compile(r"^[0-9a-f]{32}$")


def _instance_dir(instance_id: str) -> Path:
    if not _INSTANCE_ID_RE.fullmatch(str(instance_id)):
        raise ValueError("invalid dashboard MCP instance id")
    return ENDPOINT_ROOT / str(instance_id)


def _config_path(instance_id: str) -> Path:
    return _instance_dir(instance_id) / "endpoint-config.json"


def _record_path(instance_id: str) -> Path:
    return _instance_dir(instance_id) / "endpoint.json"


def process_started_at(pid: int) -> float:
    """Return the OS start time used to prevent PID-reuse confusion."""
    try:
        return float(psutil.Process(int(pid)).create_time())
    except (psutil.Error, TypeError, ValueError, OSError) as exc:
        raise RuntimeError("could not determine the dashboard process start time") from exc


def process_identity_is_live(pid: Any, started_at: Any) -> bool:
    """Check both PID and process start time, not just PID existence."""
    try:
        process = psutil.Process(int(pid))
        expected = float(started_at)
        actual = float(process.create_time())
    except (psutil.Error, TypeError, ValueError, OSError):
        return False
    return abs(actual - expected) <= 1.0


def broker_connection_is_live(connection: dict[str, Any], *, timeout: float = 0.75) -> bool:
    """Verify the exact authenticated broker named by an endpoint record."""
    try:
        request = Request(
            str(connection["url"]).rstrip("/") + "/status",
            headers={"X-CLUBB-Agent-Token": str(connection["token"])},
            method="GET",
        )
        with urlopen(request, timeout=timeout) as response:  # nosec B310: validated loopback record
            return response.status == 200
    except (KeyError, OSError, URLError, ValueError):
        return False


def _local_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind(("127.0.0.1", 0))
        return int(listener.getsockname()[1])


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _unlink_if_present(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        pass


def cleanup_instance_record(instance_id: str, *, expected_endpoint_pid: int | None = None) -> bool:
    """Remove one record/config pair without touching another instance."""
    record_path = _record_path(instance_id)
    record = _read_json(record_path)
    if expected_endpoint_pid is not None and record is not None:
        try:
            if int(record.get("endpoint_pid")) != int(expected_endpoint_pid):
                return False
        except (TypeError, ValueError):
            return False
    _unlink_if_present(record_path)
    _unlink_if_present(_config_path(instance_id))
    try:
        _instance_dir(instance_id).rmdir()
    except OSError:
        pass
    return True


def endpoint_record_is_live(record: dict[str, Any], *, check_broker: bool = True) -> bool:
    """Return whether a published endpoint still belongs to live owners."""
    try:
        version = int(record.get("version") or 0)
    except (TypeError, ValueError):
        return False
    if version != ENDPOINT_SCHEMA_VERSION:
        return False
    instance_id = str(record.get("instance_id") or "")
    if not _INSTANCE_ID_RE.fullmatch(instance_id):
        return False
    dashboard = dict(record.get("dashboard") or {})
    endpoint_pid = record.get("endpoint_pid")
    if not process_identity_is_live(dashboard.get("pid"), dashboard.get("started_at")):
        return False
    if not process_identity_is_live(endpoint_pid, record.get("endpoint_started_at")):
        return False
    broker = dict(record.get("broker") or {})
    if not broker.get("url") or not broker.get("token"):
        return False
    return not check_broker or broker_connection_is_live(broker)


def reconcile_instance_records() -> list[dict[str, Any]]:
    """Return live endpoint records and remove stale private records."""
    ENDPOINT_ROOT.mkdir(mode=0o700, parents=True, exist_ok=True)
    try:
        ENDPOINT_ROOT.chmod(0o700)
    except OSError:
        pass
    live: list[dict[str, Any]] = []
    for record_path in sorted(ENDPOINT_ROOT.glob("*/endpoint.json")):
        record = _read_json(record_path)
        instance_id = str((record or {}).get("instance_id") or record_path.parent.name)
        if record is not None and endpoint_record_is_live(record):
            live.append(record)
        else:
            try:
                cleanup_instance_record(instance_id)
            except ValueError:
                _unlink_if_present(record_path)
                _unlink_if_present(record_path.parent / "endpoint-config.json")
                try:
                    record_path.parent.rmdir()
                except OSError:
                    pass
    return live


def load_endpoint_config(path: Path | str) -> dict[str, Any]:
    config = _read_json(Path(path))
    if config is None:
        raise RuntimeError("dashboard MCP endpoint configuration is missing or invalid")
    instance_id = str(config.get("instance_id") or "")
    if not _INSTANCE_ID_RE.fullmatch(instance_id):
        raise RuntimeError("dashboard MCP endpoint configuration has an invalid instance id")
    if str(config.get("endpoint_path") or "") != ENDPOINT_PATH:
        raise RuntimeError("dashboard MCP endpoint configuration has an invalid MCP path")
    if not dict(config.get("broker") or {}).get("url") or not dict(config.get("broker") or {}).get("token"):
        raise RuntimeError("dashboard MCP endpoint configuration has no broker identity")
    return config


@dataclass
class DashboardEndpointHandle:
    """Parent-side handle for one endpoint child process."""

    instance_id: str
    config_path: Path
    record: dict[str, Any]
    process: subprocess.Popen[Any]

    @property
    def endpoint_url(self) -> str:
        return str(self.record["endpoint_url"])

    @property
    def endpoint_token(self) -> str:
        return str(self.record["endpoint_token"])

    def public_details(self) -> dict[str, Any]:
        """Return the exact values a user needs for manual MCP setup."""
        return {
            "instance_id": self.instance_id,
            "endpoint_url": self.endpoint_url,
            "endpoint_token": self.endpoint_token,
            "dashboard_pid": self.record.get("dashboard", {}).get("pid"),
            "broker_pid": self.record.get("broker", {}).get("pid"),
        }

    def stop(self) -> None:
        if self.process.poll() is None:
            try:
                os.kill(self.process.pid, signal.SIGTERM)
            except OSError:
                pass
            deadline = time.monotonic() + ENDPOINT_STOP_TIMEOUT_SECONDS
            while self.process.poll() is None and time.monotonic() < deadline:
                time.sleep(0.05)
            if self.process.poll() is None:
                try:
                    os.kill(self.process.pid, signal.SIGKILL)
                except OSError:
                    pass
        cleanup_instance_record(self.instance_id, expected_endpoint_pid=self.process.pid)


def start_dashboard_endpoint(
    broker_connection: dict[str, Any],
    *,
    dashboard_pid: int | None = None,
    dashboard_started_at: float | None = None,
) -> DashboardEndpointHandle:
    """Launch one authenticated HTTP endpoint bound to this Dash process."""
    reconcile_instance_records()
    pid = int(dashboard_pid or os.getpid())
    started_at = float(dashboard_started_at if dashboard_started_at is not None else process_started_at(pid))
    instance_id = uuid.uuid4().hex
    instance_dir = _instance_dir(instance_id)
    instance_dir.mkdir(mode=0o700, parents=True, exist_ok=False)
    instance_dir.chmod(0o700)
    config_path = _config_path(instance_id)
    endpoint_token = secrets.token_urlsafe(32)
    endpoint_port = _local_port()
    config = {
        "version": ENDPOINT_SCHEMA_VERSION,
        "instance_id": instance_id,
        "endpoint_path": ENDPOINT_PATH,
        "endpoint_host": "127.0.0.1",
        "endpoint_port": endpoint_port,
        "endpoint_token": endpoint_token,
        "dashboard": {"pid": pid, "started_at": started_at},
        "broker": {
            "url": str(broker_connection["url"]),
            "token": str(broker_connection["token"]),
            "pid": broker_connection.get("pid"),
            "version": broker_connection.get("version"),
            "repository": broker_connection.get("repository"),
        },
        "registry_path": str(_record_path(instance_id)),
    }
    atomic_write_json(config_path, config)
    log_path = instance_dir / "endpoint.log"
    with log_path.open("a", encoding="utf-8") as log_file:
        log_path.chmod(0o600)
        process = subprocess.Popen(  # noqa: S603 - fixed local module command
            [sys.executable, "-m", "dash_app.agent_integration.mcp_server", "--http-config", str(config_path)],
            cwd=str(REPO_ROOT),
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            close_fds=True,
        )

    deadline = time.monotonic() + ENDPOINT_START_TIMEOUT_SECONDS
    record: dict[str, Any] | None = None
    while time.monotonic() < deadline:
        if process.poll() is not None:
            break
        candidate = _read_json(_record_path(instance_id))
        if candidate and int(candidate.get("endpoint_pid") or 0) == int(process.pid):
            record = candidate
            break
        time.sleep(0.1)
    if record is None:
        try:
            if process.poll() is None:
                os.kill(process.pid, signal.SIGTERM)
        except OSError:
            pass
        cleanup_instance_record(instance_id, expected_endpoint_pid=process.pid)
        raise RuntimeError(f"dashboard MCP endpoint failed to start; see {log_path}")
    return DashboardEndpointHandle(instance_id, config_path, record, process)
