"""Stable MCP endpoint owned by the durable local broker.

The endpoint identity is checkout-scoped rather than Dash-process-scoped.  Its
URL and bearer credential are persisted in the private runtime directory and
the broker refreshes only the endpoint's broker credentials when it restarts.
"""

from __future__ import annotations

import json
import os
import secrets
import signal
import socket
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

import psutil

from dash_app.shared.activity import REPO_ROOT
from dash_app.shared.runtime import atomic_write_json, exclusive_file_lock, private_path


ENDPOINT_SCHEMA_VERSION = 2
ENDPOINT_PATH = "/mcp"
ENDPOINT_ROOT = private_path(REPO_ROOT, "mcp-endpoint")
IDENTITY_PATH = ENDPOINT_ROOT / "identity.json"
CONFIG_PATH = ENDPOINT_ROOT / "endpoint-config.json"
RECORD_PATH = ENDPOINT_ROOT / "endpoint.json"
LOCK_PATH = ENDPOINT_ROOT / "endpoint.lock"
ENDPOINT_START_TIMEOUT_SECONDS = 8.0


def process_started_at(pid: int) -> float:
    try:
        return float(psutil.Process(int(pid)).create_time())
    except (psutil.Error, TypeError, ValueError, OSError) as exc:
        raise RuntimeError("could not determine the endpoint process start time") from exc


def process_identity_is_live(pid: Any, started_at: Any) -> bool:
    try:
        process = psutil.Process(int(pid))
        expected = float(started_at)
        actual = float(process.create_time())
    except (psutil.Error, TypeError, ValueError, OSError):
        return False
    return abs(actual - expected) <= 1.0


def _read(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _local_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind(("127.0.0.1", 0))
        return int(listener.getsockname()[1])


def _unlink(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        pass


def _terminate_exact(record: dict[str, Any], *, timeout: float = 3.0) -> bool:
    pid = record.get("endpoint_pid")
    if not process_identity_is_live(pid, record.get("endpoint_started_at")):
        return True
    try:
        os.kill(int(pid), signal.SIGTERM)
    except (OSError, TypeError, ValueError):
        return False
    deadline = time.monotonic() + max(0.1, float(timeout))
    while process_identity_is_live(pid, record.get("endpoint_started_at")) and time.monotonic() < deadline:
        time.sleep(0.05)
    return not process_identity_is_live(pid, record.get("endpoint_started_at"))


def _reconcile_legacy_records() -> None:
    """Remove rotating endpoint records and their exact child processes."""
    legacy_root = private_path(REPO_ROOT, "mcp-instances")
    for record_path in sorted(legacy_root.glob("*/endpoint.json")):
        record = _read(record_path)
        if record:
            _terminate_exact(record)
        _unlink(record_path)
        _unlink(record_path.parent / "endpoint-config.json")
        try:
            record_path.parent.rmdir()
        except OSError:
            pass


def _identity() -> dict[str, Any]:
    current = _read(IDENTITY_PATH) or {}
    try:
        instance_id = str(current["instance_id"])
        token = str(current["endpoint_token"])
        port = int(current["endpoint_port"])
        if len(instance_id) != 32 or not token or not 1 < port < 65536:
            raise ValueError
    except (KeyError, TypeError, ValueError):
        current = {
            "version": ENDPOINT_SCHEMA_VERSION,
            "instance_id": uuid.uuid4().hex,
            "endpoint_token": secrets.token_urlsafe(32),
            "endpoint_port": _local_port(),
            "endpoint_path": ENDPOINT_PATH,
        }
        atomic_write_json(IDENTITY_PATH, current)
    return current


def stable_endpoint_record_is_live(record: dict[str, Any], *, check_broker: bool = True) -> bool:
    try:
        if int(record.get("version") or 0) != ENDPOINT_SCHEMA_VERSION:
            return False
        if not record.get("instance_id") or not process_identity_is_live(record.get("endpoint_pid"), record.get("endpoint_started_at")):
            return False
        owner = dict(record.get("owner") or {})
        if not process_identity_is_live(owner.get("pid"), owner.get("started_at")):
            return False
        broker = dict(record.get("broker") or {})
        if not broker.get("url") or not broker.get("token"):
            return False
        if not check_broker:
            return True
        request = Request(
            str(broker["url"]).rstrip("/") + "/status",
            headers={"X-CLUBB-Agent-Token": str(broker["token"])},
            method="GET",
        )
        with urlopen(request, timeout=0.75) as response:  # nosec B310: private loopback broker record
            return response.status == 200
    except (KeyError, OSError, TypeError, ValueError, URLError):
        return False


def cleanup_stable_record(*, expected_endpoint_pid: int | None = None) -> None:
    record = _read(RECORD_PATH)
    if expected_endpoint_pid is not None and record is not None:
        try:
            if int(record.get("endpoint_pid")) != int(expected_endpoint_pid):
                return
        except (TypeError, ValueError):
            return
    _unlink(RECORD_PATH)


def _write_config(identity: dict[str, Any], broker: dict[str, Any], owner: dict[str, Any]) -> None:
    config = {
        "version": ENDPOINT_SCHEMA_VERSION,
        "instance_id": identity["instance_id"],
        "endpoint_path": ENDPOINT_PATH,
        "endpoint_host": "127.0.0.1",
        "endpoint_port": int(identity["endpoint_port"]),
        "endpoint_token": identity["endpoint_token"],
        "owner": owner,
        "broker": {
            "url": str(broker["url"]),
            "token": str(broker["token"]),
            "pid": broker.get("pid"),
            "version": broker.get("version"),
            "repository": broker.get("repository"),
        },
        "registry_path": str(RECORD_PATH),
    }
    atomic_write_json(CONFIG_PATH, config)


def start_broker_endpoint(broker: dict[str, Any], *, owner_pid: int, owner_started_at: float) -> dict[str, Any]:
    """Start or reuse the one stable endpoint for this checkout."""
    ENDPOINT_ROOT.mkdir(mode=0o700, parents=True, exist_ok=True)
    ENDPOINT_ROOT.chmod(0o700)
    with exclusive_file_lock(LOCK_PATH):
        _reconcile_legacy_records()
        identity = _identity()
        old_record = _read(RECORD_PATH)
        if old_record and stable_endpoint_record_is_live(old_record, check_broker=False):
            owner = dict(old_record.get("owner") or {})
            if int(owner.get("pid") or 0) == int(owner_pid) and abs(float(owner.get("started_at")) - float(owner_started_at)) <= 1.0:
                return {
                    "instance_id": identity["instance_id"],
                    "endpoint_url": old_record.get("endpoint_url"),
                    "endpoint_token": identity["endpoint_token"],
                    "endpoint_port": identity["endpoint_port"],
                    "endpoint_log_path": str(ENDPOINT_ROOT / "endpoint.log"),
                }
            _terminate_exact(old_record)
        cleanup_stable_record()
        _write_config(identity, broker, {"pid": int(owner_pid), "started_at": float(owner_started_at)})
        log_path = ENDPOINT_ROOT / "endpoint.log"
        with log_path.open("a", encoding="utf-8") as log_file:
            log_path.chmod(0o600)
            process = subprocess.Popen(  # noqa: S603 - fixed local module command
                [sys.executable, "-m", "dash_app.agent_integration.mcp_server", "--http-config", str(CONFIG_PATH)],
                cwd=str(REPO_ROOT),
                stdin=subprocess.DEVNULL,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                close_fds=True,
            )
        deadline = time.monotonic() + ENDPOINT_START_TIMEOUT_SECONDS
        record = None
        while time.monotonic() < deadline:
            if process.poll() is not None:
                break
            candidate = _read(RECORD_PATH)
            if candidate and int(candidate.get("endpoint_pid") or 0) == int(process.pid):
                record = candidate
                break
            time.sleep(0.1)
        if record is None:
            _terminate_exact({"endpoint_pid": process.pid, "endpoint_started_at": process_started_at(process.pid) if process.poll() is None else 0})
            cleanup_stable_record(expected_endpoint_pid=process.pid)
            raise RuntimeError(f"stable broker MCP endpoint failed to start; see {log_path}")
        return {
            "instance_id": identity["instance_id"],
            "endpoint_url": record["endpoint_url"],
            "endpoint_token": identity["endpoint_token"],
            "endpoint_port": identity["endpoint_port"],
            "endpoint_log_path": str(ENDPOINT_ROOT / "endpoint.log"),
        }


def stop_broker_endpoint() -> None:
    with exclusive_file_lock(LOCK_PATH):
        record = _read(RECORD_PATH)
        if record:
            _terminate_exact(record)
        cleanup_stable_record()
