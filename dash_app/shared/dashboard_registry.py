"""Broker-owned liveness registry for the one current Dash browser."""

from __future__ import annotations

import json
import time
import uuid
from pathlib import Path
from typing import Any

import psutil

from .runtime import atomic_write_json, exclusive_file_lock, private_path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = private_path(REPO_ROOT, "dashboard.json")
REGISTRY_LOCK_PATH = private_path(REPO_ROOT, "dashboard.lock")
HEARTBEAT_TIMEOUT_SECONDS = 6.0


class DashboardAlreadyRegistered(RuntimeError):
    """Raised when another live Dash process owns this checkout."""


def process_identity_is_live(pid: Any, started_at: Any) -> bool:
    try:
        process = psutil.Process(int(pid))
        expected = float(started_at)
        actual = float(process.create_time())
    except (psutil.Error, TypeError, ValueError, OSError):
        return False
    return abs(actual - expected) <= 1.0


def _read() -> dict[str, Any] | None:
    try:
        payload = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _write(payload: dict[str, Any]) -> None:
    atomic_write_json(REGISTRY_PATH, payload)


def _matches(record: dict[str, Any], pid: int, started_at: float) -> bool:
    try:
        return int(record.get("pid") or 0) == int(pid) and abs(float(record.get("started_at")) - float(started_at)) <= 1.0
    except (TypeError, ValueError):
        return False


def _is_live(record: dict[str, Any] | None, *, now: float | None = None) -> bool:
    if not record:
        return False
    try:
        heartbeat = float(record.get("last_heartbeat"))
    except (TypeError, ValueError):
        return False
    current = time.time() if now is None else float(now)
    return process_identity_is_live(record.get("pid"), record.get("started_at")) and current - heartbeat <= HEARTBEAT_TIMEOUT_SECONDS


def dashboard_status() -> dict[str, Any]:
    """Return a bounded status object suitable for broker/API responses."""
    record = _read()
    if not _is_live(record):
        return {"status": "unavailable", "reason": "no live dashboard is registered"}
    return {
        "status": "available",
        "dashboard_id": record.get("dashboard_id"),
        "pid": record.get("pid"),
        "started_at": record.get("started_at"),
        "port": record.get("port"),
        "last_heartbeat": record.get("last_heartbeat"),
    }


def dashboard_is_live() -> bool:
    return dashboard_status().get("status") == "available"


def register_dashboard(*, pid: int, started_at: float, port: int) -> dict[str, Any]:
    """Register exactly one current Dash process for this checkout."""
    identity = {"pid": int(pid), "started_at": float(started_at)}
    with exclusive_file_lock(REGISTRY_LOCK_PATH):
        existing = _read()
        if _is_live(existing) and not _matches(existing, identity["pid"], identity["started_at"]):
            raise DashboardAlreadyRegistered(
                f"another live dashboard is already registered on port {existing.get('port')}"
            )
        record = {
            "version": 1,
            "dashboard_id": str((existing or {}).get("dashboard_id") or uuid.uuid4().hex),
            **identity,
            "port": int(port),
            "registered_at": time.time(),
            "last_heartbeat": time.time(),
        }
        _write(record)
        return dashboard_status()


def heartbeat_dashboard(*, pid: int, started_at: float) -> dict[str, Any]:
    """Refresh a registration only when it still belongs to the caller."""
    with exclusive_file_lock(REGISTRY_LOCK_PATH):
        record = _read()
        if record is None or not _matches(record, pid, started_at):
            return {"status": "unavailable", "reason": "dashboard registration is no longer current"}
        record["last_heartbeat"] = time.time()
        _write(record)
        return dashboard_status()


def unregister_dashboard(*, pid: int, started_at: float) -> bool:
    """Remove only the caller's record; never remove a replacement dashboard."""
    with exclusive_file_lock(REGISTRY_LOCK_PATH):
        record = _read()
        if record is None or not _matches(record, pid, started_at):
            return False
        try:
            REGISTRY_PATH.unlink()
        except FileNotFoundError:
            pass
        return True


def reconcile_dashboard() -> dict[str, Any]:
    """Clear stale state without touching a live dashboard registration."""
    with exclusive_file_lock(REGISTRY_LOCK_PATH):
        record = _read()
        if record is not None and not _is_live(record):
            try:
                REGISTRY_PATH.unlink()
            except FileNotFoundError:
                pass
    return dashboard_status()
