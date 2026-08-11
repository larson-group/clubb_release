"""Private heartbeat tying the dashboard broker to its foreground manager."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import psutil

from .runtime import atomic_write_json, private_path


REPO_ROOT = Path(__file__).resolve().parents[2]
MANAGER_LEASE_PATH = private_path(REPO_ROOT, "manager.json")
MANAGER_LOCK_PATH = private_path(REPO_ROOT, "manager.lock")
MANAGER_REQUIRED_ENV = "CLUBB_DASH_MANAGER_REQUIRED"
HEARTBEAT_INTERVAL_SECONDS = 2.0
LEASE_TIMEOUT_SECONDS = 30.0


def process_started_at(pid: int) -> float:
    """Return a stable process identity so a recycled PID is never adopted."""
    return float(psutil.Process(int(pid)).create_time())


def write_manager_lease(*, pid: int, started_at: float) -> dict[str, Any]:
    payload = {
        "version": 1,
        "pid": int(pid),
        "started_at": float(started_at),
        "updated_at": time.time(),
        "repository": str(REPO_ROOT),
    }
    atomic_write_json(MANAGER_LEASE_PATH, payload)
    return payload


def read_manager_lease() -> dict[str, Any] | None:
    try:
        payload = json.loads(MANAGER_LEASE_PATH.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def manager_lease_is_live(
    record: dict[str, Any] | None = None,
    *,
    now: float | None = None,
) -> bool:
    record = read_manager_lease() if record is None else record
    if not record:
        return False
    try:
        pid = int(record["pid"])
        started_at = float(record["started_at"])
        updated_at = float(record["updated_at"])
        process = psutil.Process(pid)
        same_process = abs(float(process.create_time()) - started_at) <= 1.0
    except (KeyError, OSError, TypeError, ValueError, psutil.Error):
        return False
    current = time.time() if now is None else float(now)
    return same_process and current - updated_at <= LEASE_TIMEOUT_SECONDS


def clear_manager_lease(*, pid: int, started_at: float) -> None:
    """Remove only this manager's heartbeat, never a replacement's."""
    record = read_manager_lease()
    try:
        matches = (
            record is not None
            and int(record.get("pid") or 0) == int(pid)
            and abs(float(record.get("started_at")) - float(started_at)) <= 1.0
        )
    except (TypeError, ValueError):
        matches = False
    if matches:
        try:
            MANAGER_LEASE_PATH.unlink()
        except FileNotFoundError:
            pass
