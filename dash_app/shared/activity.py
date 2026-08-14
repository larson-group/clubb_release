"""Durable local dashboard activity and broker-owned job state."""

from __future__ import annotations

import fcntl
import hashlib
import json
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from dash_app.shared.runtime import atomic_write_json, private_path, restrict_existing


REPO_ROOT = Path(__file__).resolve().parents[2]
_TOKEN = hashlib.sha256(str(REPO_ROOT).encode("utf-8")).hexdigest()[:16]
_LEGACY_ACTIVITY_PATH = Path("/tmp") / f"clubb_dash_agent_{_TOKEN}.json"
_LEGACY_LOCK_PATH = Path("/tmp") / f"clubb_dash_agent_{_TOKEN}.lock"
ACTIVITY_PATH = private_path(REPO_ROOT, "activity.json")
LOCK_PATH = private_path(REPO_ROOT, "activity.lock")
restrict_existing(_LEGACY_ACTIVITY_PATH)
restrict_existing(_LEGACY_LOCK_PATH)
MAX_EVENTS = 120
UI_HANDOFF_RETRY_SECONDS = 30.0
ACTIVE_JOB_STATES = frozenset({"queued", "submitting", "running", "stopping"})


def _initial_state() -> dict[str, Any]:
    return {
        "version": 7,
        "next_id": 1,
        "events": [],
        "plot_request": None,
        "plot_instances": None,
        "ui_request": None,
        "ui_request_queue": [],
        "ui_request_in_flight": None,
        "ui_request_ack": 0,
        "broker": {},
        "jobs": {
            "compile": None,
            "profile": None,
            "runs": {},
            "tune": None,
            "loss_runs": {},
        },
    }


def _read_state() -> dict[str, Any]:
    try:
        payload = json.loads(ACTIVITY_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        # Preserve a currently running older broker across the first upgrade.
        try:
            payload = json.loads(_LEGACY_ACTIVITY_PATH.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return _initial_state()
    if not isinstance(payload, dict):
        return _initial_state()
    defaults = _initial_state()
    for obsolete in ("next_message_id", "agents", "messages"):
        payload.pop(obsolete, None)
    payload["version"] = defaults["version"]
    for key, value in defaults.items():
        payload.setdefault(key, value)
    jobs = dict(payload.get("jobs") or {})
    for key, value in defaults["jobs"].items():
        jobs.setdefault(key, value)
    payload["jobs"] = jobs
    return payload


def read_activity() -> dict[str, Any]:
    """Return a browser-safe snapshot of local dashboard state."""
    return _read_state()


def broker_jobs() -> dict[str, Any]:
    """Return the broker-owned job snapshot used to rehydrate a Dash view."""
    jobs = read_activity().get("jobs") or {}
    return {
        "compile": dict(jobs.get("compile") or {}) or None,
        "profile": dict(jobs.get("profile") or {}) or None,
        "runs": {str(name): dict(value or {}) for name, value in (jobs.get("runs") or {}).items()},
        "tune": dict(jobs.get("tune") or {}) or None,
        "loss_runs": {str(name): dict(value or {}) for name, value in (jobs.get("loss_runs") or {}).items()},
    }


def count_active_jobs(jobs: Any) -> int:
    """Count active broker job records without depending on job categories."""
    if isinstance(jobs, dict):
        if str(jobs.get("state") or "") in ACTIVE_JOB_STATES:
            return 1
        return sum(count_active_jobs(value) for value in jobs.values())
    if isinstance(jobs, (list, tuple)):
        return sum(count_active_jobs(value) for value in jobs)
    return 0


@contextmanager
def _locked_state() -> Iterator[dict[str, Any]]:
    LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a+", encoding="utf-8") as lock_file:
        try:
            LOCK_PATH.chmod(0o600)
        except OSError:
            pass
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        state = _read_state()
        try:
            yield state
        finally:
            atomic_write_json(ACTIVITY_PATH, state)
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _append_event(state: dict[str, Any], kind: str, message: str, detail: str, status: str) -> dict[str, Any]:
    event_id = int(state.get("next_id") or 1)
    event = {
        "id": event_id,
        "time": datetime.now(timezone.utc).astimezone().strftime("%H:%M:%S"),
        "kind": str(kind),
        "message": str(message),
        "detail": str(detail or ""),
        "status": str(status),
    }
    state["next_id"] = event_id + 1
    state["events"] = [*(state.get("events") or []), event][-MAX_EVENTS:]
    return event


def publish_event(kind: str, message: str, detail: str = "", *, status: str = "info", action: dict[str, Any] | None = None) -> dict[str, Any]:
    """Append one visible activity item and optionally attach a Dash UI action."""
    with _locked_state() as state:
        event = _append_event(state, kind, message, detail, status)
        if action is not None:
            request = {"id": event["id"], **dict(action)}
            state["ui_request"] = request
            queue = [
                dict(item)
                for item in (state.get("ui_request_queue") or [])
                if isinstance(item, dict)
            ]
            queue.append(request)
            state["ui_request_queue"] = queue[-MAX_EVENTS:]
            if request.get("type") in {"profile", "budget", "plot"}:
                state["plot_request"] = request
        return event


def publish_plot_request(
    case: str,
    variables: list[str],
    *,
    output_dir: str = "output",
    time_seconds: float | None = None,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
    benchmark_sources: list[str] | None = None,
) -> dict[str, Any]:
    """Ask the active Dash browser to open a selected profile view."""
    action: dict[str, Any] = {
        "type": "profile",
        "operation": "set_view",
        "case": str(case),
        "variables": [str(value) for value in variables if str(value).strip()],
        "output_dir": str(output_dir),
    }
    # ``time_seconds`` was the original one-control API. Retain it for older
    # adapters, but use the explicit start-time field in all new requests.
    resolved_start = time_start_seconds if time_start_seconds is not None else time_seconds
    if resolved_start is not None:
        action["time_start_seconds"] = float(resolved_start)
    if average_minutes is not None:
        action["average_minutes"] = float(average_minutes)
    if window_preset:
        action["window_preset"] = str(window_preset)
    if benchmark_sources is not None:
        action["benchmark_sources"] = [str(value) for value in benchmark_sources]
    return publish_event(
        "plot",
        f"Opening {case} profile plot{'s' if len(action['variables']) != 1 else ''}",
        ", ".join(action["variables"]),
        status="running",
        action={**action, "tab": "plots"},
    )


def publish_budget_request(
    case: str,
    budget_group: str,
    *,
    output_dir: str = "output",
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
) -> dict[str, Any]:
    """Ask the active Dash browser to add one validated budget plot card."""
    action: dict[str, Any] = {
        "type": "budget",
        "operation": "add_budget",
        "case": str(case),
        "budget_group": str(budget_group),
        "output_dir": str(output_dir),
    }
    if time_start_seconds is not None:
        action["time_start_seconds"] = float(time_start_seconds)
    if average_minutes is not None:
        action["average_minutes"] = float(average_minutes)
    if window_preset:
        action["window_preset"] = str(window_preset)
    return publish_event(
        "plot",
        f"Adding {budget_group} budget plot for {case}",
        str(budget_group),
        status="running",
        action={**action, "tab": "plots"},
    )


def publish_plot_remove_request(plot_id: int, *, case: str = "") -> dict[str, Any]:
    """Ask the active Plot tab to remove one currently mounted card."""
    action: dict[str, Any] = {
        "type": "plot",
        "operation": "remove",
        "plot_id": int(plot_id),
        "tab": "plots",
    }
    if case:
        action["case"] = str(case)
    return publish_event(
        "plot",
        f"Removing Plot card {int(plot_id)}",
        str(case or ""),
        status="running",
        action=action,
    )


def set_plot_instances(snapshot: dict[str, Any] | None) -> dict[str, Any] | None:
    """Persist the sanitized current Plot-card snapshot for typed inspection."""
    with _locked_state() as state:
        value = None if snapshot is None else dict(snapshot)
        state["plot_instances"] = value
        return None if value is None else dict(value)


def claim_ui_request(last_action_id: int | None = 0) -> dict[str, Any] | None:
    """Claim the next queued UI handoff after ``last_action_id``.

    A broker can publish several typed UI actions between two Dash polling
    ticks.  Claiming one item at a time preserves those actions while making
    the poll itself the delivery acknowledgment.  Older activity files that
    only contain ``ui_request`` remain readable through the fallback below.
    """
    try:
        seen_id = int(last_action_id or 0)
    except (TypeError, ValueError):
        seen_id = 0
    with _locked_state() as state:
        try:
            acknowledged_id = int(state.get("ui_request_ack") or 0)
        except (TypeError, ValueError):
            acknowledged_id = 0
        in_flight = state.get("ui_request_in_flight")
        if isinstance(in_flight, dict):
            try:
                in_flight_id = int(in_flight.get("id") or 0)
            except (TypeError, ValueError):
                in_flight_id = 0
            if in_flight_id > acknowledged_id:
                # A new Dash page can retry a request left behind by a
                # crashed page; an existing page must wait for its render
                # acknowledgment before receiving the next request.
                try:
                    claimed_at = float(in_flight.get("_claimed_at") or 0.0)
                except (TypeError, ValueError):
                    claimed_at = 0.0
                if seen_id == 0 and claimed_at and time.time() - claimed_at >= UI_HANDOFF_RETRY_SECONDS:
                    retry_item = {key: value for key, value in in_flight.items() if key != "_claimed_at"}
                    queue = [retry_item] + [
                        dict(item)
                        for item in (state.get("ui_request_queue") or [])
                        if isinstance(item, dict)
                    ]
                    state["ui_request_queue"] = queue[-MAX_EVENTS:]
                    state["ui_request_in_flight"] = None
                else:
                    return None
        queue = []
        candidate = None
        for item in state.get("ui_request_queue") or []:
            if not isinstance(item, dict):
                continue
            try:
                request_id = int(item.get("id") or 0)
            except (TypeError, ValueError):
                continue
            if request_id <= seen_id:
                continue
            if candidate is None:
                candidate = dict(item)
            else:
                queue.append(dict(item))
        state["ui_request_queue"] = queue
        if candidate is not None:
            state["ui_request_in_flight"] = {**candidate, "_claimed_at": time.time()}
            return candidate

        # Compatibility with activity written by an older broker before the
        # queue field existed.  The current request is only delivered once
        # because the Dash-side last-action store advances to its id.
        latest = state.get("ui_request")
        if isinstance(latest, dict):
            try:
                latest_id = int(latest.get("id") or 0)
            except (TypeError, ValueError):
                latest_id = 0
            if latest_id > max(seen_id, acknowledged_id):
                state["ui_request_in_flight"] = {**latest, "_claimed_at": time.time()}
                return dict(latest)
        return None


def acknowledge_ui_request(request_id: int | None) -> bool:
    """Acknowledge one handoff after its target UI state has been applied."""
    try:
        acknowledged_id = int(request_id or 0)
    except (TypeError, ValueError):
        return False
    if acknowledged_id < 1:
        return False
    with _locked_state() as state:
        in_flight = state.get("ui_request_in_flight")
        if not isinstance(in_flight, dict):
            return False
        try:
            in_flight_id = int(in_flight.get("id") or 0)
        except (TypeError, ValueError):
            return False
        if in_flight_id != acknowledged_id:
            return False
        state["ui_request_ack"] = acknowledged_id
        state["ui_request_in_flight"] = None
        return True


def publish_tab_request(tab: str, message: str, detail: str = "", **payload: Any) -> dict[str, Any]:
    """Navigate the active browser to one named dashboard tab."""
    return publish_event(
        "navigate",
        message,
        detail,
        status="running",
        action={"type": "tab", "tab": str(tab), "operation": "navigate", **payload},
    )


def publish_run_request(
    case: str,
    proc_data: dict[str, Any],
    *,
    stats_file: str,
    config: str,
    cli_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Ask the Run tab to adopt an already-started normal SCM process."""
    return publish_event(
        "run",
        f"Showing {case} in Run",
        str(proc_data.get("log") or ""),
        status="running",
        action={
            "type": "run",
            "tab": "run",
            "operation": "start",
            "case": case,
            "proc_data": dict(proc_data),
            "stats_file": stats_file,
            "config": config,
            "cli_options": dict(cli_options or {}),
        },
    )


def publish_tune_request(
    request: dict[str, Any],
    job: dict[str, Any],
    *,
    message: str,
    detail: str = "",
    status: str = "running",
    operation: str = "launch",
) -> dict[str, Any]:
    """Ask the Tune tab to show one canonical tuning request and its worker."""
    return publish_event(
        "tune",
        message,
        detail,
        status=status,
        action={"type": "tune", "tab": "tune", "operation": str(operation), "request": dict(request), "job": dict(job)},
    )


def set_broker_metadata(**metadata: Any) -> dict[str, Any]:
    """Persist broker identity without treating a Dash reload as a new session."""
    with _locked_state() as state:
        broker = dict(state.get("broker") or {})
        broker.update(metadata)
        broker["updated_at"] = time.time()
        state["broker"] = broker
        return dict(broker)


def set_broker_job(kind: str, payload: dict[str, Any] | None) -> dict[str, Any] | None:
    """Store one JSON-safe broker job record for a browser that reconnects later."""
    if kind not in {"compile", "profile", "tune"}:
        raise ValueError("broker job kind must be compile, profile, or tune")
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        record = None if payload is None else dict(payload)
        if record is not None:
            record["updated_at"] = time.time()
        jobs[kind] = record
        jobs.setdefault("runs", {})
        jobs.setdefault("loss_runs", {})
        jobs.setdefault("profile", None)
        state["jobs"] = jobs
        return None if record is None else dict(record)


def set_broker_run(case: str, payload: dict[str, Any] | None) -> dict[str, Any] | None:
    """Store one broker-owned SCM process record keyed by its checked-in case."""
    name = str(case or "").strip()
    if not name:
        raise ValueError("broker run requires a case name")
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        runs = {str(key): dict(value or {}) for key, value in (jobs.get("runs") or {}).items()}
        if payload is None:
            runs.pop(name, None)
            record = None
        else:
            record = dict(payload)
            record["updated_at"] = time.time()
            runs[name] = record
        jobs["runs"] = runs
        jobs.setdefault("compile", None)
        jobs.setdefault("profile", None)
        jobs.setdefault("tune", None)
        jobs.setdefault("loss_runs", {})
        state["jobs"] = jobs
        return None if record is None else dict(record)


def update_broker_job(kind: str, **updates: Any) -> dict[str, Any] | None:
    """Merge status/log updates into an existing singular broker job."""
    if kind not in {"compile", "profile", "tune"}:
        raise ValueError("broker job kind must be compile, profile, or tune")
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        current = jobs.get(kind)
        if not isinstance(current, dict):
            return None
        current = dict(current)
        current.update(updates)
        current["updated_at"] = time.time()
        jobs[kind] = current
        jobs.setdefault("runs", {})
        jobs.setdefault("loss_runs", {})
        jobs.setdefault("profile", None)
        state["jobs"] = jobs
        return dict(current)


def update_broker_run(case: str, **updates: Any) -> dict[str, Any] | None:
    """Merge status/log updates into one existing broker-owned SCM run."""
    name = str(case or "").strip()
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        runs = {str(key): dict(value or {}) for key, value in (jobs.get("runs") or {}).items()}
        current = runs.get(name)
        if current is None:
            return None
        current.update(updates)
        current["updated_at"] = time.time()
        runs[name] = current
        jobs["runs"] = runs
        jobs.setdefault("compile", None)
        jobs.setdefault("profile", None)
        jobs.setdefault("tune", None)
        jobs.setdefault("loss_runs", {})
        state["jobs"] = jobs
        return dict(current)


def set_broker_loss_run(run_id: str, payload: dict[str, Any] | None) -> dict[str, Any] | None:
    """Persist a broker-owned Tune result run for status and reconnect views."""
    name = str(run_id or "").strip()
    if not name:
        raise ValueError("broker loss run requires a run id")
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        runs = {str(key): dict(value or {}) for key, value in (jobs.get("loss_runs") or {}).items()}
        if payload is None:
            runs.pop(name, None)
            record = None
        else:
            record = dict(payload)
            record["updated_at"] = time.time()
            runs[name] = record
        jobs["loss_runs"] = runs
        jobs.setdefault("compile", None)
        jobs.setdefault("profile", None)
        jobs.setdefault("runs", {})
        jobs.setdefault("tune", None)
        state["jobs"] = jobs
        return None if record is None else dict(record)


def update_broker_loss_run(run_id: str, **updates: Any) -> dict[str, Any] | None:
    """Merge a completion or log update into one Tune result run."""
    name = str(run_id or "").strip()
    with _locked_state() as state:
        jobs = dict(state.get("jobs") or {})
        runs = {str(key): dict(value or {}) for key, value in (jobs.get("loss_runs") or {}).items()}
        current = runs.get(name)
        if current is None:
            return None
        current.update(updates)
        current["updated_at"] = time.time()
        runs[name] = current
        jobs["loss_runs"] = runs
        jobs.setdefault("compile", None)
        jobs.setdefault("profile", None)
        jobs.setdefault("runs", {})
        jobs.setdefault("tune", None)
        state["jobs"] = jobs
        return dict(current)


def reset_activity() -> None:
    """Clear the local stream; useful in focused tests and manual resets."""
    with _locked_state() as state:
        state.clear()
        state.update(_initial_state())
