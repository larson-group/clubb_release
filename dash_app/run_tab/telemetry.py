"""Bounded browser telemetry projected from canonical broker SCM records."""

from __future__ import annotations

import base64
from pathlib import Path
from typing import Any

from dash_app.services.jobs import JobStore
from dash_app.shared.runtime import (
    read_file_chunk,
    read_file_tail,
    read_latest_run_progress,
)

from .runtime import build_case_command
from .state import DEFAULT_STATS_NAME


JOB_STORE = JobStore()


def normalized_scm_state(value: Any) -> str:
    state = str(value or "")
    if state in {"submitting", "starting"}:
        return "queued"
    if state in {"stopped", "rejected", "partial_failure"}:
        return "cancelled" if state == "stopped" else "error"
    return state if state in {
        "queued", "running", "stopping", "finished", "error", "cancelled"
    } else "error"


def scm_run_view(record: dict[str, Any]) -> dict[str, Any]:
    """Project one durable SCM child into the browser-safe Run contract."""
    request = dict(record.get("request") or {})
    runtime = dict(record.get("runtime") or {})
    proc = dict(runtime.get("proc_data") or {})
    display = dict(record.get("display") or {})
    error = record.get("error") or {}
    case = str(request.get("case") or display.get("case") or "")
    stats_file = str(
        display.get("stats_file") or request.get("stats_file") or DEFAULT_STATS_NAME
    )
    config = str(display.get("config") or request.get("config") or "default")
    cli_options = dict(display.get("cli_options") or runtime.get("cli_options") or {})
    raw_output_directory = (
        record.get("output_directory")
        or runtime.get("output_directory")
        or proc.get("output_directory")
        or cli_options.get("out_dir")
        or ""
    )
    output_directory = (
        str(Path(str(raw_output_directory)).expanduser().resolve())
        if raw_output_directory
        else ""
    )
    return {
        "run_id": str(record.get("run_id") or ""),
        "job_id": str(record.get("job_id") or ""),
        "case": case,
        "state": normalized_scm_state(record.get("state")),
        "origin": str(record.get("origin") or "legacy"),
        "batch_id": record.get("batch_id"),
        "batch_job_id": record.get("batch_job_id"),
        "started_at": proc.get("start_time"),
        "finished_at": record.get("finished_at_unix_seconds"),
        "returncode": record.get("returncode"),
        "message": (
            str(error.get("message") or "")
            if isinstance(error, dict)
            else str(error)
        ),
        "stats_file": stats_file,
        "config": config,
        "output_directory": output_directory,
        "command": build_case_command(case, stats_file, cli_options, config) if case else "",
        "log_available": bool(proc.get("log")),
        "updated_at": float(
            record.get("updated_at_unix_seconds")
            or record.get("created_at_unix_seconds")
            or 0
        ),
    }


def run_telemetry(
    known_revision: Any = None,
    log_cursors: dict[str, Any] | None = None,
    *,
    store: JobStore | None = None,
) -> dict[str, Any]:
    """Return one bounded lifecycle/progress/log observation for the Run tab."""
    revision, records = (store or JOB_STORE).scm_snapshot()
    visible = [
        record
        for record in records
        if str(record.get("visibility") or "user") == "user"
    ]
    by_run_id = {
        str(record.get("run_id") or ""): record
        for record in visible
        if record.get("run_id")
    }
    try:
        unchanged = int(known_revision) == revision
    except (TypeError, ValueError):
        unchanged = False

    progress: dict[str, dict[str, int]] = {}
    for run_id, record in by_run_id.items():
        if normalized_scm_state(record.get("state")) not in {
            "queued", "running", "stopping"
        }:
            continue
        runtime = dict(record.get("runtime") or {})
        log_path = str((runtime.get("proc_data") or {}).get("log") or "")
        sample = read_latest_run_progress(log_path)
        if sample is not None:
            progress[run_id] = {"iteration": sample[0], "total": sample[1]}

    logs: dict[str, dict[str, Any]] = {}
    deferred: list[str] = []
    remaining = 512 * 1024
    for run_id, raw_cursor in list(dict(log_cursors or {}).items())[:64]:
        record = by_run_id.get(str(run_id))
        if record is None:
            continue
        if remaining <= 0:
            deferred.append(str(run_id))
            continue
        runtime = dict(record.get("runtime") or {})
        log_path = str((runtime.get("proc_data") or {}).get("log") or "")
        if raw_cursor is None:
            chunk, next_cursor, eof = read_file_tail(log_path, 5000)
        else:
            try:
                cursor = max(0, int(raw_cursor))
            except (OSError, TypeError, ValueError):
                cursor = 0
            chunk, next_cursor, eof = read_file_chunk(
                log_path, cursor, min(65536, remaining)
            )
        remaining -= len(chunk)
        logs[str(run_id)] = {
            "data": base64.b64encode(chunk).decode("ascii"),
            "next_cursor": next_cursor,
            "eof": eof,
            "state": normalized_scm_state(record.get("state")),
        }

    payload: dict[str, Any] = {
        "revision": revision,
        "progress": progress,
        "logs": logs,
        "deferred_logs": deferred,
        "active": any(
            normalized_scm_state(record.get("state"))
            in {"queued", "running", "stopping"}
            for record in visible
        ),
    }
    if not unchanged:
        payload["runs"] = sorted(
            (scm_run_view(record) for record in visible),
            key=lambda item: (item["updated_at"], item["run_id"]),
        )
    return payload
