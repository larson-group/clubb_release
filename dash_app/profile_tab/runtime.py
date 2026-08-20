"""Profile command construction, process ownership, and CSV loading."""

from __future__ import annotations

import os
import shlex
import signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Sequence

import psutil

from utilities.time_clubb import REQUIRED_TIMER, parse_positive_int_list
from utilities.timing_profiles import (
    discover_profiles,
    load_profiles,
    profile_directory,
    read_profile_manifest,
    safe_profile_id,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
TIME_CLUBB = REPO_ROOT / "utilities" / "time_clubb.py"
MAX_LOG_CHARACTERS = 24_000

PROFILE_LOCK = threading.Lock()
PROFILE_PROCESSES: dict[int, subprocess.Popen[bytes]] = {}


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _integer(value: Any, name: str, *, minimum: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if parsed < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return parsed


def _positive_list(value: Any, name: str) -> tuple[int, ...]:
    try:
        return parse_positive_int_list(_clean(value))
    except Exception as exc:
        raise ValueError(f"{name} must be a comma-separated list of unique positive integers") from exc


def _resolved_output(value: Any) -> Path:
    text = _clean(value) or "output/timing"
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def profile_save_target(settings: dict[str, Any]) -> dict[str, Any]:
    """Describe the deterministic library target for a possibly partial form."""
    case_name = _clean(settings.get("case_name"))
    label = _clean(settings.get("name")) or case_name
    output = _resolved_output(settings.get("output"))
    run_id = safe_profile_id(label)
    path = profile_directory(output, run_id)
    exists = path.exists()
    manifest = read_profile_manifest(path) if exists and path.is_dir() and not path.is_symlink() else None
    return {
        "label": label,
        "run_id": run_id,
        "path": str(path),
        "exists": exists,
        "replaceable": manifest is not None,
        "case_name": str((manifest or {}).get("case_name") or ""),
    }


def overwrite_confirmation(settings: dict[str, Any]) -> str:
    """Return a confirmation message, or reject an unsafe name collision."""
    target = profile_save_target(settings)
    if not target["exists"]:
        return ""
    if not target["replaceable"]:
        raise ValueError(
            f"profile name {target['label']!r} conflicts with an unrecognized directory; "
            "choose another name or move that directory"
        )
    return (
        f"Replace the saved profile “{target['label']}”?\n\n"
        "Its CSV timing data and representative logs will be deleted before the new benchmark starts. "
        "This cannot be undone."
    )


def normalize_profile_settings(settings: dict[str, Any]) -> dict[str, Any]:
    """Validate browser settings and return one canonical JSON-safe request."""
    case_name = _clean(settings.get("case_name"))
    if not case_name:
        raise ValueError("select a case")
    case_file = REPO_ROOT / "input" / "case_setups" / f"{case_name}_model.in"
    if not case_file.is_file():
        raise ValueError(f"unknown CLUBB case: {case_name}")

    processes = _positive_list(settings.get("processes"), "processes")
    columns = _positive_list(settings.get("columns"), "batch sizes per process")
    warmups = _integer(settings.get("warmups"), "warmups", minimum=0)
    repetitions = _integer(settings.get("repetitions"), "repetitions", minimum=1)
    timeout_text = _clean(settings.get("timeout"))
    timeout = None
    if timeout_text:
        try:
            timeout = float(timeout_text)
        except ValueError as exc:
            raise ValueError("timeout must be a positive number") from exc
        if timeout <= 0:
            raise ValueError("timeout must be a positive number")

    executable = _clean(settings.get("executable"))
    implementation = _clean(settings.get("implementation")).lower() or "fortran"
    if implementation not in {"fortran", "python", "jax"}:
        raise ValueError("implementation must be Fortran, Python, or JAX")
    if executable and implementation != "fortran":
        raise ValueError("an explicit executable can only be used with Fortran")
    executable_path = Path(executable).expanduser() if executable else None
    if executable_path is not None and not executable_path.is_absolute():
        executable_path = REPO_ROOT / executable_path
    if executable_path is not None and not executable_path.is_file():
        raise ValueError(f"executable not found: {executable}")
    install_dir = _clean(settings.get("install_dir"))
    install_path = Path(install_dir).expanduser() if install_dir else None
    if install_path is not None and not install_path.is_absolute():
        install_path = REPO_ROOT / install_path
    if install_path is not None and not install_path.is_dir():
        raise ValueError(f"install directory not found: {install_dir}")
    try:
        extra_args = shlex.split(_clean(settings.get("extra_args")))
    except ValueError as exc:
        raise ValueError(f"additional run_scm.py arguments are invalid: {exc}") from exc
    managed = (
        "-multicol", "--multicol", "-batch_size", "--batch_size", "-out_dir", "--out_dir",
        "-python", "--python", "-jax", "--jax", "-exe", "--exe", "-install_dir", "--install_dir",
    )
    for token in extra_args:
        if token in managed or any(token.startswith(option + "=") for option in managed):
            raise ValueError(f"{token.split('=', 1)[0]} is managed by the Profile benchmark")

    return {
        "case_name": case_name,
        "processes": list(processes),
        "columns": list(columns),
        "warmups": warmups,
        "repetitions": repetitions,
        "output": str(_resolved_output(settings.get("output"))),
        "name": _clean(settings.get("name")) or case_name,
        "run_id": safe_profile_id(_clean(settings.get("name")) or case_name),
        "overwrite": bool(settings.get("overwrite")),
        "timeout": timeout,
        "continue_on_error": bool(settings.get("continue_on_error")),
        "config": _clean(settings.get("config")) or "default",
        "override": _clean(settings.get("override")),
        "executable": str(executable_path.resolve()) if executable_path is not None else "",
        "implementation": implementation,
        "install_dir": str(install_path.resolve()) if install_path is not None else "",
        "extra_args": extra_args,
    }


def profile_command(settings: dict[str, Any]) -> list[str]:
    """Build a shell-free command for the timing utility."""
    normalized = normalize_profile_settings(settings)
    command = [
        sys.executable,
        "-u",
        str(TIME_CLUBB),
        normalized["case_name"],
        "-processes",
        ",".join(str(value) for value in normalized["processes"]),
        "-batch_sizes",
        ",".join(str(value) for value in normalized["columns"]),
        "-warmups",
        str(normalized["warmups"]),
        "-repeats",
        str(normalized["repetitions"]),
        "-output",
        normalized["output"],
        "-name",
        normalized["name"],
    ]
    if normalized["timeout"] is not None:
        command.extend(("-timeout", str(normalized["timeout"])))
    if normalized["continue_on_error"]:
        command.append("-continue_on_error")
    if normalized["overwrite"]:
        command.append("-overwrite")

    command.extend(("-config", normalized["config"]))
    if normalized["override"]:
        command.extend(("-override", normalized["override"]))
    if normalized["executable"]:
        command.extend(("-exe", normalized["executable"]))
    else:
        if normalized["implementation"] == "python":
            command.append("-python")
        elif normalized["implementation"] == "jax":
            command.append("-jax")
    if normalized["install_dir"] and not normalized["executable"]:
        command.extend(("-install_dir", normalized["install_dir"]))
    command.extend(normalized["extra_args"])
    return command


def profile_command_display(settings: dict[str, Any]) -> str:
    return shlex.join(profile_command(settings))


def summary_run_ids(output_dir: Path) -> list[str]:
    return [
        _clean(record.get("run_id"))
        for record in reversed(discover_profiles(output_dir))
        if _clean(record.get("run_id"))
    ]


def read_profile_data(
    job: dict[str, Any] | None,
) -> tuple[str, list[dict[str, Any]], list[dict[str, Any]]]:
    """Read summary and process rows created by one broker-owned Profile job."""
    job = dict(job or {})
    output_dir = Path(_clean(job.get("output")) or REPO_ROOT / "output" / "timing")
    run_id = _clean(job.get("run_id"))
    if not run_id:
        previous = set(job.get("existing_run_ids") or [])
        candidates = [
            candidate
            for candidate in summary_run_ids(output_dir)
            if candidate not in previous
        ]
        run_id = candidates[-1] if candidates else ""
    if not run_id:
        return "", [], []
    rows, process_rows = load_profiles(output_dir, [run_id])
    return run_id, [dict(row) for row in rows], [dict(row) for row in process_rows]


def read_profile_results(job: dict[str, Any] | None) -> tuple[str, list[dict[str, Any]]]:
    """Compatibility summary-only view of one broker-owned Profile job."""
    run_id, rows, _process_rows = read_profile_data(job)
    return run_id, rows


def read_log_tail(path: str | Path | None, max_characters: int = MAX_LOG_CHARACTERS) -> str:
    try:
        with Path(str(path)).open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - max_characters))
            return handle.read().decode("utf-8", errors="replace")
    except (OSError, TypeError, ValueError):
        return ""


def _process_started_at(pid: int) -> float:
    try:
        return float(psutil.Process(pid).create_time())
    except (psutil.Error, TypeError, ValueError):
        return 0.0


def start_profile_process(settings: dict[str, Any]) -> dict[str, Any]:
    """Launch one timing wrapper and return durable broker metadata."""
    normalized = normalize_profile_settings(settings)
    command = profile_command(settings)
    output_dir = Path(normalized["output"])
    existing_run_ids = summary_run_ids(output_dir)
    log_file = tempfile.NamedTemporaryFile(
        delete=False,
        prefix="clubb_profile_",
        suffix=".log",
        dir="/tmp",
    )
    log_path = log_file.name
    try:
        process = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            env=os.environ.copy(),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    except BaseException:
        log_file.close()
        raise
    log_file.close()
    with PROFILE_LOCK:
        PROFILE_PROCESSES[process.pid] = process
    return {
        "state": "running",
        "pid": process.pid,
        "process_started_at": _process_started_at(process.pid),
        "start_time": time.time(),
        "settings": normalized,
        "command": command,
        "command_display": shlex.join(command),
        "output": normalized["output"],
        "log": log_path,
        "existing_run_ids": existing_run_ids,
        "run_id": normalized["run_id"],
        "result_rows": 0,
    }


def get_profile_process(pid: Any) -> subprocess.Popen[bytes] | None:
    try:
        numeric_pid = int(pid)
    except (TypeError, ValueError):
        return None
    with PROFILE_LOCK:
        return PROFILE_PROCESSES.get(numeric_pid)


def profile_process_status(job: dict[str, Any]) -> tuple[bool, int | None]:
    """Return ``(is_running, returncode)`` for local or recovered process state."""
    process = get_profile_process(job.get("pid"))
    if process is not None:
        return process.poll() is None, process.poll()
    try:
        process_info = psutil.Process(int(job.get("pid") or 0))
        expected_start = float(job.get("process_started_at") or 0.0)
        is_same_process = not expected_start or abs(process_info.create_time() - expected_start) < 0.01
        return bool(is_same_process and process_info.is_running()), None
    except (psutil.Error, TypeError, ValueError):
        return False, None


def release_profile_process(pid: Any) -> None:
    try:
        numeric_pid = int(pid)
    except (TypeError, ValueError):
        return
    with PROFILE_LOCK:
        PROFILE_PROCESSES.pop(numeric_pid, None)


def stop_profile_process(job: dict[str, Any]) -> None:
    """Send SIGINT so time_clubb.py can terminate all active child sessions."""
    running, _returncode = profile_process_status(job)
    if not running:
        raise ValueError("no Profile timing process is running")
    try:
        os.killpg(int(job["pid"]), signal.SIGINT)
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise ValueError("Profile timing process is no longer available") from exc


def expected_required_timer_rows(job: dict[str, Any]) -> int:
    settings = dict(job.get("settings") or {})
    return (
        len(settings.get("processes") or [])
        * len(settings.get("columns") or [])
        * int(settings.get("repetitions") or 0)
    )


def profile_results_complete(job: dict[str, Any], rows: Sequence[dict[str, Any]]) -> bool:
    settings = dict(job.get("settings") or {})
    required_rows = [row for row in rows if _clean(row.get("timer_name")) == REQUIRED_TIMER]
    return (
        len(required_rows) == expected_required_timer_rows(job)
        and all(_clean(row.get("status")) == "success" for row in required_rows)
    )
