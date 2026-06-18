"""Process launching and status polling helpers for the tuning tab."""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime

from .state import (
    LOSS_RUN_LOCK,
    LOSS_RUN_PROCS,
    OUTPUT_TUNER_DIR,
    REPO_ROOT,
    TUNE_ACTIVE,
    TUNE_LOCK,
    TUNE_STATUS_TEMPLATE,
)
from tuner.job_runtime import TunerJob, tuner_worker_env


DEFAULT_TUNABLE_PARAMS = Path(REPO_ROOT) / "input" / "tunable_parameters" / "tunable_parameters.in"
CLUBB_OUTPUT_DIR = Path(REPO_ROOT) / "output"


class ProcessGroup:
    """Small poll-compatible wrapper for multiple subprocesses launched as one UI run."""

    def __init__(self, processes):
        self.processes = list(processes)
        self.pid = self.processes[0].pid if self.processes else None

    def poll(self):
        returncodes = [proc.poll() for proc in self.processes]
        if any(returncode is None for returncode in returncodes):
            return None
        failures = [returncode for returncode in returncodes if returncode]
        return int(failures[0]) if failures else 0


def empty_status_payload():
    """Return a fresh empty tuning-status payload."""
    return dict(TUNE_STATUS_TEMPLATE)


def empty_results_payload():
    """Return a fresh empty retained-results payload."""
    return {"best_results": [], "best_results_by_case": {}}


def worker_env():
    """Build the environment for tuning subprocesses."""
    return tuner_worker_env()


def active_tuning_job(job_data):
    """Return a TunerJob handle for a JSON-safe active-job payload."""
    return TunerJob.from_dict(job_data)


def start_tuning_job(request_payload):
    """Launch the background tuning worker and return tracked job metadata."""
    with TUNE_LOCK:
        existing_proc = TUNE_ACTIVE.get("proc")
        if existing_proc is not None and existing_proc.poll() is None:
            raise RuntimeError("A tuning job is already running")
        TUNE_ACTIVE.clear()

    initial_status = empty_status_payload()
    job = TunerJob.create(
        request_payload,
        output_root=OUTPUT_TUNER_DIR,
        initial_status=initial_status,
        initial_state=initial_status.get("state", "idle"),
    )
    job.start()
    job_data = job.to_dict()
    with TUNE_LOCK:
        TUNE_ACTIVE.update(job_data)
        TUNE_ACTIVE["proc"] = job.proc
    return job_data


def _read_default_tunable_params():
    """Return scalar CLUBB tunable parameter defaults in file order."""
    params = []
    pattern = re.compile(r"^\s*([A-Za-z]\w*)\s*=\s*([^!,/]+)")
    with open(DEFAULT_TUNABLE_PARAMS, encoding="utf-8") as src:
        for raw_line in src:
            line = raw_line.split("!", 1)[0].strip()
            if not line or line.startswith("&") or line.startswith("/"):
                continue
            match = pattern.match(line)
            if not match:
                continue
            name = match.group(1)
            try:
                value = float(match.group(2).strip().replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            params.append((name, value))
    return params


def write_loss_params_file(work_dir, param_sets, filename="loss_params.in"):
    """Write a params namelist whose only changes are the supplied tuned parameters."""
    defaults = _read_default_tunable_params()
    default_names = {name for name, _value in defaults}
    override_names = sorted({name for params in param_sets for name in (params or {})})
    missing_names = [name for name in override_names if name not in default_names]
    if missing_names:
        raise RuntimeError("Tuned parameter(s) not found in default params file: " + ", ".join(missing_names))

    ngrdcol = len(param_sets)
    params_path = Path(work_dir) / filename
    with open(params_path, "w", encoding="utf-8") as dst:
        if ngrdcol > 1:
            dst.write("&multicol_def\n")
            dst.write(f"ngrdcol = {ngrdcol}\n")
            dst.write(f"batch_size = {ngrdcol}\n")
            dst.write("/\n\n")
        dst.write("&clubb_params_nl\n")
        for name, default_value in defaults:
            values = []
            for params in param_sets:
                values.append(float((params or {}).get(name, default_value)))
            dst.write(f"{name} = {', '.join(f'{value:.15g}' for value in values)}\n")
        dst.write("/\n")
    return params_path


def start_loss_run(
    case_names,
    fields,
    params,
    rank=None,
    case_configs=None,
    run_mode="window",
):
    """Launch an ad-hoc run for one or more result-table parameter sets."""
    cases = [str(case).strip() for case in (case_names or []) if str(case).strip()]
    selected_fields = [str(field).strip() for field in (fields or []) if str(field).strip()]
    run_mode = str(run_mode or "window").strip()
    if not cases:
        raise RuntimeError("Select at least one case before running a loss.")
    if run_mode not in {"window", "complete"}:
        raise RuntimeError(f"Unknown loss-run mode: {run_mode}")
    if run_mode == "window" and not selected_fields:
        raise RuntimeError("Select at least one field before running a loss.")

    if isinstance(params, list):
        param_sets = [dict(param_set or {}) for param_set in params]
    else:
        param_sets = [dict(params or {})]

    if not any(param_set for param_set in param_sets):
        raise RuntimeError("No parameters were available for this result row.")

    CLUBB_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    rank_text = f"rank_{rank}_" if rank not in (None, "") else ""
    run_id = f"{rank_text}{timestamp}_{os.getpid()}"
    work_dir = CLUBB_OUTPUT_DIR
    params_path = write_loss_params_file(work_dir, param_sets, f"{run_id}_loss_params.in")

    request_payload = {
        "cases": cases,
        "case_configs": list(case_configs or []),
        "fields": selected_fields,
        "run_mode": run_mode,
        "params": param_sets,
        "params_path": str(params_path),
        "output_dir": str(CLUBB_OUTPUT_DIR),
        "rank": rank,
    }
    request_path = work_dir / f"{run_id}_loss_request.json"
    log_name = "run_scm.log" if run_mode == "complete" else "run_scm_loss.log"
    log_path = work_dir / f"{run_id}_{log_name}"
    request_path.write_text(json.dumps(request_payload, indent=2, sort_keys=True), encoding="utf-8")

    if run_mode == "window":
        commands = [
            [
                sys.executable,
                str(Path(REPO_ROOT) / "run_scripts" / "run_scm_loss.py"),
                "-out_dir",
                str(CLUBB_OUTPUT_DIR),
                "-fields",
                ",".join(selected_fields),
                "-cases",
                ",".join(cases),
                "-params",
                str(params_path),
                "-case_config_file",
                str(request_path),
            ]
        ]
    else:
        commands = [
            [
                sys.executable,
                str(Path(REPO_ROOT) / "run_scripts" / "run_scm.py"),
                "-out_dir",
                str(CLUBB_OUTPUT_DIR),
                "-params",
                str(params_path),
                case_name,
            ]
            for case_name in cases
        ]

    processes = []
    for command in commands:
        log_handle = open(log_path, "a", encoding="utf-8")
        log_handle.write("\n=== Command: " + " ".join(command) + " ===\n")
        log_handle.flush()
        proc = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            env=worker_env(),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
        log_handle.close()
        processes.append(proc)

    proc = processes[0] if len(processes) == 1 else ProcessGroup(processes)
    with LOSS_RUN_LOCK:
        LOSS_RUN_PROCS[run_id] = proc

    return {
        "run_id": run_id,
        "pid": proc.pid,
        "rank": rank,
        "state": "running",
        "returncode": None,
        "request_path": str(request_path),
        "log_path": str(log_path),
        "work_dir": str(work_dir),
        "cmd": commands[0] if len(commands) == 1 else commands,
    }


def _pid_is_alive(pid):
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


def poll_loss_runs(loss_runs):
    """Refresh stored ad-hoc loss-run states from tracked subprocesses."""
    updated = dict(loss_runs or {})
    any_running = False

    with LOSS_RUN_LOCK:
        for key, run_data in list(updated.items()):
            if (run_data or {}).get("state") != "running":
                continue

            run_id = run_data.get("run_id")
            proc = LOSS_RUN_PROCS.get(run_id)
            if proc is None:
                if _pid_is_alive(run_data.get("pid")):
                    any_running = True
                    continue
                run_data = dict(run_data)
                run_data["state"] = "error"
                run_data["returncode"] = None
                updated[key] = run_data
                continue

            returncode = proc.poll()
            if returncode is None:
                any_running = True
                continue

            LOSS_RUN_PROCS.pop(run_id, None)
            run_data = dict(run_data)
            run_data["returncode"] = int(returncode)
            run_data["state"] = "success" if returncode == 0 else "error"
            updated[key] = run_data

    return updated, any_running


def stop_tuning_job(job_data):
    """Request termination of the active tuning worker, if one exists."""
    if not job_data:
        return
    TunerJob.from_dict(job_data).request_stop()


def clear_finished_job(job_data):
    """Forget a finished worker without deleting its retained job directory."""
    with TUNE_LOCK:
        active_pid = TUNE_ACTIVE.get("pid")
        if active_pid == (job_data or {}).get("pid"):
            TUNE_ACTIVE.clear()


def read_tuning_status(status_path):
    """Read the worker status file, tolerating partial or missing writes."""
    if not status_path:
        return empty_status_payload()
    path = Path(status_path)
    if not path.is_file():
        return empty_status_payload()
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return empty_status_payload()


def read_tuning_results(results_path):
    """Read the retained worker results file, tolerating partial or missing writes."""
    if not results_path:
        return empty_results_payload()
    path = Path(results_path)
    if not path.is_file():
        return empty_results_payload()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return empty_results_payload()
    if not isinstance(payload, dict):
        return empty_results_payload()
    if not isinstance(payload.get("best_results"), list):
        payload["best_results"] = []
    if not isinstance(payload.get("best_results_by_case"), dict):
        payload["best_results_by_case"] = {}
    return payload


def active_job_exited(job_data):
    """Return whether the tracked worker process has exited."""
    pid = (job_data or {}).get("pid")
    if not pid:
        return True

    with TUNE_LOCK:
        proc = TUNE_ACTIVE.get("proc") if TUNE_ACTIVE.get("pid") == pid else None
        if proc is not None:
            return TunerJob.from_dict(job_data, proc=proc).has_exited()
    return TunerJob.from_dict(job_data).has_exited()
