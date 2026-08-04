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
from tuner.status import atomic_write_json, write_control
from tuner.workspaces import create_workspace, replace_draft_request, workspace_display_name
from utilities.create_case_namelist import normalize_override_string, resolve_tunable_config_dir


TUNE_RESULT_OUTPUT_ROOT = Path(OUTPUT_TUNER_DIR)
_RESULT_LABEL_UNSAFE = re.compile(r"[^A-Za-z0-9_.-]+")


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


class SequentialProcessGroup:
    """Poll-compatible process group that runs parameter batches in order.

    A multi-column SCM replay writes one stats file per case.  Independent
    parameter batches must therefore not overlap: both would otherwise write
    (for example) ``arm_stats.nc`` at the same time. The first batch uses the
    stable result-run directory, and later batches use sibling batch folders.
    """

    def __init__(self, batch_specs):
        self._batch_specs = list(batch_specs)
        self._batch_index = -1
        self._current = None
        self._cancelled = False
        self._launch_next_batch()

    @property
    def pid(self):
        return self._current.pid if self._current is not None else None

    @property
    def processes(self):
        return list(getattr(self._current, "processes", []) or [])

    @property
    def active_batch(self):
        return self._batch_index + 1 if self._current is not None else len(self._batch_specs)

    @property
    def batch_count(self):
        return len(self._batch_specs)

    def snapshot(self):
        return {
            "active_batch": self.active_batch,
            "batch_count": self.batch_count,
            "pid": self.pid,
        }

    def _launch_next_batch(self):
        self._batch_index += 1
        if self._batch_index >= len(self._batch_specs):
            self._current = None
            return
        spec = self._batch_specs[self._batch_index]
        processes = []
        for command in spec["commands"]:
            with open(spec["log_path"], "a", encoding="utf-8") as log_handle:
                log_handle.write(
                    f"\n=== Batch {self._batch_index + 1}/{len(self._batch_specs)} command: "
                    + " ".join(command)
                    + " ===\n"
                )
                log_handle.flush()
                processes.append(
                    subprocess.Popen(
                        command,
                        cwd=REPO_ROOT,
                        env=worker_env(),
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                )
        self._current = ProcessGroup(processes)

    def poll(self):
        if self._cancelled:
            return 143
        if self._current is None:
            return 0
        returncode = self._current.poll()
        if returncode is None:
            return None
        if returncode:
            return int(returncode)
        self._launch_next_batch()
        return None if self._current is not None else 0

    def terminate(self):
        """Stop the current batch and prevent any subsequent launch."""
        self._cancelled = True
        for process in self.processes:
            if process.poll() is None:
                process.terminate()


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


def _track_tuning_job(job, *, workspace_id=None, revision_id=None, resume=False):
    """Start one job and retain its JSON-safe lineage metadata."""
    job.start(resume=resume)
    job_data = job.to_dict()
    if workspace_id:
        job_data["workspace_id"] = str(workspace_id)
    if revision_id:
        job_data["revision_id"] = str(revision_id)
    job_data["resume"] = bool(resume)
    with TUNE_LOCK:
        TUNE_ACTIVE.update(job_data)
        TUNE_ACTIVE["proc"] = job.proc
    return job_data


def start_tuning_job(request_payload, *, display_name=None):
    """Create a new durable workspace and launch its ``original`` execution."""
    with TUNE_LOCK:
        existing_proc = TUNE_ACTIVE.get("proc")
        if existing_proc is not None and existing_proc.poll() is None:
            raise RuntimeError("A tuning job is already running")
        TUNE_ACTIVE.clear()
    workspace_id, job = create_workspace(request_payload, display_name=display_name)
    return _track_tuning_job(job, workspace_id=workspace_id, revision_id="original")


def resume_tuning_job(job_data):
    """Resume a gracefully stopped revision in its original output directory."""
    with TUNE_LOCK:
        existing_proc = TUNE_ACTIVE.get("proc")
        if existing_proc is not None and existing_proc.poll() is None:
            raise RuntimeError("A tuning job is already running")
        TUNE_ACTIVE.clear()
    job = TunerJob.from_dict(job_data)
    status = read_tuning_status(job.status_path)
    if str(status.get("state") or "") != "stopped":
        raise RuntimeError("only a stopped tuning revision can continue")
    return _track_tuning_job(
        job,
        workspace_id=job_data.get("workspace_id"),
        revision_id=job_data.get("revision_id"),
        resume=True,
    )


def start_draft_tuning_job(workspace_id, revision_id, request_payload):
    """Start an editable draft revision after sealing its request locally."""
    with TUNE_LOCK:
        existing_proc = TUNE_ACTIVE.get("proc")
        if existing_proc is not None and existing_proc.poll() is None:
            raise RuntimeError("A tuning job is already running")
        TUNE_ACTIVE.clear()
    job = replace_draft_request(workspace_id, revision_id, request_payload)
    write_control(job.control_path, stop_requested=False)
    initial_status = empty_status_payload()
    initial_status.update({"state": "draft", "job_dir": str(job.job_dir)})
    atomic_write_json(job.status_path, initial_status)
    return _track_tuning_job(job, workspace_id=workspace_id, revision_id=revision_id)


def _read_default_tunable_params(config="default"):
    """Return scalar CLUBB tunable parameter defaults in file order for a config."""
    params = []
    pattern = re.compile(r"^\s*([A-Za-z]\w*)\s*=\s*([^!,/]+)")
    config_dir = Path(resolve_tunable_config_dir(config or "default"))
    params_path = config_dir / "tunable_parameters.in"
    with open(params_path, encoding="utf-8") as src:
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


def write_loss_params_file(work_dir, param_sets, filename="loss_params.in", config="default"):
    """Write a params namelist whose only changes are the supplied tuned parameters."""
    defaults = _read_default_tunable_params(config)
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


def _result_label(value, fallback):
    """Make a workspace display label safe and stable as one path component."""
    compact = re.sub(r"\s+", "_", str(value or "").strip())
    compact = _RESULT_LABEL_UNSAFE.sub("_", compact).strip("._-")
    return compact or fallback


def tune_result_output_dir(
    workspace_id=None,
    revision_id=None,
    run_mode="window",
    workspace_name=None,
):
    """Return the stable output location for one Tune revision result rerun."""
    workspace_id = str(workspace_id or "").strip()
    if not workspace_name and workspace_id:
        try:
            workspace_name = workspace_display_name(workspace_id)
        except (OSError, ValueError):
            workspace_name = workspace_id
    workspace_label = _result_label(workspace_name or workspace_id, "unnamed")
    revision_label = _result_label(revision_id, "original")
    run_label = "loss_window" if str(run_mode or "window").strip() == "window" else "complete"
    return TUNE_RESULT_OUTPUT_ROOT / f"{workspace_label}_{revision_label}_{run_label}"


def start_loss_run(
    case_names,
    fields,
    params,
    rank=None,
    case_configs=None,
    run_mode="window",
    config="default",
    override="",
    batch_size=8,
    workspace_id=None,
    revision_id=None,
    workspace_name=None,
):
    """Launch an ad-hoc replay, splitting parameter columns into safe batches.

    The tuning worker already evaluates at most ``batch_size`` columns per
    loss-driver call.  Result-table replays use the same limit: this avoids a
    known CLUBB stack overflow for large multi-column cases while retaining
    every selected parameter row.  Batches run in order because concurrently
    writing a case's stats file would corrupt its output.
    """
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
    try:
        batch_size = int(batch_size)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("Batch size must be a positive integer.") from exc
    if batch_size < 1:
        raise RuntimeError("Batch size must be a positive integer.")

    output_root = tune_result_output_dir(
        workspace_id=workspace_id,
        revision_id=revision_id,
        run_mode=run_mode,
        workspace_name=workspace_name,
    )
    output_root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    rank_text = f"rank_{rank}_" if rank not in (None, "") else ""
    run_id = f"{rank_text}{timestamp}_{os.getpid()}"
    work_dir = output_root / "_run_metadata"
    work_dir.mkdir(parents=True, exist_ok=True)
    config = str(config or "default").strip() or "default"
    override = normalize_override_string(override)

    param_batches = [param_sets[start : start + batch_size] for start in range(0, len(param_sets), batch_size)]
    batch_records = []
    for index, batch_params in enumerate(param_batches, start=1):
        # The first batch owns the user-facing revision/run-type directory.
        # Later batches need isolated stats files, so they occupy child folders.
        output_dir = output_root if index == 1 else output_root / f"batch_{index:03d}"
        output_dir.mkdir(parents=True, exist_ok=True)
        params_path = write_loss_params_file(
            work_dir,
            batch_params,
            f"batch_{index:03d}_loss_params.in",
            config=config,
        )
        batch_records.append(
            {
                "index": index,
                "param_count": len(batch_params),
                "params_path": str(params_path),
                "output_dir": str(output_dir),
            }
        )

    request_payload = {
        "config": config,
        "cases": cases,
        "case_configs": list(case_configs or []),
        "fields": selected_fields,
        "run_mode": run_mode,
        "params": param_sets,
        "batch_size": batch_size,
        "batch_count": len(batch_records),
        "batches": batch_records,
        "output_dir": str(output_root),
        "workspace_id": str(workspace_id or ""),
        "revision_id": str(revision_id or ""),
        "rank": rank,
        "override": override,
    }
    request_path = work_dir / f"{run_id}_loss_request.json"
    log_name = "run_scm.log" if run_mode == "complete" else "run_scm_loss.log"
    log_path = work_dir / f"{run_id}_{log_name}"
    request_path.write_text(json.dumps(request_payload, indent=2, sort_keys=True), encoding="utf-8")

    batch_specs = []
    for record in batch_records:
        if run_mode == "window":
            commands = [
                [
                    sys.executable,
                    str(Path(REPO_ROOT) / "run_scripts" / "run_scm_loss.py"),
                    "-out_dir",
                    record["output_dir"],
                    "-config",
                    config,
                    "-fields",
                    ",".join(selected_fields),
                    "-cases",
                    ",".join(cases),
                    "-params",
                    record["params_path"],
                    *(["-override", override] if override else []),
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
                    record["output_dir"],
                    "-config",
                    config,
                    "-params",
                    record["params_path"],
                    *(["-override", override] if override else []),
                    case_name,
                ]
                for case_name in cases
            ]
        batch_specs.append({"commands": commands, "log_path": log_path})

    proc = SequentialProcessGroup(batch_specs)
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
        "output_dir": str(output_root),
        "batch_size": batch_size,
        "batch_count": len(batch_specs),
        "active_batch": 1,
        "batches": batch_records,
        "cmd": batch_specs[0]["commands"][0] if len(batch_specs[0]["commands"]) == 1 else batch_specs[0]["commands"],
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
            if (run_data or {}).get("state") not in {"running", "stopping"}:
                continue

            run_id = run_data.get("run_id")
            proc = LOSS_RUN_PROCS.get(run_id)
            if proc is None:
                if _pid_is_alive(run_data.get("pid")):
                    any_running = True
                    continue
                run_data = dict(run_data)
                run_data["state"] = "stopped" if run_data.get("state") == "stopping" else "error"
                run_data["returncode"] = None
                updated[key] = run_data
                continue

            returncode = proc.poll()
            if returncode is None:
                if hasattr(proc, "snapshot"):
                    run_data = dict(run_data)
                    run_data.update(proc.snapshot())
                    updated[key] = run_data
                any_running = True
                continue

            LOSS_RUN_PROCS.pop(run_id, None)
            run_data = dict(run_data)
            run_data["returncode"] = int(returncode)
            if run_data.get("state") == "stopping":
                run_data["state"] = "stopped"
            else:
                run_data["state"] = "success" if returncode == 0 else "error"
            updated[key] = run_data

    return updated, any_running


def stop_loss_run(run_data):
    """Request a bounded Tune leaderboard rerun to stop.

    This mirrors ordinary SCM cancellation: retain the durable record while the
    process exits so the polling path can publish one final, unambiguous state.
    """
    data = dict(run_data or {})
    if data.get("state") not in {"running", "stopping"}:
        return data
    run_id = str(data.get("run_id") or "")
    with LOSS_RUN_LOCK:
        proc = LOSS_RUN_PROCS.get(run_id)
        try:
            if proc is None:
                os.kill(int(data.get("pid")), 15)
            else:
                if hasattr(proc, "terminate"):
                    proc.terminate()
                else:
                    processes = list(getattr(proc, "processes", []) or [])
                    for child in processes:
                        if child is not None and child.poll() is None:
                            child.terminate()
        except (OSError, TypeError, ValueError):
            # The normal poller will turn a missing child into an explicit
            # terminal state; a duplicate cancellation remains idempotent.
            pass
    data["state"] = "stopping"
    return data


def stop_tuning_job(job_data):
    """Request termination of the active tuning worker, if one exists."""
    if not job_data:
        return
    TunerJob.from_dict(job_data).request_stop()


def active_tuning_job_data():
    """Return the current server-owned tuner job without its live process object."""
    with TUNE_LOCK:
        if not TUNE_ACTIVE:
            return {}
        return {key: value for key, value in TUNE_ACTIVE.items() if key != "proc"}


def stop_active_tuning_job():
    """Request a graceful stop for the current server-owned tuner job."""
    job_data = active_tuning_job_data()
    if job_data:
        stop_tuning_job(job_data)
    return job_data


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
