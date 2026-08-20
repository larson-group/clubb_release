"""Dashboard application services and semantic browser operations.

This is shared runtime/application code. Agent adapters call these services,
but they do not own them.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import signal
import threading
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]

from .activity import (
    publish_event,
    publish_budget_request,
    publish_plot_remove_request,
    publish_plot_request,
    publish_tab_request,
    publish_tune_request,
    broker_jobs,
    read_activity,
    set_broker_job,
    set_broker_loss_run,
    update_broker_job,
    update_broker_loss_run,
)
from .tab_registry import (
    describe_tabs,
    invoke as invoke_tab_operation,
    register_operation,
    register_tab,
)
from dash_app.compile_tab.discovery import discover_compile_state
from dash_app.compile_tab.runtime import (
    cancel_compile_job,
    finish_compile_job,
    poll_compile_job,
    read_log_increment as read_compile_log_increment,
    start_compile_job,
    start_rebuild_job,
)
from dash_app.run_tab.namelist import cleanup_temp_files
from dash_app.run_tab.runtime import (
    get_proc,
    record_case_finish,
    start_case_process,
)
from dash_app.run_tab.telemetry import run_telemetry as collect_run_telemetry
from dash_app.run_tab.telemetry import scm_run_view as _scm_run_view
from dash_app.run_tab.state import MAX_RUN_PROCS
from dash_app.run_tab.state import DEFAULT_STATS_NAME
from dash_app.plot_tab.plot_types.profile_plot import PLOT as profile_plot
from dash_app.plot_tab.pyplotgen_runtime import (
    release_pyplotgen,
    start_pyplotgen,
    stop_pyplotgen,
)
from dash_app.profile_tab.runtime import (
    profile_process_status,
    profile_results_complete,
    read_log_tail as read_profile_log_tail,
    read_profile_results,
    release_profile_process,
    start_profile_process,
    stop_profile_process,
)
from dash_app.plot_tab.plot_types.shared import scan_output_cases
from dash_app.services import plots as plot_service
from dash_app.services.plots import validate_plot_request
from dash_app.tune_tab.discovery import (
    available_fields_for_cases,
    available_tunable_configs,
    load_case_defaults,
    load_tunable_default_ranges,
    load_tunable_names,
)
from dash_app.tune_tab.runtime import (
    active_tuning_job_data,
    active_job_exited,
    poll_loss_runs,
    read_tuning_status,
    resume_tuning_job,
    start_loss_run,
    start_draft_tuning_job,
    start_tuning_job,
    stop_loss_run,
    stop_active_tuning_job,
    stop_tuning_job,
)
from tuner.workspaces import (
    create_revision as create_tune_workspace_revision,
    delete_workspace as delete_tune_workspace,
    list_workspace_activity as list_tune_workspace_activity,
    list_workspaces as list_tune_workspaces,
    load_execution as load_tune_workspace_execution,
    rename_workspace as rename_tune_workspace,
    reset_execution as reset_tune_workspace_execution,
)
from tuner.request import (
    apply_required_parameter_links,
    evaluate_tune_settings,
    read_case_tuner_defaults,
)
from tuner.job_runtime import TunerJob
from tuner.system_defaults import available_logical_cpu_count
from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_LOSS_MODE,
    DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    LOSS_MODE_NAMES,
    TIME_WINDOW_AGGREGATION_SCOPES,
    normalize_aggregation_weights,
)
from tuner.tuning_strategy import VALID_STRATEGY_NAMES
from utilities.create_case_namelist import normalize_override_string, parse_override_pairs
from utilities.output_paths import resolve_output_dir
from utilities.clubb_settings_validation import (
    canonical_flag_name,
    canonical_parameter_name,
    is_independently_tunable,
    resolve_clubb_settings,
)
from dash_app.run_tab.discovery import list_stats_files
from dash_app.run_tab.state import NO_STATS_NAME
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name, tunable_config_file
from dash_app.services import (
    ArtifactStore,
    CompileRequest,
    JobStore,
    LeaderboardRerunRequest,
    ParameterRange,
    ProfileArtifactRequest,
    ScmRunBatchRequest,
    ScmRunRequest,
    SubmissionConflict,
    TuneRequest,
)
from dash_app.services import profiles as profile_service
from dash_app.shared.provenance import sha256_file, source_provenance
from dash_app.shared.runtime import (
    atomic_write_json,
    exclusive_file_lock,
    private_path,
    process_is_alive as _pid_is_alive,
    read_file_chunk,
)
from . import dashboard_registry


_CASE_RE = re.compile(r"^[A-Za-z0-9_]+$")
_OVERRIDE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*=[^\s;|&`$<>]+$")
_NAMELIST_KEY_RE = re.compile(r"^[A-Za-z_]\w*$")
_TUNE_PARAM_RE = re.compile(r"^[A-Za-z]\w*$")
_BATCH_ID_RE = re.compile(r"^batch_[A-Za-z0-9_-]+$")
# Native Tune presets intentionally span coupled closure families.  The
# largest current preset has thirteen logical sampled coordinates, so keep a
# modest broker-side guard without making a valid dashboard preset impossible
# to launch.  This limits logical coordinates, not their linked physical
# targets.
_TUNE_MAX_PARAMETERS = 24
_TUNE_MAX_RANDOM_SAMPLES = 1000
_TUNE_MAX_SIMANN_ITERS = 5000
_TUNE_MAX_ADAM_UPDATES = 5000
_TUNE_MAX_RESOLVE_SAMPLES = 1000
_TUNE_MAX_BATCH_SIZE = 64
_PROFILE_WINDOW_PRESETS = {"loss", "pyplotgen"}
_PROFILE_EXPORT_DIR = REPO_ROOT / "output" / "agent_exports"
_JOB_STORE = JobStore()
_ARTIFACT_STORE = ArtifactStore()
_BATCH_QUEUE_LOCK = threading.Lock()
_BATCH_RESULT_LOCK = threading.Lock()
_BATCH_WATCHERS_LOCK = threading.Lock()
_BATCH_WATCHERS: set[str] = set()
_BROKER_SHUTTING_DOWN = threading.Event()
# ``inspect_dashboard``/``invoke_dashboard`` survive only to make an agent's
# work visible in an already-open browser.  Scientific actions must use the
# typed MCP service, where they receive request-idempotency, durable job IDs,
# and immutable provenance instead of inheriting mutable Dash state.
_BROWSER_HANDOFF_OPERATIONS = {
    ("tutorial", "navigate"),
    ("tutorial", "open_lesson"),
    ("reports", "navigate"),
    ("reports", "open_report"),
    ("misc", "navigate"),
    ("misc", "open_report"),
    ("compile", "navigate"),
    ("run", "navigate"),
    ("profile", "navigate"),
    ("tune", "navigate"),
    ("plots", "navigate"),
    ("plots", "set_view"),
    ("plots", "add_budget"),
    ("plots", "list"),
    ("plots", "remove"),
}


def _tail(text: str, limit: int = 1800) -> str:
    text = str(text or "").strip()
    return text[-limit:] if len(text) > limit else text


def _read_log_tail(log_path: str | None, limit: int = 16000) -> str:
    """Read a bounded log suffix for reconnecting Dash consoles."""
    if not log_path:
        return ""
    try:
        path = Path(log_path)
        with path.open("rb") as handle:
            handle.seek(0, 2)
            size = handle.tell()
            handle.seek(max(0, size - limit))
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def _write_run_result_summary(job_id: str, output_directory: str | None, returncode: int) -> str | None:
    """Seal compact output checksums after an isolated SCM process exits."""
    if not output_directory:
        return None
    root = Path(output_directory)
    checksums: list[dict[str, Any]] = []
    try:
        candidates = sorted(path for path in root.rglob("*") if path.is_file())[:64]
    except OSError:
        candidates = []
    for path in candidates:
        try:
            checksums.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
        except OSError:
            continue
    try:
        summary = _ARTIFACT_STORE.write_summary(
            job_id,
            "execution_result.json",
            {
                "job_id": job_id,
                "returncode": int(returncode),
                "output_directory": str(root),
                "output_checksums": checksums,
                "truncated": len(candidates) >= 64,
                "completed_at_unix_seconds": time.time(),
            },
        )
        return str(summary)
    except (OSError, ValueError, FileExistsError):
        return None


def _public_scm_output_root() -> Path:
    return REPO_ROOT / "output" / "mcp_runs"


def _resolve_mcp_output_dir(value: str) -> Path:
    """Resolve one MCP output path while keeping it below the repo output root."""
    output_root = (REPO_ROOT / "output").resolve()
    requested = Path(str(value or "").strip()).expanduser()
    if not requested:
        raise ValueError("out_dir must not be empty")
    candidate = requested if requested.is_absolute() else output_root / requested
    if not requested.is_absolute() and requested.parts and requested.parts[0] == "output":
        candidate = REPO_ROOT / requested
    resolved = candidate.resolve()
    try:
        resolved.relative_to(output_root)
    except ValueError as exc:
        raise ValueError("out_dir must resolve inside the repository output directory") from exc
    return resolved


def _public_scm_batch_output_dir(batch_id: str) -> Path:
    """Return the default single flat, broker-owned directory for one batch."""
    safe_batch = str(batch_id or "").strip()
    if not _BATCH_ID_RE.fullmatch(safe_batch):
        raise ValueError("SCM batch is missing a valid batch identifier")
    return (_public_scm_output_root() / safe_batch).resolve()


def _canonical_output_directory(value: Any) -> str:
    """Return one stable filesystem identity without creating the directory."""
    text = str(value or "").strip()
    return str(Path(text).expanduser().resolve()) if text else ""


def _scm_record_output_directory(record: dict[str, Any]) -> str:
    """Read the canonical target from current and pre-migration SCM records."""
    runtime = dict(record.get("runtime") or {})
    proc_data = dict(runtime.get("proc_data") or {})
    display = dict(record.get("display") or {})
    cli_options = dict(display.get("cli_options") or {})
    for value in (
        record.get("output_directory"),
        runtime.get("output_directory"),
        proc_data.get("output_directory"),
        cli_options.get("out_dir"),
    ):
        canonical = _canonical_output_directory(value)
        if canonical:
            return canonical
    return ""


def _batch_child_summary(record: dict[str, Any]) -> dict[str, Any]:
    """Expose durable child identity/status without copying large runtime logs."""
    summary = {
        key: record.get(key)
        for key in (
            "job_id",
            "run_id",
            "state",
            "returncode",
            "error",
            "output_directory",
            "runtime",
            "finished_at_unix_seconds",
            "updated_at_unix_seconds",
            "result_summary_job_id",
            "result_summary_path",
        )
        if record.get(key) is not None
    }
    case = record.get("case") or (record.get("request") or {}).get("case")
    if case:
        summary["case"] = case
    return summary


def _batch_child_job_ids(batch: dict[str, Any]) -> list[str]:
    """Read canonical child references, with compatibility for pre-v10 batches."""
    child_ids = list(batch.get("child_job_ids") or [])
    if child_ids:
        return [str(job_id) for job_id in child_ids if job_id]
    return [
        str(child.get("job_id") or "")
        for child in batch.get("children") or []
        if child.get("job_id")
    ]


def _batch_children(batch: dict[str, Any], transaction=None) -> list[dict[str, Any]]:
    getter = transaction.get if transaction is not None else _JOB_STORE.get
    return [
        record
        for job_id in _batch_child_job_ids(batch)
        if (record := getter(job_id)) is not None
    ]


def _materialize_scm_batch(batch: dict[str, Any] | None) -> dict[str, Any] | None:
    """Attach current child summaries at public response boundaries."""
    if batch is None:
        return None
    public_batch = dict(batch)
    public_batch.pop("public_manifest_path", None)
    return {
        **public_batch,
        "children": [
            _batch_child_summary(child) for child in _batch_children(batch)
        ],
    }


def _batch_record(job_id: str) -> dict[str, Any] | None:
    record = _JOB_STORE.get(job_id)
    return record if record and record.get("kind") == "scm_batch" else None


def _batch_terminal(state: str) -> bool:
    return state in {"finished", "partial_failure", "error", "cancelled"}


def _batch_state(children: list[dict[str, Any]]) -> str:
    if not children:
        return "finished"
    states = {str(child.get("state") or "queued") for child in children}
    if states & {"starting", "running", "stopping"}:
        return "running"
    if states & {"queued", "submitting"}:
        return "queued"
    if states and all(str(child.get("state")) == "finished" and child.get("returncode") == 0 for child in children):
        return "finished"
    return "partial_failure"


def _refresh_scm_batch(batch_job_id: str, _child_job_id: str | None = None) -> dict[str, Any] | None:
    """Derive the parent lifecycle from its canonical child records."""
    def refresh(transaction):
        batch = transaction.get(batch_job_id)
        if batch is None or batch.get("kind") != "scm_batch":
            return None
        children = _batch_children(batch, transaction)
        cancellation_requested = bool(batch.get("cancellation_requested"))
        children_active = any(
            str(child.get("state") or "") in {"starting", "running", "stopping"}
            for child in children
        )
        state = (
            "stopping" if children_active else "cancelled"
        ) if cancellation_requested else _batch_state(children)
        terminal = _batch_terminal(state)
        result_path = str(_ARTIFACT_STORE.bundle_dir(batch_job_id) / "execution_result.json")
        return transaction.update(
            batch_job_id,
            state=state,
            result_summary_path=result_path if terminal else batch.get("result_summary_path"),
            finished_at_unix_seconds=time.time() if terminal else None,
        )

    updated = _materialize_scm_batch(_JOB_STORE.transaction(refresh))
    if updated is None:
        return None
    state = str(updated.get("state") or "")
    if _batch_terminal(state):
        returncode = 0 if state == "finished" else 1
        with _BATCH_RESULT_LOCK:
            result_summary_path = Path(str(updated.get("result_summary_path") or ""))
            if not result_summary_path.is_file():
                _write_run_result_summary(
                    batch_job_id, updated.get("output_directory"), returncode
                )
        _ARTIFACT_STORE.release(batch_job_id)
    return updated


def _validate_case(case: str) -> str:
    value = str(case or "").strip()
    if not _CASE_RE.fullmatch(value):
        raise ValueError("case must be a simple CLUBB case name")
    if not (REPO_ROOT / "input" / "case_setups" / f"{value}_model.in").is_file():
        raise ValueError(f"unknown or unavailable case: {value}")
    return value


def _override_args(overrides: str) -> list[str]:
    value = str(overrides or "").strip()
    if not value:
        return []
    tokens = shlex.split(value)
    if tokens and tokens[0] == "-override":
        tokens = tokens[1:]
    if not tokens or any(not _OVERRIDE_RE.fullmatch(token) for token in tokens):
        raise ValueError("overrides must be whitespace-separated key=value tokens, optionally prefixed with -override")
    # ``run_scm.py`` accepts one comma-separated override argument.  Keeping
    # the agent-facing form whitespace-separated makes individual assignments
    # easy to read, but they must be joined before passing through argparse.
    return ["-override", ",".join(tokens)]


def _watch_compile(job: dict[str, Any], job_id: str | None = None) -> None:
    offset = 0
    while True:
        chunk, offset = read_compile_log_increment(job.get("log"), offset)
        if chunk:
            publish_event("compile", "Compile output", _tail(chunk), status="running")
            update_broker_job("compile", log_tail=_read_log_tail(job.get("log")), log_offset=offset)
            if job_id:
                _JOB_STORE.update(job_id, progress={"log_offset": offset})
        returncode = poll_compile_job(job)
        if returncode is not None:
            final_job = finish_compile_job(job, returncode)
            requested_stop = str((broker_jobs().get("compile") or {}).get("state") or "") == "stopping"
            broker_state = "finished" if returncode == 0 else ("stopped" if requested_stop else "error")
            job_state = "finished" if returncode == 0 else ("cancelled" if requested_stop else "error")
            update_broker_job("compile", state=broker_state, returncode=returncode, log_offset=offset)
            publish_event(
                "compile",
                "Compile finished" if returncode == 0 else ("Compile stopped" if requested_stop else "Compile failed"),
                _tail(str(final_job.get("command") or "")),
                status="success" if returncode == 0 else ("info" if requested_stop else "error"),
            )
            if job_id:
                _JOB_STORE.update(job_id, state=job_state, returncode=returncode, log_tail=_read_log_tail(job.get("log")))
                _ARTIFACT_STORE.release(job_id)
            return
        time.sleep(0.75)


def _watch_profile(job: dict[str, Any]) -> None:
    """Mirror one timing wrapper into durable broker state and activity."""
    while True:
        run_id, rows = read_profile_results(job)
        running, returncode = profile_process_status(job)
        log_tail = read_profile_log_tail(job.get("log"))
        if run_id:
            job["run_id"] = run_id
        update_broker_job(
            "profile",
            run_id=job.get("run_id") or "",
            result_rows=len(rows),
            log_tail=log_tail,
        )
        if not running:
            record = dict(broker_jobs().get("profile") or {})
            requested_stop = str(record.get("state") or "") == "stopping"
            complete = profile_results_complete(job, rows)
            success = returncode == 0 or (returncode is None and complete)
            state = "finished" if success else ("stopped" if requested_stop else "error")
            update_broker_job(
                "profile",
                state=state,
                returncode=returncode,
                run_id=job.get("run_id") or "",
                result_rows=len(rows),
                log_tail=log_tail,
                finished_at=time.time(),
            )
            release_profile_process(job.get("pid"))
            publish_event(
                "profile",
                "Profile timing finished"
                if success
                else ("Profile timing stopped" if requested_stop else "Profile timing failed"),
                str(job.get("output") or ""),
                status="success" if success else ("info" if requested_stop else "error"),
                action={"type": "profile", "tab": "profile"},
            )
            return
        time.sleep(0.5)


def _watch_run(case: str, proc_data: dict[str, Any], job_id: str | None = None) -> None:
    pid = proc_data.get("pid")
    while True:
        proc = get_proc(pid)
        # A detached broker can restart while SCM continues.  In that case
        # there is no local Popen, but PID liveness is enough to retain the
        # log stream, stop control, and active-job guard.  Once it exits the
        # POSIX exit status is no longer observable, so report that recovery
        # limitation explicitly rather than inventing a successful result.
        returncode = proc.poll() if proc is not None else (None if _pid_is_alive(pid) else 1)
        if returncode is not None:
            requested_stop = returncode in {-signal.SIGTERM, -signal.SIGKILL}
            if job_id:
                requested_stop = requested_stop or str(
                    (_JOB_STORE.get(job_id) or {}).get("state") or ""
                ) in {"stopping", "cancelled"}
            record_case_finish(case, pid, returncode)
            cleanup_temp_files(proc_data.get("temp_files") or [])
            recovered_exit = proc is None
            job_state = "finished" if returncode == 0 else ("cancelled" if requested_stop else "error")
            batch_job_id = ""
            if job_id:
                child_record = _JOB_STORE.get(job_id) or {}
                batch_job_id = str(child_record.get("batch_job_id") or "")
                _JOB_STORE.update(
                    job_id,
                    state=job_state,
                    returncode=returncode,
                    finished_at_unix_seconds=time.time(),
                    recovery_note=(
                        "exit code unavailable after broker recovery"
                        if recovered_exit
                        else None
                    ),
                )
                if batch_job_id:
                    batch_state = _refresh_scm_batch(batch_job_id, job_id)
                    if batch_state and any(
                        str(item.get("state") or "") == "queued"
                        for item in batch_state.get("children") or []
                    ):
                        _ensure_scm_batch_watcher(batch_job_id)
            publish_event(
                "run",
                f"{case} finished" if returncode == 0 else (f"{case} stopped" if requested_stop else f"{case} failed"),
                f"exit code {returncode}",
                status="success" if returncode == 0 else ("info" if requested_stop else "error"),
            )
            if job_id:
                if not batch_job_id:
                    result_summary = _write_run_result_summary(
                        job_id, proc_data.get("output_directory"), returncode
                    )
                    _JOB_STORE.update(
                        job_id,
                        result_summary_path=result_summary,
                    )
                _ARTIFACT_STORE.release(job_id)
            return
        time.sleep(0.2)


def _background(target, *args) -> None:
    threading.Thread(target=target, args=args, daemon=True).start()


def _finite_float(value: Any, label: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if numeric != numeric or numeric in {float("inf"), float("-inf")}:
        raise ValueError(f"{label} must be finite")
    return numeric


def _profile_case_data(case_name: str, *, required: bool, output_dirs: list[str] | None = None) -> dict[str, Any] | None:
    """Load the exact Plot-tab case metadata from the standard output folder."""
    return profile_service.load_case_data(REPO_ROOT, case_name, required=required, output_dirs=output_dirs)


def _plot_output_selection(
    case_name: str,
    run_id: str | None,
    output_dir: str | None,
    output_dirs: list[str] | None = None,
) -> tuple[str, list[str] | None]:
    """Resolve the immutable run or explicit output selection for Plot actions."""
    selected_run_id = str(run_id or "").strip()
    requested_output = str(output_dir or "").strip()
    selected_output = _resolve_mcp_output_dir(requested_output) if requested_output else None

    selected_outputs = None
    if output_dirs is not None:
        if selected_run_id or requested_output:
            raise ValueError("output_dirs cannot be combined with run_id or output_dir")
        if not 1 <= len(output_dirs) <= 8:
            raise ValueError("output_dirs must contain between 1 and 8 directories")
        selected_outputs = []
        for value in output_dirs:
            requested = str(value or "").strip()
            if not requested:
                raise ValueError("output_dirs cannot contain an empty directory")
            resolved = str(_resolve_mcp_output_dir(requested))
            if resolved not in selected_outputs:
                selected_outputs.append(resolved)

    if selected_run_id:
        run = _JOB_STORE.get_run(selected_run_id)
        if run is None:
            raise ValueError("unknown run_id")
        run_case = _validate_case(str((run.get("request") or {}).get("case") or ""))
        if run_case != case_name:
            raise ValueError("run_id belongs to a different SCM case")
        run_output = str(run.get("output_directory") or "")
        if not run_output:
            raise ValueError("run has no isolated output directory")
        # A run record is broker-owned provenance and may point at private
        # artifact staging outside the public output selector root. Explicit
        # output_dir values remain constrained by _resolve_mcp_output_dir.
        resolved_run_output = Path(run_output).expanduser().resolve()
        if selected_output is not None and selected_output != resolved_run_output:
            raise ValueError("run_id and output_dir select different output directories")
        selected_output = resolved_run_output

    if selected_outputs is not None:
        return selected_run_id, selected_outputs
    return selected_run_id, [str(selected_output)] if selected_output is not None else None


def _profile_selection(
    case_name: str,
    *,
    time_seconds: float | None = None,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
    require_case_data: bool,
    output_dirs: list[str] | None = None,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    """Validate one agent plot selection against the same physical controls as Dash.

    A named preset deliberately carries its exact physical window.  Dash uses a
    nearby symbolic slider value only for display, while profile extraction and
    export retain this exact value through ``plots-time-override``.
    """
    preset = str(window_preset or "").strip().lower()
    if preset and preset not in _PROFILE_WINDOW_PRESETS:
        raise ValueError("window_preset must be 'loss' or 'pyplotgen'")
    if time_seconds is not None and time_start_seconds is not None:
        legacy_start = _finite_float(time_seconds, "time_seconds")
        explicit_start = _finite_float(time_start_seconds, "time_start_seconds")
        if abs(legacy_start - explicit_start) > 1.0e-9:
            raise ValueError("time_seconds and time_start_seconds disagree; use time_start_seconds")
    requested_start = time_start_seconds if time_start_seconds is not None else time_seconds
    if preset and (requested_start is not None or average_minutes is not None):
        raise ValueError("window_preset cannot be combined with custom start or average values")

    case_data = _profile_case_data(
        case_name,
        required=require_case_data or preset or requested_start is not None or average_minutes is not None,
        output_dirs=output_dirs,
    )
    if case_data is None:
        return None, {}

    if preset:
        start = case_data.get(f"{preset}_time_start_seconds")
        duration = case_data.get(f"{preset}_time_duration_minutes")
        if start is None or duration is None:
            raise ValueError(f"{preset} window is not available for {case_name}")
        return case_data, {
            "time_start_seconds": _finite_float(start, f"{preset} window start"),
            "average_minutes": _finite_float(duration, f"{preset} window duration"),
            "window_preset": preset,
        }

    if requested_start is None and average_minutes is None:
        return case_data, {}

    return case_data, profile_service.resolve_time_window(
        case_data,
        start_seconds=_finite_float(requested_start, "time_start_seconds") if requested_start is not None else None,
        average_minutes=_finite_float(average_minutes, "average_minutes") if average_minutes is not None else None,
    )


def _validated_profile_names(variables: list[str], case_data: dict[str, Any] | None) -> list[str]:
    return profile_service.validate_profile_names(variables, case_data)


def _profile_global_context(case_data: dict[str, Any], selection: dict[str, Any]) -> dict[str, Any]:
    """Return the non-interactive equivalent of ProfilePlotType's Dash context."""
    return profile_service.figure_context(case_data, selection)


def _write_profile_figure_png(plotly_figure: Any, path: Path) -> None:
    """Rasterize a native Plotly profile figure without requiring a browser.

    Dash's normal client-side camera button is not reachable from an external
    agent.  This small renderer consumes the same trace data produced by the
    Profile tab and keeps PNG export usable in a lightweight Dash installation
    that does not bundle Chrome/Kaleido.
    """
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_dash_matplotlib")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    source_layout = plotly_figure.layout
    source_title = str(getattr(getattr(source_layout, "title", None), "text", "") or "")
    title = re.sub(r"<[^>]+>", "", source_title).replace("<br>", "\n")
    paper_color = str(getattr(source_layout, "paper_bgcolor", None) or "#0f172a")
    plot_color = str(getattr(source_layout, "plot_bgcolor", None) or paper_color)
    font_color = str(getattr(getattr(source_layout, "font", None), "color", None) or "#e5e7eb")
    grid_color = str(getattr(getattr(source_layout, "xaxis", None), "gridcolor", None) or "#334155")
    figure = Figure(figsize=(9.0, 6.0), dpi=140, facecolor=paper_color)
    FigureCanvasAgg(figure)
    axis = figure.add_subplot(1, 1, 1, facecolor=plot_color)
    line_styles = {"solid": "-", "dash": "--", "dot": ":", "dashdot": "-.", "longdash": "--", "longdashdot": "-."}
    has_legend = False
    for trace in plotly_figure.data:
        try:
            x_values = [float(value) if value is not None else float("nan") for value in (trace.x or [])]
            y_values = [float(value) if value is not None else float("nan") for value in (trace.y or [])]
        except (TypeError, ValueError):
            continue
        if not x_values or not y_values:
            continue
        line = trace.line or {}
        color = getattr(line, "color", None) or "#38bdf8"
        dash = getattr(line, "dash", None) or "solid"
        label = str(trace.name or "") if bool(trace.showlegend) else "_nolegend_"
        has_legend = has_legend or label != "_nolegend_"
        axis.plot(x_values, y_values, color=color, linestyle=line_styles.get(str(dash), "-"), linewidth=float(getattr(line, "width", None) or 1.8), alpha=getattr(trace, "opacity", None), label=label)
    xaxis = getattr(source_layout, "xaxis", None)
    yaxis = getattr(source_layout, "yaxis", None)
    axis.set_title(title, color=font_color, pad=14)
    axis.set_xlabel(str(getattr(getattr(xaxis, "title", None), "text", "") or ""), color=font_color)
    axis.set_ylabel(str(getattr(getattr(yaxis, "title", None), "text", "") or ""), color=font_color)
    for spine in axis.spines.values():
        spine.set_color(grid_color)
    axis.tick_params(colors=font_color)
    axis.grid(True, color=grid_color, alpha=0.75)
    x_range = getattr(xaxis, "range", None)
    y_range = getattr(yaxis, "range", None)
    if x_range and len(x_range) == 2:
        axis.set_xlim(float(x_range[0]), float(x_range[1]))
    if y_range and len(y_range) == 2:
        axis.set_ylim(float(y_range[0]), float(y_range[1]))
    if has_legend:
        legend = axis.legend(facecolor=plot_color, edgecolor=grid_color)
        for text in legend.get_texts():
            text.set_color(font_color)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=140, facecolor=paper_color)


def _positive_int(value: Any, label: str, *, maximum: int | None = None) -> int:
    numeric = _finite_float(value, label)
    if int(numeric) != numeric or numeric < 1:
        raise ValueError(f"{label} must be an integer >= 1")
    result = int(numeric)
    if maximum is not None and result > maximum:
        raise ValueError(f"{label} must be <= {maximum}")
    return result


def _valid_tune_config(config: str) -> str:
    value = str(config or "default").strip() or "default"
    available = {
        str(item.get("value") or "").strip()
        for item in available_tunable_configs()
        if str(item.get("value") or "").strip()
    }
    if value not in available:
        raise ValueError(f"unknown tunable configuration: {value}")
    return value


def _normalize_tune_override(overrides: str) -> str:
    value = normalize_override_string(overrides)
    if not value:
        return ""
    if any(character in value for character in "\n\r;|&`$<>"):
        raise ValueError("tuning overrides may contain only comma-separated namelist key=value pairs")
    pairs = parse_override_pairs(value)
    if not pairs or any(not _TUNE_PARAM_RE.fullmatch(str(key).strip()) or not str(raw).strip() for key, raw in pairs):
        raise ValueError("tuning overrides must be comma-separated key=value pairs")
    return value


def _normalize_tune_cases(raw_cases: list[Any], case_data: dict[str, Any]) -> tuple[list[str], list[dict[str, Any]]]:
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("provide at least one tuning case")

    names: list[str] = []
    configurations: list[dict[str, Any]] = []
    for item in raw_cases:
        if isinstance(item, str):
            raw_item: dict[str, Any] = {"case_name": item}
        elif isinstance(item, dict):
            raw_item = dict(item)
            if "case" in raw_item and "case_name" not in raw_item:
                raw_item["case_name"] = raw_item.pop("case")
        else:
            raise ValueError("each tuning case must be a case name or an object")

        name = _validate_case(str(raw_item.get("case_name") or ""))
        if name in names:
            raise ValueError(f"duplicate tuning case: {name}")
        if name not in case_data:
            raise ValueError(f"no tuner benchmark defaults are available for case: {name}")

        allowed = {
            "case_name",
            "time_average_range",
            "altitude_comparison_range",
            "average_time_seconds",
            "num_time_windows",
        }
        unknown = sorted(set(raw_item) - allowed)
        if unknown:
            raise ValueError(f"unsupported tuning case setting(s) for {name}: {', '.join(unknown)}")

        overrides = {key: raw_item[key] for key in allowed - {"case_name"} if key in raw_item}
        try:
            defaults = read_case_tuner_defaults(name, overrides=overrides)
        except (RuntimeError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid tuning window for {name}: {exc}") from exc
        config = {
            "case_name": name,
            "time_average_range": list(defaults["time_average_range"]),
            "altitude_comparison_range": list(defaults["altitude_comparison_range"]),
            "num_time_windows": int(defaults.get("num_time_windows", 1)),
        }
        if defaults.get("average_time_seconds") is not None:
            config["average_time_seconds"] = int(defaults["average_time_seconds"])
        names.append(name)
        configurations.append(config)
    return names, configurations


def _normalize_tune_ranges(raw_ranges: list[dict[str, Any]], config: str) -> list[dict[str, Any]]:
    if not isinstance(raw_ranges, list) or not raw_ranges:
        raise ValueError("provide at least one parameter range")
    if len(raw_ranges) > _TUNE_MAX_PARAMETERS:
        raise ValueError(f"at most {_TUNE_MAX_PARAMETERS} parameter ranges may be requested through the agent tool")
    names = load_tunable_names(config)
    allowed = set(names)
    ranges: list[dict[str, Any]] = []
    seen: set[str] = set()
    owned_targets: set[str] = set()
    for item in raw_ranges:
        if not isinstance(item, dict):
            raise ValueError("each parameter range must be an object with name, min, and max")
        name = canonical_tunable_parameter_name(item.get("name"), names)
        if name not in allowed:
            raise ValueError(f"{name or 'parameter'} is not available in config {config}")
        if name in seen:
            raise ValueError(f"duplicate tuning parameter: {name}")
        # Pydantic serializes an omitted optional ``targets`` field as None;
        # native Dash requests omit the key entirely.  Both mean the ordinary
        # one-target range.
        raw_targets = item.get("targets") or [name]
        if isinstance(raw_targets, str):
            raw_targets = [raw_targets]
        targets = [
            canonical_tunable_parameter_name(target, names)
            for target in (raw_targets or [])
            if str(target).strip()
        ]
        if not targets:
            raise ValueError(f"{name} requires at least one physical target")
        if targets[0] != name:
            raise ValueError(f"{name} must be the first target of its linked range")
        if len(targets) != len(set(targets)):
            raise ValueError(f"{name} repeats a physical target")
        unknown_targets = sorted(set(targets) - allowed)
        if unknown_targets:
            raise ValueError(
                f"linked target(s) unavailable in config {config}: {', '.join(unknown_targets)}"
            )
        duplicate_targets = sorted(owned_targets.intersection(targets))
        if duplicate_targets:
            raise ValueError(
                "a physical tuning parameter may only belong to one range: "
                + ", ".join(duplicate_targets)
            )
        low = _finite_float(item.get("min"), f"{name} min")
        high = _finite_float(item.get("max"), f"{name} max")
        if low > high:
            raise ValueError(f"{name} requires min <= max")
        ranges.append({"name": name, "targets": targets, "min": low, "max": high})
        seen.add(name)
        owned_targets.update(targets)
    return ranges


def _canonical_tune_request(request: TuneRequest) -> TuneRequest:
    """Resolve a typed Tune request against its selected live configuration.

    Do this before request-id persistence, so the immutable job record,
    manifest, worker request, and later Dash rendering all use exactly the
    current namelist name rather than a legacy spelling supplied by a client.
    """
    config = _valid_tune_config(request.config)
    ranges = _normalize_tune_ranges(
        [item.model_dump() for item in request.parameter_ranges],
        config,
    )
    resolution = evaluate_tune_settings(config, _typed_override_text(request.overrides))
    errors = [str(issue.get("message") or "") for issue in resolution.get("issues", []) if issue.get("severity") == "error"]
    if errors:
        raise ValueError("; ".join(errors))
    inactive = sorted(
        target
        for item in ranges
        for target in item["targets"]
        if not is_independently_tunable(resolution.get("parameter_states", {}).get(target))
    )
    if inactive:
        details = "; ".join(
            f"{name}: {resolution['parameter_states'][name]['reason']}" for name in inactive
        )
        raise ValueError("Tune request selects inactive parameter(s): " + details)
    ranges = apply_required_parameter_links(ranges, resolution)
    aggregation_mode = str(request.aggregation_mode or DEFAULT_AGGREGATION_MODE).strip()
    if aggregation_mode not in AGGREGATION_MODE_NAMES:
        raise ValueError("unknown aggregation mode")
    aggregation_scope = str(
        request.time_window_aggregation_scope or DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE
    ).strip()
    if aggregation_scope not in TIME_WINDOW_AGGREGATION_SCOPES:
        raise ValueError("unknown time-window aggregation scope")
    aggregation_weights = normalize_aggregation_weights(request.aggregation_weights)
    return request.model_copy(
        update={
            "config": config,
            "parameter_ranges": [ParameterRange.model_validate(item) for item in ranges],
            "aggregation_mode": aggregation_mode,
            "aggregation_weights": aggregation_weights,
            "time_window_aggregation_scope": aggregation_scope,
        }
    )


def _normalize_tune_strategy(
    strategy: str,
    ranges: list[dict[str, float | str]],
    *,
    max_samples: int,
    resolve_spacing: float,
    simann_max_iters: int,
    simann_initial_temp: float,
    simann_final_temp: float,
    adam_max_updates: int,
    adam_learning_rate: float,
    adam_perturbation: float,
    adam_spsa_pairs: int,
    batch_size: int,
    max_workers: int,
    max_samples_limit: int | None = _TUNE_MAX_RANDOM_SAMPLES,
) -> dict[str, Any]:
    name = str(strategy or "random").strip().lower()
    if name not in VALID_STRATEGY_NAMES:
        raise ValueError("strategy must be random, resolve, simann, or adam")
    if name == "random":
        return {"name": name, "options": {"max_samples": _positive_int(max_samples, "max_samples", maximum=max_samples_limit)}}
    if name == "resolve":
        spacing = _finite_float(resolve_spacing, "resolve_spacing")
        if spacing <= 0.0:
            raise ValueError("resolve_spacing must be > 0")
        estimated = 1
        for item in ranges:
            width = float(item["max"]) - float(item["min"])
            estimated *= int(width / spacing + 1.0 + 1.0e-12)
            if estimated > _TUNE_MAX_RESOLVE_SAMPLES:
                raise ValueError(
                    f"the resolve grid would exceed {_TUNE_MAX_RESOLVE_SAMPLES} samples; use a coarser spacing or fewer ranges"
                )
        return {"name": name, "options": {"spacing": spacing}}

    if name == "simann":
        initial_temp = _finite_float(simann_initial_temp, "simann_initial_temp")
        final_temp = _finite_float(simann_final_temp, "simann_final_temp")
        if initial_temp <= 0.0 or final_temp <= 0.0:
            raise ValueError("SimAnn temperatures must be > 0")
        return {
            "name": name,
            "options": {
                "max_iters": _positive_int(simann_max_iters, "simann_max_iters", maximum=_TUNE_MAX_SIMANN_ITERS),
                "initial_temp": initial_temp,
                "max_final_temp": final_temp,
                "chain_count": max(1, max_workers * batch_size),
            },
        }

    learning_rate = _finite_float(adam_learning_rate, "adam_learning_rate")
    perturbation = _finite_float(adam_perturbation, "adam_perturbation")
    spsa_pairs = _positive_int(adam_spsa_pairs, "adam_spsa_pairs")
    columns_per_chain = 2 * spsa_pairs
    if learning_rate <= 0.0:
        raise ValueError("adam_learning_rate must be > 0")
    if not 0.0 < perturbation <= 0.5:
        raise ValueError("adam_perturbation must be in (0, 0.5]")
    if batch_size % columns_per_chain != 0:
        raise ValueError(
            "Adam requires batch_size to be divisible by 2 * spsa_pairs "
            f"({batch_size} is not divisible by {columns_per_chain})"
        )
    return {
        "name": name,
        "options": {
            "max_updates": _positive_int(
                adam_max_updates, "adam_max_updates", maximum=_TUNE_MAX_ADAM_UPDATES
            ),
            "learning_rate": learning_rate,
            "perturbation": perturbation,
            "spsa_pairs": spsa_pairs,
        },
    }


def _typed_override_text(values: dict[str, Any]) -> str:
    """Convert strict public key/value overrides to the legacy runner form."""

    pairs: list[str] = []
    for key, value in sorted(dict(values or {}).items()):
        name = str(key).strip()
        if not _TUNE_PARAM_RE.fullmatch(name):
            raise ValueError(f"invalid override name: {name or '<empty>'}")
        if isinstance(value, bool):
            rendered = "true" if value else "false"
        elif isinstance(value, (str, int, float)) and str(value).strip():
            rendered = str(value).strip()
        else:
            raise ValueError(f"invalid override value for {name}")
        if any(character in rendered for character in "\n\r;|&`$<>"):
            raise ValueError(f"invalid override value for {name}")
        pairs.append(f"{name}={rendered}")
    # ``create_case_namelist.parse_override_pairs`` deliberately uses commas
    # as assignment separators because values may contain spaces.  Keep the
    # typed MCP form on that exact wire format: a space-separated rendering
    # would silently turn every later ``name=value`` into the first value.
    return ",".join(pairs)


def _validated_stats_file(stats_file: str) -> str:
    value = str(stats_file or DEFAULT_STATS_NAME).strip()
    if value == NO_STATS_NAME:
        return value
    if value not in set(list_stats_files()):
        raise ValueError("stats_file must be one of the checked-in input/stats files or 'none'")
    return value


def _active_scm_processes(*, exclude_job_id: str | None = None) -> set[tuple[str, str]]:
    """Return active broker-owned SCM children, deduplicated by PID."""
    active: set[tuple[str, str]] = set()
    for item in _JOB_STORE.list_kind("scm"):
        if str(item.get("job_id") or "") == str(exclude_job_id or ""):
            continue
        if str(item.get("state") or "") not in {"starting", "running", "stopping"}:
            continue
        runtime = dict(item.get("runtime") or {})
        pid = (runtime.get("proc_data") or {}).get("pid")
        case_name = str((item.get("request") or {}).get("case") or "")
        active.add(("pid", str(pid)) if pid else ("job", str(item.get("job_id") or case_name)))
    return active


_SCM_RESERVED_STATES = {"queued", "submitting", "starting", "running", "stopping"}


def _matching_scm_records(
    case_name: str,
    output_directory: str | Path | None = None,
    *,
    states: set[str] | None = None,
    exclude_job_id: str | None = None,
    records: Any = None,
) -> list[dict[str, Any]]:
    """Return case/output matches, conservatively including old unknown paths."""
    target = _canonical_output_directory(output_directory)
    matches = []
    for raw_record in records if records is not None else _JOB_STORE.list_kind("scm"):
        if not isinstance(raw_record, dict) or raw_record.get("kind") != "scm":
            continue
        record = dict(raw_record)
        if str(record.get("job_id") or "") == str(exclude_job_id or ""):
            continue
        if states is not None and str(record.get("state") or "") not in states:
            continue
        if str((record.get("request") or {}).get("case") or "") != str(case_name):
            continue
        record_target = _scm_record_output_directory(record)
        # Pre-migration active records without a known target remain global
        # blockers rather than risking two writers for one case.
        if target and record_target and record_target != target:
            continue
        matches.append(record)
    return matches


def _scm_output_conflict(
    case_name: str,
    output_directory: str | Path,
    *,
    exclude_job_id: str | None = None,
    records: Any = None,
) -> dict[str, Any] | None:
    """Return the existing reservation for a case/output pair, if any."""
    matches = _matching_scm_records(
        case_name,
        output_directory,
        states=_SCM_RESERVED_STATES,
        exclude_job_id=exclude_job_id,
        records=records,
    )
    return max(
        matches,
        key=lambda record: float(
            record.get("updated_at_unix_seconds")
            or record.get("created_at_unix_seconds")
            or 0
        ),
        default=None,
    )


def scm_run_summary() -> dict[str, Any]:
    """Return compact authoritative counts for global dashboard labels."""
    return _JOB_STORE.scm_summary()


def run_telemetry(
    known_revision: Any = None,
    log_cursors: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return collect_run_telemetry(
        known_revision, log_cursors, store=_JOB_STORE
    )


def _assert_scm_admission(
    case_name: str,
    output_directory: str | Path,
    *,
    exclude_job_id: str | None = None,
) -> None:
    """Apply output-scoped identity and global process admission rules."""
    target = _canonical_output_directory(output_directory)
    if _scm_output_conflict(
        case_name, target, exclude_job_id=exclude_job_id
    ) is not None:
        raise ValueError(f"{case_name} is already active in {target}")
    if len(_active_scm_processes(exclude_job_id=exclude_job_id)) >= MAX_RUN_PROCS:
        raise ValueError(f"maximum active SCM runs is {MAX_RUN_PROCS}")


def _normalize_dashboard_overrides(overrides: dict[str, Any] | None) -> dict[str, dict[str, str]]:
    """Validate native Run-tab namelist deltas before the broker writes them."""
    result: dict[str, dict[str, str]] = {"flags": {}, "tunable": {}, "silhs": {}}
    raw = dict(overrides or {})
    for group in result:
        entries = raw.get(group) or {}
        if not isinstance(entries, dict):
            raise ValueError(f"{group} overrides must be a name/value mapping")
        for name, value in entries.items():
            key = str(name).strip()
            text = str(value).strip()
            if not _NAMELIST_KEY_RE.fullmatch(key) or not text or len(text) > 256:
                raise ValueError(f"invalid {group} namelist override")
            result[group][key] = text
    return result


def _canonical_scm_request(request: ScmRunRequest) -> tuple[ScmRunRequest, dict[str, Any]]:
    """Resolve typed SCM settings before sealing the immutable MCP request.

    The native Run tab and the typed MCP surface must launch with the same
    forced PDF/parameter consequences.  Keep unknown override names intact:
    the normal namelist builder remains the authority for non-settings values,
    while known flags and tunables are resolved here.
    """
    from dash_app.run_tab.namelist import is_bool_value, is_true, read_namelist_entries

    config = _valid_tune_config(request.config)
    flag_entries = read_namelist_entries(tunable_config_file(config, "configurable_model_flags.in"))
    parameter_entries = read_namelist_entries(tunable_config_file(config, "tunable_parameters.in"))
    flags = {
        entry["name"]: is_true(entry["value"]) if is_bool_value(entry["value"]) else entry["value"]
        for entry in flag_entries
    }
    parameter_names = {entry["name"] for entry in parameter_entries}
    resolved_parameters: dict[str, Any] = {}
    effective_overrides: dict[str, Any] = {}
    for raw_name, value in dict(request.overrides or {}).items():
        flag_name = canonical_flag_name(raw_name)
        parameter_name = canonical_parameter_name(raw_name)
        if flag_name in flags:
            flags[flag_name] = value
            effective_overrides[flag_name] = value
        elif parameter_name in parameter_names:
            # Pass only explicit parameter edits to the resolver.  Defaults
            # must not masquerade as user-supplied values, otherwise a
            # one-sided equality rule (C6rt -> C6thl) cannot identify its
            # follower.
            resolved_parameters[parameter_name] = value
            effective_overrides[parameter_name] = value
        else:
            effective_overrides[str(raw_name).strip()] = value

    resolution = resolve_clubb_settings(flags, resolved_parameters, auto_correct=True)
    errors = [issue.message for issue in resolution.issues if issue.severity == "error"]
    if errors:
        raise ValueError("; ".join(errors))
    effective_overrides.update(resolution.forced_flags)
    effective_overrides.update(resolution.forced_parameters)
    return request.model_copy(update={"config": config, "overrides": effective_overrides}), resolution.as_dict()


def _normalize_dashboard_cli_options(cli_options: dict[str, Any] | None) -> dict[str, Any]:
    """Keep native Run controls explicit while preserving their current inputs."""
    raw = dict(cli_options or {})
    allowed = {
        "multicol", "batch_size", "max_iters", "debug", "dt_main", "dt_rad",
        "tout", "out_dir", "extra_args", "implementation", "install_dir",
    }
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError("unsupported run option(s): " + ", ".join(unknown))
    normalized: dict[str, Any] = {}
    for key in allowed - {"extra_args"}:
        if key not in raw or raw[key] in {None, ""}:
            continue
        value = str(raw[key]).strip()
        if not value or len(value) > 512:
            raise ValueError(f"invalid {key} run option")
        normalized[key] = value
    implementation = str(normalized.get("implementation") or "fortran").lower()
    if implementation not in {"fortran", "python", "jax"}:
        raise ValueError("implementation must be fortran, python, or jax")
    if "implementation" in normalized:
        normalized["implementation"] = implementation
    if normalized.get("install_dir"):
        install_dir = Path(normalized["install_dir"]).expanduser().resolve()
        install_root = (REPO_ROOT / "install").resolve()
        if not install_dir.is_relative_to(install_root) or not install_dir.is_dir():
            raise ValueError("install_dir must name an installed CLUBB build")
        normalized["install_dir"] = str(install_dir)
    extra = raw.get("extra_args") or []
    if isinstance(extra, str):
        extra = shlex.split(extra)
    if not isinstance(extra, (list, tuple)):
        raise ValueError("extra_args must be a list of command arguments")
    tokens = [str(item).strip() for item in extra if str(item).strip()]
    if len(tokens) > 64 or any(len(item) > 512 or "\x00" in item for item in tokens):
        raise ValueError("invalid extra run arguments")
    managed = {"-python", "--python", "-jax", "--jax", "-exe", "--exe", "-install_dir", "--install_dir"}
    if any(token in managed or any(token.startswith(option + "=") for option in managed) for token in tokens):
        raise ValueError("implementation and install directory are managed by the build selector")
    if tokens:
        normalized["extra_args"] = tokens
    return normalized


def _persist_submission(kind: str, request_id: str, payload: dict[str, Any]):
    try:
        return _JOB_STORE.submit(kind, request_id, payload)
    except SubmissionConflict as exc:
        raise ValueError(f"REQUEST_ID_CONFLICT: {exc}") from exc


def _run_common_manifest_inputs(stats_file: str, config: str) -> dict[str, Any]:
    """Checksum batch-wide SCM inputs once."""
    paths = {
        "stats": None if stats_file == NO_STATS_NAME else REPO_ROOT / "input" / "stats" / stats_file,
        "clubb_params": tunable_config_file(config, "tunable_parameters.in"),
        "model_flags": tunable_config_file(config, "configurable_model_flags.in"),
        "silhs_params": tunable_config_file(config, "silhs_parameters.in"),
    }
    return {
        name: {"path": str(path), "sha256": sha256_file(Path(path))}
        for name, path in paths.items()
        if path is not None and Path(path).is_file()
    }


def _run_manifest_inputs(
    case: str,
    stats_file: str,
    config: str,
    *,
    common_inputs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Capture one case input plus batch-wide SCM inputs."""
    inputs = dict(common_inputs or _run_common_manifest_inputs(stats_file, config))
    case_path = REPO_ROOT / "input" / "case_setups" / f"{case}_model.in"
    if case_path.is_file():
        inputs["case_setup"] = {
            "path": str(case_path),
            "sha256": sha256_file(case_path),
        }
    return inputs


def _scm_build_identity(cli_options: dict[str, Any] | None = None) -> dict[str, Any]:
    """Capture the exact compiled executable selected by ``run_scm.py``.

    The regular SCM runner resolves ``install/selected`` before
    ``install/latest``.  Storing its resolved path and digest makes a later
    rebuild distinguishable without relying on mutable UI selection or a log.
    """
    options = dict(cli_options or {})
    implementation = str(options.get("implementation") or "fortran").lower()
    selected = REPO_ROOT / "install" / "selected"
    latest = REPO_ROOT / "install" / "latest"
    install = Path(options["install_dir"]) if options.get("install_dir") else (
        selected if os.path.lexists(selected) else latest
    )
    resolved_install = install.resolve(strict=False)
    executable = resolved_install / "clubb_standalone"
    return {
        "install_selector": str(install),
        "install_directory": str(resolved_install),
        "implementation": implementation,
        "executable": {
            "path": str(executable),
            "sha256": sha256_file(executable),
            "bytes": executable.stat().st_size if executable.is_file() else None,
        },
        "environment": {
            "FC": os.environ.get("FC", ""),
            "CC": os.environ.get("CC", ""),
            "LMOD_FAMILY_COMPILER": os.environ.get("LMOD_FAMILY_COMPILER", ""),
            "LOADEDMODULES": os.environ.get("LOADEDMODULES", ""),
        },
        "precision": "double",
    }


def submit_compile(request: CompileRequest) -> dict[str, Any]:
    """Idempotently submit a normal broker-owned compile job."""

    record, created = _persist_submission("compile", request.request_id, request.model_dump(mode="json"))
    if not created:
        return {"status": "existing", **record}
    try:
        if (broker_jobs().get("compile") or {}).get("state") in {"running", "stopping"}:
            raise ValueError("a CLUBB compile is already active")
        # Seal the requested source/options before a child can observe a later
        # edit.  Completion data is written as a separate compact summary.
        manifest = _ARTIFACT_STORE.create_manifest(
            record["job_id"],
            {
                "job": record,
                "requested_compile": request.model_dump(mode="json"),
                "execution": {
                    "state": "planned",
                    "runner": "compile.py",
                    "build_identity": {
                        "FC": os.environ.get("FC", ""),
                        "CC": os.environ.get("CC", ""),
                        "LMOD_FAMILY_COMPILER": os.environ.get("LMOD_FAMILY_COMPILER", ""),
                        "LOADEDMODULES": os.environ.get("LOADEDMODULES", ""),
                        "precision": "double",
                    },
                },
            },
            active=True,
        )
        _JOB_STORE.update(record["job_id"], state="starting", manifest_path=str(manifest))
        result = compile_clubb(
            debug=request.debug,
            python_bindings=request.python_bindings,
            fresh=request.fresh,
            gptl=request.gptl,
            job_id=record["job_id"],
        )
        _ARTIFACT_STORE.write_summary(record["job_id"], "submission_result.json", {"state": "started", "runtime": result})
        updated = _JOB_STORE.update(record["job_id"], state="running", runtime=result) or record
        return {"status": "started", **updated}
    except Exception as exc:
        _ARTIFACT_STORE.release(record["job_id"])
        _JOB_STORE.update(record["job_id"], state="error", error={"code": "COMPILE_SUBMISSION_FAILED", "message": str(exc)})
        raise


def _start_scm_submission(
    request: ScmRunRequest,
    record: dict[str, Any],
    stats_file: str,
    settings_resolution: dict[str, Any],
    *,
    output_dir: Path,
    batch_job_id: str,
    native_overrides: dict[str, Any] | None = None,
    native_cli_options: dict[str, Any] | None = None,
    common_inputs: dict[str, Any] | None = None,
    build_identity: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Prepare one SCM child, then atomically launch it if still active."""
    try:
        _assert_scm_admission(
            request.case,
            output_dir,
            exclude_job_id=str(record.get("job_id") or ""),
        )
    except ValueError as exc:
        code = (
            "CASE_OUTPUT_ACTIVE"
            if "already active" in str(exc)
            else "SCM_CONCURRENCY_LIMIT"
        )
        _JOB_STORE.update(record["job_id"], state="rejected", error={"code": code, "message": str(exc)})
        raise
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        manifest = _ARTIFACT_STORE.create_manifest(
            record["job_id"],
            {
                "job": record,
                "run_id": record["run_id"],
                "case": request.case,
                "batch_job_id": batch_job_id,
                "stats_file": stats_file,
                "configuration": request.config,
                "overrides": request.overrides,
                "settings_resolution": settings_resolution,
                "output_directory": str(output_dir),
                "vertical_coordinate": {"name": "zt", "units": "m"},
                "inputs": _run_manifest_inputs(
                    request.case,
                    stats_file,
                    request.config,
                    common_inputs=common_inputs,
                ),
                "execution": {
                    "state": "planned",
                    "runner": "run_scripts/run_scm.py",
                    "build_identity": build_identity or _scm_build_identity(),
                },
                "output_checksums": {
                    "state": "pending",
                    "result_resource": f"clubb-artifact://{batch_job_id}/execution_result.json",
                },
            },
            active=True,
        )
        _JOB_STORE.update(
            record["job_id"],
            output_directory=str(output_dir),
            manifest_path=str(manifest),
            batch_job_id=batch_job_id,
        )
        cli_options = dict(
            native_cli_options
            if native_cli_options is not None
            else request.run_options.model_dump(exclude_none=True)
        )
        if native_overrides is None:
            override_text = _typed_override_text(request.overrides)
            if override_text:
                cli_options["extra_args"] = [
                    *(cli_options.get("extra_args") or []),
                    "-override",
                    override_text,
                ]
            launch_overrides = {"flags": {}, "tunable": {}, "silhs": {}}
        else:
            launch_overrides = native_overrides
        cli_options["out_dir"] = str(output_dir)
        # Setup above may be relatively slow.  Serialize only the final state
        # check and process spawn so cancellation remains prompt and cannot
        # race a child from ``starting`` into an untracked running process.
        with _BATCH_QUEUE_LOCK:
            current = _JOB_STORE.get(record["job_id"]) or {}
            if str(current.get("state") or "") not in {"submitting", "starting"}:
                raise ValueError("SCM launch was cancelled")
            result = _launch_scm_process(
                request.case,
                stats_file,
                request.config,
                launch_overrides,
                cli_options,
                job_id=record["job_id"],
                run_id=record["run_id"],
                batch_job_id=batch_job_id,
            )
            updated = _JOB_STORE.update(
                record["job_id"], state="running", runtime=result
            ) or current
        _ARTIFACT_STORE.write_summary(record["job_id"], "submission_result.json", {"state": "started", "runtime": result})
        proc_data = dict(result.get("proc_data") or {})
        if proc_data:
            # Persist running before the watcher can observe a very short run
            # finishing.  Otherwise its terminal update can be overwritten by
            # this routine after the watcher has already exited.
            _background(_watch_run, request.case, proc_data, record["job_id"])
        return {"status": "started", **updated}
    except Exception as exc:
        _ARTIFACT_STORE.release(record["job_id"])
        current = _JOB_STORE.get(record["job_id"]) or {}
        if str(current.get("state") or "") not in {"cancelled", "stopping"}:
            _JOB_STORE.update(
                record["job_id"],
                state="error",
                error={"code": "SCM_SUBMISSION_FAILED", "message": str(exc)},
            )
        raise


def _canonical_scm_batch_request(request: ScmRunBatchRequest) -> tuple[ScmRunBatchRequest, dict[str, Any]]:
    """Apply the singular SCM validator once to common batch settings."""
    cases = [_validate_case(case) for case in request.cases]
    if len(cases) != len(set(cases)):
        raise ValueError("SCM batch cases must be unique")
    representative = ScmRunRequest(
        request_id="canonical-scm-batch",
        case=cases[0],
        stats_file=request.stats_file,
        config=request.config,
        overrides=request.overrides,
        run_options=request.run_options,
        out_dir=request.out_dir,
    )
    canonical, resolution = _canonical_scm_request(representative)
    return request.model_copy(
        update={
            "cases": cases,
            "stats_file": _validated_stats_file(request.stats_file),
            "config": canonical.config,
            "overrides": canonical.overrides,
        }
    ), resolution


def _batch_child_request(request: ScmRunBatchRequest, parent: dict[str, Any], case: str) -> ScmRunRequest:
    """Create a deterministic internal child request for retry-safe batches."""
    return ScmRunRequest(
        request_id=f"{parent['job_id']}-{case}",
        case=case,
        stats_file=request.stats_file,
        config=request.config,
        overrides=request.overrides,
        run_options=request.run_options,
        out_dir=request.out_dir,
    )


def _start_one_queued_scm_batch_child(
    batch_job_id: str,
) -> tuple[dict[str, Any] | None, bool]:
    """Claim and start at most one child without blocking cancellation setup."""
    child_record = None
    rejected_error = None
    with _BATCH_QUEUE_LOCK:
        batch = _batch_record(batch_job_id)
        if batch is None or batch.get("cancellation_requested") or _batch_terminal(str(batch.get("state") or "")):
            return _materialize_scm_batch(batch), False
        request = ScmRunBatchRequest.model_validate(batch.get("request") or {})
        children = _batch_children(batch)
        active = sum(
            str(child.get("state") or "") in {"starting", "running", "stopping"}
            for child in children
        )
        if active >= (request.max_workers or MAX_RUN_PROCS):
            return _materialize_scm_batch(batch), False
        for candidate in children:
            if str(candidate.get("state") or "") != "queued":
                continue
            child_job_id = str(candidate.get("job_id") or "")
            case = str((candidate.get("request") or {}).get("case") or "")
            try:
                _assert_scm_admission(
                    case,
                    str(candidate.get("output_directory") or ""),
                    exclude_job_id=child_job_id,
                )
            except ValueError as exc:
                if "maximum active SCM runs" in str(exc):
                    return _materialize_scm_batch(batch), False
                _JOB_STORE.update(
                    child_job_id,
                    state="error",
                    error={"code": "SCM_BATCH_CHILD_REJECTED", "message": str(exc)},
                )
                rejected_error = exc
                break
            child_record = _JOB_STORE.update(child_job_id, state="starting") or candidate
            break

    if rejected_error is not None:
        return _refresh_scm_batch(batch_job_id), True
    if child_record is None:
        return _refresh_scm_batch(batch_job_id), False

    child_job_id = str(child_record.get("job_id") or "")
    child_request = ScmRunRequest.model_validate(child_record.get("request") or {})
    try:
        _start_scm_submission(
            child_request,
            child_record,
            child_request.stats_file,
            dict(batch.get("settings_resolution") or {}),
            output_dir=Path(str(batch["output_directory"])),
            batch_job_id=batch_job_id,
            native_overrides=batch.get("native_overrides"),
            native_cli_options=batch.get("native_cli_options"),
            common_inputs=batch.get("common_inputs"),
            build_identity=batch.get("build_identity"),
        )
    except Exception as exc:
        current = _JOB_STORE.get(child_job_id) or {}
        if str(current.get("state") or "") not in {"cancelled", "stopping"}:
            updates: dict[str, Any] = {"state": "error"}
            if not current.get("error"):
                updates["error"] = {
                    "code": "SCM_BATCH_CHILD_FAILED",
                    "message": str(exc)[:500],
                }
            _JOB_STORE.update(child_job_id, **updates)
    return _refresh_scm_batch(batch_job_id), True


def _start_queued_scm_batch_children(batch_job_id: str) -> dict[str, Any] | None:
    """Fill available slots while yielding between launches for instant cancel."""
    if _BROKER_SHUTTING_DOWN.is_set():
        return _materialize_scm_batch(_batch_record(batch_job_id))
    latest = _materialize_scm_batch(_batch_record(batch_job_id))
    while True:
        latest, attempted = _start_one_queued_scm_batch_child(batch_job_id)
        if not attempted:
            return latest
        time.sleep(0)


def _ensure_scm_batch_watcher(batch_job_id: str) -> None:
    """Keep a durable batch queue moving even when no browser is polling."""
    with _BATCH_WATCHERS_LOCK:
        if batch_job_id in _BATCH_WATCHERS:
            return
        _BATCH_WATCHERS.add(batch_job_id)
    _background(_watch_scm_batch_queue, batch_job_id)


def _watch_scm_batch_queue(batch_job_id: str) -> None:
    try:
        while True:
            batch = _batch_record(batch_job_id)
            if batch is None or batch.get("cancellation_requested") or _batch_terminal(str(batch.get("state") or "")):
                return
            updated = _start_queued_scm_batch_children(batch_job_id) or batch
            queued = any(str(child.get("state") or "") == "queued" for child in updated.get("children") or [])
            if not queued:
                return
            time.sleep(0.2)
    finally:
        with _BATCH_WATCHERS_LOCK:
            _BATCH_WATCHERS.discard(batch_job_id)


def submit_scm_batch(
    request: ScmRunBatchRequest,
    *,
    native_overrides: dict[str, Any] | None = None,
    native_cli_options: dict[str, Any] | None = None,
    origin: str = "unknown",
    visibility: str = "user",
) -> dict[str, Any]:
    """Submit one durable multi-case SCM group into one shared output directory."""
    request, settings_resolution = _canonical_scm_batch_request(request)
    normalized_native_overrides = None
    normalized_native_cli_options = None
    native_output_directory: Path | None = None
    if native_overrides is not None or native_cli_options is not None:
        if native_overrides is not None:
            normalized_native_overrides = _normalize_dashboard_overrides(native_overrides)
        normalized_native_cli_options = _normalize_dashboard_cli_options(native_cli_options or {})
        # Native Dash owns the Run-tab output field, so preserve its selected
        # location after applying the repository's normal path rules. Public
        # MCP requests do not enter this branch and remain broker-controlled.
        requested_output = normalized_native_cli_options.get("out_dir") or "output"
        native_output_directory = resolve_output_dir(requested_output).resolve()
        normalized_native_cli_options["out_dir"] = str(native_output_directory)
    payload = request.model_dump(mode="json") | {"stats_file": request.stats_file}
    submission_origin = str(origin or "unknown")[:32]
    submission_visibility = "internal" if visibility == "internal" else "user"
    common_inputs = _run_common_manifest_inputs(request.stats_file, request.config)
    build_identity = _scm_build_identity(normalized_native_cli_options)

    def create_batch(transaction):
        parent, created = transaction.submit("scm_batch", request.request_id, payload)
        if not created:
            return parent, False
        batch_id = str(parent["batch_id"])
        output_directory = (
            native_output_directory
            or (_resolve_mcp_output_dir(request.out_dir) if request.out_dir else None)
            or _public_scm_batch_output_dir(batch_id)
        )
        output_directory = Path(output_directory).resolve()
        child_job_ids: list[str] = []
        accepted_cases: list[str] = []
        skipped_cases: list[dict[str, Any]] = []
        parent_result_summary = str(
            _ARTIFACT_STORE.bundle_dir(parent["job_id"]) / "execution_result.json"
        )
        for case in request.cases:
            conflict = _scm_output_conflict(
                case,
                output_directory,
                records=transaction.state.values(),
            )
            if conflict is not None:
                skipped_cases.append(
                    {
                        "case": case,
                        "code": "CASE_OUTPUT_ACTIVE",
                        "output_directory": str(output_directory),
                        "conflicting_job_id": str(conflict.get("job_id") or ""),
                        "conflicting_run_id": str(conflict.get("run_id") or ""),
                    }
                )
                continue
            child_request = _batch_child_request(request, parent, case)
            child_payload = child_request.model_dump(mode="json") | {
                "stats_file": request.stats_file
            }
            child, child_created = transaction.submit(
                "scm", child_request.request_id, child_payload
            )
            if not child_created:
                existing_batch = str(child.get("batch_job_id") or "")
                if existing_batch and existing_batch != parent["job_id"]:
                    raise ValueError("SCM batch child request is already owned by another batch")
            child = transaction.update(
                child["job_id"],
                state="queued",
                batch_id=batch_id,
                batch_job_id=parent["job_id"],
                output_directory=str(output_directory),
                result_summary_job_id=parent["job_id"],
                result_summary_path=parent_result_summary,
                origin=submission_origin,
                visibility=submission_visibility,
                display={
                    "case": case,
                    "stats_file": request.stats_file,
                    "config": request.config,
                    "cli_options": normalized_native_cli_options
                    or request.run_options.model_dump(exclude_none=True),
                },
            ) or child
            child_job_ids.append(str(child["job_id"]))
            accepted_cases.append(case)
        manifest_path = _ARTIFACT_STORE.bundle_dir(parent["job_id"]) / "manifest.json"
        state = "queued" if child_job_ids else "finished"
        parent = transaction.update(
            parent["job_id"],
            state=state,
            batch_id=batch_id,
            output_directory=str(output_directory),
            manifest_path=str(manifest_path),
            child_job_ids=child_job_ids,
            accepted_cases=accepted_cases,
            skipped_cases=skipped_cases,
            outcome="queued" if child_job_ids else "no_op",
            settings_resolution=settings_resolution,
            common_inputs=common_inputs,
            build_identity=build_identity,
            native_overrides=normalized_native_overrides,
            native_cli_options=normalized_native_cli_options,
            origin=submission_origin,
            visibility=submission_visibility,
            finished_at_unix_seconds=None if child_job_ids else time.time(),
        ) or parent
        return parent, True

    try:
        updated, created = _JOB_STORE.transaction(create_batch)
    except SubmissionConflict as exc:
        raise ValueError(f"REQUEST_ID_CONFLICT: {exc}") from exc
    if not created:
        existing = _materialize_scm_batch(updated) or updated
        return {"status": "existing", **existing}

    updated = _materialize_scm_batch(updated) or updated

    _ARTIFACT_STORE.create_manifest(
        updated["job_id"],
        {
            "job": updated,
            "batch_id": updated["batch_id"],
            "requested_batch": payload,
            "output_directory": updated["output_directory"],
            "settings_resolution": settings_resolution,
            "execution": {"state": "planned", "runner": "run_scripts/run_scm.py"},
        },
        active=bool(updated.get("children")),
    )
    # Return as soon as the complete durable queue exists.  The detached
    # watcher owns admission and launch, keeping the browser request fast even
    # when many SCM processes need setup work.
    if updated.get("children"):
        _ensure_scm_batch_watcher(updated["job_id"])
        return {"status": "queued", **updated}
    return {"status": "no_op", **updated}


def submit_scm_run(
    request: ScmRunRequest,
    *,
    native_cli_options: dict[str, Any] | None = None,
    origin: str = "unknown",
    visibility: str = "user",
) -> dict[str, Any]:
    """Backward-compatible one-case wrapper over the shared batch service."""
    legacy = _JOB_STORE.get_submission("scm", request.request_id)
    if legacy is not None and not legacy.get("batch_job_id"):
        return {"status": "existing", **legacy}
    batch_request = ScmRunBatchRequest(
        request_id=request.request_id,
        cases=[request.case],
        stats_file=request.stats_file,
        config=request.config,
        overrides=request.overrides,
        run_options=request.run_options,
        out_dir=request.out_dir,
    )
    batch = submit_scm_batch(
        batch_request,
        native_cli_options=native_cli_options,
        origin=origin,
        visibility=visibility,
    )
    submission_status = str(batch.get("status") or "")
    # The public one-case API historically returns the launched child.  Keep
    # that contract while native multi-case Run submissions stay asynchronous.
    if submission_status != "existing" and str(batch.get("state") or "") == "queued":
        batch = _start_queued_scm_batch_children(str(batch.get("job_id") or "")) or batch
    child = next((item for item in batch.get("children") or [] if item.get("case") == request.case), None)
    if child is None:
        return batch
    return {
        "status": (
            "existing"
            if submission_status == "existing"
            else "started"
            if child.get("state") in {"starting", "running"}
            else submission_status or batch.get("state")
        ),
        **child,
        "batch_id": batch.get("batch_id"),
        "batch_job_id": batch.get("job_id"),
    }


def submit_tune(request: TuneRequest) -> dict[str, Any]:
    """Idempotently submit a bounded Tune job through existing validation."""

    request = _canonical_tune_request(request)
    record, created = _persist_submission("tune", request.request_id, request.model_dump(mode="json"))
    if not created:
        return {"status": "existing", **record}
    try:
        recover_active_tuning_from_disk()
        if (broker_jobs().get("tune") or {}).get("state") in {"running", "stopping"}:
            raise ValueError("a Tune job is already active")
        manifest = _ARTIFACT_STORE.create_manifest(
            record["job_id"],
            {
                "job": record,
                "requested_tune": request.model_dump(mode="json"),
                "execution": {"state": "planned", "runner": "tuner", "build_identity": _scm_build_identity()},
            },
            active=True,
        )
        _JOB_STORE.update(record["job_id"], state="starting", manifest_path=str(manifest))
        result = launch_tuning(
            [item if isinstance(item, str) else item.model_dump(exclude_none=True) for item in request.cases],
            [item.model_dump() for item in request.parameter_ranges],
            request.fields,
            config=request.config,
            strategy=request.strategy,
            max_samples=request.max_samples,
            resolve_spacing=request.resolve_spacing,
            simann_max_iters=request.simann_max_iters,
            simann_initial_temp=request.simann_initial_temp,
            simann_final_temp=request.simann_final_temp,
            adam_max_updates=request.adam_max_updates,
            adam_learning_rate=request.adam_learning_rate,
            adam_perturbation=request.adam_perturbation,
            adam_spsa_pairs=request.adam_spsa_pairs,
            batch_size=request.batch_size,
            max_workers=request.max_workers,
            loss_mode=request.loss_mode or DEFAULT_LOSS_MODE,
            aggregation_mode=request.aggregation_mode or DEFAULT_AGGREGATION_MODE,
            aggregation_weights=request.aggregation_weights,
            time_window_aggregation_scope=request.time_window_aggregation_scope,
            overrides=_typed_override_text(request.overrides),
            job_id=record["job_id"],
        )
        _ARTIFACT_STORE.write_summary(record["job_id"], "submission_result.json", {"state": "started", "runtime": result})
        updated = _JOB_STORE.update(record["job_id"], state="running", runtime=result) or record
        return {"status": "started", **updated}
    except Exception as exc:
        _ARTIFACT_STORE.release(record["job_id"])
        _JOB_STORE.update(record["job_id"], state="error", error={"code": "TUNE_SUBMISSION_FAILED", "message": str(exc)})
        raise


def submit_leaderboard_rerun(request: LeaderboardRerunRequest) -> dict[str, Any]:
    """Idempotently rerun current Tune leaderboard rows with durable identity."""
    record, created = _persist_submission("leaderboard", request.request_id, request.model_dump(mode="json"))
    if not created:
        return {"status": "existing", **record}
    try:
        manifest = _ARTIFACT_STORE.create_manifest(
            record["job_id"],
            {
                "job": record,
                "requested_leaderboard_rerun": request.model_dump(mode="json"),
                "execution": {"state": "planned", "runner": "tuner loss runner", "build_identity": _scm_build_identity()},
            },
            active=True,
        )
        _JOB_STORE.update(record["job_id"], state="starting", manifest_path=str(manifest))
        result = run_tuning_loss(request.mode, request.max_results, job_id=record["job_id"])
        _ARTIFACT_STORE.write_summary(record["job_id"], "submission_result.json", {"state": "started", "runtime": result})
        updated = _JOB_STORE.update(record["job_id"], state="running", runtime=result) or record
        return {"status": "started", **updated}
    except Exception as exc:
        _ARTIFACT_STORE.release(record["job_id"])
        _JOB_STORE.update(record["job_id"], state="error", error={"code": "LEADERBOARD_SUBMISSION_FAILED", "message": str(exc)[:500]})
        raise


def recover_active_runs_from_state() -> list[dict[str, Any]]:
    """Re-adopt surviving broker-owned SCM PIDs after a broker restart.

    The original process is intentionally not recreated; the recovered watcher
    reads its log and owns cancellation/status from the durable record.  A
    terminal status after a broker gap is conservatively marked failed because
    POSIX does not preserve another process's exit code for later retrieval.
    """
    recovered: list[dict[str, Any]] = []
    for record in _JOB_STORE.list_kind("scm"):
        if str(record.get("state") or "") not in {"running", "stopping"}:
            continue
        runtime = dict(record.get("runtime") or {})
        proc_data = dict(runtime.get("proc_data") or {})
        job_id = str(record.get("job_id") or "")
        case = str((record.get("request") or {}).get("case") or "")
        if not proc_data or not _pid_is_alive(proc_data.get("pid")):
            message = "SCM process ended while its durable monitor was unavailable; exit status could not be recovered"
            if job_id:
                _JOB_STORE.update(
                    job_id,
                    state="error",
                    error={"code": "SCM_RECOVERY_EXIT_UNKNOWN", "message": message},
                    finished_at_unix_seconds=time.time(),
                )
                _ARTIFACT_STORE.release(job_id)
            continue
        _background(_watch_run, case, proc_data, job_id or None)
        recovered.append({"case": case, "pid": proc_data.get("pid"), "job_id": job_id or None})
    if recovered:
        publish_event("run", "Recovered SCM job monitoring", ", ".join(item["case"] for item in recovered), status="info")
    return recovered


def clear_terminal_scm_session() -> dict[str, Any]:
    """Clear terminal SCM control records and their temporary wrapper logs."""
    def prune(transaction):
        active_batches = {
            str(record.get("job_id") or "")
            for record in transaction.state.values()
            if isinstance(record, dict)
            and record.get("kind") == "scm_batch"
            and str(record.get("state") or "") in {"queued", "submitting", "starting", "running", "stopping"}
        }
        protected_children = {
            child_job_id
            for record in transaction.state.values()
            if isinstance(record, dict) and str(record.get("job_id") or "") in active_batches
            for child_job_id in _batch_child_job_ids(record)
        }
        terminal_records = [
            dict(record)
            for record in transaction.state.values()
            if isinstance(record, dict)
            and record.get("kind") in {"scm", "scm_batch"}
            and str(record.get("visibility") or "user") == "user"
            and str(record.get("state") or "") not in {"queued", "submitting", "starting", "running", "stopping"}
            and str(record.get("job_id") or "") not in protected_children
        ]
        terminal_ids = {
            str(record.get("job_id") or "") for record in terminal_records
        }
        transaction.delete(terminal_ids)
        return terminal_records

    removed_records = _JOB_STORE.transaction(prune)
    removed_logs = []
    for record in removed_records:
        runtime = dict(record.get("runtime") or {})
        log_path = str((runtime.get("proc_data") or {}).get("log") or "")
        if not log_path:
            continue
        try:
            Path(log_path).unlink()
            removed_logs.append(log_path)
        except FileNotFoundError:
            pass
    return {
        "status": "cleared",
        "removed_jobs": sorted(
            str(record.get("job_id") or "") for record in removed_records
        ),
        "removed_logs": removed_logs,
    }


def recover_queued_scm_batches() -> list[dict[str, Any]]:
    """Resume broker-owned batch queue advancement after a broker restart."""
    recovered: list[dict[str, Any]] = []
    for batch in _JOB_STORE.list_kind("scm_batch"):
        if str(batch.get("state") or "") not in {"queued", "running"}:
            continue
        if not any(
            str(child.get("state") or "") == "queued"
            for child in _batch_children(batch)
        ):
            continue
        updated = _start_queued_scm_batch_children(str(batch.get("job_id") or ""))
        if updated:
            _ensure_scm_batch_watcher(str(batch.get("job_id") or ""))
            recovered.append({"batch_id": updated.get("batch_id"), "job_id": updated.get("job_id"), "state": updated.get("state")})
    return recovered


def recover_active_compile_from_state() -> dict[str, Any] | None:
    """Resume compile log/status monitoring after a broker-only restart."""
    record = dict(broker_jobs().get("compile") or {})
    if str(record.get("state") or "") not in {"running", "stopping"}:
        return None
    job = dict(record.get("job") or {})
    job_id = str(record.get("job_id") or "")
    if not job or not _pid_is_alive(job.get("pid")):
        message = "compile process ended while its durable monitor was unavailable; exit status could not be recovered"
        update_broker_job(
            "compile",
            state="error",
            returncode=None,
            recovery_error=message,
        )
        if job_id:
            _JOB_STORE.update(
                job_id,
                state="error",
                error={"code": "COMPILE_RECOVERY_EXIT_UNKNOWN", "message": message},
                finished_at_unix_seconds=time.time(),
            )
            _ARTIFACT_STORE.release(job_id)
        return None
    _background(_watch_compile, job, job_id or None)
    publish_event("compile", "Recovered compile job monitoring", str(job.get("command") or job.get("pid") or ""), status="info")
    return {"pid": job.get("pid"), "job_id": job_id or None}


def get_job(job_id: str) -> dict[str, Any]:
    record = _JOB_STORE.get(str(job_id))
    if record is None:
        raise ValueError("unknown job_id")
    return record


def get_run_manifest(run_id: str) -> dict[str, Any]:
    record = _JOB_STORE.get_run(str(run_id))
    if record is None:
        raise ValueError("unknown run_id")
    manifest_path = Path(str(record.get("manifest_path") or ""))
    try:
        return json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError) as exc:
        raise ValueError("run manifest is unavailable") from exc


def get_artifact(artifact_id: str) -> dict[str, Any]:
    manifest = _ARTIFACT_STORE.get_manifest(str(artifact_id))
    if manifest is None:
        raise ValueError("unknown artifact_id")
    directory = _ARTIFACT_STORE.root / str(artifact_id)
    files = []
    try:
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.name != "manifest.json":
                files.append({"name": path.name, "uri": f"clubb-artifact://{artifact_id}/{path.name}", "bytes": path.stat().st_size})
    except OSError:
        pass
    return {"artifact_id": artifact_id, "uri": f"clubb-artifact://{artifact_id}/manifest.json", "manifest": manifest, "files": files}


def read_artifact_file(artifact_id: str, filename: str) -> bytes:
    """Read a deliberately small, whitelisted artifact resource payload."""
    safe_id = str(artifact_id or "")
    safe_name = str(filename or "")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", safe_id) or not re.fullmatch(r"[A-Za-z0-9_.-]+", safe_name):
        raise ValueError("artifact resource path is invalid")
    if Path(safe_name).name != safe_name or Path(safe_name).suffix.lower() not in {".png", ".json", ".log", ".txt"}:
        raise ValueError("artifact resource type is not available")
    path = _ARTIFACT_STORE.root / safe_id / safe_name
    if not path.is_file():
        raise ValueError("artifact resource is unavailable")
    if path.stat().st_size > 10 * 1024 * 1024:
        raise ValueError("artifact resource exceeds the 10 MiB transfer limit")
    return path.read_bytes()


def list_cases() -> dict[str, Any]:
    """List only checked-in SCM cases available to the public service."""
    cases = sorted(path.name.removesuffix("_model.in") for path in (REPO_ROOT / "input" / "case_setups").glob("*_model.in"))
    return {"cases": cases, "stats_files": sorted([NO_STATS_NAME, *list_stats_files()])}


def list_tunable_parameters(config: str = "default") -> dict[str, Any]:
    """Return the selected configuration's current tunable-parameter manifest.

    The manifest is parsed from the checked-in namelist on every call, so an
    agent never needs an embedded parameter list and a Dash source reload does
    not leave the service authoritative over stale UI state.
    """
    selected_config = _valid_tune_config(config)
    names = load_tunable_names(selected_config)
    ranges = load_tunable_default_ranges(selected_config)
    return {
        "config": selected_config,
        "parameters": [
            {
                "name": name,
                "default": (ranges.get(name) or {}).get("default"),
                "suggested_min": (ranges.get(name) or {}).get("min"),
                "suggested_max": (ranges.get(name) or {}).get("max"),
            }
            for name in names
        ],
    }


def resolve_clubb_settings_request(
    config: str = "default",
    flags: dict[str, Any] | None = None,
    parameters: dict[str, Any] | None = None,
    auto_correct: bool = False,
) -> dict[str, Any]:
    """Resolve a proposed settings subset for agent/UI preflight inspection."""
    selected_config = _valid_tune_config(config)
    from dash_app.run_tab.namelist import is_bool_value, is_true, read_namelist_entries

    base_flags = {
        entry["name"]: (
            is_true(entry["value"])
            if is_bool_value(entry["value"])
            else entry["value"]
        )
        for entry in read_namelist_entries(tunable_config_file(selected_config, "configurable_model_flags.in"))
    }
    base_flags.update(dict(flags or {}))
    return {
        "config": selected_config,
        "resolution": resolve_clubb_settings(
            base_flags,
            dict(parameters or {}),
            auto_correct=bool(auto_correct),
        ).as_dict(),
    }


def get_server_info() -> dict[str, Any]:
    """Describe the local-only service without exposing connection secrets."""
    return {
        "service": "CLUBB Dash local broker",
        "protocol_version": 3,
        "scope": "loopback-only",
        "remote_transport_supported": False,
        "artifacts_root": "output/agent_artifacts/",
        "artifacts_lifecycle": "ephemeral staging; active bundles are protected, completed bundles have bounded retention",
        "compatibility_tools": ["inspect_dashboard", "invoke_dashboard"],
        "compatibility_tools_deprecated": True,
        "tunable_parameter_discovery": "list_tunable_parameters(config)",
        "settings_resolution": "resolve_clubb_settings(config, flags, parameters)",
        "source": source_provenance(REPO_ROOT),
    }


def _signal_scm_record(record: dict[str, Any]) -> dict[str, Any]:
    """Signal the exact process recorded by one SCM child."""
    runtime = dict(record.get("runtime") or {})
    proc_data = dict(runtime.get("proc_data") or {})
    case_name = str((record.get("request") or {}).get("case") or "")
    try:
        pid = int(proc_data.get("pid"))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{case_name or 'SCM'} process is not available to stop") from exc
    try:
        os.killpg(pid, signal.SIGTERM)
    except OSError:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError as exc:
            raise ValueError(
                f"{case_name or 'SCM'} process is no longer available to stop"
            ) from exc

    def force_stop() -> None:
        time.sleep(1.0)
        if not _pid_is_alive(pid):
            return
        try:
            os.killpg(pid, signal.SIGKILL)
        except OSError:
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass

    _background(force_stop)
    publish_event(
        "run",
        f"Stop requested for {case_name}",
        str(proc_data.get("log") or ""),
        status="info",
    )
    return {
        "status": "stop_requested",
        "case": case_name,
        "job_id": str(record.get("job_id") or ""),
        "pid": pid,
    }


def _cancel_scm_job(job_id: str) -> dict[str, Any]:
    """Cancel one exact SCM child without relying on case-name uniqueness."""
    with _BATCH_QUEUE_LOCK:
        record = _JOB_STORE.get(str(job_id))
        if record is None or record.get("kind") != "scm":
            raise ValueError("SCM job was not found")
        state = str(record.get("state") or "")
        if state not in _SCM_RESERVED_STATES:
            return {"status": "already_terminal", **record}
        runtime = dict(record.get("runtime") or {})
        proc_data = dict(runtime.get("proc_data") or {})
        if proc_data.get("pid"):
            updated = _JOB_STORE.update(str(job_id), state="stopping") or record
            target = updated
        else:
            updated = _JOB_STORE.update(
                str(job_id),
                state="cancelled",
                returncode=1,
                finished_at_unix_seconds=time.time(),
            ) or record
            target = None

    batch_job_id = str(updated.get("batch_job_id") or "")
    if target is not None:
        operation = _signal_scm_record(target)
        return {"status": "stop_requested", "operation": operation, **updated}
    _ARTIFACT_STORE.release(str(job_id))
    if batch_job_id:
        _refresh_scm_batch(batch_job_id, str(job_id))
    return {"status": "cancelled", **updated}


def cancel_job(job_id: str) -> dict[str, Any]:
    """Request cancellation through the same safe lifecycle used by Dash."""
    record = get_job(job_id)
    state = str(record.get("state") or "")
    if state not in _SCM_RESERVED_STATES:
        return {"status": "already_terminal", **record}
    kind = str(record.get("kind") or "")
    if kind == "compile":
        result = stop_compile()
    elif kind == "scm":
        return _cancel_scm_job(job_id)
    elif kind == "scm_batch":
        return cancel_scm_batch(job_id)
    elif kind == "tune":
        result = stop_tuning()
    elif kind == "leaderboard":
        run_id = str(((record.get("runtime") or {}).get("run") or {}).get("run_id") or "")
        loss_run = dict((broker_jobs().get("loss_runs") or {}).get(run_id) or {})
        if not loss_run:
            raise ValueError("leaderboard run is no longer available to stop")
        update_broker_loss_run(run_id, **stop_loss_run(loss_run))
        result = {"status": "stop_requested", "run_id": run_id}
    else:
        raise ValueError("this job kind cannot be cancelled")
    updated = _JOB_STORE.update(job_id, state="stopping") or record
    return {"status": "stop_requested", "operation": result, **updated}


def cancel_scm_batch(batch_job_id: str) -> dict[str, Any]:
    """Cancel queued and active children of one SCM batch in one broker request."""
    with _BATCH_QUEUE_LOCK:
        def cancel(transaction):
            batch = transaction.get(batch_job_id)
            if batch is None or batch.get("kind") != "scm_batch":
                raise ValueError("SCM batch job was not found")
            if _batch_terminal(str(batch.get("state") or "")):
                return batch, [], True
            process_targets: list[dict[str, Any]] = []
            for child_record in _batch_children(batch, transaction):
                child_id = str(child_record.get("job_id") or "")
                child_state = str(child_record.get("state") or "")
                proc_data = dict(
                    (child_record.get("runtime") or {}).get("proc_data") or {}
                )
                if child_state not in _SCM_RESERVED_STATES:
                    continue
                if proc_data.get("pid"):
                    updated_child = transaction.update(
                        child_id, state="stopping"
                    ) or child_record
                    process_targets.append(updated_child)
                else:
                    transaction.update(
                        child_id,
                        state="cancelled",
                        returncode=1,
                        finished_at_unix_seconds=time.time(),
                    )
            batch = transaction.update(
                batch_job_id,
                state="stopping" if process_targets else "cancelled",
                cancellation_requested=True,
                finished_at_unix_seconds=(
                    None if process_targets else time.time()
                ),
            ) or batch
            return batch, process_targets, False

        batch, process_targets, already_terminal = _JOB_STORE.transaction(cancel)
        batch = _materialize_scm_batch(batch) or batch
        if already_terminal:
            return {"status": "already_terminal", **batch}

    errors = []
    stopped_cases = []
    for target in process_targets:
        case_name = str((target.get("request") or {}).get("case") or "")
        try:
            _signal_scm_record(target)
            stopped_cases.append(case_name)
        except ValueError as exc:
            errors.append(f"{case_name}: {exc}")
    for child in batch.get("children") or []:
        if str(child.get("state") or "") == "cancelled":
            _ARTIFACT_STORE.release(str(child.get("job_id") or ""))
    if not process_targets:
        _write_run_result_summary(batch_job_id, batch.get("output_directory"), 1)
        _ARTIFACT_STORE.release(batch_job_id)
    return {
        "status": "stop_requested" if process_targets else "cancelled",
        "stopped_cases": stopped_cases,
        "errors": errors,
        **batch,
    }


def cancel_all_scm_runs() -> dict[str, Any]:
    """Cancel every queued or running Run-tab SCM job in one broker action."""
    with _BATCH_QUEUE_LOCK:
        # Read process identities only after launch has left its short critical
        # section.  Every running child is then either signalled here or was
        # cancelled before it could spawn.
        targets = []
        for record in _JOB_STORE.list_kind("scm"):
            runtime = dict(record.get("runtime") or {})
            proc_data = dict(runtime.get("proc_data") or {})
            if (
                str(record.get("state") or "") not in {"running", "stopping"}
                or not proc_data.get("pid")
            ):
                continue
            targets.append(
                {
                    "case": str((record.get("request") or {}).get("case") or ""),
                    "job_id": str(record.get("job_id") or ""),
                    "pid": int(proc_data["pid"]),
                }
            )
        target_job_ids = {
            target["job_id"] for target in targets if target["job_id"]
        }
        target_cases = {target["case"] for target in targets}
        now = time.time()

        def cancel(transaction):
            cancelled_children: list[str] = []
            affected_cases: set[str] = set(target_cases)
            for record in list(transaction.state.values()):
                if not isinstance(record, dict) or record.get("kind") != "scm":
                    continue
                if str(record.get("state") or "") not in {
                    "queued",
                    "submitting",
                    "starting",
                    "running",
                    "stopping",
                }:
                    continue
                case = str((record.get("request") or {}).get("case") or "")
                affected_cases.add(case)
                has_process = str(record.get("job_id") or "") in target_job_ids
                updates = {"state": "stopping" if has_process else "cancelled"}
                if not has_process:
                    updates.update(returncode=1, finished_at_unix_seconds=now)
                    cancelled_children.append(str(record.get("job_id") or ""))
                transaction.update(str(record.get("job_id") or ""), **updates)

            batches: dict[str, dict[str, Any]] = {}
            for record in list(transaction.state.values()):
                if (
                    not isinstance(record, dict)
                    or record.get("kind") != "scm_batch"
                    or _batch_terminal(str(record.get("state") or ""))
                ):
                    continue
                children = _batch_children(record, transaction)
                has_process = any(
                    str(child.get("state") or "") == "stopping"
                    for child in children
                )
                updated = transaction.update(
                    str(record.get("job_id") or ""),
                    state="stopping" if has_process else "cancelled",
                    cancellation_requested=True,
                    finished_at_unix_seconds=None if has_process else now,
                )
                if updated:
                    batches[str(updated["job_id"])] = updated
            return batches, affected_cases, cancelled_children

        batches, affected_cases, cancelled_children = _JOB_STORE.transaction(cancel)

    errors = []
    signalled_pids = []
    for target in targets:
        pid = target["pid"]
        try:
            os.killpg(pid, signal.SIGTERM)
        except OSError:
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError as exc:
                errors.append(f"{target['case']}: {exc}")
                continue
        signalled_pids.append(pid)

    def force_stop() -> None:
        time.sleep(1.0)
        for pid in signalled_pids:
            if not _pid_is_alive(pid):
                continue
            try:
                os.killpg(pid, signal.SIGKILL)
            except OSError:
                try:
                    os.kill(pid, signal.SIGKILL)
                except OSError:
                    pass

    if signalled_pids:
        _background(force_stop)

    def finish_cancellation() -> None:
        for batch in batches.values():
            if str(batch.get("state") or "") == "cancelled":
                _write_run_result_summary(
                    batch["job_id"], batch.get("output_directory"), 1
                )
                _ARTIFACT_STORE.release(batch["job_id"])
        for job_id in cancelled_children:
            if job_id:
                _ARTIFACT_STORE.release(job_id)

    if batches or cancelled_children:
        _background(finish_cancellation)
    return {
        "status": (
            "stop_requested" if targets else "cancelled" if affected_cases else "idle"
        ),
        "stopped_cases": sorted(target_cases),
        "errors": errors,
        "summary": scm_run_summary(),
    }


def read_job_log(job_id: str, cursor: int = 0, max_bytes: int = 8192) -> dict[str, Any]:
    record = get_job(job_id)
    runtime = dict(record.get("runtime") or {})
    run = dict(runtime.get("run") or {})
    path = Path(str(runtime.get("log") or runtime.get("log_path") or run.get("log") or run.get("log_path") or ""))
    offset = max(0, int(cursor))
    limit = min(max(1, int(max_bytes)), 16384)
    raw, next_cursor, _eof = read_file_chunk(path, offset, limit)
    return {
        "job_id": job_id,
        "cursor": offset,
        "next_cursor": next_cursor,
        "text": raw.decode("utf-8", errors="replace"),
    }


def create_profile_artifact(request: ProfileArtifactRequest) -> dict[str, Any]:
    """Export native Plot-tab data as immutable PNG artifacts for one run.

    The selection code is shared with Dash's profile figure path.  Unlike the
    legacy ``save_profile_png`` action it resolves output from a durable run
    record rather than whatever later run happens to occupy ``output/<case>``.
    """
    run = _JOB_STORE.get_run(request.run_id)
    if run is None:
        raise ValueError("unknown run_id")
    output_directory = str(run.get("output_directory") or "")
    if not output_directory:
        raise ValueError("run has no isolated output directory")
    case_name = _validate_case(str((run.get("request") or {}).get("case") or ""))
    payload = request.model_dump(mode="json")
    record, created = _persist_submission("profile", request.request_id, payload)
    if not created:
        return {"status": "existing", **record}
    try:
        duration_minutes = request.time_window.duration_seconds / 60.0
        case_data, selection = _profile_selection(
            case_name,
            time_start_seconds=request.time_window.start_seconds,
            average_minutes=duration_minutes,
            require_case_data=True,
            output_dirs=[output_directory],
        )
        assert case_data is not None
        names = _validated_profile_names(request.variables, case_data)
        context = _profile_global_context(case_data, selection)
        artifact_dir = _ARTIFACT_STORE.bundle_dir(record["job_id"])
        files: list[dict[str, str]] = []
        for index, name in enumerate(names, start=1):
            figure = profile_plot.build_figure({"var": name, "size": "normal"}, context)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "profile"
            image = artifact_dir / f"{index:02d}_{safe_name}.png"
            _write_profile_figure_png(figure, image)
            files.append({"name": image.name, "uri": f"clubb-artifact://{record['job_id']}/{image.name}", "sha256": sha256_file(image)})
        actual = {
            "start_seconds": float(selection.get("time_start_seconds", request.time_window.start_seconds)),
            "duration_seconds": float(selection.get("average_minutes", duration_minutes)) * 60.0,
        }
        manifest = _ARTIFACT_STORE.create_manifest(
            record["job_id"],
            {
                "job": record,
                "source_run_id": request.run_id,
                "case": case_name,
                "requested_time_window": request.time_window.model_dump(),
                "actual_time_window": actual,
                "vertical_coordinate": {"name": "zt", "units": "m", "selection": request.vertical_coordinate},
                "variables": names,
                "files": files,
            },
        )
        updated = _JOB_STORE.update(record["job_id"], state="finished", artifact_id=record["job_id"], manifest_path=str(manifest)) or record
        publish_event(
            "plot",
            "Created immutable profile artifact",
            f"artifact {record['job_id']} · source run {request.run_id} · {manifest}",
            status="success",
        )
        return {"status": "finished", **updated, "uri": f"clubb-artifact://{record['job_id']}/manifest.json", "files": files}
    except Exception as exc:
        _JOB_STORE.update(record["job_id"], state="error", error={"code": "PROFILE_ARTIFACT_FAILED", "message": str(exc)[:500]})
        raise


def launch_compile_request(options: dict[str, Any], *, env_id: str = "current", job_id: str | None = None) -> dict[str, Any]:
    """Start one compile through the broker-owned lifecycle used by Dash and MCP."""
    normalized = dict(options or {})
    if (broker_jobs().get("compile") or {}).get("state") in {"running", "stopping"}:
        raise ValueError("a CLUBB compile is already active")
    publish_event("compile", "Starting CLUBB compile", "Resolving the current compiler environment.", status="running")
    job = start_compile_job(discover_compile_state(), str(env_id or "current"), normalized)
    set_broker_job(
        "compile",
        {
            "state": "running",
            "job": dict(job),
            "log_tail": "",
            "log_offset": 0,
            "returncode": None,
            "broker_managed": True,
            "job_id": job_id,
        },
    )
    publish_event(
        "compile",
        "Compile command started",
        str(job.get("command") or ""),
        status="running",
        action={"type": "compile", "tab": "compile", "operation": "start", "job": dict(job)},
    )
    _background(_watch_compile, job, job_id)
    return {"status": "started", "pid": job.get("pid"), "command": job.get("command"), "job": dict(job)}


def launch_profile_request(settings: dict[str, Any]) -> dict[str, Any]:
    """Start one process-based CLUBB timing sweep under broker ownership."""
    existing = dict(broker_jobs().get("profile") or {})
    if str(existing.get("state") or "") in {"running", "stopping"}:
        raise ValueError("a Profile timing sweep is already active")
    job = start_profile_process(dict(settings or {}))
    record = {
        **job,
        "broker_managed": True,
        "log_tail": "",
        "returncode": None,
    }
    set_broker_job("profile", record)
    publish_event(
        "profile",
        "Profile timing started",
        str(job.get("command_display") or ""),
        status="running",
        action={"type": "profile", "tab": "profile"},
    )
    _background(_watch_profile, dict(job))
    return {"status": "started", "job": record}


def _watch_pyplotgen(job: dict[str, Any], process) -> None:
    returncode = process.wait()
    success = returncode == 0 and Path(str(job.get("html_path") or "")).is_file()
    current = dict(broker_jobs().get("pyplotgen") or {})
    requested_stop = str(current.get("state") or "") == "stopping"
    state = "finished" if success else ("stopped" if requested_stop else "error")
    error = None
    if not success and not requested_stop:
        error = f"PyPlotGen exited with status {returncode}; see {job.get('log_path')}"
    update_broker_job(
        "pyplotgen",
        state=state,
        returncode=returncode,
        error=error,
        finished_at=time.time(),
    )
    release_pyplotgen(process)
    publish_event(
        "plot",
        "PyPlotGen gallery finished" if success else "PyPlotGen gallery failed",
        str(job.get("html_path") if success else job.get("log_path") or ""),
        status="success" if success else ("info" if requested_stop else "error"),
    )


def launch_pyplotgen_request(output_dirs: list[str]) -> dict[str, Any]:
    """Start one fixed-destination Plot-tab PyPlotGen export."""
    existing = dict(broker_jobs().get("pyplotgen") or {})
    if str(existing.get("state") or "") in {"queued", "submitting", "running", "stopping"}:
        raise ValueError("a PyPlotGen export is already active")
    job, process = start_pyplotgen(list(output_dirs or []))
    set_broker_job("pyplotgen", job)
    publish_event(
        "plot",
        "PyPlotGen gallery started",
        str(job.get("command_display") or ""),
        status="running",
    )
    _background(_watch_pyplotgen, dict(job), process)
    return {"status": "started", "job": job}


def stop_pyplotgen_request() -> dict[str, Any]:
    """Immediately stop the active PyPlotGen process group."""
    record = dict(broker_jobs().get("pyplotgen") or {})
    state = str(record.get("state") or "")
    if state == "stopping":
        return {"status": "stop_requested", "job": record}
    if state not in {"queued", "submitting", "running"}:
        raise ValueError("no PyPlotGen export is active")
    pid = int(record.get("pid"))
    updated = update_broker_job("pyplotgen", state="stopping") or {**record, "state": "stopping"}
    stop_pyplotgen(pid)

    def force_stop() -> None:
        time.sleep(1.0)
        try:
            os.killpg(pid, signal.SIGKILL)
        except OSError:
            pass

    _background(force_stop)
    publish_event("plot", "Stopping PyPlotGen gallery", str(record.get("output_dir") or ""), status="info")
    return {"status": "stop_requested", "job": updated}


def launch_rebuild_request(
    builds: list[dict[str, Any]],
    discovery: dict[str, Any],
    label: str,
    *,
    job_id: str | None = None,
) -> dict[str, Any]:
    """Start a selected-build rebuild through the durable compile lifecycle.

    Rebuilds are still a normal Dash-only operation, but they must share the
    compile admission guard and broker-owned process handle with typed compile
    requests.  Otherwise a reload could let a rebuild and a compile overlap.
    """
    if (broker_jobs().get("compile") or {}).get("state") in {"running", "stopping"}:
        raise ValueError("a CLUBB compile is already active")
    selected = [dict(build) for build in (builds or []) if isinstance(build, dict) and build.get("path")]
    if not selected:
        raise ValueError("select at least one existing build to rebuild")
    publish_event("compile", "Starting CLUBB rebuild", str(label or "selected builds"), status="running")
    job = start_rebuild_job(selected, dict(discovery or {}), str(label or "selected builds"))
    set_broker_job(
        "compile",
        {
            "state": "running",
            "job": dict(job),
            "log_tail": "",
            "log_offset": 0,
            "returncode": None,
            "broker_managed": True,
            "job_id": job_id,
            "operation": "rebuild",
        },
    )
    publish_event(
        "compile",
        "Rebuild command started",
        str(job.get("command") or ""),
        status="running",
        action={
            "type": "compile",
            "tab": "compile",
            "operation": "start",
            "job": dict(job),
            "preserve_tab": True,
        },
    )
    _background(_watch_compile, job, job_id)
    return {"status": "started", "pid": job.get("pid"), "command": job.get("command"), "job": dict(job)}


def compile_clubb(
    *,
    debug: bool = True,
    python_bindings: bool = False,
    fresh: bool = False,
    gptl: bool = False,
    job_id: str | None = None,
) -> dict[str, Any]:
    """Start one normal CLUBB compile and stream real progress into Dash."""
    options = {"debug": bool(debug), "run_tests": False, "python": bool(python_bindings), "fresh": bool(fresh), "openmp": False, "tuning": False, "gptl": bool(gptl), "precision": "double", "gpu": "none", "toolchain": "auto", "extra_args": "", "module_stack": []}
    return launch_compile_request(options, job_id=job_id)


def run_scm(
    case: str,
    overrides: str = "",
    config: str = "default",
    stats_file: str = DEFAULT_STATS_NAME,
    *,
    cli_options: dict[str, Any] | None = None,
    output_dir: str | None = None,
) -> dict[str, Any]:
    """Route the compatibility SCM surface through the typed submission path."""
    case_name = _validate_case(case)
    stats_name = _validated_stats_file(stats_file)
    config_name = _valid_tune_config(config)
    normalized_options = _normalize_dashboard_cli_options(cli_options)
    override_args = _override_args(overrides)
    assignments = override_args[1].split(",") if override_args else []
    typed_overrides = {
        name: value for name, value in (item.split("=", 1) for item in assignments)
    }
    typed_options = {
        name: normalized_options[name]
        for name in ("max_iters", "dt_main", "dt_rad", "tout")
        if name in normalized_options
    }
    request = ScmRunRequest(
        request_id=f"compat-scm-{time.time_ns()}",
        case=case_name,
        stats_file=stats_name,
        config=config_name,
        overrides=typed_overrides,
        run_options=typed_options,
        out_dir=str(output_dir or normalized_options.get("out_dir") or "") or None,
    )
    return submit_scm_run(request, native_cli_options=normalized_options)


def _launch_scm_process(
    case: str,
    stats_file: str,
    config: str,
    overrides: dict[str, Any] | None,
    cli_options: dict[str, Any] | None,
    *,
    job_id: str | None = None,
    run_id: str | None = None,
    batch_job_id: str | None = None,
) -> dict[str, Any]:
    """Launch the process owned by one already-persisted typed SCM child."""
    case_name = str(case)
    stats_name = str(stats_file)
    config_name = str(config)
    normalized_overrides = dict(overrides or {})
    normalized_options = dict(cli_options or {})
    with exclusive_file_lock(private_path(REPO_ROOT, "scm-submission.lock")):
        _assert_scm_admission(
            case_name,
            normalized_options.get("out_dir") or "",
            exclude_job_id=job_id,
        )
        detail = " ".join(normalized_options.get("extra_args") or []) or "dashboard configuration"
        publish_event("run", f"Starting {case_name}", detail, status="running")
        proc_data = start_case_process(case_name, stats_name, normalized_overrides, normalized_options, config_name)
        if normalized_options.get("out_dir"):
            proc_data["output_directory"] = str(normalized_options["out_dir"])
    return {
        "status": "started",
        "case": case_name,
        "pid": proc_data.get("pid"),
        "log": proc_data.get("log"),
        "proc_data": dict(proc_data),
        "stats_file": stats_name,
        "config": config_name,
        "cli_options": normalized_options,
        "output_directory": normalized_options.get("out_dir"),
    }


def _watch_tuning(job: dict[str, Any], request: dict[str, Any], job_id: str | None = None) -> None:
    """Keep the file-backed tuner alive independently of Dash and mirror its status."""
    last_heartbeat = 0.0
    while True:
        now = time.monotonic()
        if now - last_heartbeat >= 5.0:
            TunerJob.from_dict(job).heartbeat()
            last_heartbeat = now
        status = read_tuning_status(job.get("status_path"))
        state = str(status.get("state") or "")
        requested_stop = str((broker_jobs().get("tune") or {}).get("state") or "") == "stopping"
        update_broker_job(
            "tune",
            status=status,
            state="stopping" if requested_stop else "running",
            log_tail=_tail(_read_log_tail(job.get("log_path")), 16000),
        )
        if job_id:
            _JOB_STORE.update(job_id, progress={"samples_evaluated": status.get("samples_evaluated", 0), "best_total_loss": status.get("best_total_loss")})
        # A manual Tune job recovered after a Dash reload has a durable job
        # directory but no in-process Popen/pid. Its status file remains the
        # authority; treating a missing pid as "exited" would immediately
        # discard a worker that is still running.
        exited = bool(job.get("pid")) and active_job_exited(job)
        if state in {"finished", "stopped", "error"} or exited:
            if state == "finished":
                message, event_status = "Tuning finished", "success"
            elif state == "stopped":
                message, event_status = "Tuning stopped", "info"
            else:
                message, event_status = "Tuning failed", "error"
            detail = (
                f"samples: {status.get('samples_evaluated', 0)} | "
                f"best loss: {status.get('best_total_loss', '--')} | "
                f"log: {job.get('log_path') or '--'}"
            )
            update_broker_job("tune", status=status, state=state if state in {"finished", "stopped", "error"} else "error", finished_at=time.time())
            publish_event("tune", message, detail, status=event_status)
            if job_id:
                _JOB_STORE.update(job_id, state=state if state in {"finished", "stopped", "error"} else "error", finished_at_unix_seconds=time.time())
                _ARTIFACT_STORE.release(job_id)
            return
        time.sleep(1.0)


def recover_active_tuning_from_disk() -> dict[str, Any] | None:
    """Adopt the newest active file-backed Tune worker after a Dash-only reload.

    A worker created by an older/manual Tune-tab callback may outlive the Dash
    process that started it. The job directory is enough to renew its lease and
    restore its UI; process identity is deliberately optional here because the
    status/control files are the tuner's durable control contract.
    """
    existing = broker_jobs().get("tune") or {}
    if str(existing.get("state") or "") in {"running", "stopping"}:
        existing_job = dict(existing.get("job") or {})
        # A recovered worker intentionally has no PID, so its status file is
        # its liveness authority.  A broker-owned job has a PID: retain it only
        # while that PID is alive.  This admission-time check covers a watcher
        # that was interrupted while cancellation was in progress.
        existing_status = read_tuning_status(existing_job.get("status_path"))
        persisted_state = str(existing_status.get("state") or "").lower()
        pid = existing_job.get("pid")
        worker_exited = bool(pid) and active_job_exited(existing_job)
        recovered_without_live_status = not pid and persisted_state not in {"running", "stopping"}
        if worker_exited or recovered_without_live_status:
            final_state = persisted_state if persisted_state in {"finished", "stopped", "error"} else "error"
            if worker_exited:
                detail = (
                    "Discarded stale Tune record: recorded worker PID "
                    f"{pid} has exited; status is {persisted_state or 'missing'}"
                )
            else:
                detail = (
                    "Discarded stale recovered Tune record: no worker PID and "
                    f"status is {persisted_state or 'missing'}"
                )
            update_broker_job(
                "tune",
                state=final_state,
                status=existing_status,
                recovery_error=detail,
                finished_at=time.time(),
            )
            publish_event("tune", "Stale Tune recovery cleared", detail, status="error")
        else:
            return dict(existing)

    output_roots = [REPO_ROOT / "output" / "tuner", REPO_ROOT / "output_tuner"]
    candidates: list[tuple[float, Path, dict[str, Any]]] = []
    status_paths = []
    for output_root in output_roots:
        try:
            status_paths.extend(output_root.glob("*/status.json"))
            status_paths.extend(output_root.glob("*/*/status.json"))
        except OSError:
            continue
    for status_path in status_paths:
        try:
            status = read_tuning_status(str(status_path))
            if str(status.get("state") or "") not in {"running", "stopping"}:
                continue
            candidates.append((status_path.stat().st_mtime, status_path.parent, status))
        except OSError:
            continue
    if not candidates:
        return None

    _, job_dir, status = max(candidates, key=lambda item: item[0])
    try:
        request = json.loads((job_dir / "request.json").read_text(encoding="utf-8"))
        if not isinstance(request, dict):
            raise ValueError("request is not an object")
        job = TunerJob.from_dir(job_dir).to_dict()
    except (OSError, RuntimeError, ValueError, TypeError) as exc:
        publish_event("tune", "Could not recover Tune worker", str(exc), status="error")
        return None

    record = {
        "state": str(status.get("state") or "running"),
        "job": job,
        "request": request,
        "status": status,
        "log_tail": _read_log_tail(job.get("log_path")),
        "broker_managed": True,
        "recovered": True,
    }
    set_broker_job("tune", record)
    publish_event("tune", "Recovered running Tune job", str(job_dir), status="info")
    _background(_watch_tuning, job, request)
    return record


def launch_tuning(
    cases: list[Any],
    parameter_ranges: list[dict[str, Any]],
    fields: list[str] | None = None,
    *,
    config: str = "default",
    strategy: str = "random",
    max_samples: int = 12,
    resolve_spacing: float = 0.1,
    simann_max_iters: int = 200,
    simann_initial_temp: float = 1.0,
    simann_final_temp: float = 1.0e-12,
    adam_max_updates: int = 100,
    adam_learning_rate: float = 0.01,
    adam_perturbation: float = 0.05,
    adam_spsa_pairs: int = 2,
    batch_size: int = 1,
    max_workers: int = 1,
    loss_mode: str = DEFAULT_LOSS_MODE,
    aggregation_mode: str = DEFAULT_AGGREGATION_MODE,
    aggregation_weights: list[float] | tuple[float, ...] = DEFAULT_AGGREGATION_WEIGHTS,
    time_window_aggregation_scope: str = DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    overrides: str = "",
    job_id: str | None = None,
    workspace_id: str | None = None,
    revision_id: str | None = None,
    workspace_display_name: str | None = None,
    _max_samples_limit: int | None = _TUNE_MAX_RANDOM_SAMPLES,
) -> dict[str, Any]:
    """Start a bounded, fully visible tuning job from a structured request.

    The agent API deliberately accepts the same physics choices the Tune tab
    displays, but not arbitrary executable text.  Dash receives the canonical
    request and worker metadata immediately, so the user can inspect the exact
    cases, windows, fields, ranges, flags, and stop control in the normal UI.
    """
    # The same guard applies to typed MCP, legacy bridge actions, and native
    # Dash handoff.  Recovering first matters after a Dash/broker reload: the
    # file-backed worker may still be running even when no in-memory Popen is.
    recover_active_tuning_from_disk()
    active_tune = broker_jobs().get("tune") or {}
    if str(active_tune.get("state") or "") in {"running", "stopping"}:
        raise ValueError("a Tune job is already active")
    selected_config = _valid_tune_config(config)
    case_data = load_case_defaults()
    case_names, case_configs = _normalize_tune_cases(cases, case_data)
    common_fields = set(available_fields_for_cases(case_names, case_data))
    requested_fields = [str(field).strip() for field in (fields or []) if str(field).strip()]
    if not requested_fields:
        requested_fields = [
            field
            for field in (case_data.get(case_names[0], {}) or {}).get("default_clubb_fields", [])
            if field in common_fields
        ]
    if not requested_fields:
        raise ValueError("no common tunable benchmark fields are available for the selected cases")
    invalid_fields = sorted(set(requested_fields) - common_fields)
    if invalid_fields:
        raise ValueError("field(s) unavailable for every selected case: " + ", ".join(invalid_fields))

    ranges = _normalize_tune_ranges(parameter_ranges, selected_config)
    normalized_override = _normalize_tune_override(overrides)
    resolution = evaluate_tune_settings(selected_config, normalized_override)
    errors = [str(issue.get("message") or "") for issue in resolution.get("issues", []) if issue.get("severity") == "error"]
    if errors:
        raise ValueError("; ".join(errors))
    inactive = sorted(
        target
        for item in ranges
        for target in item["targets"]
        if not is_independently_tunable(resolution.get("parameter_states", {}).get(target))
    )
    if inactive:
        details = "; ".join(
            f"{name}: {resolution['parameter_states'][name]['reason']}" for name in inactive
        )
        raise ValueError("Tune request selects inactive parameter(s): " + details)
    ranges = apply_required_parameter_links(ranges, resolution)
    # The UI default intentionally uses roughly physical-core parallelism,
    # but a user-selected maximum must not be silently capped at that
    # conservative default.  The process affinity is the real machine limit
    # exposed to this local dashboard/broker.
    worker_cap = max(1, available_logical_cpu_count())
    requested_workers = _positive_int(max_workers, "max_workers", maximum=worker_cap)
    requested_batch_size = _positive_int(batch_size, "batch_size", maximum=_TUNE_MAX_BATCH_SIZE)
    selected_loss_mode = str(loss_mode or DEFAULT_LOSS_MODE).strip()
    selected_aggregation_mode = str(aggregation_mode or DEFAULT_AGGREGATION_MODE).strip()
    if selected_loss_mode not in LOSS_MODE_NAMES:
        raise ValueError("unknown loss mode")
    if selected_aggregation_mode not in AGGREGATION_MODE_NAMES:
        raise ValueError("unknown aggregation mode")
    try:
        selected_aggregation_weights = normalize_aggregation_weights(aggregation_weights)
    except ValueError as exc:
        raise ValueError(str(exc)) from exc
    selected_aggregation_scope = str(time_window_aggregation_scope or DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE).strip()
    if selected_aggregation_scope not in TIME_WINDOW_AGGREGATION_SCOPES:
        raise ValueError("unknown time-window aggregation scope")
    request = {
        "config": selected_config,
        "override": normalized_override,
        "cases": case_names,
        "case_configs": case_configs,
        "selected_fields": requested_fields,
        "batch_size": requested_batch_size,
        "max_workers": requested_workers,
        "loss_mode": selected_loss_mode,
        "aggregation_mode": selected_aggregation_mode,
        "time_window_aggregation_mode": selected_aggregation_mode,
        "aggregation_weights": selected_aggregation_weights,
        "time_window_aggregation_scope": selected_aggregation_scope,
        "parameter_ranges": ranges,
        "settings_resolution": resolution,
        "strategy": _normalize_tune_strategy(
            strategy,
            ranges,
            max_samples=max_samples,
            resolve_spacing=resolve_spacing,
            simann_max_iters=simann_max_iters,
            simann_initial_temp=simann_initial_temp,
            simann_final_temp=simann_final_temp,
            adam_max_updates=adam_max_updates,
            adam_learning_rate=adam_learning_rate,
            adam_perturbation=adam_perturbation,
            adam_spsa_pairs=adam_spsa_pairs,
            batch_size=requested_batch_size,
            max_workers=requested_workers,
            max_samples_limit=_max_samples_limit,
        ),
    }
    publish_event(
        "tune",
        "Starting tuning job",
        f"cases: {', '.join(case_names)} | ranges: {', '.join(str(item['name']) for item in ranges)}",
        status="running",
    )
    if bool(workspace_id) != bool(revision_id):
        raise ValueError("workspace_id and revision_id must be supplied together")
    if workspace_id:
        job = start_draft_tuning_job(workspace_id, revision_id, request)
    else:
        job = start_tuning_job(request, display_name=workspace_display_name)
    set_broker_job(
        "tune",
        {
            "state": "running",
            "job": dict(job),
            "request": dict(request),
            "status": read_tuning_status(job.get("status_path")),
            "log_tail": "",
            "broker_managed": True,
            "job_id": job_id,
        },
    )
    event = publish_tune_request(
        request,
        job,
        message="Showing agent tuning job in Tune",
        detail=str(job.get("log_path") or ""),
    )
    _background(_watch_tuning, job, request, job_id)
    return {
        "status": "started",
        "activity_id": event["id"],
        "pid": job.get("pid"),
        "job_dir": job.get("job_dir"),
        "log": job.get("log_path"),
        "request": request,
        "job": dict(job),
    }


def launch_tuning_request(request: dict[str, Any]) -> dict[str, Any]:
    """Broker-launch a validated native Tune-tab request without shrinking its UI budget."""
    payload = dict(request or {})
    strategy = dict(payload.get("strategy") or {})
    options = dict(strategy.get("options") or {})
    return launch_tuning(
        list(payload.get("case_configs") or []),
        list(payload.get("parameter_ranges") or []),
        list(payload.get("selected_fields") or []),
        config=payload.get("config", "default"),
        strategy=strategy.get("name", "random"),
        max_samples=options.get("max_samples", 12),
        resolve_spacing=options.get("spacing", 0.1),
        simann_max_iters=options.get("max_iters", 200),
        simann_initial_temp=options.get("initial_temp", 1.0),
        simann_final_temp=options.get("max_final_temp", 1.0e-12),
        adam_max_updates=options.get("max_updates", 100),
        adam_learning_rate=options.get("learning_rate", 0.01),
        adam_perturbation=options.get("perturbation", 0.05),
        adam_spsa_pairs=options.get("spsa_pairs", 2),
        batch_size=payload.get("batch_size", 1),
        max_workers=payload.get("max_workers", 1),
        loss_mode=payload.get("loss_mode", DEFAULT_LOSS_MODE),
        aggregation_mode=payload.get("aggregation_mode", DEFAULT_AGGREGATION_MODE),
        aggregation_weights=payload.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS),
        time_window_aggregation_scope=payload.get("time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE),
        overrides=payload.get("override", ""),
        workspace_id=payload.get("workspace_id"),
        revision_id=payload.get("revision_id"),
        workspace_display_name=payload.get("workspace_display_name"),
        _max_samples_limit=None,
    )


def list_tuning_workspaces() -> dict[str, Any]:
    """Return all durable Tune workspaces and their execution revisions."""
    return {"workspaces": list_tune_workspaces()}


def list_tuning_workspace_activity() -> dict[str, Any]:
    """Return cheap revision-state data for the Tune activity indicator."""
    return {"activity": list_tune_workspace_activity()}


def load_tuning_workspace(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Read one saved Tune execution for a readonly UI/MCP handoff."""
    return load_tune_workspace_execution(workspace_id, revision_id)


def create_tuning_revision(workspace_id: str, revision_id: str, *, restart: bool = False) -> dict[str, Any]:
    """Create an unstarted revision or an exact restart attempt."""
    job = create_tune_workspace_revision(workspace_id, revision_id, restart=restart)
    return {
        "workspace_id": workspace_id,
        "revision_id": job.job_dir.name,
        "job": job.to_dict(),
        "kind": "restart" if restart else "revision",
    }


def start_tuning_draft_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Start an unchanged saved draft through the ordinary validated Tune path."""
    loaded = load_tune_workspace_execution(workspace_id, revision_id)
    if str((loaded.get("execution") or {}).get("state") or "draft") not in {"draft", "idle"}:
        raise ValueError("only an unstarted draft revision can be started")
    request = dict(loaded.get("request") or {})
    if not request:
        raise ValueError("saved draft has no readable request")
    request.update({"workspace_id": workspace_id, "revision_id": revision_id})
    return launch_tuning_request(request)


def restart_tuning_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Create and immediately start a fresh exact restart attempt."""
    created = create_tuning_revision(workspace_id, revision_id, restart=True)
    started = start_tuning_draft_revision(workspace_id, str(created["revision_id"]))
    return {**created, "status": started.get("status"), "job": started.get("job"), "request": started.get("request")}


def reset_tuning_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Destructively erase one inactive revision's generated data into a draft."""
    job = reset_tune_workspace_execution(workspace_id, revision_id)
    return {"workspace_id": workspace_id, "revision_id": revision_id, "job": job.to_dict(), "status": "draft"}


def resume_tuning_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Continue a stopped revision from its durable scheduler checkpoint."""
    loaded = load_tune_workspace_execution(workspace_id, revision_id)
    job_data = dict(loaded["job"])
    job_data.update({"workspace_id": workspace_id, "revision_id": revision_id})
    job = resume_tuning_job(job_data)
    request = dict(loaded["request"])
    set_broker_job(
        "tune",
        {
            "state": "running",
            "job": dict(job),
            "request": request,
            "status": read_tuning_status(job.get("status_path")),
            "log_tail": "",
            "broker_managed": True,
        },
    )
    _background(_watch_tuning, job, request)
    return {"status": "resumed", "job": job, "request": request}


def rename_tuning_workspace(workspace_id: str, display_name: str) -> dict[str, Any]:
    """Rename a saved Tune workspace without moving its immutable files."""
    return rename_tune_workspace(workspace_id, display_name)


def delete_tuning_workspace(workspace_id: str) -> dict[str, Any]:
    """Delete an inactive Tune workspace after an explicit UI/MCP request."""
    delete_tune_workspace(workspace_id)
    return {"deleted": workspace_id}


def stop_tuning() -> dict[str, Any]:
    """Request a graceful stop for the active Dash-owned tuning job."""
    job = stop_active_tuning_job()
    if not job:
        raise ValueError("no active tuning job is owned by this dashboard")
    try:
        request = json.loads(Path(str(job.get("request_path") or "")).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        request = {}
    event = publish_tune_request(
        request,
        job,
        message="Stop requested for tuning job",
        detail=str(job.get("log_path") or ""),
        status="info",
        operation="stop",
    )
    update_broker_job("tune", state="stopping", request=request)
    return {"status": "stop_requested", "activity_id": event["id"], "pid": job.get("pid"), "job_dir": job.get("job_dir")}


def stop_compile() -> dict[str, Any]:
    """Ask the broker-owned compile subprocess to stop, including after Dash reload."""
    record = broker_jobs().get("compile") or {}
    job = dict(record.get("job") or {})
    if str(record.get("state") or "") not in {"running", "stopping"} or not job:
        raise ValueError("no broker-owned compile is running")
    cancelled = cancel_compile_job(job)
    if not cancelled:
        try:
            os.killpg(os.getpgid(int(job["pid"])), signal.SIGTERM)
            cancelled = True
        except (KeyError, OSError, TypeError, ValueError):
            cancelled = False
    if not cancelled:
        raise ValueError("compile process is no longer available to stop")
    update_broker_job("compile", state="stopping")
    publish_event("compile", "Compile stop requested", str(job.get("command") or ""), status="info")
    return {"status": "stop_requested", "pid": job.get("pid")}


def stop_profile() -> dict[str, Any]:
    """Gracefully stop the broker-owned timing wrapper and its child groups."""
    job = dict(broker_jobs().get("profile") or {})
    if str(job.get("state") or "") not in {"running", "stopping"}:
        raise ValueError("no broker-owned Profile timing sweep is running")
    stop_profile_process(job)
    update_broker_job("profile", state="stopping")
    publish_event(
        "profile",
        "Profile timing stop requested",
        str(job.get("log") or ""),
        status="info",
    )
    return {"status": "stop_requested", "pid": job.get("pid")}


def stop_run(case: str, output_dir: str | None = None) -> dict[str, Any]:
    """Stop one unambiguous broker-owned SCM case process."""
    case_name = _validate_case(case)
    target = (
        _canonical_output_directory(resolve_output_dir(output_dir))
        if output_dir is not None
        else None
    )
    candidates = [
        record
        for record in _matching_scm_records(
            case_name,
            target,
            states={"running", "stopping"},
        )
        if ((record.get("runtime") or {}).get("proc_data") or {}).get("pid")
    ]
    if not candidates:
        raise ValueError(f"no broker-owned {case_name} run is running")
    if len(candidates) > 1:
        raise ValueError(
            f"multiple {case_name} runs are active; stop one by job_id or output_dir"
        )
    result = _cancel_scm_job(str(candidates[0].get("job_id") or ""))
    return dict(result.get("operation") or result)


def stop_all_broker_work(*, reason: str = "dashboard manager is stopping") -> dict[str, Any]:
    """Best-effort, idempotent shutdown of every broker-owned worker."""
    _BROKER_SHUTTING_DOWN.set()
    snapshot = broker_jobs()
    requested: list[str] = []
    errors: list[str] = []

    compile_record = dict(snapshot.get("compile") or {})
    if str(compile_record.get("state") or "") in {"running", "stopping"}:
        try:
            stop_compile()
            requested.append("compile")
            job_id = str(compile_record.get("job_id") or "")
            if job_id:
                _JOB_STORE.update(job_id, state="stopping")
        except ValueError as exc:
            errors.append(f"compile: {exc}")

    profile_record = dict(snapshot.get("profile") or {})
    if str(profile_record.get("state") or "") in {"running", "stopping"}:
        try:
            stop_profile()
            requested.append("Profile")
        except ValueError as exc:
            errors.append(f"Profile: {exc}")

    pyplotgen_record = dict(snapshot.get("pyplotgen") or {})
    if str(pyplotgen_record.get("state") or "") in {"running", "stopping"}:
        try:
            stop_pyplotgen(pyplotgen_record.get("pid"))
            update_broker_job("pyplotgen", state="stopping")
            requested.append("PyPlotGen")
        except (OSError, TypeError, ValueError) as exc:
            errors.append(f"PyPlotGen: {exc}")

    if _JOB_STORE.list_kind("scm") or _JOB_STORE.list_kind("scm_batch"):
        try:
            scm_result = cancel_all_scm_runs()
            requested.extend(
                f"SCM:{case}" for case in scm_result.get("stopped_cases") or []
            )
            errors.extend(f"SCM:{message}" for message in scm_result.get("errors") or [])
        except ValueError as exc:
            errors.append(f"SCM: {exc}")

    tune_record = dict(snapshot.get("tune") or {})
    if str(tune_record.get("state") or "") in {"running", "stopping"}:
        job = dict(tune_record.get("job") or {})
        try:
            if not job:
                raise ValueError("durable Tune process metadata is missing")
            stop_tuning_job(job)
            update_broker_job("tune", state="stopping")
            requested.append("Tune")
            job_id = str(tune_record.get("job_id") or "")
            if job_id:
                _JOB_STORE.update(job_id, state="stopping")
        except (KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
            errors.append(f"Tune: {exc}")

    for run_id, raw_record in (snapshot.get("loss_runs") or {}).items():
        record = dict(raw_record or {})
        if str(record.get("state") or "") not in {"running", "stopping"}:
            continue
        try:
            update_broker_loss_run(str(run_id), **stop_loss_run(record))
            requested.append(f"Tune-result:{run_id}")
            job_id = str(record.get("job_id") or "")
            if job_id:
                _JOB_STORE.update(job_id, state="stopping")
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            errors.append(f"Tune-result:{run_id}: {exc}")

    publish_event(
        "broker",
        "Stopping broker-owned work",
        f"{reason} | requested: {', '.join(requested) if requested else 'none'}"
        + (f" | errors: {'; '.join(errors)}" if errors else ""),
        status="error" if errors else "info",
    )
    return {"status": "stop_requested", "requested": requested, "errors": errors}


def inspect_tuning() -> dict[str, Any]:
    """Return the current Tune worker's exact request and file-backed status."""
    job = active_tuning_job_data()
    if not job:
        record = broker_jobs().get("tune") or {}
        job = dict(record.get("job") or {})
        request = dict(record.get("request") or {})
        if job:
            return {"status": read_tuning_status(job.get("status_path")), "job": job, "request": request}
    if not job:
        return {"status": "idle", "job": None, "request": None}
    try:
        request = json.loads(Path(str(job.get("request_path") or "")).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        request = {}
    return {"status": read_tuning_status(job.get("status_path")), "job": job, "request": request}


def _watch_tuning_loss_run(run_data: dict[str, Any], job_id: str | None = None) -> None:
    """Mirror one broker-owned Tune result run into the durable activity stream."""
    run_id = str(run_data.get("run_id") or "")
    last_log = ""
    while run_id:
        updated, _any_running = poll_loss_runs({run_id: run_data})
        current = dict(updated.get(run_id) or run_data)
        log_tail = _tail(_read_log_tail(current.get("log_path")), 4000)
        if log_tail and log_tail != last_log:
            publish_event("tune", "Tune result run output", log_tail, status="running")
            last_log = log_tail
        update_broker_loss_run(run_id, **current, log_tail=log_tail)
        if job_id:
            _JOB_STORE.update(job_id, progress={"run_id": run_id, "log_tail": log_tail[-4000:]})
        if current.get("state") != "running":
            success = current.get("state") == "success"
            stopped = current.get("state") == "stopped"
            publish_event(
                "tune",
                "Tune result run finished" if success else ("Tune result run stopped" if stopped else "Tune result run failed"),
                f"pid: {current.get('pid')} | log: {current.get('log_path')}",
                status="success" if success else ("info" if stopped else "error"),
            )
            if job_id:
                _JOB_STORE.update(job_id, state="finished" if success else ("cancelled" if stopped else "error"), result=current)
                _ARTIFACT_STORE.release(job_id)
            return
        run_data = current
        time.sleep(1.0)


def run_tuning_loss(mode: str = "window", max_results: int = 16, *, job_id: str | None = None) -> dict[str, Any]:
    """Run the current Tune leaderboard through its loss window or full SCM path."""
    selected_mode = str(mode or "window").strip().lower()
    if selected_mode not in {"window", "complete"}:
        raise ValueError("mode must be 'window' or 'complete'")
    result_limit = _positive_int(max_results, "max_results", maximum=16)
    record = broker_jobs().get("tune") or {}
    job = dict(record.get("job") or active_tuning_job_data() or {})
    request = dict(record.get("request") or {})
    if not job:
        raise ValueError("no Tune job is available; run a tuning job before requesting a result run")
    if not request:
        try:
            request = json.loads(Path(str(job.get("request_path") or "")).read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            request = {}
    status = read_tuning_status(job.get("status_path"))
    top_results = list(status.get("top_results") or [])[:result_limit]
    param_sets = [dict(result.get("params") or {}) for result in top_results if isinstance(result, dict)]
    param_sets = [params for params in param_sets if params]
    if not param_sets:
        raise ValueError("the Tune leaderboard has no parameter rows to run yet")

    case_configs = [dict(item or {}) for item in (request.get("case_configs") or [])]
    case_names = [str(item.get("case_name") or "").strip() for item in case_configs]
    if not case_names:
        case_names = [str(item).strip() for item in (request.get("cases") or []) if str(item).strip()]
    if not case_names:
        raise ValueError("the Tune request does not contain any cases")
    fields = [str(item).strip() for item in (request.get("selected_fields") or []) if str(item).strip()]
    if selected_mode == "window" and (not case_configs or not fields):
        raise ValueError("the Tune request lacks its editable case windows or selected fields")

    run_data = start_loss_run(
        case_names,
        fields if selected_mode == "window" else [],
        param_sets,
        rank=f"agent_{selected_mode}",
        case_configs=case_configs if selected_mode == "window" else None,
        run_mode=selected_mode,
        config=str(request.get("config") or "default"),
        override=str(request.get("override") or ""),
        workspace_id=job.get("workspace_id"),
        revision_id=job.get("revision_id"),
    )
    run_id = str(run_data["run_id"])
    set_broker_loss_run(run_id, {**run_data, "mode": selected_mode, "result_count": len(param_sets), "broker_managed": True, "job_id": job_id})
    event = publish_event(
        "tune",
        f"Started {selected_mode} result run",
        f"{len(param_sets)} leaderboard row(s) | pid: {run_data.get('pid')} | log: {run_data.get('log_path')}",
        status="running",
        action={"type": "tab", "tab": "tune"},
    )
    _background(_watch_tuning_loss_run, run_data, job_id)
    return {
        "status": "started",
        "activity_id": event["id"],
        "mode": selected_mode,
        "result_count": len(param_sets),
        "column_count": run_data.get("column_count"),
        "run": run_data,
    }


def plot_profiles(
    case: str,
    variables: list[str],
    time_seconds: float | None = None,
    *,
    run_id: str | None = None,
    output_dir: str | None = None,
    output_dirs: list[str] | None = None,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
    benchmark_sources: list[str] | None = None,
) -> dict[str, Any]:
    """Open profile variables with optional validated benchmark overlays."""
    case_name = _validate_case(case)
    selected_run_id, selected_output_dirs = _plot_output_selection(
        case_name,
        run_id,
        output_dir,
        output_dirs,
    )
    case_data, selection = _profile_selection(
        case_name,
        time_seconds=time_seconds,
        time_start_seconds=time_start_seconds,
        average_minutes=average_minutes,
        window_preset=window_preset,
        require_case_data=False,
        output_dirs=selected_output_dirs,
    )
    request = {
        "operation": "set_view",
        "case": case_name,
        "variables": variables,
        "_normalized_common": selection,
    }
    if benchmark_sources is not None:
        request["benchmark_sources"] = benchmark_sources
    validated = validate_plot_request(case_data or {}, request)
    names = validated["variables"]
    selection = {
        key: value
        for key, value in validated.items()
        if key in {"time_start_seconds", "average_minutes", "window_preset", "benchmark_sources"}
    }
    event = publish_plot_request(
        case_name,
        names,
        output_dir=selected_output_dirs[0] if selected_output_dirs else "output",
        output_dirs=(
            selected_output_dirs
            if selected_output_dirs is not None and len(selected_output_dirs) > 1
            else None
        ),
        **selection,
    )
    return {
        "status": "requested",
        "activity_id": event["id"],
        "case": case_name,
        "variables": names,
        "run_id": selected_run_id or None,
        "output_directory": selected_output_dirs[0] if selected_output_dirs else "output",
        "output_directories": selected_output_dirs or ["output"],
        **selection,
    }


def add_budget_plot(
    case: str,
    budget_group: str = "wp2",
    *,
    run_id: str | None = None,
    output_dir: str | None = None,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
) -> dict[str, Any]:
    """Add a validated budget plot request to the active Plot tab."""
    case_name = _validate_case(case)
    selected_run_id, output_dirs = _plot_output_selection(case_name, run_id, output_dir)

    case_data = _profile_case_data(case_name, required=True, output_dirs=output_dirs)
    validated = validate_plot_request(
        case_data or {},
        {
            "operation": "add_budget",
            "case": case_name,
            "budget_group": budget_group,
            "time_start_seconds": time_start_seconds,
            "average_minutes": average_minutes,
            "window_preset": window_preset,
        },
    )
    event = publish_budget_request(
        case_name,
        validated["budget_group"],
        output_dir=output_dirs[0] if output_dirs else "output",
        time_start_seconds=validated.get("time_start_seconds"),
        average_minutes=validated.get("average_minutes"),
        window_preset=validated.get("window_preset"),
    )
    return {
        "status": "requested",
        "activity_id": event["id"],
        "case": case_name,
        "plot_type": "budget",
        "budget_group": validated["budget_group"],
        "run_id": selected_run_id or None,
        "output_directory": output_dirs[0] if output_dirs else "output",
        **{
            key: value
            for key, value in validated.items()
            if key in {"time_start_seconds", "average_minutes", "window_preset"}
        },
    }


def list_plot_instances() -> dict[str, Any]:
    """List the currently mounted Plot cards from the live dashboard snapshot."""
    snapshot = read_activity().get("plot_instances") or {}
    if not isinstance(snapshot, dict) or snapshot.get("plots") is None:
        raise ValueError("the active dashboard has not published its Plot state yet")
    return {
        "status": "available",
        "case": str(snapshot.get("case") or ""),
        "output_dirs": [str(path) for path in snapshot.get("output_dirs") or []],
        "plots": [dict(item) for item in snapshot.get("plots") or [] if isinstance(item, dict)],
        "next_id": int(snapshot.get("next_id") or 0),
    }


def remove_plot_instance(plot_id: int) -> dict[str, Any]:
    """Publish a validated request to remove one live Plot card by stable ID."""
    snapshot = list_plot_instances()
    normalized_id = plot_service.normalize_plot_id(plot_id)
    owned_ids = {
        plot_service.normalize_plot_id(item.get("id"))
        for item in snapshot["plots"]
        if isinstance(item, dict) and item.get("id") is not None
    }
    if normalized_id not in owned_ids:
        raise ValueError(f"plot_id {normalized_id} is not an active Plot card owned by this dashboard")
    event = publish_plot_remove_request(normalized_id, case=snapshot.get("case", ""))
    return {
        "status": "requested",
        "activity_id": event["id"],
        "tab": "plots",
        "operation": "remove",
        "plot_id": normalized_id,
        "case": snapshot.get("case") or None,
    }


def save_profile_png(
    case: str,
    variables: list[str],
    *,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
) -> dict[str, Any]:
    """Render selected native profile views to durable PNG files under output/."""
    case_name = _validate_case(case)
    case_data, selection = _profile_selection(
        case_name,
        time_start_seconds=time_start_seconds,
        average_minutes=average_minutes,
        window_preset=window_preset,
        require_case_data=True,
    )
    names = _validated_profile_names(variables, case_data)
    context = _profile_global_context(case_data, selection)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    saved_paths: list[str] = []
    for index, var_name in enumerate(names, start=1):
        figure = profile_plot.build_figure({"var": var_name, "size": "normal"}, context)
        safe_var = re.sub(r"[^A-Za-z0-9_.-]+", "_", var_name).strip("._") or "profile"
        destination = _PROFILE_EXPORT_DIR / f"{case_name}_{safe_var}_{stamp}_{index}.png"
        _write_profile_figure_png(figure, destination)
        saved_paths.append(str(destination))
    event = publish_event(
        "plot",
        f"Saved {len(saved_paths)} profile PNG{'s' if len(saved_paths) != 1 else ''}",
        "\n".join(saved_paths),
        status="success",
        action={"type": "profile", "tab": "plots", "operation": "set_view", "case": case_name, "variables": names, **selection},
    )
    return {
        "status": "saved",
        "activity_id": event["id"],
        "case": case_name,
        "variables": names,
        "paths": saved_paths,
        **selection,
    }


def open_dashboard(tab: str) -> dict[str, Any]:
    """Open a safe top-level dashboard tab without issuing a model command."""
    aliases = {
        "notes": "misc",
    }
    value = aliases.get(str(tab or "").strip(), str(tab or "").strip())
    allowed = {"tutorial", "compile", "run", "profile", "tune", "plots", "reports", "misc"}
    if value not in allowed:
        raise ValueError(f"tab must be one of: {', '.join(sorted(allowed))}")
    event = publish_tab_request(value, f"Opening {value.replace('_', ' ').title()}")
    return {"status": "requested", "activity_id": event["id"], "tab": value}


def open_note(slug: str) -> dict[str, Any]:
    """Open one registered investigation report in Misc by stable slug."""
    from dash_app.misc_tab.registry import discover_subtabs
    value = str(slug or "").strip()
    valid = {subtab.slug for subtab in discover_subtabs()}
    if value not in valid:
        raise ValueError(f"unknown notes report: {value}")
    event = publish_tab_request(
        "misc",
        f"Opening Misc report: {value}",
        report_slug=value,
        type="misc",
        operation="open_report",
    )
    return {"status": "requested", "activity_id": event["id"], "slug": value}


def open_static_report(report_id: str) -> dict[str, Any]:
    """Open one immutable report bundle in the dashboard's Reports tab."""
    from dash_app.reports_tab.catalog import report_by_id

    report = report_by_id(str(report_id or "").strip())
    if report is None:
        raise ValueError(f"unknown static report: {report_id}")
    event = publish_tab_request(
        "reports",
        f"Opening static report: {report.title}",
        type="reports",
        operation="open_report",
        report_id=report.report_id,
    )
    return {"status": "requested", "activity_id": event["id"], "report_id": report.report_id}


def open_tutorial_lesson(lesson: str) -> dict[str, Any]:
    """Open one known Tutorial page in the active Dash browser."""
    value = str(lesson or "").strip()
    allowed = {"tutorial-welcome", "tutorial-equations", "tutorial-adg1-explorer"}
    if value not in allowed:
        raise ValueError("unknown tutorial lesson")
    event = publish_tab_request(
        "tutorial",
        f"Opening tutorial lesson: {value.removeprefix('tutorial-').replace('-', ' ')}",
        type="tutorial",
        operation="open_lesson",
        lesson=value,
    )
    return {"status": "requested", "activity_id": event["id"], "lesson": value}


def inspect_compile() -> dict[str, Any]:
    """Return the broker's current compile record without requiring a live browser."""
    record = broker_jobs().get("compile")
    return {"status": "idle" if not record else str(record.get("state") or "unknown"), "job": record}


def inspect_runs(case: str | None = None) -> dict[str, Any]:
    """Return broker-owned SCM records, optionally narrowed to one case."""
    runs = {
        str(record.get("run_id") or record.get("job_id") or ""): _scm_run_view(record)
        for record in _JOB_STORE.list_kind("scm")
        if str(record.get("visibility") or "user") == "user"
    }
    if case is None or not str(case).strip():
        return {"runs": runs}
    case_name = _validate_case(str(case))
    matching = [record for record in runs.values() if record.get("case") == case_name]
    latest = max(matching, key=lambda record: record.get("updated_at") or 0, default=None)
    return {"case": case_name, "run": latest}


def inspect_dashboard(tab: str | None = None) -> dict[str, Any]:
    """Describe agent-ready tabs together with lightweight live choices/status."""
    manifest = describe_tabs(tab)
    manifest["dashboard"] = dashboard_registry.dashboard_status()
    for entry in manifest["tabs"]:
        tab_name = entry["tab"]
        entry["operations"] = [
            operation
            for operation in entry["operations"]
            if (tab_name, str(operation.get("name") or "")) in _BROWSER_HANDOFF_OPERATIONS
        ]
        if tab_name == "tutorial":
            entry["context"] = {
                "lessons": ["tutorial-welcome", "tutorial-equations", "tutorial-adg1-explorer"]
            }
        elif tab_name == "misc":
            from dash_app.misc_tab.registry import discover_subtabs

            entry["context"] = {
                "subtabs": [
                    {"slug": subtab.slug, "title": subtab.title}
                    for subtab in discover_subtabs()
                ],
            }
        elif tab_name == "reports":
            from dash_app.reports_tab.catalog import discover_reports

            entry["context"] = {
                "root": "doc/reports",
                "reports": [
                    {
                        "id": report.report_id,
                        "title": report.title,
                        "summary": report.summary,
                        "tags": list(report.tags),
                        "created_at": report.created_at,
                    }
                    for report in discover_reports()
                ],
            }
        elif tab_name == "compile":
            entry["context"] = inspect_compile()
        elif tab_name == "run":
            entry["context"] = {
                "available_cases": sorted(
                    path.name.removesuffix("_model.in")
                    for path in (REPO_ROOT / "input" / "case_setups").glob("*_model.in")
                ),
                **inspect_runs(),
            }
        elif tab_name == "tune":
            entry["context"] = {
                "configs": [
                    str(item.get("value") or "")
                    for item in available_tunable_configs()
                    if str(item.get("value") or "")
                ],
                "active": inspect_tuning(),
            }
        elif tab_name == "plots":
            entry["context"] = {
                "output_cases": sorted(scan_output_cases([str(REPO_ROOT / "output")])),
                "time_presets": ["loss", "pyplotgen"],
            }
    return manifest


def invoke_dashboard(tab: str, operation: str, arguments: dict[str, Any] | None = None) -> dict[str, Any]:
    """Request a visible browser handoff through the deprecated generic API."""
    key = (str(tab or "").strip(), str(operation or "").strip())
    if key not in _BROWSER_HANDOFF_OPERATIONS:
        raise ValueError(
            "invoke_dashboard is compatibility-only for visible browser handoff; "
            "use a typed MCP domain operation for scientific work or job state"
        )
    dashboard = dashboard_registry.dashboard_status()
    if dashboard.get("status") != "available":
        return {
            "tab": str(tab),
            "operation": str(operation),
            "error": {
                "code": "DASHBOARD_UNAVAILABLE",
                "message": "no live dashboard is registered for this checkout",
            },
            "dashboard": dashboard,
        }
    result = invoke_tab_operation(tab, operation, arguments)
    return {"tab": str(tab), "operation": str(operation), **dict(result or {})}


def _register_dashboard_operations() -> None:
    """Register top-level tabs once; handlers retain their normal tab-owned validation."""
    tabs = {
        "tutorial": ("Tutorial", "Interactive teaching demos and CLUBB-equation walkthroughs."),
        "reports": ("Reports", "Immutable static investigation bundles published under doc/reports."),
        "misc": ("Misc", "Interactive investigations and focused diagnostics."),
        "compile": ("Compile", "Configured CLUBB builds, logs, and build status."),
        "run": ("Run", "Checked-in SCM case execution and live output."),
        "profile": ("Profile", "Concurrent CLUBB timing sweeps and performance plots."),
        "tune": ("Tune", "Structured parameter tuning, status, and result reruns."),
        "plots": ("Plots", "Profile and diagnostic plots backed by CLUBB NetCDF output."),
    }
    for tab_name, (title, description) in tabs.items():
        register_tab(tab_name, title, description)
        register_operation(
            tab_name,
            "navigate",
            f"Open the {title} tab in the active Dash browser.",
            {"type": "object", "properties": {}},
            lambda _arguments, target=tab_name: open_dashboard(target),
        )

    register_operation(
        "reports",
        "open_report",
        "Open one published static report by immutable report id.",
        {"type": "object", "required": ["report_id"], "properties": {"report_id": {"type": "string"}}},
        lambda arguments: open_static_report(arguments.get("report_id", "")),
    )
    register_operation(
        "misc",
        "open_report",
        "Open one registered investigation report in Misc by stable slug.",
        {"type": "object", "required": ["slug"], "properties": {"slug": {"type": "string"}}},
        lambda arguments: open_note(arguments.get("slug", "")),
    )
    register_operation(
        "tutorial",
        "open_lesson",
        "Open a Tutorial page without relying on a browser click.",
        {
            "type": "object",
            "required": ["lesson"],
            "properties": {
                "lesson": {
                    "enum": ["tutorial-welcome", "tutorial-equations", "tutorial-adg1-explorer"]
                }
            },
        },
        lambda arguments: open_tutorial_lesson(arguments.get("lesson", "")),
    )
    register_operation(
        "compile",
        "start",
        "Start a normal CLUBB compile and show its native console.",
        {
            "type": "object",
            "properties": {
                "debug": {"type": "boolean", "default": True},
                "python_bindings": {"type": "boolean", "default": False},
                "fresh": {"type": "boolean", "default": False},
                "gptl": {"type": "boolean", "default": False},
            },
        },
        lambda arguments: compile_clubb(
            debug=arguments.get("debug", True),
            python_bindings=arguments.get("python_bindings", False),
            fresh=arguments.get("fresh", False),
            gptl=arguments.get("gptl", False),
        ),
    )
    register_operation(
        "compile",
        "stop",
        "Request a stop for the current broker-owned compile.",
        {"type": "object", "properties": {}},
        lambda _arguments: stop_compile(),
    )
    register_operation(
        "compile",
        "inspect",
        "Read the current broker-owned compile state and log metadata.",
        {"type": "object", "properties": {}},
        lambda _arguments: inspect_compile(),
    )
    register_operation(
        "run",
        "start",
        "Run one checked-in SCM case with safe key=value overrides.",
        {
            "type": "object",
            "required": ["case"],
            "properties": {
                "case": {"type": "string"},
                "overrides": {"type": "string", "default": ""},
                "config": {"type": "string", "default": "default"},
                "stats_file": {"type": "string", "default": DEFAULT_STATS_NAME},
            },
        },
        lambda arguments: run_scm(
            arguments.get("case", ""),
            arguments.get("overrides", ""),
            arguments.get("config", "default"),
            arguments.get("stats_file", DEFAULT_STATS_NAME),
        ),
    )
    register_operation(
        "run",
        "stop",
        "Stop one broker-owned SCM case, narrowed by output directory if needed.",
        {
            "type": "object",
            "required": ["case"],
            "properties": {
                "case": {"type": "string"},
                "output_dir": {"type": "string"},
            },
        },
        lambda arguments: stop_run(
            arguments.get("case", ""), arguments.get("output_dir")
        ),
    )
    register_operation(
        "run",
        "inspect",
        "Read broker-owned SCM case records, optionally for one case.",
        {"type": "object", "properties": {"case": {"type": "string"}}},
        lambda arguments: inspect_runs(arguments.get("case")),
    )
    register_operation(
        "plots",
        "set_view",
        "Open selected profile variables with custom physical time controls and optional validated benchmark overlays.",
        {
            "type": "object",
            "required": ["case", "variables"],
            "properties": {
                "case": {"type": "string"},
                "variables": {"type": "array", "items": {"type": "string"}},
                "run_id": {
                    "type": "string",
                    "description": "Optional immutable SCM run selection; preferred over mutable output/.",
                },
                "output_dir": {
                    "type": "string",
                    "description": "Optional output directory below repository output/; use when run_id is unavailable.",
                },
                "output_dirs": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                    "maxItems": 8,
                    "description": "Ordered output directories for native Plot-tab comparison mode.",
                },
                "time_start_seconds": {"type": "number"},
                "average_minutes": {"type": "number"},
                "window_preset": {"enum": ["loss", "pyplotgen"]},
                "benchmark_sources": {
                    "type": "array",
                    "items": {"type": "string", "enum": ["sam", "coamps"]},
                    "description": "Optional exact overlay selection; omitted preserves the current/default UI selection.",
                },
            },
        },
        lambda arguments: plot_profiles(
            arguments.get("case", ""),
            arguments.get("variables") or [],
            run_id=arguments.get("run_id"),
            output_dir=arguments.get("output_dir"),
            output_dirs=arguments.get("output_dirs") if "output_dirs" in arguments else None,
            time_start_seconds=arguments.get("time_start_seconds"),
            average_minutes=arguments.get("average_minutes"),
            window_preset=arguments.get("window_preset"),
            benchmark_sources=arguments.get("benchmark_sources") if "benchmark_sources" in arguments else None,
        ),
    )
    register_operation(
        "plots",
        "add_budget",
        "Add one validated budget plot card, defaulting to the WP2 budget group.",
        {
            "type": "object",
            "required": ["case"],
            "properties": {
                "case": {"type": "string"},
                "budget_group": {
                    "type": "string",
                    "default": "wp2",
                    "description": "Validated budget group from the selected output; wp2 is the default.",
                },
                "run_id": {
                    "type": "string",
                    "description": "Optional immutable SCM run selection; preferred over mutable output/.",
                },
                "output_dir": {
                    "type": "string",
                    "description": "Optional output directory below repository output/; use when run_id is unavailable.",
                },
                "time_start_seconds": {"type": "number"},
                "average_minutes": {"type": "number"},
                "window_preset": {"enum": ["loss", "pyplotgen"]},
            },
        },
        lambda arguments: add_budget_plot(
            arguments.get("case", ""),
            arguments.get("budget_group", "wp2"),
            run_id=arguments.get("run_id"),
            output_dir=arguments.get("output_dir"),
            time_start_seconds=arguments.get("time_start_seconds"),
            average_minutes=arguments.get("average_minutes"),
            window_preset=arguments.get("window_preset"),
        ),
    )
    register_operation(
        "plots",
        "list",
        "List the active Plot cards with stable IDs and validated family selections.",
        {"type": "object", "properties": {}},
        lambda _arguments: list_plot_instances(),
    )
    register_operation(
        "plots",
        "remove",
        "Remove one active Plot card by its stable ID.",
        {
            "type": "object",
            "required": ["plot_id"],
            "properties": {
                "plot_id": {
                    "type": "integer",
                    "minimum": 0,
                    "description": "Stable ID returned by plots.list.",
                },
            },
        },
        lambda arguments: remove_plot_instance(arguments.get("plot_id")),
    )
    register_operation(
        "plots",
        "export_profiles",
        "Save selected profile figures as PNG artifacts under output/agent_exports.",
        {
            "type": "object",
            "required": ["case", "variables"],
            "properties": {
                "case": {"type": "string"},
                "variables": {"type": "array", "items": {"type": "string"}},
                "time_start_seconds": {"type": "number"},
                "average_minutes": {"type": "number"},
                "window_preset": {"enum": ["loss", "pyplotgen"]},
            },
        },
        lambda arguments: save_profile_png(
            arguments.get("case", ""),
            arguments.get("variables") or [],
            time_start_seconds=arguments.get("time_start_seconds"),
            average_minutes=arguments.get("average_minutes"),
            window_preset=arguments.get("window_preset"),
        ),
    )
    register_operation(
        "tune",
        "launch",
        "Launch a bounded structured tuning request and populate the native Tune controls.",
        {
            "type": "object",
            "required": ["cases", "parameter_ranges"],
            "properties": {
                "cases": {"type": "array"},
                "parameter_ranges": {"type": "array"},
                "fields": {"type": "array"},
                "config": {"type": "string", "default": "default"},
                "strategy": {"type": "string", "default": "random"},
                "max_samples": {"type": "integer", "default": 12},
                "resolve_spacing": {"type": "number", "default": 0.1},
                "simann_max_iters": {"type": "integer", "default": 200},
                "simann_initial_temp": {"type": "number", "default": 1.0},
                "simann_final_temp": {"type": "number", "default": 1.0e-12},
                "adam_max_updates": {"type": "integer", "default": 100},
                "adam_learning_rate": {"type": "number", "default": 0.01},
                "adam_perturbation": {"type": "number", "default": 0.05},
                "adam_spsa_pairs": {"type": "integer", "default": 2},
                "batch_size": {"type": "integer", "default": 1},
                "max_workers": {"type": "integer", "default": 1},
                "loss_mode": {"type": "string"},
                "aggregation_mode": {"type": "string"},
                "aggregation_weights": {"type": "array", "items": {"type": "number"}, "minItems": 4, "maxItems": 4},
                "time_window_aggregation_scope": {"type": "string", "enum": ["overall", "by_case"]},
                "overrides": {"type": "string", "default": ""},
            },
        },
        lambda arguments: launch_tuning(
            arguments.get("cases") or [],
            arguments.get("parameter_ranges") or [],
            arguments.get("fields"),
            config=arguments.get("config", "default"),
            strategy=arguments.get("strategy", "random"),
            max_samples=arguments.get("max_samples", 12),
            resolve_spacing=arguments.get("resolve_spacing", 0.1),
            simann_max_iters=arguments.get("simann_max_iters", 200),
            simann_initial_temp=arguments.get("simann_initial_temp", 1.0),
            simann_final_temp=arguments.get("simann_final_temp", 1.0e-12),
            adam_max_updates=arguments.get("adam_max_updates", 100),
            adam_learning_rate=arguments.get("adam_learning_rate", 0.01),
            adam_perturbation=arguments.get("adam_perturbation", 0.05),
            adam_spsa_pairs=arguments.get("adam_spsa_pairs", 2),
            batch_size=arguments.get("batch_size", 1),
            max_workers=arguments.get("max_workers", 1),
            loss_mode=arguments.get("loss_mode", DEFAULT_LOSS_MODE),
            aggregation_mode=arguments.get("aggregation_mode", DEFAULT_AGGREGATION_MODE),
            aggregation_weights=arguments.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS),
            time_window_aggregation_scope=arguments.get("time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE),
            overrides=arguments.get("overrides", ""),
        ),
    )
    register_operation(
        "tune",
        "stop",
        "Request a graceful stop for the active tuning worker.",
        {"type": "object", "properties": {}},
        lambda _arguments: stop_tuning(),
    )
    register_operation(
        "tune",
        "inspect",
        "Read the active Tune request, status, results path, and log path.",
        {"type": "object", "properties": {}},
        lambda _arguments: inspect_tuning(),
    )
    register_operation(
        "tune",
        "list_workspaces",
        "List saved Tune workspaces, revisions, status, timestamps, and disk use.",
        {"type": "object", "properties": {}},
        lambda _arguments: list_tuning_workspaces(),
    )
    register_operation(
        "tune",
        "load_workspace",
        "Load one immutable Tune workspace revision and its request, status, and results.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: load_tuning_workspace(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "new_revision",
        "Clone a stopped/finished revision into an editable unstarted revN draft.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: create_tuning_revision(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "restart_revision",
        "Create and start a fresh exact <revision>_restartN attempt.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: restart_tuning_revision(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "reset_revision",
        "Delete one inactive revision's generated data and restore its editable draft.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: reset_tuning_revision(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "start_draft_revision",
        "Start an unchanged saved draft revision.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: start_tuning_draft_revision(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "resume_revision",
        "Continue one gracefully stopped Tune revision from its saved scheduler checkpoint.",
        {
            "type": "object",
            "required": ["workspace_id", "revision_id"],
            "properties": {"workspace_id": {"type": "string"}, "revision_id": {"type": "string"}},
        },
        lambda arguments: resume_tuning_revision(arguments.get("workspace_id", ""), arguments.get("revision_id", "")),
    )
    register_operation(
        "tune",
        "run_leaderboard",
        "Rerun up to 16 current leaderboard rows through window-loss or complete SCM mode.",
        {
            "type": "object",
            "properties": {"mode": {"enum": ["window", "complete"], "default": "window"}, "max_results": {"type": "integer", "default": 16}},
        },
        lambda arguments: run_tuning_loss(arguments.get("mode", "window"), arguments.get("max_results", 16)),
    )


_register_dashboard_operations()


def allowed_action_names() -> list[str]:
    """Return the external browser-handoff actions accepted by the broker."""
    return [
        "inspect_dashboard",
        "invoke_dashboard",
    ]


def dispatch(action: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Dispatch only registered semantic operations; never execute raw commands."""
    if action == "stop_all_broker_work":
        return stop_all_broker_work(reason=str(payload.get("reason") or "dashboard manager is stopping"))
    if _BROKER_SHUTTING_DOWN.is_set():
        raise ValueError("the dashboard runtime broker is shutting down")
    # These typed domain mutations are accepted only through the gateway's
    # private internal-action path.  Keeping their validation and launch here
    # makes the detached broker own every watcher thread, so an MCP adapter
    # restart cannot orphan a Compile, SCM, Tune, or artifact job.
    if action == "domain_submit_compile":
        return submit_compile(CompileRequest.model_validate(payload.get("request") or {}))
    if action == "domain_submit_scm_run":
        return submit_scm_run(
            ScmRunRequest.model_validate(payload.get("request") or {}),
            origin=str(payload.get("submission_origin") or "unknown"),
        )
    if action == "domain_submit_scm_batch":
        return submit_scm_batch(
            ScmRunBatchRequest.model_validate(payload.get("request") or {}),
            native_overrides=payload.get("native_overrides") if "native_overrides" in payload else None,
            native_cli_options=payload.get("native_cli_options") if "native_cli_options" in payload else None,
            origin=str(payload.get("submission_origin") or "unknown"),
        )
    if action == "domain_submit_tune":
        return submit_tune(TuneRequest.model_validate(payload.get("request") or {}))
    if action == "domain_submit_leaderboard_rerun":
        return submit_leaderboard_rerun(
            LeaderboardRerunRequest.model_validate(payload.get("request") or {})
        )
    if action == "domain_create_profile_artifact":
        return create_profile_artifact(
            ProfileArtifactRequest.model_validate(payload.get("request") or {})
        )
    if action == "domain_cancel_job":
        return cancel_job(str(payload.get("job_id") or ""))
    if action == "domain_cancel_all_scm":
        return cancel_all_scm_runs()
    if action == "clear_terminal_scm_session":
        return clear_terminal_scm_session()
    if action == "domain_create_tune_revision":
        return create_tuning_revision(
            str(payload.get("workspace_id") or ""),
            str(payload.get("revision_id") or ""),
        )
    if action == "domain_restart_tune_revision":
        return restart_tuning_revision(
            str(payload.get("workspace_id") or ""),
            str(payload.get("revision_id") or ""),
        )
    if action == "domain_reset_tune_revision":
        return reset_tuning_revision(
            str(payload.get("workspace_id") or ""),
            str(payload.get("revision_id") or ""),
        )
    if action == "domain_start_tune_draft_revision":
        return start_tuning_draft_revision(
            str(payload.get("workspace_id") or ""),
            str(payload.get("revision_id") or ""),
        )
    if action == "domain_resume_tune_revision":
        return resume_tuning_revision(
            str(payload.get("workspace_id") or ""),
            str(payload.get("revision_id") or ""),
        )
    if action == "domain_rename_tune_workspace":
        return rename_tuning_workspace(
            str(payload.get("workspace_id") or ""),
            str(payload.get("display_name") or ""),
        )
    if action == "domain_delete_tune_workspace":
        return delete_tuning_workspace(str(payload.get("workspace_id") or ""))
    if action == "inspect_dashboard":
        return inspect_dashboard(payload.get("tab"))
    if action == "invoke_dashboard":
        return invoke_dashboard(
            payload.get("tab", ""),
            payload.get("operation", ""),
            payload.get("arguments") or {},
        )
    if action == "compile_clubb":
        return compile_clubb(
            debug=payload.get("debug", True),
            python_bindings=payload.get("python_bindings", False),
            fresh=payload.get("fresh", False),
            gptl=payload.get("gptl", False),
        )
    if action == "launch_compile_request":
        return launch_compile_request(dict(payload.get("options") or {}), env_id=payload.get("env_id", "current"))
    if action == "launch_profile_request":
        return launch_profile_request(dict(payload.get("settings") or {}))
    if action == "launch_pyplotgen_request":
        return launch_pyplotgen_request(list(payload.get("output_dirs") or []))
    if action == "stop_pyplotgen_request":
        return stop_pyplotgen_request()
    if action == "stop_profile":
        return stop_profile()
    if action == "launch_rebuild_request":
        return launch_rebuild_request(
            list(payload.get("builds") or []),
            dict(payload.get("discovery") or {}),
            str(payload.get("label") or "selected builds"),
        )
    if action == "run_scm":
        return run_scm(
            payload.get("case", ""),
            payload.get("overrides", ""),
            payload.get("config", "default"),
            payload.get("stats_file", DEFAULT_STATS_NAME),
        )
    if action == "plot_profiles":
        return plot_profiles(
            payload.get("case", ""),
            payload.get("variables") or [],
            payload.get("time_seconds"),
            time_start_seconds=payload.get("time_start_seconds"),
            average_minutes=payload.get("average_minutes"),
            window_preset=payload.get("window_preset"),
            benchmark_sources=payload.get("benchmark_sources") if "benchmark_sources" in payload else None,
        )
    if action == "save_profile_png":
        return save_profile_png(
            payload.get("case", ""),
            payload.get("variables") or [],
            time_start_seconds=payload.get("time_start_seconds"),
            average_minutes=payload.get("average_minutes"),
            window_preset=payload.get("window_preset"),
        )
    if action == "open_dashboard":
        return open_dashboard(payload.get("tab", ""))
    if action == "open_note":
        return open_note(payload.get("slug", ""))
    if action == "launch_tuning":
        return launch_tuning(
            payload.get("cases") or [],
            payload.get("parameter_ranges") or [],
            payload.get("fields"),
            config=payload.get("config", "default"),
            strategy=payload.get("strategy", "random"),
            max_samples=payload.get("max_samples", 12),
            resolve_spacing=payload.get("resolve_spacing", 0.1),
            simann_max_iters=payload.get("simann_max_iters", 200),
            simann_initial_temp=payload.get("simann_initial_temp", 1.0),
            simann_final_temp=payload.get("simann_final_temp", 1.0e-12),
            adam_max_updates=payload.get("adam_max_updates", 100),
            adam_learning_rate=payload.get("adam_learning_rate", 0.01),
            adam_perturbation=payload.get("adam_perturbation", 0.05),
            adam_spsa_pairs=payload.get("adam_spsa_pairs", 2),
            batch_size=payload.get("batch_size", 1),
            max_workers=payload.get("max_workers", 1),
            loss_mode=payload.get("loss_mode", DEFAULT_LOSS_MODE),
            aggregation_mode=payload.get("aggregation_mode", DEFAULT_AGGREGATION_MODE),
            aggregation_weights=payload.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS),
            time_window_aggregation_scope=payload.get("time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE),
            overrides=payload.get("overrides", ""),
        )
    if action == "launch_tuning_request":
        return launch_tuning_request(payload.get("request") or {})
    if action == "list_tuning_workspaces":
        return list_tuning_workspaces()
    if action == "list_tuning_workspace_activity":
        return list_tuning_workspace_activity()
    if action == "load_tuning_workspace":
        return load_tuning_workspace(payload.get("workspace_id", ""), payload.get("revision_id", ""))
    if action == "create_tuning_revision":
        return create_tuning_revision(
            payload.get("workspace_id", ""),
            payload.get("revision_id", ""),
            restart=bool(payload.get("restart", False)),
        )
    if action == "start_tuning_draft_revision":
        return start_tuning_draft_revision(payload.get("workspace_id", ""), payload.get("revision_id", ""))
    if action == "restart_tuning_revision":
        return restart_tuning_revision(payload.get("workspace_id", ""), payload.get("revision_id", ""))
    if action == "reset_tuning_revision":
        return reset_tuning_revision(payload.get("workspace_id", ""), payload.get("revision_id", ""))
    if action == "resume_tuning_revision":
        return resume_tuning_revision(payload.get("workspace_id", ""), payload.get("revision_id", ""))
    if action == "rename_tuning_workspace":
        return rename_tuning_workspace(payload.get("workspace_id", ""), payload.get("display_name", ""))
    if action == "delete_tuning_workspace":
        return delete_tuning_workspace(payload.get("workspace_id", ""))
    if action == "stop_tuning":
        return stop_tuning()
    if action == "stop_compile":
        return stop_compile()
    if action == "stop_run":
        return stop_run(payload.get("case", ""), payload.get("output_dir"))
    if action == "inspect_tuning":
        return inspect_tuning()
    if action == "run_tuning_loss":
        return run_tuning_loss(payload.get("mode", "window"), payload.get("max_results", 16))
    raise ValueError("unknown action; allowed actions are " + ", ".join(allowed_action_names()))
