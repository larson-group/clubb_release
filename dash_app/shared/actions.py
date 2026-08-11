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
    publish_run_request,
    publish_tab_request,
    publish_tune_request,
    broker_jobs,
    read_activity,
    set_broker_job,
    set_broker_loss_run,
    set_broker_run,
    update_broker_job,
    update_broker_loss_run,
    update_broker_run,
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
from dash_app.run_tab.runtime import get_proc, is_case_active, read_log_increment as read_run_log_increment, record_case_finish, snapshot_active_cases, start_case_process
from dash_app.run_tab.state import MAX_RUN_PROCS
from dash_app.run_tab.state import DEFAULT_STATS_NAME
from dash_app.plot_tab.plot_types.profile_plot import PLOT as profile_plot
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
from dash_app.shared.runtime import atomic_write_json, exclusive_file_lock, private_path
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
_TUNE_MAX_RESOLVE_SAMPLES = 1000
_TUNE_MAX_BATCH_SIZE = 64
_PROFILE_WINDOW_PRESETS = {"loss", "pyplotgen"}
_PROFILE_EXPORT_DIR = REPO_ROOT / "output" / "agent_exports"
_JOB_STORE = JobStore()
_ARTIFACT_STORE = ArtifactStore()
_BATCH_QUEUE_LOCK = threading.Lock()
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


def _pid_is_alive(pid: Any) -> bool:
    """Check a recovered child without relying on an in-memory ``Popen``."""
    try:
        numeric_pid = int(pid)
        if numeric_pid < 1:
            return False
        os.kill(numeric_pid, 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


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


def _public_scm_output_dir(case: str, run_id: str) -> Path:
    """Return the default durable output location for a typed MCP SCM run.

    A run ID is already minted by the durable job store, so it gives each
    request a stable, discoverable default directory. Public requests may
    override this with a validated path below the repository output root.
    """
    safe_case = str(case or "").strip()
    if not _CASE_RE.fullmatch(safe_case):
        raise ValueError("typed SCM run is missing a valid case name")
    safe_run_id = str(run_id or "").strip()
    if not re.fullmatch(r"run_[A-Za-z0-9_-]+", safe_run_id):
        raise ValueError("typed SCM run is missing a valid run identifier")
    return _public_scm_output_root() / safe_case / safe_run_id


def _public_scm_batch_output_dir(batch_id: str) -> Path:
    """Return the default single flat, broker-owned directory for one batch."""
    safe_batch = str(batch_id or "").strip()
    if not _BATCH_ID_RE.fullmatch(safe_batch):
        raise ValueError("SCM batch is missing a valid batch identifier")
    return _public_scm_output_root() / safe_batch


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
        )
        if record.get(key) is not None
    }
    case = record.get("case") or (record.get("request") or {}).get("case")
    if case:
        summary["case"] = case
    return summary


def _write_scm_batch_manifest(batch: dict[str, Any]) -> str:
    """Atomically publish the compact parent/child batch manifest."""
    output_directory = Path(str(batch.get("output_directory") or ""))
    output_directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "batch_id": batch.get("batch_id"),
        "job_id": batch.get("job_id"),
        "request_id": batch.get("request_id"),
        "state": batch.get("state"),
        "request": batch.get("request") or {},
        "output_directory": str(output_directory),
        "children": batch.get("children") or [],
        "created_at_unix_seconds": batch.get("created_at_unix_seconds"),
        "updated_at_unix_seconds": batch.get("updated_at_unix_seconds"),
        "finished_at_unix_seconds": batch.get("finished_at_unix_seconds"),
    }
    path = output_directory / "batch_manifest.json"
    atomic_write_json(path, payload)
    return str(path)


def _batch_record(job_id: str) -> dict[str, Any] | None:
    record = _JOB_STORE.get(job_id)
    return record if record and record.get("kind") == "scm_batch" else None


def _batch_terminal(state: str) -> bool:
    return state in {"finished", "partial_failure", "error", "cancelled"}


def _batch_state(children: list[dict[str, Any]]) -> str:
    states = {str(child.get("state") or "queued") for child in children}
    if states & {"starting", "running", "stopping"}:
        return "running"
    if states & {"queued", "submitting"}:
        return "queued"
    if states and all(str(child.get("state")) == "finished" and child.get("returncode") == 0 for child in children):
        return "finished"
    return "partial_failure"


def _refresh_scm_batch(batch_job_id: str, child_job_id: str | None = None) -> dict[str, Any] | None:
    """Mirror one child transition into the durable parent and public manifest."""
    batch = _batch_record(batch_job_id)
    if batch is None:
        return None
    children = []
    for item in batch.get("children") or []:
        child = dict(item or {})
        job_id = str(child.get("job_id") or "")
        current = _JOB_STORE.get(job_id)
        if current and (child_job_id is None or job_id == child_job_id or current.get("state") != child.get("state")):
            child = _batch_child_summary(current)
        children.append(child)
    state = _batch_state(children)
    updates: dict[str, Any] = {"children": children, "state": state}
    if _batch_terminal(state):
        updates["finished_at_unix_seconds"] = time.time()
    updated = _JOB_STORE.update(batch_job_id, **updates) or batch | updates
    updated["children"] = children
    updated["state"] = state
    manifest_path = _write_scm_batch_manifest(updated)
    if manifest_path != updated.get("public_manifest_path"):
        updated = _JOB_STORE.update(batch_job_id, public_manifest_path=manifest_path) or updated | {"public_manifest_path": manifest_path}
    if _batch_terminal(state):
        _ARTIFACT_STORE.release(batch_job_id)
    return updated


def _write_public_run_manifest(output_directory: Path, payload: dict[str, Any]) -> Path:
    """Write compact durable provenance beside a typed MCP run's NetCDF data."""
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = output_directory / "run_manifest.json"
    atomic_write_json(manifest_path, payload)
    return manifest_path


def _seal_public_run_manifest(output_directory: str | None, returncode: int) -> str | None:
    """Record terminal status without making private artifact retention public."""
    if not output_directory:
        return None
    manifest_path = Path(output_directory) / "run_manifest.json"
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            return None
        execution = dict(payload.get("execution") or {})
        execution.update(
            {
                "state": "finished" if returncode == 0 else "error",
                "returncode": int(returncode),
                "finished_at_unix_seconds": time.time(),
            }
        )
        payload["execution"] = execution
        _write_public_run_manifest(Path(output_directory), payload)
        return str(manifest_path)
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


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


def _watch_run(case: str, proc_data: dict[str, Any], job_id: str | None = None) -> None:
    offset = 0
    pid = proc_data.get("pid")
    while True:
        chunk, offset = read_run_log_increment(proc_data.get("log"), offset)
        if chunk:
            # The Run tab owns live SCM log presentation. Repeating those
            # chunks in the global activity drawer turns ordinary iteration
            # progress into noisy snapshots, while the durable log offset
            # below is still needed for console reconnect/recovery.
            update_broker_run(case, log_tail=_read_log_tail(proc_data.get("log")), log_offset=offset)
            if job_id:
                _JOB_STORE.update(job_id, progress={"log_offset": offset})
        proc = get_proc(pid)
        # A detached broker can restart while SCM continues.  In that case
        # there is no local Popen, but PID liveness is enough to retain the
        # log stream, stop control, and active-job guard.  Once it exits the
        # POSIX exit status is no longer observable, so report that recovery
        # limitation explicitly rather than inventing a successful result.
        returncode = proc.poll() if proc is not None else (None if _pid_is_alive(pid) else 1)
        if returncode is not None:
            requested_stop = str(
                ((broker_jobs().get("runs") or {}).get(case) or {}).get("state") or ""
            ) == "stopping"
            record_case_finish(case, pid, returncode)
            cleanup_temp_files(proc_data.get("temp_files") or [])
            recovered_exit = proc is None
            broker_state = "finished" if returncode == 0 else ("stopped" if requested_stop else "error")
            job_state = "finished" if returncode == 0 else ("cancelled" if requested_stop else "error")
            update_broker_run(
                case,
                state=broker_state,
                returncode=returncode,
                finished_at=time.time(),
                log_offset=offset,
            )
            publish_event(
                "run",
                f"{case} finished" if returncode == 0 else (f"{case} stopped" if requested_stop else f"{case} failed"),
                f"exit code {returncode}",
                status="success" if returncode == 0 else ("info" if requested_stop else "error"),
            )
            if job_id:
                child_record = _JOB_STORE.get(job_id) or {}
                batch_job_id = str(child_record.get("batch_job_id") or "")
                result_summary = _write_run_result_summary(job_id, proc_data.get("output_directory"), returncode)
                public_manifest = _seal_public_run_manifest(proc_data.get("output_directory"), returncode)
                _JOB_STORE.update(
                    job_id,
                    state=job_state,
                    returncode=returncode,
                    finished_at_unix_seconds=time.time(),
                    log_tail=_read_log_tail(proc_data.get("log")),
                    recovery_note="exit code unavailable after broker recovery" if recovered_exit else None,
                    result_summary_path=result_summary,
                    public_manifest_path=public_manifest,
                )
                _ARTIFACT_STORE.release(job_id)
                if batch_job_id:
                    batch_state = _refresh_scm_batch(batch_job_id, job_id)
                    if batch_state and any(str(item.get("state") or "") == "queued" for item in batch_state.get("children") or []):
                        _ensure_scm_batch_watcher(batch_job_id)
                        _start_queued_scm_batch_children(batch_job_id)
            return
        time.sleep(0.75)


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


def _plot_output_selection(case_name: str, run_id: str | None, output_dir: str | None) -> tuple[str, list[str] | None]:
    """Resolve the immutable run or explicit output selection for Plot actions."""
    selected_run_id = str(run_id or "").strip()
    requested_output = str(output_dir or "").strip()
    selected_output = _resolve_mcp_output_dir(requested_output) if requested_output else None

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
    batch_size: int,
    max_workers: int,
    max_samples_limit: int | None = _TUNE_MAX_RANDOM_SAMPLES,
) -> dict[str, Any]:
    name = str(strategy or "random").strip().lower()
    if name not in VALID_STRATEGY_NAMES:
        raise ValueError("strategy must be random, resolve, or simann")
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


def _active_scm_processes() -> set[tuple[str, str]]:
    """Return the union of broker and legacy in-memory SCM children by PID."""
    active: set[tuple[str, str]] = set()
    for case_name, item in (broker_jobs().get("runs") or {}).items():
        item = dict(item or {})
        if str(item.get("state") or "") not in {"running", "stopping"}:
            continue
        pid = (item.get("proc_data") or {}).get("pid")
        active.add(("pid", str(pid)) if pid else ("broker", str(case_name)))
    for case_name, item in snapshot_active_cases().items():
        pid = (item or {}).get("pid")
        active.add(("pid", str(pid)) if pid else ("native", str(case_name)))
    return active


def _assert_scm_admission(case_name: str) -> None:
    """Apply the single cross-client case and global-concurrency rule."""
    broker_record = dict((broker_jobs().get("runs") or {}).get(case_name) or {})
    if str(broker_record.get("state") or "") in {"running", "stopping"} or is_case_active(case_name):
        raise ValueError(f"{case_name} is already running; wait for it to finish before starting another run")
    if len(_active_scm_processes()) >= MAX_RUN_PROCS:
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
    allowed = {"multicol", "batch_size", "max_iters", "debug", "dt_main", "dt_rad", "tout", "out_dir", "extra_args"}
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
    extra = raw.get("extra_args") or []
    if isinstance(extra, str):
        extra = shlex.split(extra)
    if not isinstance(extra, (list, tuple)):
        raise ValueError("extra_args must be a list of command arguments")
    tokens = [str(item).strip() for item in extra if str(item).strip()]
    if len(tokens) > 64 or any(len(item) > 512 or "\x00" in item for item in tokens):
        raise ValueError("invalid extra run arguments")
    if tokens:
        normalized["extra_args"] = tokens
    return normalized


def _persist_submission(kind: str, request_id: str, payload: dict[str, Any]):
    try:
        return _JOB_STORE.submit(kind, request_id, payload)
    except SubmissionConflict as exc:
        raise ValueError(f"REQUEST_ID_CONFLICT: {exc}") from exc


def _run_manifest_inputs(case: str, stats_file: str, config: str) -> dict[str, Any]:
    """Capture compact checksums for every checked-in input used by an SCM run."""
    paths = {
        "case_setup": REPO_ROOT / "input" / "case_setups" / f"{case}_model.in",
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


def _scm_build_identity() -> dict[str, Any]:
    """Capture the exact compiled executable selected by ``run_scm.py``.

    The regular SCM runner resolves ``install/selected`` before
    ``install/latest``.  Storing its resolved path and digest makes a later
    rebuild distinguishable without relying on mutable UI selection or a log.
    """
    selected = REPO_ROOT / "install" / "selected"
    latest = REPO_ROOT / "install" / "latest"
    install = selected if os.path.lexists(selected) else latest
    resolved_install = install.resolve(strict=False)
    executable = resolved_install / "clubb_standalone"
    return {
        "install_selector": str(install),
        "install_directory": str(resolved_install),
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
        )
        _ARTIFACT_STORE.activate(record["job_id"])
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
    output_dir: Path | None = None,
    batch_job_id: str | None = None,
    native_overrides: dict[str, Any] | None = None,
    native_cli_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Start one SCM child after the cross-client admission lock is held."""
    try:
        _assert_scm_admission(request.case)
    except ValueError as exc:
        code = "CASE_ALREADY_RUNNING" if "already running" in str(exc) else "SCM_CONCURRENCY_LIMIT"
        _JOB_STORE.update(record["job_id"], state="rejected", error={"code": code, "message": str(exc)})
        raise
    output_dir = output_dir or _public_scm_output_dir(request.case, record["run_id"])
    try:
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
                "inputs": _run_manifest_inputs(request.case, stats_file, request.config),
                "execution": {
                    "state": "planned",
                    "runner": "run_scripts/run_scm.py",
                    "build_identity": _scm_build_identity(),
                },
                "output_checksums": {"state": "pending", "result_resource": f"clubb-artifact://{record['job_id']}/execution_result.json"},
            },
        )
        _ARTIFACT_STORE.activate(record["job_id"])
        public_manifest = None
        if batch_job_id is None:
            public_manifest = _write_public_run_manifest(
                output_dir,
                {
                    "schema_version": 1,
                    "run_id": record["run_id"],
                    "job_id": record["job_id"],
                    "case": request.case,
                    "stats_file": stats_file,
                    "configuration": request.config,
                    "overrides": request.overrides,
                    "run_options": request.run_options.model_dump(exclude_none=True),
                    "output_directory": str(output_dir),
                    "private_artifact_manifest": f"clubb-artifact://{record['job_id']}/manifest.json",
                    "execution": {"state": "planned", "runner": "run_scripts/run_scm.py"},
                },
            )
        _JOB_STORE.update(
            record["job_id"],
            state="starting",
            output_directory=str(output_dir),
            manifest_path=str(manifest),
            public_manifest_path=str(public_manifest) if public_manifest else None,
            batch_job_id=batch_job_id,
        )
        if native_overrides is not None or native_cli_options is not None:
            cli_options = _normalize_dashboard_cli_options(native_cli_options or {})
            cli_options["out_dir"] = str(output_dir)
            result = launch_dashboard_scm_request(
                request.case,
                stats_file,
                request.config,
                native_overrides or {},
                cli_options,
                job_id=record["job_id"],
                run_id=record["run_id"],
            )
        else:
            result = run_scm(
                request.case,
                _typed_override_text(request.overrides),
                request.config,
                stats_file,
                cli_options=request.run_options.model_dump(exclude_none=True),
                output_dir=str(output_dir),
                job_id=record["job_id"],
                run_id=record["run_id"],
            )
        _ARTIFACT_STORE.write_summary(record["job_id"], "submission_result.json", {"state": "started", "runtime": result})
        updated = _JOB_STORE.update(
            record["job_id"], state="running", runtime=result
        ) or record
        return {"status": "started", **updated}
    except Exception as exc:
        _ARTIFACT_STORE.release(record["job_id"])
        _JOB_STORE.update(record["job_id"], state="error", error={"code": "SCM_SUBMISSION_FAILED", "message": str(exc)})
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


def _start_queued_scm_batch_children_locked(batch_job_id: str) -> dict[str, Any] | None:
    """Start as many queued children as admission and the batch limit allow."""
    batch = _batch_record(batch_job_id)
    if batch is None or _batch_terminal(str(batch.get("state") or "")):
        return batch
    request = ScmRunBatchRequest.model_validate(batch.get("request") or {})
    active = sum(
        1
        for child in (batch.get("children") or [])
        if str(child.get("state") or "") in {"starting", "running", "stopping"}
    )
    max_workers = request.max_workers or MAX_RUN_PROCS
    for child_summary in list(batch.get("children") or []):
        if active >= max_workers:
            break
        if str(child_summary.get("state") or "") != "queued":
            continue
        child_job_id = str(child_summary.get("job_id") or "")
        child_record = _JOB_STORE.get(child_job_id)
        if child_record is None:
            continue
        case = str(child_record.get("request", {}).get("case") or child_summary.get("case") or "")
        try:
            _assert_scm_admission(case)
        except ValueError as exc:
            if "maximum active SCM runs" in str(exc):
                break
            _JOB_STORE.update(child_job_id, state="error", error={"code": "SCM_BATCH_CHILD_REJECTED", "message": str(exc)})
            continue
        child_request = ScmRunRequest.model_validate(child_record.get("request") or {})
        _JOB_STORE.update(child_job_id, state="starting")
        try:
            _start_scm_submission(
                child_request,
                child_record,
                child_request.stats_file,
                _canonical_scm_request(child_request)[1],
                output_dir=Path(str(batch["output_directory"])),
                batch_job_id=batch_job_id,
                native_overrides=batch.get("native_overrides"),
                native_cli_options=batch.get("native_cli_options"),
            )
            active += 1
        except ValueError as exc:
            _JOB_STORE.update(child_job_id, state="error", error={"code": "SCM_BATCH_CHILD_FAILED", "message": str(exc)[:500]})
        batch = _batch_record(batch_job_id) or batch
    return _refresh_scm_batch(batch_job_id)


def _start_queued_scm_batch_children(batch_job_id: str) -> dict[str, Any] | None:
    """Serialize queue advancement so completion/recovery cannot double-start a child."""
    if _BROKER_SHUTTING_DOWN.is_set():
        return _batch_record(batch_job_id)
    with _BATCH_QUEUE_LOCK:
        return _start_queued_scm_batch_children_locked(batch_job_id)


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
            if batch is None or _batch_terminal(str(batch.get("state") or "")):
                return
            updated = _start_queued_scm_batch_children(batch_job_id) or batch
            queued = any(str(child.get("state") or "") == "queued" for child in updated.get("children") or [])
            if not queued:
                return
            time.sleep(0.75)
    finally:
        with _BATCH_WATCHERS_LOCK:
            _BATCH_WATCHERS.discard(batch_job_id)


def submit_scm_batch(
    request: ScmRunBatchRequest,
    *,
    native_overrides: dict[str, Any] | None = None,
    native_cli_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Submit one durable multi-case SCM group into one flat public folder."""
    request, settings_resolution = _canonical_scm_batch_request(request)
    normalized_native_overrides = None
    normalized_native_cli_options = None
    native_output_directory: Path | None = None
    if native_overrides is not None or native_cli_options is not None:
        normalized_native_overrides = _normalize_dashboard_overrides(native_overrides or {})
        normalized_native_cli_options = _normalize_dashboard_cli_options(native_cli_options or {})
        # Native Dash owns the Run-tab output field, so preserve its selected
        # location after applying the repository's normal path rules. Public
        # MCP requests do not enter this branch and remain broker-controlled.
        requested_output = normalized_native_cli_options.get("out_dir")
        if requested_output:
            native_output_directory = resolve_output_dir(requested_output)
            normalized_native_cli_options["out_dir"] = str(native_output_directory)
    payload = request.model_dump(mode="json") | {"stats_file": request.stats_file}
    parent, created = _persist_submission("scm_batch", request.request_id, payload)
    if not created:
        refreshed = _refresh_scm_batch(parent["job_id"])
        return {"status": "existing", **(refreshed or parent)}

    batch_id = str(parent["batch_id"])
    output_directory = (
        native_output_directory
        or (_resolve_mcp_output_dir(request.out_dir) if request.out_dir else None)
        or _public_scm_batch_output_dir(batch_id)
    )
    private_manifest = _ARTIFACT_STORE.create_manifest(
        parent["job_id"],
        {
            "job": parent,
            "batch_id": batch_id,
            "requested_batch": payload,
            "output_directory": str(output_directory),
            "settings_resolution": settings_resolution,
            "execution": {"state": "planned", "runner": "run_scripts/run_scm.py"},
        },
    )
    _ARTIFACT_STORE.activate(parent["job_id"])
    children: list[dict[str, Any]] = []
    for case in request.cases:
        child_request = _batch_child_request(request, parent, case)
        child_payload = child_request.model_dump(mode="json") | {"stats_file": request.stats_file}
        child, child_created = _persist_submission("scm", child_request.request_id, child_payload)
        if not child_created:
            existing_batch = str(child.get("batch_job_id") or "")
            if existing_batch and existing_batch != parent["job_id"]:
                raise ValueError("SCM batch child request is already owned by another batch")
        child = _JOB_STORE.update(
            child["job_id"],
            state="queued",
            batch_id=batch_id,
            batch_job_id=parent["job_id"],
            output_directory=str(output_directory),
        ) or child
        children.append(_batch_child_summary(child))

    updated = _JOB_STORE.update(
        parent["job_id"],
        state="queued",
        batch_id=batch_id,
        output_directory=str(output_directory),
        manifest_path=str(private_manifest),
        public_manifest_path=str(output_directory / "batch_manifest.json"),
        children=children,
        native_overrides=normalized_native_overrides,
        native_cli_options=normalized_native_cli_options,
    ) or parent
    _write_scm_batch_manifest(updated)
    started = _start_queued_scm_batch_children(parent["job_id"])
    _ensure_scm_batch_watcher(parent["job_id"])
    return {"status": "started" if started and started.get("state") == "running" else "queued", **(started or updated)}


def submit_scm_run(request: ScmRunRequest) -> dict[str, Any]:
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
    batch = submit_scm_batch(batch_request)
    child = next((item for item in batch.get("children") or [] if item.get("case") == request.case), None)
    if child is None:
        return batch
    return {
        "status": batch.get("status") or batch.get("state"),
        **child,
        "batch_id": batch.get("batch_id"),
        "batch_job_id": batch.get("job_id"),
        "batch_manifest_path": batch.get("public_manifest_path"),
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
        )
        _ARTIFACT_STORE.activate(record["job_id"])
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
        )
        _ARTIFACT_STORE.activate(record["job_id"])
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
    for case, raw_record in (broker_jobs().get("runs") or {}).items():
        record = dict(raw_record or {})
        if str(record.get("state") or "") not in {"running", "stopping"}:
            continue
        proc_data = dict(record.get("proc_data") or {})
        job_id = str(record.get("job_id") or "")
        if not proc_data or not _pid_is_alive(proc_data.get("pid")):
            message = "SCM process ended while its durable monitor was unavailable; exit status could not be recovered"
            update_broker_run(
                str(record.get("case") or case),
                state="error",
                returncode=None,
                recovery_error=message,
                finished_at=time.time(),
            )
            if job_id:
                _JOB_STORE.update(
                    job_id,
                    state="error",
                    error={"code": "SCM_RECOVERY_EXIT_UNKNOWN", "message": message},
                    finished_at_unix_seconds=time.time(),
                )
                _ARTIFACT_STORE.release(job_id)
            continue
        _background(_watch_run, str(record.get("case") or case), proc_data, job_id or None)
        recovered.append({"case": str(record.get("case") or case), "pid": proc_data.get("pid"), "job_id": job_id or None})
    if recovered:
        publish_event("run", "Recovered SCM job monitoring", ", ".join(item["case"] for item in recovered), status="info")
    return recovered


def recover_queued_scm_batches() -> list[dict[str, Any]]:
    """Resume broker-owned batch queue advancement after a broker restart."""
    recovered: list[dict[str, Any]] = []
    for batch in _JOB_STORE.list_kind("scm_batch"):
        if str(batch.get("state") or "") not in {"queued", "running"}:
            continue
        if not any(str(child.get("state") or "") == "queued" for child in batch.get("children") or []):
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


def cancel_job(job_id: str) -> dict[str, Any]:
    """Request cancellation through the same safe lifecycle used by Dash."""
    record = get_job(job_id)
    state = str(record.get("state") or "")
    if state not in {"submitting", "running", "stopping"}:
        return {"status": "already_terminal", **record}
    kind = str(record.get("kind") or "")
    if kind == "compile":
        result = stop_compile()
    elif kind == "scm":
        result = stop_run(str((record.get("request") or {}).get("case") or ""))
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


def read_job_log(job_id: str, cursor: int = 0, max_bytes: int = 8192) -> dict[str, Any]:
    record = get_job(job_id)
    runtime = dict(record.get("runtime") or {})
    run = dict(runtime.get("run") or {})
    path = Path(str(runtime.get("log") or runtime.get("log_path") or run.get("log") or run.get("log_path") or ""))
    if not path.is_file():
        return {"job_id": job_id, "cursor": max(0, int(cursor)), "next_cursor": max(0, int(cursor)), "text": ""}
    offset = max(0, int(cursor))
    limit = min(max(1, int(max_bytes)), 16384)
    with path.open("rb") as handle:
        handle.seek(offset)
        raw = handle.read(limit)
    return {"job_id": job_id, "cursor": offset, "next_cursor": offset + len(raw), "text": raw.decode("utf-8", errors="replace")}


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
        action={"type": "compile", "tab": "compile", "operation": "start", "job": dict(job)},
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
    job_id: str | None = None,
    run_id: str | None = None,
) -> dict[str, Any]:
    """Run one checked-in SCM case with limited, structured overrides."""
    case_name = _validate_case(case)
    stats_name = _validated_stats_file(stats_file)
    config_name = _valid_tune_config(config)
    override_args = _override_args(overrides)
    normalized_options = _normalize_dashboard_cli_options(cli_options)
    extra_args = list(normalized_options.get("extra_args") or [])
    extra_args.extend(override_args)
    if extra_args:
        normalized_options["extra_args"] = extra_args
    if output_dir:
        isolated = Path(output_dir).resolve()
        artifact_root = _ARTIFACT_STORE.root.resolve()
        public_root = (REPO_ROOT / "output").resolve()
        if artifact_root not in isolated.parents and public_root not in isolated.parents:
            raise ValueError("MCP output directory must be in the repository output directory")
        isolated.mkdir(parents=True, exist_ok=True)
        normalized_options["out_dir"] = str(isolated)
    return launch_dashboard_scm_request(
        case_name,
        stats_name,
        config_name,
        {"flags": {}, "tunable": {}, "silhs": {}},
        normalized_options,
        job_id=job_id,
        run_id=run_id,
        check_admission=True,
    )


def launch_dashboard_scm_request(
    case: str,
    stats_file: str,
    config: str,
    overrides: dict[str, Any] | None,
    cli_options: dict[str, Any] | None,
    *,
    job_id: str | None = None,
    run_id: str | None = None,
    check_admission: bool = True,
) -> dict[str, Any]:
    """Launch one Run-tab request through the durable broker lifecycle.

    The entry point is internal to Dash, retaining its richer control surface
    while making process ownership, cancellation, recovery, and admission the
    same as for typed MCP SCM jobs.
    """
    case_name = _validate_case(case)
    stats_name = _validated_stats_file(stats_file)
    config_name = _valid_tune_config(config)
    normalized_overrides = _normalize_dashboard_overrides(overrides)
    normalized_options = _normalize_dashboard_cli_options(cli_options)
    with exclusive_file_lock(private_path(REPO_ROOT, "scm-submission.lock")):
        if check_admission:
            _assert_scm_admission(case_name)
        detail = " ".join(normalized_options.get("extra_args") or []) or "dashboard configuration"
        publish_event("run", f"Starting {case_name}", detail, status="running")
        proc_data = start_case_process(case_name, stats_name, normalized_overrides, normalized_options, config_name)
        if normalized_options.get("out_dir"):
            proc_data["output_directory"] = str(normalized_options["out_dir"])
        set_broker_run(
            case_name,
            {
                "state": "running",
                "case": case_name,
                "proc_data": dict(proc_data),
                "stats_file": stats_name,
                "config": config_name,
                "cli_options": dict(normalized_options),
                "log_tail": "",
                "log_offset": 0,
                "returncode": None,
                "broker_managed": True,
                "job_id": job_id,
                "run_id": run_id,
            },
        )
        publish_run_request(case_name, proc_data, stats_file=stats_name, config=config_name, cli_options=normalized_options)
        _background(_watch_run, case_name, proc_data, job_id)
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


def stop_run(case: str) -> dict[str, Any]:
    """Ask the broker-owned SCM case process to stop without requiring Dash state."""
    case_name = _validate_case(case)
    record = dict((broker_jobs().get("runs") or {}).get(case_name) or {})
    proc_data = dict(record.get("proc_data") or {})
    if str(record.get("state") or "") not in {"running", "stopping"} or not proc_data:
        raise ValueError(f"no broker-owned {case_name} run is running")
    proc = get_proc(proc_data.get("pid"))
    try:
        if proc is not None and proc.poll() is None:
            proc.terminate()
        elif proc is None:
            os.kill(int(proc_data["pid"]), signal.SIGTERM)
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise ValueError(f"{case_name} process is no longer available to stop") from exc
    update_broker_run(case_name, state="stopping")
    publish_event("run", f"Stop requested for {case_name}", str(proc_data.get("log") or ""), status="info")
    return {"status": "stop_requested", "case": case_name, "pid": proc_data.get("pid")}


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

    for case, raw_record in (snapshot.get("runs") or {}).items():
        record = dict(raw_record or {})
        if str(record.get("state") or "") not in {"running", "stopping"}:
            continue
        try:
            stop_run(str(case))
            requested.append(f"SCM:{case}")
            job_id = str(record.get("job_id") or "")
            if job_id:
                _JOB_STORE.update(job_id, state="stopping")
        except ValueError as exc:
            errors.append(f"SCM:{case}: {exc}")

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

    # A batch may have children that have not started yet and therefore do not
    # appear in broker_jobs(). Mark them terminal so a later broker cannot
    # recover and launch them after the manager has deliberately stopped.
    for batch in _JOB_STORE.list_kind("scm_batch"):
        if str(batch.get("state") or "") not in {"queued", "running"}:
            continue
        children: list[dict[str, Any]] = []
        for raw_child in batch.get("children") or []:
            child = dict(raw_child or {})
            if str(child.get("state") or "") == "queued":
                child["state"] = "cancelled"
                child_id = str(child.get("job_id") or "")
                if child_id:
                    _JOB_STORE.update(child_id, state="cancelled")
                    _ARTIFACT_STORE.release(child_id)
            children.append(child)
        updated = _JOB_STORE.update(
            str(batch.get("job_id") or ""),
            state="cancelled",
            children=children,
            finished_at_unix_seconds=time.time(),
            cancellation_reason=str(reason)[:500],
        )
        if updated:
            try:
                _write_scm_batch_manifest(updated)
            except OSError as exc:
                errors.append(f"SCM batch manifest: {exc}")
        requested.append(f"SCM-batch:{batch.get('batch_id') or batch.get('job_id')}")

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
        batch_size=request.get("batch_size", 8),
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
        "batch_count": run_data.get("batch_count"),
        "batch_size": run_data.get("batch_size"),
        "run": run_data,
    }


def plot_profiles(
    case: str,
    variables: list[str],
    time_seconds: float | None = None,
    *,
    run_id: str | None = None,
    output_dir: str | None = None,
    time_start_seconds: float | None = None,
    average_minutes: float | None = None,
    window_preset: str | None = None,
    benchmark_sources: list[str] | None = None,
) -> dict[str, Any]:
    """Open profile variables with optional validated benchmark overlays."""
    case_name = _validate_case(case)
    selected_run_id, output_dirs = _plot_output_selection(case_name, run_id, output_dir)
    case_data, selection = _profile_selection(
        case_name,
        time_seconds=time_seconds,
        time_start_seconds=time_start_seconds,
        average_minutes=average_minutes,
        window_preset=window_preset,
        require_case_data=False,
        output_dirs=output_dirs,
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
        output_dir=output_dirs[0] if output_dirs else "output",
        **selection,
    )
    return {
        "status": "requested",
        "activity_id": event["id"],
        "case": case_name,
        "variables": names,
        "run_id": selected_run_id or None,
        "output_directory": output_dirs[0] if output_dirs else "output",
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
    allowed = {"tutorial", "compile", "run", "tune", "plots", "reports", "misc"}
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
    runs = broker_jobs().get("runs") or {}
    if case is None or not str(case).strip():
        return {"runs": runs}
    case_name = _validate_case(str(case))
    return {"case": case_name, "run": runs.get(case_name)}


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
        "Stop one broker-owned SCM case.",
        {"type": "object", "required": ["case"], "properties": {"case": {"type": "string"}}},
        lambda arguments: stop_run(arguments.get("case", "")),
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
        return submit_scm_run(ScmRunRequest.model_validate(payload.get("request") or {}))
    if action == "domain_submit_scm_batch":
        return submit_scm_batch(
            ScmRunBatchRequest.model_validate(payload.get("request") or {}),
            native_overrides=payload.get("native_overrides") if "native_overrides" in payload else None,
            native_cli_options=payload.get("native_cli_options") if "native_cli_options" in payload else None,
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
    if action == "launch_rebuild_request":
        return launch_rebuild_request(
            list(payload.get("builds") or []),
            dict(payload.get("discovery") or {}),
            str(payload.get("label") or "selected builds"),
        )
    if action == "launch_dashboard_scm_request":
        return launch_dashboard_scm_request(
            payload.get("case", ""),
            payload.get("stats_file", DEFAULT_STATS_NAME),
            payload.get("config", "default"),
            dict(payload.get("overrides") or {}),
            dict(payload.get("cli_options") or {}),
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
        return stop_run(payload.get("case", ""))
    if action == "inspect_tuning":
        return inspect_tuning()
    if action == "run_tuning_loss":
        return run_tuning_loss(payload.get("mode", "window"), payload.get("max_results", 16))
    raise ValueError("unknown action; allowed actions are " + ", ".join(allowed_action_names()))
