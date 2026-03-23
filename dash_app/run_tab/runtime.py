"""Process launching, status tracking, and log streaming helpers for the run tab."""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

from .namelist import cleanup_temp_files, write_temp_namelist
from .state import (
    CUDA_MPS_LOG_DIR,
    CUDA_MPS_PIPE_DIR,
    DEFAULT_STATS_NAME,
    ENABLE_CUDA_MPS,
    FLAGS_FILE,
    ITERATION_RE,
    MAX_RUN_PROCS,
    MAX_UI_LOG_LINES,
    MPS_LOCK,
    MPS_STATE,
    NO_STATS_NAME,
    REPO_ROOT,
    RUN_ACTIVE_CASES,
    RUN_FINALIZED,
    RUN_LOCK,
    RUN_PROCS,
    RUN_STATUS,
    SILHS_FILE,
    STATS_DIR,
    TUNABLE_FILE,
    set_child_stack_limit,
)


def ensure_cuda_mps():
    """Start CUDA MPS once and return child env overrides if successful."""
    if not ENABLE_CUDA_MPS:
        return {}
    with MPS_LOCK:
        if MPS_STATE["attempted"]:
            return dict(MPS_STATE["env"])
        MPS_STATE["attempted"] = True
        mps_ctl = shutil.which("nvidia-cuda-mps-control")
        if not mps_ctl:
            MPS_STATE["message"] = "nvidia-cuda-mps-control not found; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}
        try:
            os.makedirs(CUDA_MPS_PIPE_DIR, exist_ok=True)
            os.makedirs(CUDA_MPS_LOG_DIR, exist_ok=True)
        except OSError as exc:
            MPS_STATE["message"] = f"failed to create MPS directories: {exc}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        mps_env = os.environ.copy()
        mps_env["CUDA_MPS_PIPE_DIRECTORY"] = CUDA_MPS_PIPE_DIR
        mps_env["CUDA_MPS_LOG_DIRECTORY"] = CUDA_MPS_LOG_DIR
        try:
            start = subprocess.run(
                [mps_ctl, "-d"],
                env=mps_env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        except Exception as exc:
            MPS_STATE["message"] = f"failed to launch MPS daemon: {exc}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        if start.returncode != 0:
            err = (start.stderr or "").strip()
            suffix = f" ({err})" if err else ""
            MPS_STATE["message"] = f"MPS daemon start failed with code {start.returncode}{suffix}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        MPS_STATE["enabled"] = True
        MPS_STATE["env"] = {
            "CUDA_MPS_PIPE_DIRECTORY": CUDA_MPS_PIPE_DIR,
            "CUDA_MPS_LOG_DIRECTORY": CUDA_MPS_LOG_DIR,
        }
        MPS_STATE["message"] = "CUDA MPS enabled for dashboard-launched runs."
        print(f"[run-tab] {MPS_STATE['message']}")
        return dict(MPS_STATE["env"])


def run_child_env():
    """Build the environment inherited by dashboard-launched SCM children."""
    child_env = os.environ.copy()
    child_env.update(ensure_cuda_mps())
    return child_env


def read_log_increment(log_path, offset):
    """Read newly appended text from a log file starting at a byte offset."""
    if not log_path or not os.path.exists(log_path):
        return "", 0
    try:
        with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
            handle.seek(0, os.SEEK_END)
            end_pos = handle.tell()
            if offset < 0 or offset > end_pos:
                offset = 0
            if offset == end_pos:
                return "", end_pos
            handle.seek(offset)
            return handle.read(), handle.tell()
    except Exception:
        return "", offset


def cleanup_log_file(log_path):
    """Delete a temporary run log if it still exists."""
    if not log_path:
        return
    try:
        if os.path.exists(log_path):
            os.remove(log_path)
    except Exception:
        # Log cleanup should never make the UI brittle.
        pass


def append_log_tail(existing, chunk, max_lines=MAX_UI_LOG_LINES):
    """Append new log text and trim the retained UI buffer to the last N lines."""
    text = f"{existing or ''}{chunk or ''}"
    if not text:
        return ""
    lines = text.splitlines(keepends=True)
    if len(lines) <= max_lines:
        return text
    return "".join(lines[-max_lines:])


def mark_case_started(case_name, proc, proc_data):
    """Record a newly launched process in the shared run-tab state."""
    with RUN_LOCK:
        RUN_PROCS[proc.pid] = proc
        RUN_ACTIVE_CASES[case_name] = dict(proc_data)
        RUN_STATUS.pop(case_name, None)
        RUN_FINALIZED.discard(case_name)


def get_proc(pid):
    """Return the tracked subprocess for a pid, if it still exists."""
    if not pid:
        return None
    with RUN_LOCK:
        return RUN_PROCS.get(pid)


def record_case_finish(case_name, pid, status):
    """Finalize one case status and return whether it had already been finalized."""
    with RUN_LOCK:
        already_finalized = case_name in RUN_FINALIZED
        if pid:
            RUN_PROCS.pop(pid, None)
        RUN_ACTIVE_CASES.pop(case_name, None)
        RUN_STATUS[case_name] = int(status)
        RUN_FINALIZED.add(case_name)
    return already_finalized


def get_cached_status(case_name):
    """Return the cached terminal status for a case, if any."""
    with RUN_LOCK:
        return RUN_STATUS.get(case_name)


def clear_case_status(case_name):
    """Clear cached terminal status for a case."""
    with RUN_LOCK:
        RUN_STATUS.pop(case_name, None)
        RUN_FINALIZED.discard(case_name)


def snapshot_status_lists():
    """Return sorted completed and failed case-name lists from cached status."""
    with RUN_LOCK:
        completed = sorted([name for name, status in RUN_STATUS.items() if status == 0])
        failed = sorted([name for name, status in RUN_STATUS.items() if status != 0])
    return completed, failed


def is_case_active(case_name):
    """Return whether a case is currently tracked as running."""
    with RUN_LOCK:
        return case_name in RUN_ACTIVE_CASES


def snapshot_active_cases():
    """Return a shallow copy of the active-case metadata map."""
    with RUN_LOCK:
        return {name: dict(data) for name, data in RUN_ACTIVE_CASES.items()}


def clean_cli_option(value):
    """Normalize optional CLI values by stripping whitespace and None."""
    if value is None:
        return ""
    return str(value).strip()


def normalize_task_limit(value):
    """Parse a user-specified task limit, falling back to the dashboard default."""
    cleaned = clean_cli_option(value)
    if not cleaned:
        return MAX_RUN_PROCS
    try:
        return max(1, int(cleaned))
    except (TypeError, ValueError):
        return MAX_RUN_PROCS


def build_case_command(case_name, stats_name, cli_options=None):
    """Build the exact run_scm.py command shown in the UI copy button."""
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join("input", "stats", stats_value)
    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg]
    cli_options = cli_options or {}
    for flag, key in (("-ngrdcol", "ngrdcol"), ("-hr", "hr"), ("-max_iters", "max_iters"), ("-debug", "debug"), ("-dt_main", "dt_main"), ("-dt_rad", "dt_rad"), ("-tout", "tout"), ("-out_dir", "out_dir")):
        value = clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    cmd.append(case_name)
    return " ".join(shlex.quote(str(part)) for part in cmd)


def start_case_process(case_name, stats_name, overrides, cli_options=None):
    """Launch one SCM case with temporary override files and return runtime metadata."""
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join(STATS_DIR, stats_value)

    params_path = write_temp_namelist(TUNABLE_FILE, overrides.get("tunable"), "clubb_params_")
    flags_path = write_temp_namelist(FLAGS_FILE, overrides.get("flags"), "clubb_flags_")
    silhs_path = write_temp_namelist(SILHS_FILE, overrides.get("silhs"), "clubb_silhs_")

    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg]
    cli_options = cli_options or {}
    for flag, key in (("-ngrdcol", "ngrdcol"), ("-hr", "hr"), ("-max_iters", "max_iters"), ("-debug", "debug"), ("-dt_main", "dt_main"), ("-dt_rad", "dt_rad"), ("-tout", "tout"), ("-out_dir", "out_dir")):
        value = clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    if params_path:
        cmd.extend(["-params", params_path])
    if flags_path:
        cmd.extend(["-flags", flags_path])
    if silhs_path:
        cmd.extend(["-silhs_params", silhs_path])
    cmd.append(case_name)

    log_file = tempfile.NamedTemporaryFile(delete=False, prefix="clubb_run_", suffix=".log", dir="/tmp")
    log_path = log_file.name
    proc = subprocess.Popen(
        cmd,
        cwd=REPO_ROOT,
        env=run_child_env(),
        stdout=log_file,
        stderr=subprocess.STDOUT,
        text=True,
        preexec_fn=set_child_stack_limit,
    )
    log_file.close()

    temp_files = [path for path in (params_path, flags_path, silhs_path) if path]
    proc_data = {
        "pid": proc.pid,
        "log": log_path,
        "start_time": time.time(),
        "temp_files": temp_files,
    }
    mark_case_started(case_name, proc, proc_data)
    return proc_data


def launch_from_queue(running, queued, logs, max_run_procs=None):
    """Launch queued cases until the run concurrency limit is reached."""
    queue = list(queued or [])
    launched = False
    limit = normalize_task_limit(max_run_procs)
    while queue and len(running) < limit:
        item = queue.pop(0)
        case_name = item.get("case")
        stats_name = item.get("stats") or DEFAULT_STATS_NAME
        overrides = item.get("overrides") or {"flags": {}, "tunable": {}, "silhs": {}}
        cli_options = item.get("cli_options") or {}
        if not case_name or case_name in running or is_case_active(case_name) or get_cached_status(case_name) is not None:
            continue
        running[case_name] = start_case_process(case_name, stats_name, overrides, cli_options)
        logs[case_name] = f"--- Running {case_name} ({stats_name}) ---\n"
        launched = True
    return queue, launched


def extract_progress(log_text):
    """Extract fractional iteration progress from run_scm.py log output."""
    if not log_text:
        return None
    last_match = None
    for match in ITERATION_RE.finditer(log_text):
        last_match = match
    if not last_match:
        return None
    try:
        current = int(last_match.group(1))
        total = int(last_match.group(2))
    except (TypeError, ValueError):
        return None
    if total <= 0:
        return None
    return max(0.0, min(1.0, current / total))


def format_runtime(seconds):
    """Format elapsed runtime for display in the console headers."""
    if seconds is None:
        return ""
    total = max(0, int(round(float(seconds))))
    hours, rem = divmod(total, 3600)
    mins, secs = divmod(rem, 60)
    if hours > 0:
        return f"{hours}h {mins:02d}m {secs:02d}s"
    if mins > 0:
        return f"{mins}m {secs:02d}s"
    return f"{secs}s"


def format_eta(seconds):
    """Format remaining runtime as m:ss for progress display."""
    if seconds is None:
        return ""
    total = max(0, int(round(float(seconds))))
    mins, secs = divmod(total, 60)
    return f"{mins}m {secs:02d}s"
