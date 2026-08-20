"""Process launching helpers for broker-owned Run-tab SCM jobs."""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

from .namelist import write_temp_namelist
from .state import (
    CUDA_MPS_LOG_DIR,
    CUDA_MPS_PIPE_DIR,
    DEFAULT_STATS_NAME,
    ENABLE_CUDA_MPS,
    MAX_RUN_PROCS,
    MPS_LOCK,
    MPS_STATE,
    NO_STATS_NAME,
    REPO_ROOT,
    RUN_LOCK,
    RUN_PROCS,
    STATS_DIR,
    set_child_stack_limit,
)
from dash_app.shared.tunable_configs import tunable_config_file


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
    preserve_ld_path = child_env.get("CLUBB_DASH_PRESERVE_LD_LIBRARY_PATH", "").strip().lower()
    if preserve_ld_path not in {"1", "true", "yes", "on"}:
        child_env.pop("LD_LIBRARY_PATH", None)
    child_env.update(ensure_cuda_mps())
    return child_env


def mark_case_started(proc):
    """Retain the one Popen handle needed to poll and signal a broker child."""
    with RUN_LOCK:
        RUN_PROCS[proc.pid] = proc


def get_proc(pid):
    """Return the tracked subprocess for a pid, if it still exists."""
    if not pid:
        return None
    with RUN_LOCK:
        return RUN_PROCS.get(pid)


def record_case_finish(case_name, pid, status):
    """Release a terminal broker child's Popen handle."""
    with RUN_LOCK:
        if pid:
            RUN_PROCS.pop(pid, None)


def clean_cli_option(value):
    """Normalize optional CLI values by stripping whitespace and None."""
    if value is None:
        return ""
    return str(value).strip()


def split_extra_cli_args(value):
    """Split free-form run_scm.py arguments using shell-like quoting rules."""
    cleaned = clean_cli_option(value)
    if not cleaned:
        return []
    return shlex.split(cleaned)


def extra_cli_args(cli_options):
    """Return normalized free-form CLI tokens stored in run-tab options."""
    value = (cli_options or {}).get("extra_args")
    if isinstance(value, (list, tuple)):
        return [cleaned for cleaned in (clean_cli_option(item) for item in value) if cleaned]
    return split_extra_cli_args(value)


def normalize_task_limit(value):
    """Parse a user-specified task limit, falling back to the dashboard default."""
    cleaned = clean_cli_option(value)
    if not cleaned:
        return MAX_RUN_PROCS
    try:
        numeric = float(cleaned)
    except (TypeError, ValueError):
        return MAX_RUN_PROCS
    return int(numeric) if numeric >= 1 and numeric.is_integer() else MAX_RUN_PROCS


def append_launch_target(cmd, cli_options):
    """Append the explicit implementation and supporting install snapshot."""
    implementation = clean_cli_option((cli_options or {}).get("implementation")).lower()
    if implementation == "python":
        cmd.append("-python")
    elif implementation == "jax":
        cmd.append("-jax")
    install_dir = clean_cli_option((cli_options or {}).get("install_dir"))
    if install_dir:
        cmd.extend(["-install_dir", install_dir])


def build_case_command(case_name, stats_name, cli_options=None, config_name=None):
    """Build the exact run_scm.py command shown in the UI copy button."""
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join("input", "stats", stats_value)
    config_value = clean_cli_option(config_name) or "default"
    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg, "-config", config_value]
    cli_options = cli_options or {}
    append_launch_target(cmd, cli_options)
    for flag, key in (("-multicol", "multicol"), ("-batch_size", "batch_size"), ("-max_iters", "max_iters"), ("-debug", "debug"), ("-dt_main", "dt_main"), ("-dt_rad", "dt_rad"), ("-tout", "tout"), ("-out_dir", "out_dir")):
        value = clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    cmd.extend(extra_cli_args(cli_options))
    cmd.append(case_name)
    return " ".join(shlex.quote(str(part)) for part in cmd)


def start_case_process(case_name, stats_name, overrides, cli_options=None, config_name=None):
    """Launch one SCM case with temporary override files and return runtime metadata."""
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join(STATS_DIR, stats_value)

    config_value = clean_cli_option(config_name) or "default"
    params_path = write_temp_namelist(
        tunable_config_file(config_value, "tunable_parameters.in"),
        overrides.get("tunable"),
        "clubb_params_",
    )
    flags_path = write_temp_namelist(
        tunable_config_file(config_value, "configurable_model_flags.in"),
        overrides.get("flags"),
        "clubb_flags_",
    )
    silhs_path = write_temp_namelist(
        tunable_config_file(config_value, "silhs_parameters.in"),
        overrides.get("silhs"),
        "clubb_silhs_",
    )

    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg, "-config", config_value]
    cli_options = cli_options or {}
    append_launch_target(cmd, cli_options)
    for flag, key in (("-multicol", "multicol"), ("-batch_size", "batch_size"), ("-max_iters", "max_iters"), ("-debug", "debug"), ("-dt_main", "dt_main"), ("-dt_rad", "dt_rad"), ("-tout", "tout"), ("-out_dir", "out_dir")):
        value = clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    if params_path:
        cmd.extend(["-params", params_path])
    if flags_path:
        cmd.extend(["-flags", flags_path])
    if silhs_path:
        cmd.extend(["-silhs_params", silhs_path])
    cmd.extend(extra_cli_args(cli_options))
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
        start_new_session=True,
    )
    log_file.close()

    temp_files = [path for path in (params_path, flags_path, silhs_path) if path]
    proc_data = {
        "pid": proc.pid,
        "log": log_path,
        "start_time": time.time(),
        "temp_files": temp_files,
        "config": config_value,
    }
    mark_case_started(proc)
    return proc_data
