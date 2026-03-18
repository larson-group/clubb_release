"""Shared constants and mutable runtime state for the run tab."""

import os
import re
import resource
import threading

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CASE_DIR = os.path.join(REPO_ROOT, "input", "case_setups")
STATS_DIR = os.path.join(REPO_ROOT, "input", "stats")
DEFAULT_STATS_NAME = "standard_stats.in"
NO_STATS_NAME = "none"
FLAGS_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "configurable_model_flags.in")
TUNABLE_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "tunable_parameters.in")
SILHS_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "silhs_parameters.in")
RUN_SCM_ALL = os.path.join(REPO_ROOT, "run_scripts", "run_scm_all.py")
RUN_PROCS = {}
RUN_ACTIVE_CASES = {}
RUN_STATUS = {}
RUN_FINALIZED = set()
RUN_LOCK = threading.Lock()
RUN_STREAM_LOCK = threading.Lock()
ITERATION_RE = re.compile(
    r"iteration:\s*(\d+)\s*/\s*(\d+)\s*--\s*time\s*=\s*([-+0-9.eE]+)\s*/\s*([-+0-9.eE]+)",
    re.IGNORECASE,
)
MAX_UI_LOG_LINES = 2000
RUN_STACK_KB = int(os.environ.get("CLUBB_RUN_STACK_KB", "8192000"))
ENABLE_CUDA_MPS = os.environ.get("CLUBB_ENABLE_CUDA_MPS", "1").strip().lower() not in {"0", "false", "no"}
_uid = os.getuid() if hasattr(os, "getuid") else "user"
CUDA_MPS_PIPE_DIR = os.path.join("/tmp", f"clubb_mps_pipe_{_uid}")
CUDA_MPS_LOG_DIR = os.path.join("/tmp", f"clubb_mps_log_{_uid}")
MPS_LOCK = threading.Lock()
MPS_STATE = {
    "attempted": False,
    "enabled": False,
    "message": "",
    "env": {},
}


def physical_core_count():
    """Return the best-effort physical core count to avoid oversubscribing runs."""
    if os.name == "posix" and os.path.exists("/proc/cpuinfo"):
        try:
            pairs = set()
            phys_id = None
            core_id = None
            with open("/proc/cpuinfo", "r", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line:
                        if phys_id is not None and core_id is not None:
                            pairs.add((phys_id, core_id))
                        phys_id = None
                        core_id = None
                        continue
                    if line.startswith("physical id"):
                        phys_id = line.split(":", 1)[1].strip()
                    elif line.startswith("core id"):
                        core_id = line.split(":", 1)[1].strip()
            if phys_id is not None and core_id is not None:
                pairs.add((phys_id, core_id))
            if pairs:
                return len(pairs)
        except Exception:
            # Fall back to logical cores if physical core probing fails.
            pass

    if os.name == "posix":
        try:
            import subprocess

            out = subprocess.check_output(["sysctl", "-n", "hw.physicalcpu"], text=True).strip()
            cores = int(out)
            if cores > 0:
                return cores
        except Exception:
            # sysctl is only available on some platforms; ignore failures.
            pass

    return max(1, os.cpu_count() or 1)


MAX_RUN_PROCS = max(1, physical_core_count() - 2)


def set_child_stack_limit():
    """Raise child stack size before exec when the host permits it."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        target = max(soft, RUN_STACK_KB * 1024)
        if hard != resource.RLIM_INFINITY:
            target = min(target, hard)
        resource.setrlimit(resource.RLIMIT_STACK, (target, hard))
    except Exception:
        # Keep launches robust if changing the stack limit is not allowed.
        pass
