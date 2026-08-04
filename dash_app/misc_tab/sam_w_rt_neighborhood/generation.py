"""Small asynchronous job runner for interactive Misc-laboratory generation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
MAX_LOG_CHARS = 4_000


@dataclass
class GenerationJob:
    process: subprocess.Popen
    log_path: Path
    started_at: float
    label: str


_LOCK = threading.Lock()
_JOBS: dict[str, GenerationJob] = {}


def _tail(path: Path) -> str:
    try:
        with path.open("rb") as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - MAX_LOG_CHARS))
            return stream.read().decode("utf-8", errors="replace").strip()
    except OSError:
        return ""


def _start(kind: str, label: str, command: list[str]) -> str:
    with _LOCK:
        existing = _JOBS.get(kind)
        if existing is not None and existing.process.poll() is None:
            return f"{existing.label} is already running."
        log_file = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            delete=False,
            prefix=f"clubb_misc_{kind}_",
            suffix=".log",
            dir="/tmp",
        )
        process = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        log_file.close()
        _JOBS[kind] = GenerationJob(process, Path(log_file.name), time.time(), label)
    return f"Started {label}. This may take several minutes."


def start_atlas_generation() -> str:
    return _start(
        "sam_wrt_atlas",
        "SAM w–rₜ atlas generation",
        [
            sys.executable,
            "-u",
            "-m",
            "dash_app.misc_tab.sam_w_rt_neighborhood.atlas",
            "--force",
        ],
    )


def job_status(kind: str) -> tuple[bool, bool, str]:
    """Return ``running, succeeded, display_message`` for one job."""

    with _LOCK:
        job = _JOBS.get(kind)
    if job is None:
        return False, False, ""
    returncode = job.process.poll()
    log = _tail(job.log_path)
    if returncode is None:
        suffix = f"\n\n{log}" if log else ""
        return True, False, f"{job.label} is running…{suffix}"
    if returncode == 0:
        suffix = f"\n\n{log}" if log else ""
        return False, True, f"{job.label} finished successfully.{suffix}"
    suffix = f"\n\n{log}" if log else ""
    return False, False, f"{job.label} failed (exit {returncode}).{suffix}"
