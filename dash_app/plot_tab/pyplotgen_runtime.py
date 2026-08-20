"""Launch and track Plot-tab PyPlotGen exports."""

from __future__ import annotations

import os
import re
import shlex
import signal
import subprocess
import sys
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
PYPLOTGEN_ROOT = REPO_ROOT / "postprocessing" / "pyplotgen"
PYPLOTGEN_SCRIPT = PYPLOTGEN_ROOT / "pyplotgen.py"
PYPLOTGEN_OUTPUT_ROOT = REPO_ROOT / "output" / "pyplotgen"
_PROCESSES: dict[int, subprocess.Popen[bytes]] = {}
_PROGRESS_RE = re.compile(rb"Progress:\s*(\d+)\s+of\s+(\d+)\s+total")


def _selected_output_dirs(values: list[str]) -> list[str]:
    selected = []
    seen = set()
    for value in values or []:
        path = Path(str(value or "").strip()).expanduser().resolve()
        if not path.is_dir():
            raise ValueError(f"PyPlotGen input directory does not exist: {path}")
        text = str(path)
        if text not in seen:
            selected.append(text)
            seen.add(text)
    if not selected:
        raise ValueError("select at least one active output before generating PyPlotGen plots")
    return selected


def read_pyplotgen_progress(log_path: str, max_bytes: int = 32768) -> tuple[int, int] | None:
    """Read the most recent PyPlotGen panel counter from a bounded log tail."""
    try:
        path = Path(log_path)
        size = path.stat().st_size
        with path.open("rb") as handle:
            handle.seek(max(0, size - max_bytes))
            matches = list(_PROGRESS_RE.finditer(handle.read()))
    except OSError:
        return None
    if not matches:
        return None
    return int(matches[-1].group(1)), int(matches[-1].group(2))


def start_pyplotgen(output_dirs: list[str]) -> tuple[dict[str, Any], subprocess.Popen[bytes]]:
    """Start one export in a fixed, unique directory below output/pyplotgen."""
    selected = _selected_output_dirs(output_dirs)
    run_id = uuid.uuid4().hex
    folder = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + run_id[:8]
    output_dir = PYPLOTGEN_OUTPUT_ROOT / folder
    log_dir = PYPLOTGEN_OUTPUT_ROOT / ".logs"
    matplotlib_dir = PYPLOTGEN_OUTPUT_ROOT / ".matplotlib"
    log_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{run_id}.log"
    html_path = output_dir / "index.html"
    command = [sys.executable, "-u", str(PYPLOTGEN_SCRIPT), "-l", "-c", *selected, "-o", str(output_dir)]

    log_handle = log_path.open("wb")
    try:
        process = subprocess.Popen(
            command,
            cwd=str(PYPLOTGEN_ROOT),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env={
                **os.environ,
                "MPLBACKEND": "Agg",
                "MPLCONFIGDIR": str(matplotlib_dir),
                "PYTHONUNBUFFERED": "1",
            },
        )
    except Exception:
        log_handle.close()
        raise
    process._clubb_log_handle = log_handle  # type: ignore[attr-defined]
    _PROCESSES[process.pid] = process
    job = {
        "state": "running",
        "run_id": run_id,
        "pid": process.pid,
        "output_dirs": selected,
        "output_dir": str(output_dir),
        "html_path": str(html_path),
        "html_url": f"/_clubb-pyplotgen/{folder}/index.html",
        "log_path": str(log_path),
        "command": command,
        "command_display": shlex.join(command),
        "started_at": datetime.now().timestamp(),
        "returncode": None,
        "broker_managed": True,
    }
    return job, process


def release_pyplotgen(process: subprocess.Popen[bytes]) -> None:
    """Release the broker's process and log handles after completion."""
    _PROCESSES.pop(process.pid, None)
    handle = getattr(process, "_clubb_log_handle", None)
    if handle is not None:
        handle.close()


def stop_pyplotgen(pid: Any) -> None:
    """Stop one broker-owned PyPlotGen process group."""
    process_id = int(pid)
    try:
        os.killpg(process_id, signal.SIGTERM)
    except OSError:
        os.kill(process_id, signal.SIGTERM)
