"""Prepare help-driven options and launch Plot-tab PyPlotGen galleries.

pyplotgen_help() runs the real CLI with -h, returning its verbatim text and
selectable option/value hints. pyplotgen_arguments() turns selected options
and quoted value text into literal argument tokens. start_pyplotgen() starts
the broker-owned export in a unique output/pyplotgen directory; the dashboard
supplies its active folders and export destination. No shell is used.

read_pyplotgen_progress() reads the latest panel counter from the export log.
stop_pyplotgen() signals the process group; release_pyplotgen() closes the log
handle after it finishes. Help is cached until the CLI script changes, with
Matplotlib's cache kept under the same ignored output root as export state.
"""

from __future__ import annotations

import os
import re
import shlex
import signal
import subprocess
import sys
import uuid
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
PYPLOTGEN_ROOT = REPO_ROOT / "pyplotgen"
PYPLOTGEN_SCRIPT = PYPLOTGEN_ROOT / "pyplotgen.py"
PYPLOTGEN_OUTPUT_ROOT = REPO_ROOT / "output" / "pyplotgen"
_PROGRESS_RE = re.compile(rb"Progress:\s*(\d+)\s+of\s+(\d+)\s+total")
_MANAGED_OPTIONS = {"-h", "-help", "-clubb", "-output_dir"}


def _parse_help_options(text):
    """Read option names, value hints and descriptions from argparse's help section."""
    options = []
    current = None
    for line in text.splitlines():
        if re.match(r"^  --?\w", line):
            parts = re.split(r"\s{2,}", line.strip(), maxsplit=1)
            match = re.fullmatch(r"(--?[\w-]+(?:,\s*--?[\w-]+)*)(?:\s+(.*))?", parts[0])
            if match is None:
                raise ValueError(f"Could not read PyPlotGen option: {parts[0]}")
            aliases = [value.strip() for value in match[1].split(',')]
            current = {"flag": aliases[0], "aliases": aliases, "metavar": match[2] or "",
                       "description": parts[1] if len(parts) > 1 else ""}
            options.append(current)
        elif current is not None and line.startswith('    '):
            current['description'] = (current['description'] + ' ' + line.strip()).strip()
    if not options:
        raise ValueError("PyPlotGen -h did not list any options")
    return [option for option in options if not _MANAGED_OPTIONS.intersection(option['aliases'])]


@lru_cache(maxsize=1)
def _read_help(script, modified_ns, size):
    """Cache actual CLI help until the PyPlotGen entry point changes."""
    matplotlib_dir = PYPLOTGEN_OUTPUT_ROOT / '.matplotlib'
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [sys.executable, str(script), '-h'], cwd=str(PYPLOTGEN_ROOT),
        env={**os.environ, 'MPLBACKEND': 'Agg', 'MPLCONFIGDIR': str(matplotlib_dir)},
        capture_output=True, text=True, timeout=30, check=True,
    )
    return {"text": result.stdout, "options": _parse_help_options(result.stdout)}


def pyplotgen_help():
    """Return verbatim PyPlotGen -h output and selectable options derived from that text."""
    stat = PYPLOTGEN_SCRIPT.stat()
    return _read_help(str(PYPLOTGEN_SCRIPT), stat.st_mtime_ns, stat.st_size)


def pyplotgen_arguments(options, option_values=None):
    """Build literal CLI arguments from selected help-listed flags and their value fields."""
    if options is None:
        options = ['-les']
    if not isinstance(options, list) or any(not isinstance(flag, str) for flag in options):
        raise ValueError('PyPlotGen options must be a list of option names')
    values = {} if option_values is None else option_values
    if not isinstance(values, dict):
        raise ValueError('PyPlotGen option values must be an option-to-text mapping')
    catalog = {option['flag']: option for option in pyplotgen_help()['options']}
    arguments = []
    for flag in dict.fromkeys(options):
        if flag not in catalog:
            raise ValueError(f'PyPlotGen option is unavailable or controlled by Dash: {flag}')
        metavar = catalog[flag]['metavar']
        arguments.append(flag)
        if not metavar:
            continue
        tokens = shlex.split(str(values.get(flag) or ''))
        if not tokens and not metavar.startswith('['):
            raise ValueError(f'{flag} requires a value: {metavar}')
        if any(re.match(r'^--?[A-Za-z_]', token) for token in tokens):
            raise ValueError(f'Use the option selector for flags; enter only values for {flag}')
        if '...' not in metavar and len(tokens) > 1:
            raise ValueError(f'{flag} takes one value; quote a path containing spaces')
        arguments.extend(tokens)
    return arguments


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


def start_pyplotgen(output_dirs: list[str], options=None, option_values=None) -> tuple[dict[str, Any], subprocess.Popen[bytes]]:
    """Start one export in a fixed, unique directory below output/pyplotgen."""
    selected = _selected_output_dirs(output_dirs)
    extra_arguments = pyplotgen_arguments(options, option_values)
    run_id = uuid.uuid4().hex
    folder = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + run_id[:8]
    output_dir = PYPLOTGEN_OUTPUT_ROOT / folder
    log_dir = PYPLOTGEN_OUTPUT_ROOT / ".logs"
    matplotlib_dir = PYPLOTGEN_OUTPUT_ROOT / ".matplotlib"
    log_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{run_id}.log"
    html_path = output_dir / "index.html"
    command = [sys.executable, "-u", str(PYPLOTGEN_SCRIPT), *extra_arguments, "-clubb", *selected, "-output_dir", str(output_dir)]

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
