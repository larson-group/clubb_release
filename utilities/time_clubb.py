#!/usr/bin/env python3
"""Benchmark CLUBB with groups of concurrent ``run_scm.py`` processes.

Timing options are parsed here.  Every other option is forwarded verbatim to
``run_scripts/run_scm.py`` so executable, configuration, and namelist choices
stay under the existing runner's control.  For example::

    utilities/time_clubb.py arm -processes 1,4 -batch_sizes 1,8 \
        -max_iters 20 -config my_config -override C2=1.2

``-multicol``, ``-batch_size``, and ``-out_dir`` are reserved because this
utility assigns them independently to every child process.  Each invocation
writes one compact, self-contained profile directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

try:
    from utilities.timing_profiles import (
        GROUP_WALL_TIMER,
        TIMING_FIELDS,
        TIMINGS_FILE,
        create_profile_manifest,
        reserve_profile_directory,
        update_profile_manifest,
        write_batches,
        write_profile_manifest,
        write_profile_readme,
    )
except ModuleNotFoundError:  # Support direct execution from utilities/.
    from timing_profiles import (  # type: ignore[no-redef]
        GROUP_WALL_TIMER,
        TIMING_FIELDS,
        TIMINGS_FILE,
        create_profile_manifest,
        reserve_profile_directory,
        update_profile_manifest,
        write_batches,
        write_profile_manifest,
        write_profile_readme,
    )


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_SCM = REPO_ROOT / "run_scripts" / "run_scm.py"
REQUIRED_TIMER = "advance_clubb_to_end"
_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_CPU_ROW = re.compile(
    rf"^\s*(.*?)\s+(\d+)\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*$"
)
_GPTL_THREAD = re.compile(r"^\s*Stats for thread\s+(\d+)\s*:", re.IGNORECASE)
_GPTL_ROW = re.compile(
    rf"^(?P<marker>[*!]?)(?P<indent>\s*)(?P<name>\S.*?)\s+"
    rf"(?P<calls>\d+)\s+(?:-|\d+)\s+(?P<wall>{_FLOAT})(?:\s|$)"
)


def _restore_sigint_handler() -> None:
    """Restore Python interrupt handling after launch from an ignoring parent."""
    signal.signal(signal.SIGINT, signal.default_int_handler)


class TimerParseError(ValueError):
    """Raised when a CLUBB timing file is missing or unrecognized."""


@dataclass(frozen=True)
class TimerRecord:
    backend: str
    time_basis: str
    thread: int
    timer_name: str
    calls: int
    inclusive_seconds: float
    exclusive_seconds: float | None = None


@dataclass
class ProcessResult:
    process_index: int
    return_code: int | None
    status: str
    timers: list[TimerRecord]
    timing_file: Path | None
    log_file: Path
    message: str = ""


@dataclass(frozen=True)
class BenchmarkOptions:
    case_name: str
    processes: tuple[int, ...]
    batch_sizes: tuple[int, ...]
    warmups: int
    repetitions: int
    output: Path
    name: str
    timeout: float | None
    continue_on_error: bool
    overwrite: bool = False


class CsvSink:
    """Append rows to a CSV while checking and preserving its schema."""

    def __init__(self, path: Path, fields: Sequence[str]):
        self.path = path
        self.fields = tuple(fields)
        self.row_count = 0
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.stat().st_size:
            with path.open(newline="", encoding="utf-8") as existing:
                header = next(csv.reader(existing), [])
            if tuple(header) != self.fields:
                raise ValueError(f"Existing CSV has an incompatible header: {path}")
        self._file = path.open("a", newline="", encoding="utf-8")
        self._writer = csv.DictWriter(self._file, fieldnames=self.fields)
        if self._file.tell() == 0:
            self._writer.writeheader()
            self._file.flush()

    def write(self, row: dict[str, object]) -> None:
        self._writer.writerow(row)
        self._file.flush()
        self.row_count += 1

    def close(self) -> None:
        self._file.close()

    def __enter__(self) -> "CsvSink":
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        self.close()


def parse_positive_int_list(value: str) -> tuple[int, ...]:
    """Parse a comma-separated, duplicate-free list of positive integers."""
    items = value.split(",")
    if any(not item.strip() for item in items):
        raise argparse.ArgumentTypeError("must be a comma-separated list of integers")
    try:
        values = tuple(int(item.strip()) for item in items)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a comma-separated list of integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("must contain at least one integer")
    if any(item < 1 for item in values):
        raise argparse.ArgumentTypeError("values must be positive")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("values must not be repeated")
    return values


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Time CLUBB across concurrent process counts and per-process batch sizes. "
            "Unrecognized arguments are passed to run_scm.py."
        ),
        allow_abbrev=False,
    )
    parser.add_argument("case_name", help="SCM case name; place it before forwarded options")
    parser.add_argument(
        "-processes",
        "--processes",
        type=parse_positive_int_list,
        default=(1,),
        metavar="N[,N...]",
        help="Concurrent CLUBB process counts (default: 1)",
    )
    parser.add_argument(
        "-batch_sizes",
        "--batch-sizes",
        "-columns",
        "--columns",
        dest="batch_sizes",
        type=parse_positive_int_list,
        default=(1,),
        metavar="N[,N...]",
        help="CLUBB runtime batch size assigned to each process (default: 1)",
    )
    parser.add_argument("-warmups", "--warmups", type=nonnegative_int, default=1)
    parser.add_argument("-repeats", "--repeats", dest="repetitions", type=positive_int, default=3)
    parser.add_argument(
        "-output",
        "--output",
        type=Path,
        default=REPO_ROOT / "output" / "timing",
        help="Directory in which compact profile folders are stored (default: output/timing)",
    )
    parser.add_argument("-name", "--name", default=None, help="Benchmark label (default: case name)")
    parser.add_argument(
        "-overwrite",
        "--overwrite",
        action="store_true",
        help="Replace the existing compact profile with the same normalized name",
    )
    parser.add_argument("-timeout", "--timeout", type=positive_float, default=None, metavar="SECONDS")
    parser.add_argument(
        "-continue_on_error",
        "--continue-on-error",
        action="store_true",
        help="Continue the sweep after a failed or timed-out process group",
    )
    return parser


def _option_occurrences(arguments: Sequence[str], names: Iterable[str]) -> list[str | None]:
    names = tuple(names)
    occurrences: list[str | None] = []
    for index, token in enumerate(arguments):
        if token in names:
            occurrences.append(arguments[index + 1] if index + 1 < len(arguments) else None)
            continue
        for name in names:
            if token.startswith(name + "="):
                occurrences.append(token.split("=", 1)[1])
                break
    return occurrences


def validate_passthrough(arguments: Sequence[str], parser: argparse.ArgumentParser) -> list[str]:
    forwarded = list(arguments)
    if forwarded[:1] == ["--"]:
        forwarded.pop(0)

    managed_options = (
        "-multicol", "--multicol", "-batch_size", "--batch_size",
        "-out_dir", "--out_dir",
    )
    for option in managed_options:
        if _option_occurrences(forwarded, (option,)):
            parser.error(f"{option} is managed by time_clubb.py and may not be forwarded")
    if _option_occurrences(forwarded, ("-primary_timer", "--primary-timer")):
        parser.error("the required run timer is internal and is not configurable")

    debug_values = _option_occurrences(forwarded, ("-debug", "--debug"))
    for value in debug_values:
        if value is None:
            parser.error("-debug requires a value")
        try:
            debug_level = int(value)
        except ValueError:
            continue  # Let run_scm.py report its own type error.
        if debug_level < 0:
            parser.error("negative -debug values disable timers and cannot be benchmarked")

    if not _option_occurrences(forwarded, ("-stats", "--stats")):
        forwarded.extend(("-stats", "none"))
    if not debug_values:
        forwarded.extend(("-debug", "0"))
    return forwarded


def parse_arguments(argv: Sequence[str] | None = None) -> tuple[BenchmarkOptions, list[str]]:
    parser = build_argument_parser()
    arguments = list(argv) if argv is not None else sys.argv[1:]
    # ``argparse`` treats ``-batch_size`` as a prefix of ``-batch_sizes`` even
    # with long-option abbreviation disabled, so reject managed child options
    # before parsing the benchmark's own options.
    for option in ("-multicol", "--multicol", "-batch_size", "--batch_size", "-out_dir", "--out_dir"):
        if _option_occurrences(arguments, (option,)):
            parser.error(f"{option} is managed by time_clubb.py and may not be forwarded")
    namespace, passthrough = parser.parse_known_args(arguments)
    forwarded = validate_passthrough(passthrough, parser)
    output = namespace.output.expanduser()
    if not output.is_absolute():
        output = Path.cwd() / output
    options = BenchmarkOptions(
        case_name=namespace.case_name,
        processes=namespace.processes,
        batch_sizes=namespace.batch_sizes,
        warmups=namespace.warmups,
        repetitions=namespace.repetitions,
        output=output.resolve(),
        name=namespace.name or namespace.case_name,
        timeout=namespace.timeout,
        continue_on_error=namespace.continue_on_error,
        overwrite=namespace.overwrite,
    )
    return options, forwarded


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_timer_text(text: str) -> list[TimerRecord]:
    """Normalize the built-in CLUBB or GPTL timer text."""
    lines = text.splitlines()
    if any("CLUBB TIMER SUMMARY" in line for line in lines):
        backend = "cpu_time"
        for line in lines:
            if line.strip().lower().startswith("backend:"):
                backend = line.split(":", 1)[1].strip()
                break
        records: list[TimerRecord] = []
        in_table = False
        for line in lines:
            if "Inclusive(s)" in line and "Exclusive(s)" in line:
                in_table = True
                continue
            if not in_table or not line.strip() or set(line.strip()) <= {"-", "="}:
                continue
            if line.lstrip().startswith("Total exclusive"):
                break
            match = _CPU_ROW.match(line)
            if match:
                records.append(
                    TimerRecord(
                        backend=backend,
                        time_basis="cpu_time",
                        thread=0,
                        timer_name=match.group(1).strip(),
                        calls=int(match.group(2)),
                        inclusive_seconds=_as_float(match.group(3)),
                        exclusive_seconds=_as_float(match.group(4)),
                    )
                )
        if not records:
            raise TimerParseError("CLUBB cpu_time summary contains no timer rows")
        return records

    gptl_rows: dict[int, list[tuple[TimerRecord, int]]] = {}
    current_thread: int | None = None
    in_table = False
    for line in lines:
        thread_match = _GPTL_THREAD.match(line)
        if thread_match:
            current_thread = int(thread_match.group(1))
            in_table = False
            continue
        if current_thread is not None and "Called" in line and "Recurse" in line and "Wall" in line:
            in_table = True
            continue
        if not in_table or current_thread is None:
            continue
        if line.lstrip().startswith(("Overhead sum", "Total calls")):
            in_table = False
            continue
        match = _GPTL_ROW.match(line)
        if match:
            marker = match.group("marker")
            # GPTL places the multiple-parent marker in column one.  Count it
            # as indentation so marked and unmarked rows retain the same tree
            # depth, while keeping the marker out of the timer name.
            visual_indent = len(match.group("indent")) + int(bool(marker))
            gptl_rows.setdefault(current_thread, []).append(
                (
                    TimerRecord(
                        backend="gptl",
                        time_basis="wallclock",
                        thread=current_thread,
                        timer_name=match.group("name").strip(),
                        calls=int(match.group("calls")),
                        inclusive_seconds=_as_float(match.group("wall")),
                    ),
                    visual_indent,
                )
            )
    if not gptl_rows:
        raise TimerParseError("timing file is neither a populated CLUBB cpu_time summary nor GPTL output")

    records = []
    for thread_rows in gptl_rows.values():
        minimum_indent = min(indent for _record, indent in thread_rows)
        ordered: dict[str, TimerRecord] = {}
        children: dict[str, list[str]] = {}
        assigned_parent: dict[str, str] = {}
        stack: list[str] = []
        for record, indent in thread_rows:
            depth = max(0, (indent - minimum_indent) // 2)
            if len(stack) > depth:
                del stack[depth:]
            parent = stack[depth - 1] if depth > 0 and len(stack) >= depth else ""
            if record.timer_name not in ordered:
                ordered[record.timer_name] = record
                if parent and parent != record.timer_name and record.timer_name not in assigned_parent:
                    assigned_parent[record.timer_name] = parent
                    children.setdefault(parent, []).append(record.timer_name)
            if len(stack) == depth:
                stack.append(record.timer_name)
            else:
                stack[depth] = record.timer_name

        for name, record in ordered.items():
            child_seconds = sum(
                ordered[child].inclusive_seconds
                for child in children.get(name, [])
                if child in ordered
            )
            records.append(
                TimerRecord(
                    backend=record.backend,
                    time_basis=record.time_basis,
                    thread=record.thread,
                    timer_name=record.timer_name,
                    calls=record.calls,
                    inclusive_seconds=record.inclusive_seconds,
                    exclusive_seconds=max(0.0, record.inclusive_seconds - child_seconds),
                )
            )
    return records


def parse_timer_file(path: Path) -> list[TimerRecord]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise TimerParseError(f"unable to read {path}: {exc}") from exc
    return parse_timer_text(text)


def _terminate_process(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return


def _kill_process(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return


def _stop_processes(processes: Sequence[subprocess.Popen[bytes]]) -> None:
    for process in processes:
        _terminate_process(process)
    deadline = time.monotonic() + 2.0
    for process in processes:
        if process.poll() is None:
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
    for process in processes:
        _kill_process(process)
    for process in processes:
        try:
            process.wait(timeout=1.0)
        except subprocess.TimeoutExpired:
            pass


def _find_timing_file(output_dir: Path) -> Path | None:
    files = sorted(output_dir.glob("*.timing"))
    return files[0] if len(files) == 1 else None


def run_process_group(
    *,
    options: BenchmarkOptions,
    forwarded: Sequence[str],
    process_count: int,
    batch_size: int,
    group_dir: Path,
    run_scm_path: Path = RUN_SCM,
    invocation_cwd: Path | None = None,
) -> tuple[list[ProcessResult], float]:
    """Launch and collect one equal-work group of ``run_scm.py`` processes."""
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    cwd = invocation_cwd or Path.cwd()
    running: list[tuple[int, subprocess.Popen[bytes], object, Path, Path]] = []
    started = time.perf_counter()
    timed_out: set[int] = set()

    try:
        for process_index in range(process_count):
            process_dir = group_dir / f"process_{process_index:03d}"
            process_dir.mkdir(parents=True, exist_ok=True)
            log_path = process_dir / "run_scm.log"
            log_file = log_path.open("wb")
            command = [
                sys.executable,
                str(run_scm_path),
                options.case_name,
                *forwarded,
                "-multicol",
                str(batch_size),
                "-batch_size",
                str(batch_size),
                "-out_dir",
                str(process_dir),
            ]
            try:
                process = subprocess.Popen(
                    command,
                    cwd=cwd,
                    env=environment,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except BaseException:
                log_file.close()
                raise
            running.append((process_index, process, log_file, process_dir, log_path))

        deadline = started + options.timeout if options.timeout is not None else None
        while any(process.poll() is None for _, process, _, _, _ in running):
            if deadline is not None and time.perf_counter() >= deadline:
                timed_out = {
                    index for index, process, _, _, _ in running if process.poll() is None
                }
                _stop_processes([process for _, process, _, _, _ in running])
                break
            time.sleep(0.02)
    except BaseException:
        _stop_processes([process for _, process, _, _, _ in running])
        raise
    finally:
        for _, _, log_file, _, _ in running:
            log_file.close()

    group_wall_seconds = time.perf_counter() - started
    results: list[ProcessResult] = []
    for process_index, process, _, process_dir, log_path in running:
        timing_file = _find_timing_file(process_dir)
        timers: list[TimerRecord] = []
        message = ""
        if process_index in timed_out:
            status = "timeout"
            message = f"process group exceeded {options.timeout:g} seconds"
        elif process.returncode != 0:
            status = "failed"
            message = f"run_scm.py exited with status {process.returncode}"
        elif timing_file is None:
            status = "missing_timing"
            message = "expected exactly one .timing file"
        else:
            try:
                timers = parse_timer_file(timing_file)
            except TimerParseError as exc:
                status = "timer_parse_error"
                message = str(exc)
            else:
                if any(timer.timer_name == REQUIRED_TIMER for timer in timers):
                    status = "success"
                else:
                    status = "missing_required_timer"
                    message = f"required timer {REQUIRED_TIMER!r} was not found"
        results.append(
            ProcessResult(
                process_index=process_index,
                return_code=process.returncode,
                status=status,
                timers=timers,
                timing_file=timing_file,
                log_file=log_path,
                message=message,
            )
        )
    return results, group_wall_seconds


def write_timing_rows(
    sink: CsvSink,
    results: Sequence[ProcessResult],
    *,
    batch_id: str,
    phase: str,
    repetition: int,
    group_wall_seconds: float,
) -> None:
    """Persist raw process/thread/timer observations for one execution."""
    for result in results:
        timer_rows: Sequence[TimerRecord | None] = result.timers or (None,)
        for timer in timer_rows:
            sink.write(
                {
                    "batch_id": batch_id,
                    "phase": phase,
                    "repetition": repetition,
                    "process": result.process_index,
                    "thread": timer.thread if timer else "",
                    "timer": timer.timer_name if timer else "",
                    "calls": timer.calls if timer else "",
                    "inclusive_s": timer.inclusive_seconds if timer else "",
                    "exclusive_s": (
                        timer.exclusive_seconds
                        if timer is not None and timer.exclusive_seconds is not None
                        else ""
                    ),
                    "status": result.status,
                    "return_code": result.return_code,
                    "message": result.message,
                }
            )
    group_status = "success" if all(result.status == "success" for result in results) else "failed"
    messages = "; ".join(
        f"process {result.process_index}: {result.message}"
        for result in results
        if result.message
    )
    sink.write(
        {
            "batch_id": batch_id,
            "phase": phase,
            "repetition": repetition,
            "process": "",
            "thread": "",
            "timer": GROUP_WALL_TIMER,
            "calls": 1,
            "inclusive_s": group_wall_seconds,
            "exclusive_s": "",
            "status": group_status,
            "return_code": "",
            "message": messages,
        }
    )


def choose_representative(
    candidates: Sequence[tuple[str, int, ProcessResult]],
) -> tuple[str, int, ProcessResult] | None:
    """Prefer the first successful measured process, then the best fallback."""
    for phase, require_success in (
        ("measured", True),
        ("warmup", True),
        ("measured", False),
        ("warmup", False),
    ):
        for candidate in candidates:
            if candidate[0] == phase and (not require_success or candidate[2].status == "success"):
                return candidate
    return None


def archive_representative(
    profile_dir: Path,
    batch_id: str,
    case_name: str,
    representative: tuple[str, int, ProcessResult] | None,
) -> tuple[str, int | str, int | str]:
    """Copy the compact audit set from one representative child process."""
    if representative is None:
        return "", "", ""
    phase, repetition, result = representative
    source_dir = result.log_file.parent
    destination = Path(profile_dir) / "logs" / batch_id
    destination.mkdir(parents=True, exist_ok=True)
    model_log = source_dir / f"{case_name}_log"
    if not model_log.is_file():
        model_log = result.log_file
    candidates = (
        (source_dir / f"{case_name}.in", destination / f"{case_name}.in"),
        (source_dir / f"{case_name}_setup.txt", destination / f"{case_name}_setup.txt"),
        (model_log, destination / f"{case_name}.log"),
        (result.timing_file, destination / f"{case_name}.timing"),
    )
    for source, target in candidates:
        if source is not None and source.is_file():
            shutil.copy2(source, target)
    setup_target = destination / f"{case_name}_setup.txt"
    if not setup_target.is_file():
        setup_target.write_text(
            "CLUBB did not emit a native setup file for this run.\n"
            "This normally occurs at debug level 0; see the preserved input namelist "
            "and profile.json for effective settings and provenance.\n",
            encoding="utf-8",
        )
    return phase, repetition, result.process_index


def _git_revision() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def _git_dirty() -> bool | None:
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return bool(result.stdout.strip())


def _option_value(arguments: Sequence[str], names: Sequence[str]) -> str:
    values = _option_occurrences(arguments, names)
    return str(values[-1] or "") if values else ""


def _file_sha256(path: Path) -> str:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return ""


def _cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name") and ":" in line:
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def _profile_build_metadata(forwarded: Sequence[str]) -> dict[str, object]:
    requested = _option_value(forwarded, ("-exe", "--exe"))
    executable = Path(requested).expanduser() if requested else None
    if executable is not None and not executable.is_absolute():
        executable = (REPO_ROOT / executable).resolve()
    install_dir = _option_value(forwarded, ("-install_dir", "--install_dir"))
    implementation = (
        "python" if _option_occurrences(forwarded, ("-python", "--python"))
        else "jax" if _option_occurrences(forwarded, ("-jax", "--jax"))
        else "fortran"
    )
    return {
        "implementation": implementation,
        "requested_executable": str(executable or ""),
        "executable_sha256": _file_sha256(executable) if executable is not None else "",
        "requested_install_dir": install_dir,
    }


def _profile_environment() -> dict[str, object]:
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_model": _cpu_model(),
        "logical_cpu_count": os.cpu_count(),
        "python": platform.python_version(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "1"),
    }


def _reported_executable(results: Sequence[ProcessResult]) -> str:
    for result in results:
        try:
            text = result.log_file.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        matches = re.findall(r"(?mi)^\s*-\s*using executable:\s*(.+?)\s*$", text)
        if matches:
            return matches[-1].strip()
    return ""


def _resolved_build_metadata(reported_executable: str) -> dict[str, str]:
    if not reported_executable:
        return {}
    if any(character.isspace() for character in reported_executable):
        return {"reported_executable": reported_executable}
    candidate = Path(reported_executable).expanduser()
    if not candidate.is_absolute():
        candidate = REPO_ROOT / candidate
    resolved = candidate.resolve()
    return {
        "reported_executable": reported_executable,
        "resolved_executable": str(resolved),
        "executable_sha256": _file_sha256(resolved),
    }


def _namelist_value(text: str, name: str) -> str:
    matches = re.findall(
        rf"(?mi)^\s*{re.escape(name)}\s*=\s*([^!,\n]+)",
        text,
    )
    return matches[-1].strip().strip("'\"") if matches else ""


def _namelist_float(text: str, name: str) -> float | None:
    value = _namelist_value(text, name)
    try:
        return _as_float(value)
    except ValueError:
        return None


def _observed_profile_metadata(profile_root: Path, case_name: str) -> dict[str, dict[str, object]]:
    """Recover the executable and effective case grid from generated artifacts."""
    logs = sorted((profile_root / "logs").glob(f"*/{case_name}.log"))
    reported_executable = ""
    for log_path in logs:
        try:
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        matches = re.findall(r"(?mi)^\s*-\s*using executable:\s*(.+?)\s*$", log_text)
        if matches:
            reported_executable = matches[-1].strip()
            break
    benchmark: dict[str, object] = {}
    namelists = sorted((profile_root / "logs").glob(f"*/{case_name}.in"))
    if namelists:
        try:
            text = namelists[0].read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        grid_type_value = _namelist_float(text, "grid_type")
        nzmax_value = _namelist_float(text, "nzmax")
        grid_type = int(grid_type_value) if grid_type_value is not None else None
        nzmax = int(nzmax_value) if nzmax_value is not None else None
        vertical_levels = nzmax
        deltaz = _namelist_float(text, "deltaz_nl")
        zm_init = _namelist_float(text, "zm_init_nl")
        zm_top = _namelist_float(text, "zm_top_nl")
        if grid_type == 1 and None not in (deltaz, zm_init, zm_top) and float(deltaz) != 0:
            vertical_levels = int((float(zm_top) - float(zm_init) + abs(float(deltaz))) // abs(float(deltaz)))
        benchmark.update(
            {
                "vertical_levels": vertical_levels,
                "grid_type": grid_type,
                "nzmax": nzmax,
                "deltaz": deltaz,
                "zm_init": zm_init,
                "zm_top": zm_top,
                "zt_grid": _namelist_value(text, "zt_grid_fname"),
                "zm_grid": _namelist_value(text, "zm_grid_fname"),
                "time_initial": _namelist_float(text, "time_initial"),
                "time_final": _namelist_float(text, "time_final"),
                "dt_main": _namelist_float(text, "dt_main"),
            }
        )

    return {
        "benchmark": {key: value for key, value in benchmark.items() if value is not None},
        "build": _resolved_build_metadata(reported_executable),
    }


def run_benchmark(
    options: BenchmarkOptions,
    forwarded: Sequence[str],
    *,
    run_scm_path: Path = RUN_SCM,
    invocation_cwd: Path | None = None,
) -> int:
    if not run_scm_path.is_file():
        raise FileNotFoundError(f"run_scm.py not found: {run_scm_path}")

    options.output.mkdir(parents=True, exist_ok=True)
    run_id, profile_root = reserve_profile_directory(
        options.output,
        options.name,
        overwrite=options.overwrite,
    )
    started_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    git_revision = _git_revision()
    forwarded_args = shlex.join(forwarded)
    failed = False
    completed_batches = 0
    failed_batches = 0
    observed_model_steps: set[int] = set()
    observed_backends: set[str] = set()
    observed_time_bases: set[str] = set()
    reported_executable = ""
    terminal_status = "error"
    batch_rows: list[dict[str, object]] = []

    manifest = create_profile_manifest(
        run_id=run_id,
        label=options.name,
        case_name=options.case_name,
        started_utc=started_utc,
        git_revision=git_revision,
        git_dirty=_git_dirty(),
        forwarded_args=forwarded_args,
        omp_num_threads=os.environ.get("OMP_NUM_THREADS", "1"),
        benchmark={
            "processes": list(options.processes),
            "batch_sizes": list(options.batch_sizes),
            "warmups": options.warmups,
            "repetitions": options.repetitions,
            "timeout_seconds": options.timeout,
            "continue_on_error": options.continue_on_error,
        },
        environment=_profile_environment(),
        build=_profile_build_metadata(forwarded),
    )
    write_profile_manifest(profile_root, manifest)
    write_profile_readme(profile_root, manifest)
    write_batches(profile_root, batch_rows)

    print(f"Profile library: {options.output}")
    print(f"Run ID: {run_id}")
    print(f"Portable profile: {profile_root}")
    print(f"Forwarding to run_scm.py: {forwarded_args}")

    try:
        with CsvSink(profile_root / TIMINGS_FILE, TIMING_FIELDS) as timing_sink:
            for process_count in options.processes:
                for batch_size in options.batch_sizes:
                    batch_id = f"p{process_count:04d}_b{batch_size:06d}"
                    candidates: list[tuple[str, int, ProcessResult]] = []
                    batch_failed_runs = 0
                    warmups_completed = 0
                    measurements_completed = 0
                    stop_after_batch = False
                    batch_exception: BaseException | None = None
                    batch_row: dict[str, object] = {
                        "batch_id": batch_id,
                        "process_count": process_count,
                        "batch_size": batch_size,
                        "total_batch_size": process_count * batch_size,
                        "status": "running",
                        "warmups_completed": 0,
                        "measurements_completed": 0,
                        "failed_runs": 0,
                        "representative_phase": "",
                        "representative_repetition": "",
                        "representative_process": "",
                    }
                    # Publish the workload coordinates before any timing rows.
                    # Dash can then join each flushed observation to meaningful
                    # process and batch-size metadata while this point is active.
                    batch_rows.append(batch_row)
                    write_batches(profile_root, batch_rows)
                    with tempfile.TemporaryDirectory(prefix=f"clubb_profile_{batch_id}_") as temporary:
                        temporary_root = Path(temporary)
                        try:
                            for phase, count in (("warmup", options.warmups), ("measured", options.repetitions)):
                                for repetition in range(1, count + 1):
                                    label = "Warm-up" if phase == "warmup" else "Repetition"
                                    print(
                                        f"{label} {repetition}/{count}: {process_count} process(es) "
                                        f"x batch size {batch_size}"
                                    )
                                    results, group_wall = run_process_group(
                                        options=options,
                                        forwarded=forwarded,
                                        process_count=process_count,
                                        batch_size=batch_size,
                                        group_dir=temporary_root / f"{phase}_{repetition:03d}",
                                        run_scm_path=run_scm_path,
                                        invocation_cwd=invocation_cwd,
                                    )
                                    write_timing_rows(
                                        timing_sink,
                                        results,
                                        batch_id=batch_id,
                                        phase=phase,
                                        repetition=repetition,
                                        group_wall_seconds=group_wall,
                                    )
                                    candidates.extend((phase, repetition, result) for result in results)
                                    if not reported_executable:
                                        reported_executable = _reported_executable(results)
                                    for result in results:
                                        for timer in result.timers:
                                            observed_backends.add(timer.backend)
                                            observed_time_bases.add(timer.time_basis)
                                            if timer.timer_name == "advance_clubb_core_api":
                                                observed_model_steps.add(timer.calls)
                                    run_failed = any(result.status != "success" for result in results)
                                    batch_failed_runs += int(run_failed)
                                    if phase == "warmup":
                                        warmups_completed += 1
                                    else:
                                        measurements_completed += 1
                                    batch_row.update(
                                        {
                                            "warmups_completed": warmups_completed,
                                            "measurements_completed": measurements_completed,
                                            "failed_runs": batch_failed_runs,
                                        }
                                    )
                                    write_batches(profile_root, batch_rows)
                                    if run_failed:
                                        failed = True
                                        print(f"{label} failed; compact timing rows were preserved.", file=sys.stderr)
                                        if not options.continue_on_error:
                                            stop_after_batch = True
                                            break
                                if stop_after_batch:
                                    break
                        except BaseException as exc:
                            batch_exception = exc

                        representative = choose_representative(candidates)
                        representative_phase, representative_repetition, representative_process = archive_representative(
                            profile_root,
                            batch_id,
                            options.case_name,
                            representative,
                        )

                    batch_status = (
                        "complete"
                        if batch_failed_runs == 0 and measurements_completed == options.repetitions
                        else "partial"
                        if measurements_completed or warmups_completed
                        else "failed"
                    )
                    batch_row.update(
                        {
                            "status": batch_status,
                            "warmups_completed": warmups_completed,
                            "measurements_completed": measurements_completed,
                            "failed_runs": batch_failed_runs,
                            "representative_phase": representative_phase,
                            "representative_repetition": representative_repetition,
                            "representative_process": representative_process,
                        }
                    )
                    write_batches(profile_root, batch_rows)
                    completed_batches += 1
                    failed_batches += int(batch_status != "complete")
                    manifest_results = dict(manifest.get("results") or {})
                    manifest_results.update(
                        {
                            "observed_model_steps": sorted(observed_model_steps),
                            "timer_backends": sorted(observed_backends),
                            "time_bases": sorted(observed_time_bases),
                        }
                    )
                    build = dict(manifest.get("build") or {})
                    build.update(_resolved_build_metadata(reported_executable))
                    manifest = {**manifest, "build": build, "results": manifest_results}
                    manifest = update_profile_manifest(
                        profile_root,
                        manifest,
                        timing_rows=timing_sink.row_count,
                        batches_completed=completed_batches,
                        failed_batches=failed_batches,
                    )
                    if batch_exception is not None:
                        raise batch_exception
                    if stop_after_batch:
                        terminal_status = "partial" if completed_batches else "failed"
                        return 1

            terminal_status = "partial" if failed else "complete"
            return 1 if failed else 0
    except KeyboardInterrupt:
        terminal_status = "cancelled"
        raise
    finally:
        try:
            observed = _observed_profile_metadata(profile_root, options.case_name)
            benchmark = dict(manifest.get("benchmark") or {})
            benchmark.update(observed.get("benchmark") or {})
            build = dict(manifest.get("build") or {})
            for key, value in (observed.get("build") or {}).items():
                if value not in (None, ""):
                    build[key] = value
            for key, value in _resolved_build_metadata(reported_executable).items():
                if value not in (None, ""):
                    build[key] = value
            manifest_results = dict(manifest.get("results") or {})
            manifest_results.update(
                {
                    "observed_model_steps": sorted(observed_model_steps),
                    "timer_backends": sorted(observed_backends),
                    "time_bases": sorted(observed_time_bases),
                }
            )
            manifest = {
                **manifest,
                "benchmark": benchmark,
                "build": build,
                "results": manifest_results,
            }
            update_profile_manifest(
                profile_root,
                manifest,
                status=terminal_status,
                timing_rows=(timing_sink.row_count if "timing_sink" in locals() else 0),
                batches_completed=completed_batches,
                failed_batches=failed_batches,
                completed=True,
            )
        except OSError as exc:
            print(f"Warning: unable to finalize profile metadata: {exc}", file=sys.stderr)


def main(argv: Sequence[str] | None = None) -> int:
    # The managed Dash broker ignores SIGINT so Ctrl-C remains owned by the
    # foreground manager. Signal dispositions survive exec, so explicitly
    # restore Python's handler before doing any work. This lets Profile cancel
    # raise KeyboardInterrupt and unwind every active run_scm.py process group.
    _restore_sigint_handler()
    options, forwarded = parse_arguments(argv)
    try:
        return run_benchmark(options, forwarded)
    except KeyboardInterrupt:
        print("Timing run interrupted; completed CSV rows were preserved.", file=sys.stderr)
        return 130
    except (OSError, ValueError) as exc:
        print(f"time_clubb.py: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
