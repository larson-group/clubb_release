#!/usr/bin/env python3
"""Measure complete CLUBB runs with groups of reusable case processes.

For each process count and per-process column count, warm up a group of workers
and measure full-case repetitions in those same processes. Every run includes
model initialization, advancement and cleanup; process reuse lets JAX keep its
compiled programs between runs. Runtime selection belongs to clubb_case_runner
and model input preparation belongs to create_case_namelist.

Usage:
    From the repository root, with a current native build:

        ./utilities/time_clubb.py bomex
        ./utilities/time_clubb.py bomex -process_counts 1,2 -batch_sizes 1,4

    CLI and execution:
    - main(argv), parse_arguments(argv) and build_argument_parser(): run the CLI
      or parse BenchmarkOptions and forwarded SCM options.
    - validate_passthrough(): add stats/debug defaults and reject options whose
      values the timing sweep assigns itself. parse_positive_int_list,
      nonnegative_int, positive_int and positive_float validate CLI values.
    - run_benchmark(options, forwarded): write a sweep and return its exit code.
    - ProcessGroup: start a worker group on context entry, run(group_dir) to
      execute one complete case per peer, and close() to release every runner.
      run_process_group() is the convenience entry point for a single execution.
    - BenchmarkOptions, ProcessResult and TimerRecord: workload settings,
      per-process results and normalized timer values.

    Profile output:
    - parse_timer_text(), parse_timer_file(): read built-in or GPTL summaries.
    - CsvSink: append rows; write_timing_rows() writes process and group values.
    - choose_representative(), archive_representative(): preserve one example's
      inputs, log and timing file. report_failed_group() preserves diagnostics.

Understanding:
    Put the case before forwarded options, which use the same parsers as
    run_scm.py. Defaults are one process, one column, one full warmup and three
    repetitions, with stats disabled and debug 0. -multicol, -batch_size and
    -output_dir are reserved: this script assigns columns and distinct output
    directories to peers. A workload is one process-count/batch-size pair.

    All workers finish a warmup/repetition before the next is submitted. Inputs
    are prepared before sequential submissions to all peers; no peer is waited
    on until all are submitted. Group wall time covers submission through the
    last completion, excluding startup, input preparation and inter-run waiting.
    Named .timing regions retain their backend's clock and nesting semantics.

    Each profile contains profile.json, batches.csv, timings.csv and selected
    logs/inputs. Raw timing phases are startup, warmup and measured. Startup
    includes input validation and runtime/environment setup. Temporary per-run
    outputs are removed after representative/failure artifacts are preserved.
    Failed or timed-out workers stop that workload; -continue_on_error permits
    the next workload, without retrying the failed one. JAX summaries use host
    wall time; native summaries use CPU time or GPTL wall time. The Python SCM
    does not emit timing summaries.
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

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

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


from utilities.clubb_case_runner import clubb_case_runner, parse_case_arguments
from input.scripts.create_case_namelist import create_case_namelist_from_args


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
    """One named timer's calls and seconds for a backend/thread.

    Inclusive time contains children; exclusive time omits them. GPTL exclusive
    time is inferred from its printed call tree. time_basis identifies the clock.
    """

    backend: str
    time_basis: str
    thread: int
    timer_name: str
    calls: int
    inclusive_seconds: float
    exclusive_seconds: float | None = None


@dataclass
class ProcessResult:
    """A peer's status, normalized timers and paths to its run artifacts."""

    process_index: int
    return_code: int | None
    status: str
    timers: list[TimerRecord]
    timing_file: Path | None
    log_file: Path
    message: str = ""


@dataclass(frozen=True)
class BenchmarkOptions:
    """Sweep settings: process counts, columns per process, repeats and output policy."""

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
        """Open the CSV at path, checking its existing header against fields."""
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
        """Append and flush one schema-compatible row, updating row_count."""
        self._writer.writerow(row)
        self._file.flush()
        self.row_count += 1

    def close(self) -> None:
        """Close the CSV output file."""
        self._file.close()

    def __enter__(self) -> "CsvSink":
        """Return this open sink for context-managed row writing."""
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        """Close the CSV output file when leaving the context."""
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
    """Parse a CLI integer and reject values below zero."""
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def positive_int(value: str) -> int:
    """Parse a CLI integer and reject values below one."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def positive_float(value: str) -> float:
    """Parse a CLI floating-point value and reject values at or below zero."""
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def build_argument_parser() -> argparse.ArgumentParser:
    """Return the timing-option parser; model/runtime options are parsed separately."""
    parser = argparse.ArgumentParser(
        description=(
            "Time CLUBB across concurrent process counts and per-process batch sizes. "
            "Unrecognized arguments are passed to run_scm.py."
        ),
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument("case_name", help="SCM case name; place it before forwarded options")
    parser.add_argument(
        '-process_counts', dest='processes',
        type=parse_positive_int_list,
        default=(1,),
        metavar="N[,N...]",
        help="Concurrent CLUBB process counts (default: 1)",
    )
    parser.add_argument(
        '-batch_sizes',
        dest="batch_sizes",
        type=parse_positive_int_list,
        default=(1,),
        metavar="N[,N...]",
        help="CLUBB runtime batch size assigned to each process (default: 1)",
    )
    parser.add_argument('-warmups', dest='warmups', type=nonnegative_int, default=1)
    parser.add_argument('-repeats', dest="repetitions", type=positive_int, default=3)
    parser.add_argument(
        '-output_root', dest='output',
        type=Path,
        default=REPO_ROOT / "output" / "timing",
        help="Directory in which compact profile folders are stored (default: output/timing)",
        metavar='DIR',
    )
    parser.add_argument('-name', dest='name', default=None, help="Benchmark label (default: case name)")
    parser.add_argument(
        '-overwrite', dest='overwrite',
        action="store_true",
        help="Replace the existing compact profile with the same normalized name",
    )
    parser.add_argument('-timeout', dest='timeout', type=positive_float, default=None, metavar="SECONDS")
    parser.add_argument(
        '-continue_on_error', dest='continue_on_error',
        action="store_true",
        help="Continue the sweep after a failed or timed-out process group",
    )
    return parser


def _option_occurrences(arguments: Sequence[str], names: Iterable[str]) -> list[str | None]:
    """Return values for matching separate/attached options, using None if missing."""
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
    """Validate forwarded SCM options and add timing defaults for statistics and debug level."""
    forwarded = list(arguments)
    if forwarded[:1] == ["--"]:
        forwarded.pop(0)

    managed_options = (
        "-multicol", "--multicol", "-batch_size", "--batch_size",
        "-output_dir", "--out_dir", "-benchmark_worker",
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

    if not _option_occurrences(forwarded, ("-stats",)):
        forwarded.extend(("-stats", "none"))
    if not debug_values:
        forwarded.extend(("-debug", "0"))
    return forwarded


def parse_arguments(argv: Sequence[str] | None = None) -> tuple[BenchmarkOptions, list[str]]:
    """Return resolved sweep options and forwarded SCM arguments from argv or sys.argv."""
    parser = build_argument_parser()
    arguments = list(argv) if argv is not None else sys.argv[1:]
    # ``argparse`` treats ``-batch_size`` as a prefix of ``-batch_sizes`` even
    # with long-option abbreviation disabled, so reject managed child options
    # before parsing the benchmark's own options.
    for option in ("-multicol", "--multicol", "-batch_size", "--batch_size", "-output_dir", "--out_dir"):
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
    """Convert a timer number, including Fortran D-exponent notation, to float."""
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
                        time_basis="wallclock" if backend == "wall_clock" else "cpu_time",
                        thread=0,
                        timer_name=match.group(1).strip(),
                        calls=int(match.group(2)),
                        inclusive_seconds=_as_float(match.group(3)),
                        exclusive_seconds=_as_float(match.group(4)),
                    )
                )
        if not records:
            raise TimerParseError("CLUBB timer summary contains no timer rows")
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
    """Return normalized timers from path, raising TimerParseError for unreadable output."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise TimerParseError(f"unable to read {path}: {exc}") from exc
    return parse_timer_text(text)


def _find_timing_file(output_dir: Path) -> Path | None:
    """Return the sole .timing file directly in output_dir, or None if ambiguous/missing."""
    files = sorted(output_dir.glob("*.timing"))
    return files[0] if len(files) == 1 else None


class ProcessGroup:
    """Keep one worker per peer through a workload's warmups and repetitions.

    Context entry prepares/validates inputs and initializes every runner. run()
    prepares new outputs, submits to all peers and waits for all completions.
    The case runner owns subprocesses/output; this group owns their scheduling.
    Use close() or a context manager to release the entire group after failure.
    """

    def __init__(self, *, options, forwarded, process_count, batch_size, group_dir,
                 runner_factory=clubb_case_runner,
                 prepare_case=create_case_namelist_from_args):
        """Store one process_count/batch_size workload and parse its case options.

        forwarded contains ordinary SCM options; options supplies timeout and
        sweep settings. group_dir holds startup files. runner_factory and
        prepare_case default to the shared runtime and namelist owners; focused
        checks substitute synthetic implementations. No worker starts here.
        """
        self.options = options
        self.forwarded = forwarded
        self.process_count = process_count
        self.batch_size = batch_size
        self.group_dir = Path(group_dir)
        self.runner_factory = runner_factory
        self.prepare_case = prepare_case
        self.workers = []
        self.broken = False
        self.timed_out = set()
        self.startup_wall_seconds = 0.0
        self.startup_error = ''
        self.startup_return_code = None
        self.args = parse_case_arguments([
            options.case_name, *forwarded, '-multicol', str(batch_size),
            '-batch_size', str(batch_size),
        ])

    def __enter__(self):
        """Validate initial inputs and start all workers, recording startup wall time."""
        started = time.perf_counter()
        environment = os.environ.copy()
        environment.setdefault('OMP_NUM_THREADS', '1')
        try:
            for index in range(self.process_count):
                directory = self.group_dir / f'process_{index:03d}'
                directory.mkdir(parents=True, exist_ok=True)
                log_path = directory / 'startup.log'
                runner = self.runner_factory(
                    self.args, env=environment, persistent=True,
                    log_file=log_path, timeout=self.options.timeout,
                )
                worker = {'index': index, 'runner': runner, 'log_path': log_path,
                          'startup_path': log_path}
                self.workers.append(worker)
                # Validate generated inputs before starting the model runtime.
                # Each execution later gets its own input and output directory.
                with log_path.open('a', encoding='utf-8') as log:
                    self.prepare_case(
                        self.args, directory, log_file=log, timeout=self.options.timeout,
                    )
                runner.init()
            self.startup_wall_seconds = time.perf_counter() - started
            return self
        except Exception as exc:
            if not self.workers:
                raise
            self.broken = True
            self.startup_error = str(exc)
            self.startup_return_code = getattr(exc, 'returncode', 1)
            if isinstance(exc, (TimeoutError, subprocess.TimeoutExpired)):
                self.timed_out.add(self.workers[-1]['index'])
            with self.workers[-1]['startup_path'].open('a', encoding='utf-8') as log:
                log.write(str(exc) + '\n')
            self.close()
            self.startup_wall_seconds = time.perf_counter() - started
            return self
        except BaseException:
            self.close()
            raise

    def _wait(self):
        """Wait for every submitted peer; close the cohort on failure or timeout."""
        pending = {worker['index'] for worker in self.workers}
        deadline = (time.perf_counter() + self.options.timeout
                    if self.options.timeout is not None else None)
        while pending:
            for worker in self.workers:
                status = worker['runner'].status()
                if status['state'] == 'completed':
                    pending.discard(worker['index'])
                elif status['state'] in ('failed', 'closed'):
                    self.broken = True
            if deadline is not None and time.perf_counter() >= deadline and pending:
                self.timed_out.update(pending)
                self.broken = True
            if self.broken:
                self.close()
                break
            if pending:
                time.sleep(0.001)

    def run(self, group_dir):
        """Return (peer results, group wall seconds) for one complete execution.

        Write each peer's inputs/output under group_dir, preparing all inputs
        before starting the wall clock for submission and model execution. Each
        peer must finish and produce a readable .timing file containing
        advance_clubb_to_end. A process/preparation failure closes all workers;
        result status also identifies missing or malformed timing output.
        """
        group_dir = Path(group_dir)
        started = time.perf_counter()
        deadline = (started + self.options.timeout
                    if self.options.timeout is not None else None)
        for worker in self.workers:
            worker['submitted'] = False
            worker.pop('namelist', None)
            directory = group_dir / f"process_{worker['index']:03d}"
            directory.mkdir(parents=True, exist_ok=True)
            worker['log_path'] = directory / 'run_scm.log'
            worker['log_path'].write_bytes(worker['startup_path'].read_bytes())
        for worker in self.workers if not self.broken else []:
            directory = worker['log_path'].parent
            try:
                remaining = None if deadline is None else max(0., deadline-time.perf_counter())
                with worker['log_path'].open('a', encoding='utf-8') as log:
                    worker['namelist'] = self.prepare_case(
                        self.args, directory, log_file=log, timeout=remaining,
                    )
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                self.broken = True
                if isinstance(exc, (TimeoutError, subprocess.TimeoutExpired)):
                    self.timed_out.add(worker['index'])
                with worker['log_path'].open('a', encoding='utf-8') as log:
                    log.write(str(exc) + '\n')
                break
        # Exclude input generation; submit every peer before polling completion.
        if not self.broken:
            started = time.perf_counter()
            for worker in self.workers:
                try:
                    worker['runner'].run(worker['namelist'], log_file=worker['log_path'])
                    worker['submitted'] = True
                except (OSError, RuntimeError) as exc:
                    with worker['log_path'].open('a', encoding='utf-8') as log:
                        log.write(str(exc) + '\n')
                    self.broken = True
                    break
            if not self.broken:
                self._wait()
        if self.broken:
            self.close()
        elapsed = time.perf_counter() - started
        results = []
        for worker in self.workers:
            index = worker['index']
            current = worker['runner'].status()
            snapshot = worker.get('closed_status', current)
            timing_file = _find_timing_file(group_dir / f'process_{index:03d}')
            timers = []
            code = snapshot['return_code'] if snapshot['return_code'] is not None else current['return_code']
            if code is None and self.startup_error:
                code = self.startup_return_code
            message = ''
            if index in self.timed_out:
                status = 'timeout'
                message = f'process group exceeded {self.options.timeout:g} seconds'
            elif not worker['submitted'] or snapshot['state'] != 'completed':
                status = 'failed'
                message = snapshot['error'] or self.startup_error or f'CLUBB runner exited with status {code} before completing the case'
            elif timing_file is None:
                status, message = 'missing_timing', 'expected exactly one .timing file'
            else:
                try:
                    timers = parse_timer_file(timing_file)
                except TimerParseError as exc:
                    status, message = 'timer_parse_error', str(exc)
                else:
                    status = ('success' if any(t.timer_name == REQUIRED_TIMER for t in timers)
                              else 'missing_required_timer')
                    if status != 'success':
                        message = f'required timer {REQUIRED_TIMER!r} was not found'
            results.append(ProcessResult(index, code, status, timers, timing_file,
                                         worker['log_path'], message))
        return results, elapsed

    def close(self):
        """Save final status and close every runner owned by this cohort."""
        for worker in self.workers:
            worker.setdefault('closed_status', worker['runner'].status())
            worker['runner'].close()

    def __exit__(self, exc_type, exc_value, traceback):
        """Close every worker when leaving the cohort context."""
        self.close()


def run_process_group(
    *, options: BenchmarkOptions, forwarded: Sequence[str], process_count: int,
    batch_size: int, group_dir: Path, runner_factory=clubb_case_runner,
    prepare_case=create_case_namelist_from_args,
):
    """Run a full case per peer and return (peer results, group wall seconds).

    Unlike run_benchmark, this helper does not reuse a group across warmups or
    repetitions. Arguments select one workload and its input/output directory.
    """
    with ProcessGroup(
        options=options, forwarded=forwarded, process_count=process_count,
        batch_size=batch_size, group_dir=group_dir, runner_factory=runner_factory,
        prepare_case=prepare_case,
    ) as group:
        return group.run(group_dir)


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


def report_failed_group(
    results: Sequence[ProcessResult],
    profile_dir: Path,
    batch_id: str,
    case_name: str,
    phase: str,
    repetition: int,
) -> None:
    """Show one failed child's diagnostics and retain its full log."""
    failed = [result for result in results if result.status != "success"]
    result = failed[0]
    destination = profile_dir / "logs" / batch_id
    destination.mkdir(parents=True, exist_ok=True)
    log_path = destination / f"{phase}_{repetition:03d}_process_{result.process_index:03d}_failed.log"
    sources = [result.log_file]
    model_log = result.log_file.parent / f"{case_name}_log"
    if model_log.is_file():
        sources.append(model_log)
    with log_path.open("wb") as output:
        for source in sources:
            output.write(f"--- {source.name} ---\n".encode())
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, output)
            output.write(b"\n")
    # Keep the main console bounded even when a child wrote a very large log.
    with log_path.open("rb") as handle:
        handle.seek(max(0, log_path.stat().st_size - 8192))
        tail = handle.read().decode("utf-8", errors="replace")
    print(
        f"{len(failed)}/{len(results)} processes failed. "
        f"Process {result.process_index}: {result.message}",
        file=sys.stderr,
    )
    print("\n".join(tail.splitlines()[-40:]), file=sys.stderr)
    print(f"Full failure log: {log_path}", file=sys.stderr, flush=True)


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
    """Return the checkout commit hash, or an empty string when Git is unavailable."""
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
    """Report tracked working-tree changes, or None when Git status is unavailable."""
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
    """Return the last supplied option value, or an empty string when absent."""
    values = _option_occurrences(arguments, names)
    return str(values[-1] or "") if values else ""


def _file_sha256(path: Path) -> str:
    """Return the file content digest, or an empty string when it cannot be read."""
    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return ""


def _cpu_model() -> str:
    """Return a CPU description from /proc/cpuinfo or the platform fallback."""
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name") and ":" in line:
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def _profile_build_metadata(forwarded: Sequence[str]) -> dict[str, object]:
    """Record requested backend/install options and the explicit executable digest."""
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
    """Return host, Python and OpenMP settings for the saved profile manifest."""
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
    """Read the first available runtime selection message from the peer logs."""
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
    """Record a reported command and resolve/hash it when it is a single file path."""
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
    """Read the last simple scalar assignment for name from the generated namelist text."""
    matches = re.findall(
        rf"(?mi)^\s*{re.escape(name)}\s*=\s*([^!,\n]+)",
        text,
    )
    return matches[-1].strip().strip("'\"") if matches else ""


def _namelist_float(text: str, name: str) -> float | None:
    """Read a numeric namelist assignment, returning None if absent or unparseable."""
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
    runner_factory=clubb_case_runner,
    prepare_case=create_case_namelist_from_args,
) -> int:
    """Write a profile for options' sweep and return 0 if every execution succeeds.

    forwarded supplies runtime/model options. Each workload owns its worker
    group through all warmups and measured repetitions, with fresh model setup
    per execution. Flush raw rows and workload progress for live Dash readers;
    retain representative/failure artifacts before deleting temporary outputs.
    Failed executions return 1; unexpected errors propagate after finalizing
    the manifest. runner_factory/prepare_case allow isolated protocol checks.
    """
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
            "execution_model": "persistent_workers",
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
    print(f"Case/runtime options: {forwarded_args}")

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
                    # Keep workers through every full warmup and repetition for this shape.
                    with tempfile.TemporaryDirectory(prefix=f".workers_{batch_id}_", dir=profile_root) as temporary:
                        temporary_root = Path(temporary)
                        try:
                            with ProcessGroup(
                                options=options, forwarded=forwarded,
                                process_count=process_count, batch_size=batch_size,
                                group_dir=temporary_root / "startup",
                                runner_factory=runner_factory, prepare_case=prepare_case,
                            ) as group:
                                if not group.broken:
                                    write_timing_rows(
                                        timing_sink, [], batch_id=batch_id, phase="startup",
                                        repetition=0, group_wall_seconds=group.startup_wall_seconds,
                                    )
                                # Each run waits for all peers, including the warmup boundary.
                                abort_workload = False
                                for phase, count in (("warmup", options.warmups), ("measured", options.repetitions)):
                                    for repetition in range(1, count + 1):
                                        label = "Warm-up" if phase == "warmup" else "Repetition"
                                        print(
                                            f"{label} {repetition}/{count}: {process_count} process(es) "
                                            f"x batch size {batch_size}"
                                        )
                                        results, group_wall = group.run(
                                            temporary_root / f"{phase}_{repetition:03d}",
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
                                            report_failed_group(
                                                results, profile_root, batch_id, options.case_name,
                                                phase, repetition,
                                            )
                                            abort_workload = True
                                            stop_after_batch = not options.continue_on_error
                                            break
                                    if abort_workload:
                                        break
                        except BaseException as exc:
                            batch_exception = exc

                        # Preserve usable diagnostics before temporary run outputs disappear.
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
    """Run the timing CLI; return 0/1 for sweep success/failure, 130 for interruption.

    argv defaults to sys.argv. File/value errors return 2; argparse errors exit
    during parsing. Restore interrupt handling so cancellation closes workers.
    """
    # The managed Dash broker ignores SIGINT so Ctrl-C remains owned by the
    # foreground manager. Signal dispositions survive exec, so explicitly
    # restore Python's handler before doing any work. This lets Profile cancel
    # raise KeyboardInterrupt and release every active case-runner process group.
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
