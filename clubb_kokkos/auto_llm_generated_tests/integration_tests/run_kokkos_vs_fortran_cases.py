#!/usr/bin/env python3
"""Compare the Kokkos standalone with native Fortran on the JAX case set.

Adapted from human_reviewed_tests/run_jax_vs_fortran_cases.py. Both executables
come from the separate clubb_kokkos CMake build and run through run_scm.py -exe.
The enabled case names and warm-rain/SILHS variants match that harness; default
runs preserve native case durations/timesteps and use standard stats. Exact
comparison is the default (absolute and percentage limits both zero). Zero
limits require byte-identical saved NetCDF files; explicit nonzero limits use
the shared comparator's numeric criteria.

From the repository root, after building both standalone targets:
    python3 clubb_kokkos/auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py -build_dir build/clubb_kokkos
    python3 clubb_kokkos/auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py -build_dir build/clubb_kokkos -cases bomex -max_iters 5 -stats multi_col_stats

parse_args reads comparison controls and forwards shared model options to both
runs. main runs four C8 parameter columns per case/flag set, using concurrent
CPU case pairs or one GPU case pair at a time. CMake records the configured
Kokkos backends and requested/effective KOKKOS_BFB in kokkos_backends.json;
both modes are recorded in the summary. CPU workers default to half the available
logical CPUs. main writes logs, per-case bindiff JSON and a combined summary to
output/tests/kokkos_driver_test_results. A failed run, numeric/schema difference
or final directory audit returns a nonzero exit status. Limits may be overridden
explicitly for diagnosis; the summary records them. This checks saved outputs,
not every internal temporary or unexercised option combination. No JAX/F2PY
installation is needed. Only this harness's previous result tree is replaced.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path

CLUBB_ROOT = Path(__file__).resolve().parents[3]
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from input.scripts.flag_sets import (
    build_override_arg,
    format_override_value,
    get_flag_sets,
    read_flag_settings,
)  # noqa: E402


from input.scripts.create_case_namelist import parse_forwarded_args
from input.scripts.case_json_to_namelist import CASE_NAME
from run_scripts.run_scm_all import positive_int, split_values
from tuner.system_defaults import default_max_workers


@dataclass(frozen=True)
class CaseConfig:
    """Case, step limit, optional main timestep, and namelist overrides.

    None preserves the case namelist's duration or timestep, respectively.
    Radiation and statistics intervals remain at their case defaults.
    """

    case: str
    max_iters: int | None = None
    dt_main: int | None = None
    overrides: dict | None = None
    # A named variant can reuse a native case without duplicating its inputs.
    input_case: str | None = None
    # Percentage points, not fractional relative error (1e-3 means 0.001%).
    percent_threshold: float = 0.0

    @property
    def source_case(self) -> str:
        """Select the native case supplying this variant's input data."""
        return self.input_case or self.case


# Same enabled cases and physics variants as the JAX comparison. Use native
# duration/dt and exact comparison limits; JAX roundoff caps do not apply here.
DEFAULT_CASES = (
    CaseConfig("arm"),
    CaseConfig("atex"),
    CaseConfig("bomex"),
    CaseConfig("cobra"),
    CaseConfig("dycoms2_rf01"),
    CaseConfig("dycoms2_rf01_fixed_sst"),
    CaseConfig("dycoms2_rf02_nd"),
    CaseConfig("fire"),
    CaseConfig("gabls2"),
    CaseConfig("gabls3_night"),
    CaseConfig("jun25_altocu"),
    CaseConfig("neutral"),
    CaseConfig("wangara"),
    CaseConfig("dycoms2_rf02_do"),
    CaseConfig("dycoms2_rf02_ds"),
    CaseConfig(
        "lba_kk",
        input_case="lba",
        overrides={
            "microphysics_setting.microphys_scheme": '"khairoutdinov_kogan"',
            "microphysics_setting.lh_microphys_type": '"disabled"',
            "microphysics_setting.l_ice_microphys": False,
            "microphysics_setting.l_graupel": False,
        },
    ),
    CaseConfig("rico"),
    CaseConfig(
        "lba_kk_silhs",
        input_case="lba",
        overrides={
            "microphysics_setting.microphys_scheme": '"khairoutdinov_kogan"',
            "microphysics_setting.lh_microphys_type": '"interactive"',
            "microphysics_setting.l_local_kk": True,
            "microphysics_setting.l_ice_microphys": False,
            "microphysics_setting.l_graupel": False,
            "configurable_silhs_flags_nl.l_lh_deterministic_test": True,
            "configurable_silhs_flags_nl.l_lh_importance_sampling": False,
            "configurable_silhs_flags_nl.l_random_k_lh_start": False,
        },
    ),
    CaseConfig(
        "rico_silhs",
        overrides={
            "configurable_silhs_flags_nl.l_lh_deterministic_test": True,
            "configurable_silhs_flags_nl.l_lh_importance_sampling": False,
            "configurable_silhs_flags_nl.l_random_k_lh_start": False,
        },
    ),
)

RESULTS_DIRNAME = Path("output") / "tests" / "kokkos_driver_test_results"
KOKKOS_OUTPUT_DIRNAME = "kokkos_output"
FORTRAN_OUTPUT_DIRNAME = "fortran_output"
# Run logs are kept out of both output roots: run_bindiff_all.py -flag_sets treats
# every immediate child directory of a root as a flag set to compare.
LOGS_DIRNAME = "logs"
SUMMARY_FILENAME = "case_compare_summary.json"
FINAL_BINDIFF_LOG_FILENAME = "final_bindiff.log"
HR_SPEC = "C8/0.2:0.8/4"


@dataclass
class FlagData:
    """One flag set: its name, the config file it came from, and its overrides."""

    flag_name: str
    flag_dict: dict | None
    # Stored as a plain string (not a Path) so asdict() output stays JSON serializable.
    flag_file: str | None = None

    @property
    def dir_name(self) -> str:
        """Directory name identifying this flag set under both output roots."""
        if not self.flag_dict or self.flag_file is None:
            return self.flag_name
        return f"{Path(self.flag_file).stem}_{self.flag_name}"


@dataclass
class CaseResult:
    """Outcome, completed steps and comparison evidence for one case/flag pair."""

    case: str
    flag_data: FlagData
    status: str
    kokkos_rc: int
    fortran_rc: int
    bindiff_rc: int
    kokkos_elapsed_s: float
    fortran_elapsed_s: float
    elapsed_s: float
    note: str = ""
    first_failing_timestep: int | None = None
    kokkos_timesteps: int | None = None
    fortran_timesteps: int | None = None
    max_iters: int | None = None
    dt_main: int | None = None
    namelist_overrides: dict | None = None
    input_case: str | None = None
    bindiff_threshold: float = 0.0
    bindiff_percent_threshold: float = 0.0
    byte_identical: bool | None = None


def _progress_text(log_path: Path) -> str:
    """Rejoin progress text split by the native stderr completion message."""
    if not log_path.exists():
        return ""
    text = log_path.read_text(encoding="utf-8", errors="replace")
    # Native stdout/stderr share a pipe. Its final stderr write can interrupt
    # a buffered stdout record, e.g. "3 Program exited normally\n60 / 360".
    # Remove only that known message; the remaining text still has to contain
    # an actual completed-iteration record rather than a configured step cap.
    return re.sub(r"[ \t]?Program exited normally\r?\n", "", text)


def _parse_completed_timesteps(log_path: Path) -> int | None:
    """Read completed model steps, not the configured cap or saved-output count."""
    text = _progress_text(log_path)
    # Both standalones print an iteration after completing that model step.
    pattern = r"^iteration:\s*(\d+)\s*/\s*\d+\s*-- time"
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    return int(matches[-1]) if matches else None


def _parse_model_timing(log_path: Path) -> tuple[float, float] | None:
    """Read the model start time and main timestep from progress lines."""
    text = _progress_text(log_path)
    lines = re.findall(
        r"^iteration:\s*(\d+)\s*/\s*\d+\s*-- time =\s*([-+]?\d+(?:\.\d+)?)",
        text,
        flags=re.MULTILINE,
    )
    if len(lines) < 2:
        return None
    step_1, time_1 = int(lines[0][0]), float(lines[0][1])
    step_2, time_2 = int(lines[1][0]), float(lines[1][1])
    if step_2 <= step_1 or time_2 <= time_1:
        return None
    dt = (time_2 - time_1) / (step_2 - step_1)
    return time_1 - step_1 * dt, dt


def _first_failing_timestep(
    report_path: Path,
    case: str,
    total: int,
    model_timing: tuple[float, float] | None,
) -> int | None:
    """Map bindiff's first failing saved records to completed model steps."""
    try:
        files = json.loads(report_path.read_text(encoding="utf-8"))["cases"][case][
            "files"
        ]
    except (OSError, ValueError, KeyError, TypeError):
        return None

    first_failure = None
    for file in files:
        if not isinstance(file, dict):
            continue
        entry = file.get("first_failing_prefix")
        if not isinstance(entry, dict):
            continue
        try:
            record = entry["record"]
            n_saved = entry["saved_records"]
            if not (
                isinstance(record, int)
                and isinstance(n_saved, int)
                and 0 <= record < n_saved <= total
            ):
                continue
            if n_saved == total:
                step = record + 1
            else:
                if model_timing is None or entry["output_time"] is None:
                    continue
                start_time, dt = model_timing
                step_value = (entry["output_time"] - start_time) / dt
                if not math.isfinite(step_value):
                    continue
                step = round(step_value)
                if abs(step_value - step) > 1e-3:
                    continue
            if 1 <= step <= total:
                first_failure = (
                    step if first_failure is None else min(first_failure, step)
                )
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            continue
    return first_failure


class ComparisonCancelled(Exception):
    """Stop the current command and do not launch the rest of its case pair."""


class RunSupervisor:
    """Start each command, wait for it, and stop any programs it starts.

    For a model run, this script starts run_scm.py, which starts the model.
    The supervisor keeps them in one group so that interrupting this script
    stops both, rather than leaving the model running in the background.
    """

    def __init__(self):
        """Keep cancellation state shared by the signal handler and workers."""
        self.stopping = threading.Event()
        self.signal_number = None

    def handle_signal(self, signum, _frame):
        """Record an interruption without raising during child-process creation."""
        # Do not raise inside Popen: a signal can arrive between creating a
        # child and recording its handle. Polling lets that launch finish and
        # then cleans up its entire group, even on repeated Ctrl-C.
        self.signal_number = signum

    def check_cancelled(self):
        """Stop before launching another command when cancellation was requested."""
        if self.signal_number is not None or self.stopping.is_set():
            raise ComparisonCancelled()

    @staticmethod
    def _stop_group(proc):
        """Terminate and reap the launcher and its model descendants."""
        # The launcher may already have exited while its model is still alive.
        # Always address the group rather than checking only proc.poll().
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            proc.wait()
            return
        deadline = time.monotonic() + 2.0
        while time.monotonic() < deadline:
            proc.poll()  # Reap the direct child while descendants shut down.
            try:
                os.killpg(proc.pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.05)
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()

    def run_and_log(self, cmd: list[str], cwd: Path, log_path: Path) -> int:
        """Run one process group, saving its exact command, output and exit code."""
        self.check_cancelled()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log:
            log.write("$ " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")
            log.flush()
            proc = subprocess.Popen(
                cmd,
                cwd=str(cwd),
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                while True:
                    self.check_cancelled()
                    try:
                        returncode = proc.wait(timeout=0.2)
                        break
                    except subprocess.TimeoutExpired:
                        continue
                self.check_cancelled()
                log.write(f"\n[exit_code] {returncode}\n")
                return returncode
            finally:
                self._stop_group(proc)


def _acquire_run_lock(repo_root: Path):
    """Prevent concurrent harnesses from sharing and deleting result paths."""
    lock_path = repo_root / "output" / "tests" / ".run_kokkos_vs_fortran_cases.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_file = lock_path.open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock_file.close()
        print(
            "ERROR: another run_kokkos_vs_fortran_cases.py process is already "
            f"using {repo_root / RESULTS_DIRNAME}."
        )
        return None

    lock_file.write(f"pid={os.getpid()}\n")
    lock_file.flush()
    return lock_file


def _tail(path: Path, n: int = 25) -> str:
    """Return a short diagnostic excerpt from a failed model's log."""
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-n:])


@dataclass
class TaskCtx:
    """Executables, shared options and output ownership for one comparison."""

    config: CaseConfig
    flag_data: FlagData
    repo_root: Path
    run_scm_args: list[str]
    bindiff_threshold: float
    bindiff_percent_threshold: float
    run_output_root: Path
    supervisor: RunSupervisor
    kokkos_exe: Path
    fortran_exe: Path

    @property
    def output_group(self) -> str:
        """Isolate variants that reuse the same native case and NetCDF filename."""
        # Variants share the native case's namelist and NetCDF filenames. Give
        # them separate immediate directories, also understood by -flag_sets.
        suffix = f"__{self.config.case}" if self.config.input_case else ""
        return self.flag_data.dir_name + suffix


def _run_case_w_flags(task: TaskCtx) -> CaseResult:
    """Run both executables with identical inputs and compare saved outputs."""
    run_scm = task.repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = task.repo_root / "run_scripts" / "run_bindiff_all.py"

    start = time.time()

    common_args = [str(run_scm), "-multicol", HR_SPEC]
    if (
        task.config.max_iters is not None
        and _forwarded_value(task.run_scm_args, "-max_iters") is None
    ):
        common_args += ["-max_iters", str(task.config.max_iters)]
    if (
        task.config.dt_main is not None
        and _forwarded_value(task.run_scm_args, "-dt_main") is None
    ):
        common_args += ["-dt_main", str(task.config.dt_main)]

    # Separate flag sets/variants on both sides for the final directory audit.
    flag_dir_name = task.output_group
    kokkos_run_out_dir = task.run_output_root / KOKKOS_OUTPUT_DIRNAME / flag_dir_name
    f90_run_out_dir = task.run_output_root / FORTRAN_OUTPUT_DIRNAME / flag_dir_name
    kokkos_run_out_dir.mkdir(parents=True, exist_ok=True)
    f90_run_out_dir.mkdir(parents=True, exist_ok=True)

    # Explicit flag-set overrides take precedence over the curated case setup.
    overrides = {**(task.config.overrides or {}), **(task.flag_data.flag_dict or {})}
    override_arg = build_override_arg(overrides)
    if override_arg is not None:
        common_args += ["-override", override_arg]

    # User overrides apply last, after the curated case and flag-set settings.
    common_args += task.run_scm_args

    kokkos_cmd = [
        sys.executable,
        *common_args,
        "-exe",
        str(task.kokkos_exe),
        "-output_dir",
        str(kokkos_run_out_dir),
    ]
    f90_cmd = [
        sys.executable,
        *common_args,
        "-exe",
        str(task.fortran_exe),
        "-output_dir",
        str(f90_run_out_dir),
    ]

    kokkos_cmd.append(task.config.source_case)
    f90_cmd.append(task.config.source_case)

    log_dir = task.run_output_root / LOGS_DIRNAME / flag_dir_name
    kokkos_log = log_dir / f"{task.config.case}_run_kokkos.log"
    f90_log = log_dir / f"{task.config.case}_run_fortran.log"
    diff_log = log_dir / f"{task.config.case}_bindiff.log"
    diff_report = log_dir / f"{task.config.case}_bindiff.json"

    kokkos_start = time.time()
    kokkos_rc = task.supervisor.run_and_log(kokkos_cmd, task.repo_root, kokkos_log)
    kokkos_elapsed = time.time() - kokkos_start

    # Fortran runs even when Kokkos failed: a flag set that breaks both drivers is a
    # configuration problem, while one that breaks only Kokkos is a missing Kokkos feature.
    f90_start = time.time()
    f90_rc = task.supervisor.run_and_log(f90_cmd, task.repo_root, f90_log)
    f90_elapsed = time.time() - f90_start
    kokkos_timesteps = (
        _parse_completed_timesteps(kokkos_log) if kokkos_rc == 0 else None
    )
    fortran_timesteps = _parse_completed_timesteps(f90_log) if f90_rc == 0 else None

    if kokkos_rc != 0 or f90_rc != 0:
        if kokkos_rc != 0 and f90_rc != 0:
            status = "both_failed"
            note = f"{_tail(kokkos_log)}\n--- fortran ---\n{_tail(f90_log)}"
        elif kokkos_rc != 0:
            status = "kokkos_failed"
            note = _tail(kokkos_log)
        else:
            status = "fortran_failed"
            note = _tail(f90_log)
        return CaseResult(
            case=task.config.case,
            input_case=task.config.input_case,
            bindiff_threshold=task.bindiff_threshold,
            bindiff_percent_threshold=task.bindiff_percent_threshold,
            max_iters=task.config.max_iters,
            dt_main=task.config.dt_main,
            namelist_overrides=overrides or None,
            status=status,
            kokkos_rc=kokkos_rc,
            fortran_rc=f90_rc,
            bindiff_rc=-1,
            kokkos_elapsed_s=kokkos_elapsed,
            fortran_elapsed_s=f90_elapsed,
            elapsed_s=time.time() - start,
            flag_data=task.flag_data,
            note=note,
            kokkos_timesteps=kokkos_timesteps,
            fortran_timesteps=fortran_timesteps,
        )

    diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-verbose",
        "2",
        "-strict",
        "-case",
        task.config.source_case,
        "-threshold",
        str(task.bindiff_threshold),
        "-percent_threshold",
        str(task.bindiff_percent_threshold),
        "-result_json",
        str(diff_report),
        str(kokkos_run_out_dir),
        str(f90_run_out_dir),
    ]
    diff_report.unlink(missing_ok=True)
    diff_rc = task.supervisor.run_and_log(diff_cmd, task.repo_root, diff_log)

    # Matched runs show their completed steps directly. For a numerical diff,
    # map bindiff's saved-prefix result to the model step shown in the table.
    first_failing_timestep = None
    if (
        diff_rc != 0
        and kokkos_timesteps is not None
        and kokkos_timesteps == fortran_timesteps
    ):
        first_failing_timestep = _first_failing_timestep(
            diff_report,
            task.config.source_case,
            kokkos_timesteps,
            _parse_model_timing(f90_log),
        )

    status = "match" if diff_rc == 0 else "diff"
    note = f"Bindiff report: {diff_report}"
    if diff_rc == 0 and (
        kokkos_timesteps is None or kokkos_timesteps != fortran_timesteps
    ):
        status = "step_mismatch"
        note += (
            f"; completed steps disagree or are missing: Kokkos={kokkos_timesteps}, "
            f"Fortran={fortran_timesteps}"
        )
    # Consume the comparison owner's structured evidence, rather than parsing
    # prose or equating a tolerance match with byte-identical output.
    byte_identical = None
    try:
        files = json.loads(diff_report.read_text())["cases"][task.config.source_case][
            "files"
        ]
        if not files:
            raise ValueError("No compared files in the structured report")
        byte_identical = all(file["comparison"] == "byte_identical" for file in files)
        # The shared comparator clips small fields when computing percentages.
        # A numeric match alone therefore cannot establish the zero-limit BFB gate.
        if (
            status == "match"
            and task.bindiff_threshold == 0
            and task.bindiff_percent_threshold == 0
            and not byte_identical
        ):
            status = "diff"
            note += "; zero-limit BFB check failed: saved files are not byte-identical"
    except (OSError, ValueError, KeyError, TypeError):
        if diff_rc == 0:
            status = "diff"
            note += "; missing or invalid structured comparison report"
    return CaseResult(
        case=task.config.case,
        input_case=task.config.input_case,
        bindiff_threshold=task.bindiff_threshold,
        bindiff_percent_threshold=task.bindiff_percent_threshold,
        max_iters=task.config.max_iters,
        dt_main=task.config.dt_main,
        namelist_overrides=overrides or None,
        flag_data=task.flag_data,
        status=status,
        kokkos_rc=kokkos_rc,
        fortran_rc=f90_rc,
        bindiff_rc=diff_rc,
        kokkos_elapsed_s=kokkos_elapsed,
        fortran_elapsed_s=f90_elapsed,
        elapsed_s=time.time() - start,
        note=note,
        byte_identical=byte_identical,
        first_failing_timestep=first_failing_timestep,
        kokkos_timesteps=kokkos_timesteps,
        fortran_timesteps=fortran_timesteps,
    )


def _forwarded_value(arguments: list[str], option: str) -> str | None:
    """Read the last value of a forwarded run_scm option."""
    value = None
    for i, argument in enumerate(arguments):
        if argument == option and i + 1 < len(arguments):
            value = arguments[i + 1]
        elif argument.startswith(option + "="):
            value = argument.split("=", 1)[1]
    return value


def parse_args(argv=None) -> argparse.Namespace:
    """Parse harness controls and retain the shared model options for both runs."""
    parser = argparse.ArgumentParser(
        description=(
            "Run SCM cases with both the Kokkos driver and Fortran standalone, "
            "then compare outputs using run_bindiff_all.py. Other single-dash "
            "run_scm.py options are forwarded to both runs."
        ),
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        "-build_dir",
        type=Path,
        default=CLUBB_ROOT / "build" / "clubb_kokkos",
        help="Separate CMake build containing clubb_standalone and fortran/clubb_standalone (default: build/clubb_kokkos).",
    )
    parser.add_argument(
        "-workers",
        dest="jobs",
        type=positive_int,
        default=None,
        help="Concurrent case pairs; CPU default is half the available logical CPUs, GPU requires 1.",
        metavar="N",
    )
    parser.add_argument(
        "-bindiff_verbose",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="Verbosity level for the final combined run_bindiff_all.py run.",
    )
    parser.add_argument(
        "-bindiff_threshold",
        type=float,
        default=None,
        help="Override both absolute and percentage thresholds (default: zero for both; exact numeric agreement).",
    )
    parser.add_argument(
        "-bindiff_percent_threshold",
        type=float,
        default=None,
        help="Override only the percentage threshold, in percentage points; takes precedence over -bindiff_threshold.",
    )
    parser.add_argument(
        "-cases",
        nargs="+",
        default=None,
        help="Case names to run (default is the curated supported set).",
    )
    parser.add_argument(
        "-flag_config_file",
        type=str,
        default=None,
        help="JSON file describing alternate flag settings.",
    )
    args, forwarded = parse_forwarded_args(parser, argv)
    args.cases = split_values(args.cases) if args.cases is not None else None
    if args.cases == []:
        parser.error("-cases must contain at least one case name")
    if args.cases and len(args.cases) != len(set(args.cases)):
        parser.error("-cases must not repeat names that share output paths")
    if args.cases and any(not CASE_NAME.fullmatch(case) for case in args.cases):
        parser.error("-cases accepts case names, not paths")
    for option in ("bindiff_threshold", "bindiff_percent_threshold"):
        value = getattr(args, option)
        if value is not None and (not math.isfinite(value) or value < 0):
            parser.error(f"-{option} must be finite and nonnegative")
    # These options belong to the comparison, not to either individual model.
    reserved = {
        "-jax",
        "-python",
        "-exe",
        "-driver_test",
        "-gdb",
        "-output_dir",
        "-multicol",
        "-install_dir",
    }
    for argument in forwarded:
        option = argument.split("=", 1)[0]
        if option.startswith("--"):
            parser.error(f"Use single-dash run_scm.py options: {option}")
        if option in reserved:
            parser.error(f"{option} is controlled by the comparison harness")
    for option in ("-max_iters", "-dt_main"):
        value = _forwarded_value(forwarded, option)
        if value is not None:
            try:
                positive = int(value) > 0
            except ValueError:
                positive = False
            if not positive:
                parser.error(f"{option} must be a positive integer")
    backend_path = args.build_dir.expanduser() / "kokkos_backends.json"
    try:
        build_config = json.loads(backend_path.read_text(encoding="utf-8"))
        backends = build_config["devices"]
        if (
            not isinstance(backends, list)
            or not backends
            or any(not isinstance(device, str) or not device for device in backends)
        ):
            raise ValueError("devices must be a nonempty list of backend names")
        for key in ("kokkos_bfb_requested", "kokkos_bfb"):
            if not isinstance(build_config[key], bool):
                raise ValueError(f"{key} must be a boolean")
        if build_config["kokkos_bfb"] and not build_config["kokkos_bfb_requested"]:
            raise ValueError("effective KOKKOS_BFB cannot be true unless requested")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(
            f"Cannot read configured Kokkos backends from {backend_path}: {error}. "
            "Reconfigure the selected clubb_kokkos CMake build."
        )
    cpu_backends = {"SERIAL", "OPENMP", "THREADS", "HPX"}
    cpu_comparison = all(device in cpu_backends for device in backends)
    if args.jobs is None:
        args.jobs = default_max_workers() if cpu_comparison else 1
    if not cpu_comparison and args.jobs != 1:
        parser.error("GPU or unknown Kokkos backends require -workers 1")
    args.kokkos_backends = backends
    args.kokkos_bfb_requested = build_config["kokkos_bfb_requested"]
    args.kokkos_bfb = build_config["kokkos_bfb"]
    args.run_scm_args = forwarded
    return args


def main() -> int:
    """Run the comparison matrix and stop child process groups on interruption."""
    args = parse_args()

    supervisor = RunSupervisor()
    previous_handlers = {
        sig: signal.signal(sig, supervisor.handle_signal)
        for sig in (signal.SIGINT, signal.SIGTERM)
    }
    try:
        return _run_comparisons(args, supervisor)
    except ComparisonCancelled:
        print(
            "Comparison interrupted; stopped all active case process groups.",
            flush=True,
        )
        return 128 + (supervisor.signal_number or signal.SIGINT)
    finally:
        supervisor.stopping.set()
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)


def _run_tasks(tasks, jobs, emit, supervisor):
    """Supervise external case pairs and report each completed result."""
    # Workers only supervise external programs; processes add no compute
    # parallelism here and can orphan grandchildren when Pool terminates them.
    executor = ThreadPoolExecutor(max_workers=min(jobs, len(tasks)))
    try:
        futures = [executor.submit(_run_case_w_flags, task) for task in tasks]
        for future in as_completed(futures):
            supervisor.check_cancelled()
            emit(future.result())
    except BaseException:
        supervisor.stopping.set()
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def _run_comparisons(args, supervisor) -> int:
    """Validate inputs, replace this harness's results and write the final audit."""
    defaults = {config.case: config for config in DEFAULT_CASES}
    cases = list(args.cases) if args.cases else list(defaults)
    max_iters = _forwarded_value(args.run_scm_args, "-max_iters")
    dt_main = _forwarded_value(args.run_scm_args, "-dt_main")
    configs = []
    for case in cases:
        config = defaults.get(case, CaseConfig(case))
        if max_iters is not None:
            config = replace(config, max_iters=int(max_iters))
        if dt_main is not None:
            config = replace(config, dt_main=int(dt_main))
        configs.append(config)

    repo_root = CLUBB_ROOT
    build_dir = args.build_dir.expanduser().resolve()
    kokkos_exe = build_dir / "clubb_standalone"
    fortran_exe = build_dir / "fortran" / "clubb_standalone"
    for executable in (kokkos_exe, fortran_exe):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            print(f"ERROR: executable not found or not executable: {executable}")
            print(
                "Build both with: "
                + shlex.join(
                    [
                        "cmake",
                        "--build",
                        str(build_dir),
                        "--target",
                        "clubb_kokkos_standalone",
                        "clubb_standalone",
                    ]
                )
            )
            return 2
    run_lock = _acquire_run_lock(repo_root)
    if run_lock is None:
        return 2
    results_root = (repo_root / RESULTS_DIRNAME).resolve()
    kokkos_output_root = results_root / KOKKOS_OUTPUT_DIRNAME
    f90_output_root = results_root / FORTRAN_OUTPUT_DIRNAME
    run_scm = repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = repo_root / "run_scripts" / "run_bindiff_all.py"
    if not run_scm.exists() or not run_bindiff.exists():
        print("ERROR: could not locate run_scripts/run_scm.py or run_bindiff_all.py")
        return 2

    # get_flag_sets injects the unmodified "default" run and rejects a JSON flag set
    # that would collide with that name. Validate before the previous results are
    # deleted, so a bad config file doesn't cost the last run's output.
    try:
        flag_config = (
            read_flag_settings(args.flag_config_file) if args.flag_config_file else {}
        )
        flag_sets = get_flag_sets(False, flag_config)
    except (OSError, ValueError) as exc:
        print(f"ERROR: could not load flag sets: {exc}")
        return 2
    flag_data_list = [
        FlagData(flag_name=name, flag_dict=overrides, flag_file=args.flag_config_file)
        for name, overrides in flag_sets.items()
    ]

    tasks = []
    for flag_data in flag_data_list:
        for config in configs:
            tasks.append(
                TaskCtx(
                    config=config,
                    flag_data=flag_data,
                    repo_root=repo_root,
                    run_scm_args=args.run_scm_args,
                    bindiff_threshold=(
                        args.bindiff_threshold
                        if args.bindiff_threshold is not None
                        else 0.0
                    ),
                    bindiff_percent_threshold=(
                        args.bindiff_percent_threshold
                        if args.bindiff_percent_threshold is not None
                        else (
                            args.bindiff_threshold
                            if args.bindiff_threshold is not None
                            else config.percent_threshold
                        )
                    ),
                    run_output_root=results_root,
                    supervisor=supervisor,
                    kokkos_exe=kokkos_exe,
                    fortran_exe=fortran_exe,
                )
            )

    output_owners = set()
    for task in tasks:
        if Path(task.output_group).name != task.output_group or task.output_group in {
            ".",
            "..",
        }:
            print(
                f"ERROR: flag set must name one output directory: {task.output_group}"
            )
            return 2
        key = (task.output_group, task.config.source_case)
        if key in output_owners:
            print(f"ERROR: duplicate output destination: {key}")
            return 2
        output_owners.add(key)
    if results_root.exists():
        shutil.rmtree(results_root)

    print("\nRun settings:")
    print(f"  Kokkos executable: {kokkos_exe}")
    print(f"  Fortran executable: {fortran_exe}")
    print(f"  Output: {results_root}")
    print(
        f"  Statistics: {_forwarded_value(args.run_scm_args, '-stats') or 'standard'}"
    )
    print(f"  Debug: {_forwarded_value(args.run_scm_args, '-debug') or 'case default'}")
    if args.run_scm_args:
        print(f"  Forwarded to both runs: {shlex.join(args.run_scm_args)}")
    print("  Columns: 4")
    print(f"  Workers: {min(args.jobs, len(tasks))}")
    default_percent = (
        args.bindiff_percent_threshold
        if args.bindiff_percent_threshold is not None
        else args.bindiff_threshold if args.bindiff_threshold is not None else 0.0
    )
    print(
        f"  Comparison limits: absolute {tasks[0].bindiff_threshold:g}; "
        f"percentage {default_percent:g}%"
    )

    print("  Kokkos backends: " + ", ".join(args.kokkos_backends))
    print(
        f"  KOKKOS_BFB: requested {args.kokkos_bfb_requested}, "
        f"effective {args.kokkos_bfb}; Cholesky "
        + ("Fortran/f2k" if args.kokkos_bfb else "custom")
    )
    print(f"  Concurrent case pairs: {min(args.jobs, len(tasks))}")
    print("\nCases and flag sets:")
    case_label = "case" if len(cases) == 1 else "cases"
    flag_label = "flag set" if len(flag_data_list) == 1 else "flag sets"
    pair_label = "pair" if len(tasks) == 1 else "pairs"
    print(
        f"  {len(cases)} {case_label} x {len(flag_data_list)} {flag_label}"
        f" = {len(tasks)} Kokkos/Fortran {pair_label}"
    )
    print("  Flag sets:")
    for flag_data in flag_data_list:
        print(f"    {flag_data.dir_name}")
        for setting, value in (flag_data.flag_dict or {}).items():
            print(f"      - {setting} = {format_override_value(value)}")
    print("  Cases:")
    percent_by_case = {
        task.config.case: task.bindiff_percent_threshold for task in tasks
    }
    for config in configs:
        iterations = (
            f"{config.max_iters} iterations"
            if config.max_iters is not None
            else "native duration"
        )
        print(f"    {config.case} ({iterations})")
        if percent_by_case[config.case] != default_percent:
            print(f"      - percentage limit: {percent_by_case[config.case]:g}%")
        if config.input_case:
            print(f"      - input case: {config.input_case}")
        if config.dt_main is not None:
            print(f"      - dt_main: {config.dt_main} s")
        for setting, value in (config.overrides or {}).items():
            print(f"      - {setting} = {format_override_value(value)}")

    start = time.time()
    results: list[CaseResult] = []
    flag_width = max(24, max(len(fd.dir_name) for fd in flag_data_list))
    case_width = max(22, max(len(config.case) for config in configs))
    table_header = (
        f"[{'Status':14}] {'Flag set':{flag_width}} {'Case':{case_width}} "
        f"{'Timesteps':>13} {'Kokkos (s)':>10} {'Fortran (s)':>12} {'Total (s)':>10}"
    )
    print(f"\n{table_header}", flush=True)
    print("-" * len(table_header), flush=True)

    def _emit(result: CaseResult) -> None:
        """Print the completed steps or first failing step with each result."""
        results.append(result)
        total_steps = result.fortran_timesteps or result.kokkos_timesteps
        if (
            result.status == "match"
            and total_steps is not None
            and result.kokkos_timesteps == result.fortran_timesteps
        ):
            steps = f"{total_steps} / {total_steps}"
        elif result.status == "diff" and result.first_failing_timestep is not None:
            steps = f"{result.first_failing_timestep} / {total_steps}"
        else:
            steps = "-"
        print(
            f"[{result.status:14}] {result.flag_data.dir_name:{flag_width}} {result.case:{case_width}} "
            f"{steps:>13} {result.kokkos_elapsed_s:10.1f} {result.fortran_elapsed_s:12.1f} "
            f"{result.elapsed_s:10.1f}",
            flush=True,
        )
        if result.status != "match":
            print(f"  {result.note}", flush=True)

    _run_tasks(tasks, args.jobs, _emit, supervisor)

    results.sort(key=lambda r: (r.flag_data.dir_name, r.case))
    statuses = (
        "match",
        "diff",
        "step_mismatch",
        "kokkos_failed",
        "fortran_failed",
        "both_failed",
    )

    # The combined audit also detects unmatched flag-set directories/files.
    # It cannot override a failed per-case comparison.
    final_bindiff_log = results_root / FINAL_BINDIFF_LOG_FILENAME
    final_diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-flag_sets",
        "-strict",
        "-verbose",
        str(args.bindiff_verbose),
        "-threshold",
        str(tasks[0].bindiff_threshold),
        "-percent_threshold",
        str(max(task.bindiff_percent_threshold for task in tasks)),
        str(kokkos_output_root),
        str(f90_output_root),
    ]
    final_bindiff_rc = supervisor.run_and_log(
        final_diff_cmd, repo_root, final_bindiff_log
    )

    summary = {
        "total_cases": len(results),
        "workers": min(args.jobs, len(tasks)),
        "kokkos_backends": args.kokkos_backends,
        "kokkos_bfb_requested": args.kokkos_bfb_requested,
        "kokkos_bfb": args.kokkos_bfb,
        "executables": {"kokkos": str(kokkos_exe), "fortran": str(fortran_exe)},
        "byte_identical": sum(r.byte_identical is True for r in results),
        **{status: sum(r.status == status for r in results) for status in statuses},
        "elapsed_s": time.time() - start,
        "run_scm_args": args.run_scm_args,
        "final_bindiff_rc": final_bindiff_rc,
        "flag_sets": {
            flag_data.dir_name: {
                status: sum(
                    r.status == status and r.flag_data.dir_name == flag_data.dir_name
                    for r in results
                )
                for status in statuses
            }
            for flag_data in flag_data_list
        },
        "cases": [asdict(r) for r in results],
    }

    summary_path = results_root / SUMMARY_FILENAME
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def _print_count_table(title: str, counts: dict) -> None:
        """Summarize numeric differences separately from model/step failures."""
        failed = sum(
            counts[status]
            for status in (
                "step_mismatch",
                "kokkos_failed",
                "fortran_failed",
                "both_failed",
            )
        )
        print(f"\n{title}:")
        print(f"  {'Match':>5} {'Diff':>4} {'Failed':>6}")
        print(f"  {'-----':>5} {'----':>4} {'------':>6}")
        print(f"  {counts['match']:>5} {counts['diff']:>4} {failed:>6}")

    if len(flag_data_list) > 1:
        _print_count_table("Overall cases", summary)
    for flag_data in flag_data_list:
        _print_count_table(
            f"Flag set {flag_data.dir_name}", summary["flag_sets"][flag_data.dir_name]
        )

    case_failure = any(summary[status] > 0 for status in statuses if status != "match")
    if final_bindiff_rc != 0 and not case_failure:
        print("\nFinal bindiff audit failed despite all case comparisons matching.")

    hours, remainder = divmod(summary["elapsed_s"], 3600)
    minutes, seconds = divmod(remainder, 60)
    duration = (
        f"{int(hours)}h {int(minutes)}m {seconds:.1f}s"
        if hours >= 1
        else f"{int(minutes)}m {seconds:.1f}s" if minutes >= 1 else f"{seconds:.1f}s"
    )
    audit_status = (
        "passed" if final_bindiff_rc == 0 else f"failed (exit {final_bindiff_rc})"
    )
    print(f"\nTotal time: {duration}")
    print(f"Byte-identical case pairs: {summary['byte_identical']} / {len(results)}")
    print(f"Results JSON: {summary_path}")
    print(f"Bindiff log: {final_bindiff_log} [{audit_status}]")

    return 1 if case_failure or final_bindiff_rc != 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
