#!/usr/bin/env python3
"""Run selected SCM cases with both JAX and Fortran drivers and diff outputs."""

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

CLUBB_ROOT = Path(__file__).resolve().parents[1]
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from clubb_jax.run_jax import ensure_environment  # noqa: E402
from utilities.flag_sets import build_override_arg, format_override_value, get_flag_sets, read_flag_settings  # noqa: E402


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
    percent_threshold: float = 1.0e-7

    @property
    def source_case(self) -> str:
        return self.input_case or self.case


# Curated JAX coverage, including supported microphysics. Unsupported driver
# features (e.g. SILHS) remain excluded.
DEFAULT_CASES = (
    # --- Cases without microphysics ---
    # Turbulence, forcing and radiation coverage at the standard tolerances.
    CaseConfig("arm", 360),       # diffs expected after ~600 60s-timesteps
    CaseConfig("atex", 360),      # diffs expected after ~400 60s-timesteps
    CaseConfig("bomex"),
    CaseConfig("cobra", 360),     # very stable, limited for speed
    CaseConfig("dycoms2_rf01"),
    CaseConfig("dycoms2_rf01_fixed_sst", 300),  # diffs after ~360 steps with l_diag_Lscale_from_tau=.false.
    CaseConfig("dycoms2_rf02_nd"),
    CaseConfig("fire"),
    CaseConfig("gabls2", 360),    # very stable, limited for speed
    CaseConfig("gabls3_night", 360),  # very stable, limited for speed
    CaseConfig("jun25_altocu", 180),  # diffs expected after ~200 60s-timesteps
    CaseConfig("neutral"),
    CaseConfig("wangara"),

    # --- KK microphysics ---
    # Warm-rain coverage with the standard absolute/percentage tolerances.
    CaseConfig("dycoms2_rf02_do"),  # Full native 360 steps at 60s, including every saved prefix.
    CaseConfig("dycoms2_rf02_ds", 240),  # First cumulative failure at step 242; late roundoff amplification.
    # Warm-rain KK variant of the native Morrison/SILHS case; native dt is 60s.
    CaseConfig("lba_kk", 360, input_case="lba", overrides={
        "microphysics_setting.microphys_scheme": '"khairoutdinov_kogan"',
        "microphysics_setting.lh_microphys_type": '"disabled"',
        "microphysics_setting.l_ice_microphys": False,
        "microphysics_setting.l_graupel": False,
    }),
    CaseConfig("rico", 380, 60),  # First failing saved prefix: 385 steps; Fortran O0/O2: 390.

    # --- Morrison microphysics ---
    # The float32 core is sensitive to FMA/intermediate rounding: Fortran debug
    # versus release builds reproduce these passing limits. Keep absolute 1e-7,
    # but allow 1e-3 percent (0.001%). Short runs cover cloud processes and rain/ice
    # startup; they do not establish long-term or graupel accuracy.
    CaseConfig("lba", 93, percent_threshold=1.0e-3, overrides={
        "microphysics_setting.lh_microphys_type": '"disabled"',
    }),
    CaseConfig("clex9_oct14", 73, percent_threshold=1.0e-3),  # 13 active rain/ice steps.
    CaseConfig("nov11_altocu", 62, percent_threshold=1.0e-3),  # 2 active rain/ice steps.
)

RESULTS_DIRNAME = Path("output") / "tests" / "jax_driver_test_results"
JAX_OUTPUT_DIRNAME = "jax_output"
FORTRAN_OUTPUT_DIRNAME = "fortran_output"
# Run logs are kept out of both output roots: run_bindiff_all.py --flag-sets treats
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
    case: str
    flag_data: FlagData
    status: str
    jax_rc: int
    fortran_rc: int
    bindiff_rc: int
    jax_elapsed_s: float
    fortran_elapsed_s: float
    elapsed_s: float
    note: str = ""
    first_failing_timestep: int | None = None
    jax_timesteps: int | None = None
    fortran_timesteps: int | None = None
    max_iters: int | None = None
    dt_main: int | None = None
    namelist_overrides: dict | None = None
    input_case: str | None = None
    bindiff_threshold: float = 1.0e-7
    bindiff_percent_threshold: float = 1.0e-7


def _parse_completed_timesteps(log_path: Path, *, jax: bool) -> int | None:
    """Read completed model steps, not the configured cap or saved-output count."""
    if not log_path.exists():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    # JAX reports completion explicitly. Fortran prints each iteration after
    # advancing it, so its last progress line gives the completed count.
    pattern = (r"^Completed (\d+) timesteps\b" if jax else
               r"^iteration:\s*(\d+)\s*/\s*\d+\s*-- time")
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    return int(matches[-1]) if matches else None


def _parse_model_timing(log_path: Path) -> tuple[float, float] | None:
    """Read the model start time and main timestep from progress lines."""
    if not log_path.exists():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    lines = re.findall(
        r"^iteration:\s*(\d+)\s*/\s*\d+\s*-- time =\s*([-+]?\d+(?:\.\d+)?)",
        text, flags=re.MULTILINE,
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
    report_path: Path, case: str, total: int,
    model_timing: tuple[float, float] | None,
) -> int | None:
    """Map bindiff's first failing saved records to completed model steps."""
    try:
        files = json.loads(report_path.read_text(encoding="utf-8"))["cases"][case]["files"]
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
            if not (isinstance(record, int) and isinstance(n_saved, int)
                    and 0 <= record < n_saved <= total):
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
                first_failure = step if first_failure is None else min(first_failure, step)
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
        self.stopping = threading.Event()
        self.signal_number = None

    def handle_signal(self, signum, _frame):
        # Do not raise inside Popen: a signal can arrive between creating a
        # child and recording its handle. Polling lets that launch finish and
        # then cleans up its entire group, even on repeated Ctrl-C.
        self.signal_number = signum

    def check_cancelled(self):
        if self.signal_number is not None or self.stopping.is_set():
            raise ComparisonCancelled()

    @staticmethod
    def _stop_group(proc):
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
        self.check_cancelled()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log:
            log.write("$ " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")
            log.flush()
            proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True)
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
    lock_path = repo_root / "output" / "tests" / ".run_jax_vs_fortran_cases.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_file = lock_path.open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock_file.close()
        print(
            "ERROR: another run_jax_vs_fortran_cases.py process is already "
            f"using {repo_root / RESULTS_DIRNAME}."
        )
        return None

    lock_file.write(f"pid={os.getpid()}\n")
    lock_file.flush()
    return lock_file


def _tail(path: Path, n: int = 25) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-n:])

@dataclass
class TaskCtx:
    config: CaseConfig
    flag_data: FlagData
    repo_root: Path
    run_scm_args: list[str]
    bindiff_threshold: float
    bindiff_percent_threshold: float
    run_output_root: Path
    supervisor: RunSupervisor

    @property
    def output_group(self) -> str:
        # Variants share the native case's namelist and NetCDF filenames. Give
        # them separate immediate directories, also understood by --flag-sets.
        suffix = f"__{self.config.case}" if self.config.input_case else ""
        return self.flag_data.dir_name + suffix


def _run_case_w_flags(task: TaskCtx) -> CaseResult:
    run_scm = task.repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = task.repo_root / "run_scripts" / "run_bindiff_all.py"

    start = time.time()

    common_args = [str(run_scm), "-multicol", HR_SPEC]
    if task.config.max_iters is not None and _forwarded_value(task.run_scm_args, "-max_iters") is None:
        common_args += ["-max_iters", str(task.config.max_iters)]
    if task.config.dt_main is not None and _forwarded_value(task.run_scm_args, "-dt_main") is None:
        common_args += ["-dt_main", str(task.config.dt_main)]
    common_args += task.run_scm_args

    #A subdirectory for each flagset is created to store the .nc output from that flagset.
    # The output directory will look like this:
    # output/tests/jax_driver_test_results/jax_output/<file_name + flag_name> and output/tests/jax_driver_test_results/fortran_output/<file_name + flag_name>
    flag_dir_name = task.output_group
    jax_run_out_dir = task.run_output_root / JAX_OUTPUT_DIRNAME / flag_dir_name
    f90_run_out_dir = task.run_output_root / FORTRAN_OUTPUT_DIRNAME / flag_dir_name
    jax_run_out_dir.mkdir(parents=True, exist_ok=True)
    f90_run_out_dir.mkdir(parents=True, exist_ok=True)

    jax_cmd = [sys.executable, *common_args, "-jax", "-out_dir", str(jax_run_out_dir)]
    f90_cmd = [sys.executable, *common_args, "-out_dir", str(f90_run_out_dir)]

    # Explicit flag-set overrides take precedence over the curated case setup.
    overrides = {**(task.config.overrides or {}), **(task.flag_data.flag_dict or {})}
    override_arg = build_override_arg(overrides)
    if override_arg is not None:
        jax_cmd += ["-override", override_arg]
        f90_cmd += ["-override", override_arg]

    jax_cmd.append(task.config.source_case)
    f90_cmd.append(task.config.source_case)

    log_dir = task.run_output_root / LOGS_DIRNAME / flag_dir_name
    jax_log = log_dir / f"{task.config.case}_run_jax.log"
    f90_log = log_dir / f"{task.config.case}_run_fortran.log"
    diff_log = log_dir / f"{task.config.case}_bindiff.log"
    diff_report = log_dir / f"{task.config.case}_bindiff.json"

    jax_start = time.time()
    jax_rc = task.supervisor.run_and_log(jax_cmd, task.repo_root, jax_log)
    jax_elapsed = time.time() - jax_start

    # Fortran runs even when JAX failed: a flag set that breaks both drivers is a
    # configuration problem, while one that breaks only JAX is a missing JAX feature.
    f90_start = time.time()
    f90_rc = task.supervisor.run_and_log(f90_cmd, task.repo_root, f90_log)
    f90_elapsed = time.time() - f90_start
    jax_timesteps = _parse_completed_timesteps(jax_log, jax=True) if jax_rc == 0 else None
    fortran_timesteps = _parse_completed_timesteps(f90_log, jax=False) if f90_rc == 0 else None

    if jax_rc != 0 or f90_rc != 0:
        if jax_rc != 0 and f90_rc != 0:
            status = "both_failed"
            note = f"{_tail(jax_log)}\n--- fortran ---\n{_tail(f90_log)}"
        elif jax_rc != 0:
            status = "jax_failed"
            note = _tail(jax_log)
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
            jax_rc=jax_rc,
            fortran_rc=f90_rc,
            bindiff_rc=-1,
            jax_elapsed_s=jax_elapsed,
            fortran_elapsed_s=f90_elapsed,
            elapsed_s=time.time() - start,
            flag_data=task.flag_data,
            note=note,
            jax_timesteps=jax_timesteps,
            fortran_timesteps=fortran_timesteps,
        )

    diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-v", "2",
        "-strict",
        "-case", task.config.source_case,
        "-t", str(task.bindiff_threshold),
        "-pt", str(task.bindiff_percent_threshold),
        "--result-json", str(diff_report),
        str(jax_run_out_dir),
        str(f90_run_out_dir),
    ]
    diff_report.unlink(missing_ok=True)
    diff_rc = task.supervisor.run_and_log(diff_cmd, task.repo_root, diff_log)

    # Matched runs show their completed steps directly. For a numerical diff,
    # map bindiff's saved-prefix result to the model step shown in the table.
    first_failing_timestep = None
    if diff_rc != 0 and jax_timesteps is not None and jax_timesteps == fortran_timesteps:
        first_failing_timestep = _first_failing_timestep(
            diff_report, task.config.source_case, jax_timesteps,
            _parse_model_timing(f90_log),
        )

    status = "match" if diff_rc == 0 else "diff"
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
        jax_rc=jax_rc,
        fortran_rc=f90_rc,
        bindiff_rc=diff_rc,
        jax_elapsed_s=jax_elapsed,
        fortran_elapsed_s=f90_elapsed,
        elapsed_s=time.time() - start,
        note=f"Bindiff report: {diff_report}",
        first_failing_timestep=first_failing_timestep,
        jax_timesteps=jax_timesteps,
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run SCM cases with both the JAX driver and Fortran standalone, "
            "then compare outputs using run_bindiff_all.py. Other single-dash "
            "run_scm.py options are forwarded to both runs."
        ),
        allow_abbrev=False,
    )
    parser.add_argument("-j", "-jobs", dest="jobs", type=int, default=1,
                        help="Number of concurrent case pairs (default: 1; each JAX compilation can use several GiB).")
    parser.add_argument(
        "-bindiff_verbose", type=int, default=2, choices=[0, 1, 2],
        help="Verbosity level for the final combined run_bindiff_all.py run.",
    )
    parser.add_argument(
        "-bindiff_threshold", type=float, default=None,
        help="Override both absolute and percentage thresholds (default: absolute 1e-7, per-case percentage).",
    )
    parser.add_argument(
        "-bindiff_percent_threshold", type=float, default=None,
        help="Override only the percentage threshold, in percentage points; takes precedence over -bindiff_threshold.",
    )
    parser.add_argument(
        "-cases", nargs="+", default=None,
        help="Case names to run (default is the curated supported set).",
    )
    parser.add_argument(
        "-flag_config_file", type=str, default=None,
        help="JSON file describing alternate flag settings.",
    )
    args, forwarded = parser.parse_known_args()
    if args.jobs <= 0:
        parser.error("-jobs must be positive")
    # These options belong to the comparison, not to either individual model.
    reserved = {"-jax", "-python", "-exe", "-driver_test", "-gdb",
                "-out_dir", "-multicol", "-override", "-install_dir"}
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
    args.run_scm_args = forwarded
    return args


def main() -> int:
    args = parse_args()
    ensure_environment()
    accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").lower()
    if accelerator in {"cuda13", "metal"} and args.jobs != 1:
        print(
            "WARNING: GPU comparison currently runs one case process at a time to avoid "
            "multiple JAX workers contending for the same device; forcing -j 1."
        )
        args.jobs = 1

    supervisor = RunSupervisor()
    previous_handlers = {sig: signal.signal(sig, supervisor.handle_signal)
                         for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        return _run_comparisons(args, supervisor)
    except ComparisonCancelled:
        print("Comparison interrupted; stopped all active case process groups.", flush=True)
        return 128 + (supervisor.signal_number or signal.SIGINT)
    finally:
        supervisor.stopping.set()
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)


def _run_tasks(tasks, jobs, emit, supervisor):
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

    repo_root = Path(__file__).resolve().parents[1]
    run_lock = _acquire_run_lock(repo_root)
    if run_lock is None:
        return 2
    results_root = (repo_root / RESULTS_DIRNAME).resolve()
    jax_output_root = results_root / JAX_OUTPUT_DIRNAME
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
        flag_config = read_flag_settings(args.flag_config_file) if args.flag_config_file else {}
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
            tasks.append(TaskCtx(
                config=config,
                flag_data=flag_data,
                repo_root=repo_root,
                run_scm_args=args.run_scm_args,
                bindiff_threshold=args.bindiff_threshold if args.bindiff_threshold is not None else 1.0e-7,
                bindiff_percent_threshold=(
                    args.bindiff_percent_threshold if args.bindiff_percent_threshold is not None
                    else args.bindiff_threshold if args.bindiff_threshold is not None
                    else config.percent_threshold
                ),
                run_output_root=results_root,
                supervisor=supervisor,
            ))

    output_owners = set()
    for task in tasks:
        key = (task.output_group, task.config.source_case)
        if key in output_owners:
            print(f"ERROR: duplicate output destination: {key}")
            return 2
        output_owners.add(key)
    if results_root.exists():
        shutil.rmtree(results_root)

    print("\nRun settings:")
    print(f"  Output: {results_root}")
    print(f"  Statistics: {_forwarded_value(args.run_scm_args, '-stats') or 'input/stats/standard_stats.in'}")
    print(f"  Debug: {_forwarded_value(args.run_scm_args, '-debug') or 'case default'}")
    if args.run_scm_args:
        print(f"  Forwarded to both runs: {shlex.join(args.run_scm_args)}")
    print("  Columns: 4")
    print(f"  Workers: {min(args.jobs, len(tasks))}")
    default_percent = (
        args.bindiff_percent_threshold if args.bindiff_percent_threshold is not None
        else args.bindiff_threshold if args.bindiff_threshold is not None
        else 1.0e-7
    )
    print(
        f"  Comparison limits: absolute {tasks[0].bindiff_threshold:g}; "
        f"percentage {default_percent:g}%"
    )

    print("\nCases and flag sets:")
    case_label = "case" if len(cases) == 1 else "cases"
    flag_label = "flag set" if len(flag_data_list) == 1 else "flag sets"
    pair_label = "pair" if len(tasks) == 1 else "pairs"
    print(
        f"  {len(cases)} {case_label} x {len(flag_data_list)} {flag_label}"
        f" = {len(tasks)} JAX/Fortran {pair_label}"
    )
    print("  Flag sets:")
    for flag_data in flag_data_list:
        print(f"    {flag_data.dir_name}")
        for setting, value in (flag_data.flag_dict or {}).items():
            print(f"      - {setting} = {format_override_value(value)}")
    print("  Cases:")
    percent_by_case = {task.config.case: task.bindiff_percent_threshold for task in tasks}
    for config in configs:
        iterations = f"{config.max_iters} iterations" if config.max_iters is not None else "native duration"
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
        f"{'Timesteps':>13} {'JAX (s)':>10} {'Fortran (s)':>12} {'Total (s)':>10}"
    )
    print(f"\n{table_header}", flush=True)
    print("-" * len(table_header), flush=True)

    def _emit(result: CaseResult) -> None:
        results.append(result)
        total_steps = result.fortran_timesteps or result.jax_timesteps
        if (result.status == "match" and total_steps is not None
            and result.jax_timesteps == result.fortran_timesteps):
            steps = f"{total_steps} / {total_steps}"
        elif result.status == "diff" and result.first_failing_timestep is not None:
            steps = f"{result.first_failing_timestep} / {total_steps}"
        else:
            steps = "-"
        print(
            f"[{result.status:14}] {result.flag_data.dir_name:{flag_width}} {result.case:{case_width}} "
            f"{steps:>13} {result.jax_elapsed_s:10.1f} {result.fortran_elapsed_s:12.1f} "
            f"{result.elapsed_s:10.1f}",
            flush=True,
        )

    _run_tasks(tasks, args.jobs, _emit, supervisor)

    results.sort(key=lambda r: (r.flag_data.dir_name, r.case))
    statuses = ("match", "diff", "jax_failed", "fortran_failed", "both_failed")

    # Per-case comparisons enforce each case's tolerances. This additional audit
    # uses the loosest selected percentage threshold so allowed Morrison differences
    # do not fail it; it cannot override any per-case failure. --flag-sets also
    # detects output directories/files present on only one side.
    final_bindiff_log = results_root / FINAL_BINDIFF_LOG_FILENAME
    final_diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "--flag-sets",
        "-strict",
        "-v", str(args.bindiff_verbose),
        "-t", str(tasks[0].bindiff_threshold),
        "-pt", str(max(task.bindiff_percent_threshold for task in tasks)),
        str(jax_output_root),
        str(f90_output_root),
    ]
    final_bindiff_rc = supervisor.run_and_log(final_diff_cmd, repo_root, final_bindiff_log)

    summary = {
        "total_cases": len(results),
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
        failed = sum(counts[status] for status in ("jax_failed", "fortran_failed", "both_failed"))
        print(f"\n{title}:")
        print(f"  {'Match':>5} {'Diff':>4} {'Failed':>6}")
        print(f"  {'-----':>5} {'----':>4} {'------':>6}")
        print(f"  {counts['match']:>5} {counts['diff']:>4} {failed:>6}")

    if len(flag_data_list) > 1:
        _print_count_table("Overall cases", summary)
    for flag_data in flag_data_list:
        _print_count_table(f"Flag set {flag_data.dir_name}", summary["flag_sets"][flag_data.dir_name])

    case_failure = any(summary[status] > 0 for status in statuses if status != "match")
    if final_bindiff_rc != 0 and not case_failure:
        print("\nFinal bindiff audit failed despite all case comparisons matching.")

    hours, remainder = divmod(summary["elapsed_s"], 3600)
    minutes, seconds = divmod(remainder, 60)
    duration = (
        f"{int(hours)}h {int(minutes)}m {seconds:.1f}s" if hours >= 1
        else f"{int(minutes)}m {seconds:.1f}s" if minutes >= 1
        else f"{seconds:.1f}s"
    )
    audit_status = "passed" if final_bindiff_rc == 0 else f"failed (exit {final_bindiff_rc})"
    print(f"\nTotal time: {duration}")
    print(f"Results JSON: {summary_path}")
    print(f"Bindiff log: {final_bindiff_log} [{audit_status}]")

    return 1 if case_failure or final_bindiff_rc != 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
