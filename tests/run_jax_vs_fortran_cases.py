#!/usr/bin/env python3
"""Run selected SCM cases with both JAX and Fortran drivers and diff outputs."""

from __future__ import annotations

import argparse
import fcntl
import json
import multiprocessing as mp
import os
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

def _reexec_with_repo_jax_python() -> None:
    """Initialize and use the repository-local JAX environment."""
    repo_root = Path(__file__).resolve().parents[1]
    launcher = repo_root / "clubb_jax" / "run_jax_wrapper.sh"
    if not launcher.is_file():
        return

    initialized_env_var = "_CLUBB_JAX_HARNESS_ENV_INITIALIZED"
    if os.environ.get(initialized_env_var) != "1":
        init_env = subprocess.run([str(launcher), "--init_env"], check=False)
        if init_env.returncode != 0:
            raise SystemExit(init_env.returncode)

    accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").lower()
    default_venv = ".venv-jax-cuda13" if accelerator == "cuda13" else ".venv-jax"
    venv_dir = Path(os.environ.get("CLUBB_JAX_VENV", repo_root / default_venv))
    if not venv_dir.is_absolute():
        venv_dir = repo_root / venv_dir
    venv_python = venv_dir / "bin" / "python"
    if not venv_python.is_file():
        return

    if Path(sys.executable).absolute() == venv_python.absolute():
        return

    exec_env = os.environ.copy()
    exec_env[initialized_env_var] = "1"
    os.execve(
        str(venv_python),
        [str(venv_python), str(Path(__file__).resolve()), *sys.argv[1:]],
        exec_env,
    )


_reexec_with_repo_jax_python()

CLUBB_ROOT = Path(__file__).resolve().parents[1]
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from utilities.flag_sets import build_override_arg, get_flag_sets, read_flag_settings  # noqa: E402


# The python driver is much simpler and doesn't support all
# features used by some cases (e.g. microphysics, BUGS, sponge layer, SILHS), so we
# run a curated set of cases that avoid those features.
# Values are per-case max_iters (number of timesteps to run).
# None means run the full case (don't pass -max_iters to run_scm.py).
DEFAULT_CASES = {
    "arm":                    360,      # stable, diffs expected after ~600 60s-timesteps
    "atex":                   360,      # stable, diffs expected after ~400 60s-timesteps
    "bomex":                  None,
    "cobra":                  360,      # very stable, limited for speed
    "dycoms2_rf01":           None,
    "dycoms2_rf01_fixed_sst": 300,      # stablish, switching to l_diag_Lscale_from_tau=.false.
                                        # starting causing diffs after ~360 timesteps
    "dycoms2_rf02_nd":        None,
    "fire":                   None,
    "gabls2":                 360,      # very stable, limited for speed
    "gabls3_night":           360,      # very stable, limited for speed
    "jun25_altocu":           180,      # stablish, diffs expected after ~200 60s-timesteps
    "neutral":                None,
    "wangara":                None,
}

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
    avg_diff_timestep: float = -1.0


def _parse_earliest_timesteps(log_path: Path) -> list[int]:
    if not log_path.exists():
        return []
    text = log_path.read_text(encoding="utf-8", errors="replace")
    timesteps: list[int] = []
    in_table = False
    for line in text.splitlines():
        if "Earliest Timestep" in line:
            in_table = True
            continue
        if in_table:
            stripped = line.strip()
            if not stripped or stripped.startswith("="):
                in_table = False
                continue
            if set(stripped.replace(" ", "")) <= {"-"}:
                continue
            for tok in reversed(stripped.split()):
                try:
                    timesteps.append(int(tok))
                    break
                except ValueError:
                    continue
    return timesteps


def _run_and_log(cmd: list[str], cwd: Path, log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")
        # The child writes straight to the file descriptor, so flush the header
        # first or it only lands, after the child's output, when the file closes.
        log.flush()
        proc = subprocess.run(cmd, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT)
        log.write(f"\n[exit_code] {proc.returncode}\n")
    return proc.returncode


def _check_jax_runtime() -> str | None:
    """Validate the interpreter before launching every case with it."""
    accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").lower()
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import jax, jaxlib, netCDF4, sys, tabulate; "
                "backend = jax.default_backend(); "
                "expected = 'gpu' if sys.argv[1] == 'cuda13' else 'cpu'; "
                "assert backend == expected, "
                "f'requested {sys.argv[1]} but JAX initialized {backend}: {jax.devices()}'; "
                "print(f'jax={jax.__version__} jaxlib={jaxlib.__version__} "
                "backend={backend} devices={jax.devices()}')"
            ),
            accelerator,
        ],
        text=True,
        capture_output=True,
    )
    if probe.returncode == 0:
        return probe.stdout.strip()

    detail = (probe.stderr or probe.stdout).strip()
    print(f"ERROR: JAX runtime check failed with {sys.executable}:\n{detail}")
    print(
        "Create or repair the selected repository JAX environment, or invoke this script with a Python "
        "environment containing compatible jax and jaxlib packages."
    )
    return None


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
    case: str
    flag_data: FlagData
    repo_root: Path
    stats: str
    debug: str | None
    max_iters: int | None
    bindiff_threshold: float
    run_output_root: Path


def _run_case_w_flags(task: TaskCtx) -> CaseResult:
    run_scm = task.repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = task.repo_root / "run_scripts" / "run_bindiff_all.py"

    start = time.time()

    common_args = [
        str(run_scm),
        "-stats", task.stats,
        "-multicol", HR_SPEC,
    ]
    if task.debug is not None:
        common_args += ["-debug", task.debug]
    if task.max_iters is not None:
        common_args += ["-max_iters", str(task.max_iters)]


    #A subdirectory for each flagset is created to store the .nc output from that flagset.
    # The output directory will look like this:
    # output/tests/jax_driver_test_results/jax_output/<file_name + flag_name> and output/tests/jax_driver_test_results/fortran_output/<file_name + flag_name>
    flag_dir_name = task.flag_data.dir_name
    jax_run_out_dir = task.run_output_root / JAX_OUTPUT_DIRNAME / flag_dir_name
    f90_run_out_dir = task.run_output_root / FORTRAN_OUTPUT_DIRNAME / flag_dir_name
    jax_run_out_dir.mkdir(parents=True, exist_ok=True)
    f90_run_out_dir.mkdir(parents=True, exist_ok=True)

    jax_cmd = [sys.executable, *common_args, "-jax", "-out_dir", str(jax_run_out_dir)]
    f90_cmd = [sys.executable, *common_args, "-out_dir", str(f90_run_out_dir)]

    override_arg = build_override_arg(task.flag_data.flag_dict)
    if override_arg is not None:
        jax_cmd += ["-override", override_arg]
        f90_cmd += ["-override", override_arg]

    jax_cmd.append(task.case)
    f90_cmd.append(task.case)

    log_dir = task.run_output_root / LOGS_DIRNAME / flag_dir_name
    jax_log = log_dir / f"{task.case}_run_jax.log"
    f90_log = log_dir / f"{task.case}_run_fortran.log"
    diff_log = log_dir / f"{task.case}_bindiff.log"

    jax_start = time.time()
    jax_rc = _run_and_log(jax_cmd, task.repo_root, jax_log)
    jax_elapsed = time.time() - jax_start

    # Fortran runs even when JAX failed: a flag set that breaks both drivers is a
    # configuration problem, while one that breaks only JAX is a missing JAX feature.
    f90_start = time.time()
    f90_rc = _run_and_log(f90_cmd, task.repo_root, f90_log)
    f90_elapsed = time.time() - f90_start

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
            case=task.case,
            status=status,
            jax_rc=jax_rc,
            fortran_rc=f90_rc,
            bindiff_rc=-1,
            jax_elapsed_s=jax_elapsed,
            fortran_elapsed_s=f90_elapsed,
            elapsed_s=time.time() - start,
            flag_data=task.flag_data,
            note=note,
        )

    diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-v", "2",
        "-case", task.case,
        "-t", str(task.bindiff_threshold),
        "-pt", str(task.bindiff_threshold),
        str(jax_run_out_dir),
        str(f90_run_out_dir),
    ]
    diff_rc = _run_and_log(diff_cmd, task.repo_root, diff_log)

    ts_list = _parse_earliest_timesteps(diff_log) if diff_rc != 0 else []
    avg_ts = sum(ts_list) / len(ts_list) if ts_list else -1.0

    status = "match" if diff_rc == 0 else "diff"
    return CaseResult(
        case=task.case,
        flag_data=task.flag_data,
        status=status,
        jax_rc=jax_rc,
        fortran_rc=f90_rc,
        bindiff_rc=diff_rc,
        jax_elapsed_s=jax_elapsed,
        fortran_elapsed_s=f90_elapsed,
        elapsed_s=time.time() - start,
        note=_tail(diff_log),
        avg_diff_timestep=avg_ts,
    )

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run SCM cases with both the JAX driver and Fortran standalone, "
            "then compare outputs using run_bindiff_all.py."
        )
    )
    parser.add_argument("-j", "--jobs", type=int, default=8, help="Number of parallel case workers.")
    parser.add_argument(
        "-stats",
        "--stats",
        default="input/stats/standard_stats.in",
        help="Stats setting forwarded to run_scm.py (use 'none' to disable stats output).",
    )
    parser.add_argument(
        "-debug",
        "--debug",
        default=None,
        help="Debug level forwarded to run_scm.py (0-3).",
    )
    parser.add_argument(
        "--max-iters",
        type=int,
        default=None,
        help="Override max_iters for all cases (default: use per-case values from DEFAULT_CASES).",
    )
    parser.add_argument(
        "--bindiff-verbose",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="Verbosity level for the final combined run_bindiff_all.py run.",
    )
    parser.add_argument(
        "--bindiff-threshold",
        type=float,
        default=1.0e-7,
        help="Difference threshold passed to run_bindiff_all.py via -t.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Do not delete existing output dirs before rerun.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        default=None,
        help="Case names to run (default is the curated supported set).",
    )

    parser.add_argument(
        "--flag-config-file", type=str, default=None,
        help="JSON file describing alternate flag settings."
    )
    parser.add_argument(
        "--skip-default-flags",
        action="store_true",
        help="Do not run the unmodified default flag configuration.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").lower()
    if accelerator == "cuda13" and args.jobs != 1:
        print(
            "WARNING: GPU comparison currently runs one case process at a time to avoid "
            "multiple JAX workers contending for the same device; forcing -j 1."
        )
        args.jobs = 1

    jax_runtime = _check_jax_runtime()
    if jax_runtime is None:
        return 2
    print(f"Python runtime: {sys.executable} ({jax_runtime})")

    cases = list(args.cases) if args.cases else list(DEFAULT_CASES.keys())
    case_iters = {
        case: args.max_iters if args.max_iters is not None else DEFAULT_CASES.get(case)
        for case in cases
    }

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
        flag_sets = get_flag_sets(args.skip_default_flags, flag_config)
    except (OSError, ValueError) as exc:
        print(f"ERROR: could not load flag sets: {exc}")
        return 2
    if not flag_sets:
        print("ERROR: no flag sets to run; --skip-default-flags needs --flag-config-file.")
        return 2

    if results_root.exists() and not args.keep_existing:
        shutil.rmtree(results_root)

    flag_data_list = [
        FlagData(flag_name=name, flag_dict=overrides, flag_file=args.flag_config_file)
        for name, overrides in flag_sets.items()
    ]

    for flag_data in flag_data_list:
        (jax_output_root / flag_data.dir_name).mkdir(parents=True, exist_ok=True)
        (f90_output_root / flag_data.dir_name).mkdir(parents=True, exist_ok=True)

    tasks = []
    for flag_data in flag_data_list:
        for case in cases:
            tasks.append(TaskCtx(
                case=case,
                flag_data=flag_data,
                repo_root=repo_root,
                stats=args.stats,
                debug=args.debug,
                max_iters=case_iters[case],
                bindiff_threshold=args.bindiff_threshold,
                run_output_root=results_root,
            ))

    print(
        f"Running {len(cases)} case(s) x {len(flag_data_list)} flag set(s) "
        f"= {len(tasks)} run pair(s) with {args.jobs} worker(s)"
    )
    print(f"Results root: {results_root}")
    print(f"JAX output: {jax_output_root}")
    print(f"Fortran output: {f90_output_root}")
    print(f"Flag sets: {', '.join(fd.dir_name for fd in flag_data_list)}")
    print(f"Stats: {args.stats}")
    print(f"Debug: {args.debug if args.debug is not None else '(case default)'}")
    print(f"JAX accelerator: {accelerator}")

    start = time.time()
    results: list[CaseResult] = []

    def _emit(result: CaseResult) -> None:
        results.append(result)
        ts_info = ""
        if result.avg_diff_timestep >= 0:
            ts_info = f" avg_diff_ts={result.avg_diff_timestep:.1f}"
        print(
            f"[{result.status:14}] {result.flag_data.dir_name:24} {result.case:22} "
            f"jax={result.jax_elapsed_s:.1f}s f90={result.fortran_elapsed_s:.1f}s "
            f"total={result.elapsed_s:.1f}s{ts_info}"
        )

    if args.jobs == 1:
        for task in tasks:
            _emit(_run_case_w_flags(task))
    else:
        try:
            with mp.Pool(processes=args.jobs) as pool:
                for result in pool.imap_unordered(_run_case_w_flags, tasks):
                    _emit(result)
        except (PermissionError, OSError) as exc:
            print(f"WARNING: multiprocessing unavailable ({exc}); falling back to serial execution.")
            for task in tasks:
                _emit(_run_case_w_flags(task))

    results.sort(key=lambda r: (r.flag_data.dir_name, r.case))
    statuses = ("match", "diff", "jax_failed", "fortran_failed", "both_failed")

    # The final run re-diffs the same pairs the per-case runs already covered, but in
    # --flag-sets mode it also reports flag sets or NetCDF files that exist on only one
    # side -- which the per-case comparisons cannot see, since they only ever look at
    # files present in both directories.
    final_bindiff_log = results_root / FINAL_BINDIFF_LOG_FILENAME
    final_diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "--flag-sets",
        "-v", str(args.bindiff_verbose),
        "-t", str(args.bindiff_threshold),
        "-pt", str(args.bindiff_threshold),
        str(jax_output_root),
        str(f90_output_root),
    ]
    final_bindiff_rc = _run_and_log(final_diff_cmd, repo_root, final_bindiff_log)

    # Aggregate from the per-case results rather than re-parsing the combined log, so
    # the number survives a structurally failed final bindiff.
    diff_timesteps = [r.avg_diff_timestep for r in results if r.avg_diff_timestep >= 0]
    final_avg_ts = sum(diff_timesteps) / len(diff_timesteps) if diff_timesteps else -1.0

    summary = {
        "total_cases": len(results),
        **{status: sum(r.status == status for r in results) for status in statuses},
        "elapsed_s": time.time() - start,
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

    print("\nSummary:")
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=2))
    if final_avg_ts >= 0:
        print(f"Average earliest diff timestep (across all cases): {final_avg_ts:.1f}")
    print(f"Detailed results: {summary_path}")
    print(f"Final bindiff log: {final_bindiff_log} (rc={final_bindiff_rc})")

    if any(summary[status] > 0 for status in statuses if status != "match"):
        return 1
    if final_bindiff_rc != 0:
        print(
            "ERROR: every case matched but the final flag-set bindiff still reported a "
            f"problem (rc={final_bindiff_rc}); see {final_bindiff_log}."
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
