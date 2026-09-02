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
SUMMARY_FILENAME = "case_compare_summary.json"
FINAL_BINDIFF_LOG_FILENAME = "final_bindiff.log"
HR_SPEC = "C8/0.2:0.8/4"


@dataclass
class CaseResult:
    case: str
    status: str
    jax_rc: int
    fortran_rc: int
    bindiff_rc: int
    jax_elapsed_s: float
    fortran_elapsed_s: float
    elapsed_s: float
    case_dir: str
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


def _run_case(
    case: str,
    repo_root: Path,
    stats: str,
    debug: str | None,
    max_iters: int | None,
    bindiff_threshold: float,
    jax_out: Path,
    f90_out: Path,
    results_root: Path,
) -> CaseResult:
    run_scm = repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = repo_root / "run_scripts" / "run_bindiff_all.py"

    start = time.time()

    common_args = [
        str(run_scm),
        "-stats", stats,
        "-multicol", HR_SPEC,
    ]
    if debug is not None:
        common_args += ["-debug", debug]
    if max_iters is not None:
        common_args += ["-max_iters", str(max_iters)]

    jax_cmd = [sys.executable, *common_args, "-jax", "-out_dir", str(jax_out), case]
    f90_cmd = [sys.executable, *common_args, "-out_dir", str(f90_out), case]

    jax_log = results_root / f"{case}_run_jax.log"
    f90_log = results_root / f"{case}_run_fortran.log"
    diff_log = results_root / f"{case}_bindiff.log"

    jax_start = time.time()
    jax_rc = _run_and_log(jax_cmd, repo_root, jax_log)
    jax_elapsed = time.time() - jax_start
    if jax_rc != 0:
        return CaseResult(
            case=case,
            status="jax_failed",
            jax_rc=jax_rc,
            fortran_rc=-1,
            bindiff_rc=-1,
            jax_elapsed_s=jax_elapsed,
            fortran_elapsed_s=0.0,
            elapsed_s=time.time() - start,
            case_dir=str(results_root),
            note=_tail(jax_log),
        )

    f90_start = time.time()
    f90_rc = _run_and_log(f90_cmd, repo_root, f90_log)
    f90_elapsed = time.time() - f90_start
    if f90_rc != 0:
        return CaseResult(
            case=case,
            status="fortran_failed",
            jax_rc=jax_rc,
            fortran_rc=f90_rc,
            bindiff_rc=-1,
            jax_elapsed_s=jax_elapsed,
            fortran_elapsed_s=f90_elapsed,
            elapsed_s=time.time() - start,
            case_dir=str(results_root),
            note=_tail(f90_log),
        )

    diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-v", "2",
        "-case", case,
        "-t", str(bindiff_threshold),
        "-pt", str(bindiff_threshold),
        str(jax_out),
        str(f90_out),
    ]
    diff_rc = _run_and_log(diff_cmd, repo_root, diff_log)

    ts_list = _parse_earliest_timesteps(diff_log) if diff_rc != 0 else []
    avg_ts = sum(ts_list) / len(ts_list) if ts_list else -1.0

    status = "match" if diff_rc == 0 else "diff"
    return CaseResult(
        case=case,
        status=status,
        jax_rc=jax_rc,
        fortran_rc=f90_rc,
        bindiff_rc=diff_rc,
        jax_elapsed_s=jax_elapsed,
        fortran_elapsed_s=f90_elapsed,
        elapsed_s=time.time() - start,
        case_dir=str(results_root),
        note=_tail(diff_log),
        avg_diff_timestep=avg_ts,
    )


def _worker(task: tuple) -> CaseResult:
    (
        case,
        repo_root,
        stats,
        debug,
        max_iters,
        bindiff_threshold,
        jax_output_root,
        f90_output_root,
        results_root,
    ) = task

    return _run_case(
        case=case,
        repo_root=Path(repo_root),
        stats=stats,
        debug=debug,
        max_iters=max_iters,
        bindiff_threshold=bindiff_threshold,
        jax_out=Path(jax_output_root),
        f90_out=Path(f90_output_root),
        results_root=Path(results_root),
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

    if results_root.exists() and not args.keep_existing:
        shutil.rmtree(results_root)
    jax_output_root.mkdir(parents=True, exist_ok=True)
    f90_output_root.mkdir(parents=True, exist_ok=True)

    tasks = [
        (
            case,
            str(repo_root),
            args.stats,
            args.debug,
            case_iters[case],
            args.bindiff_threshold,
            str(jax_output_root),
            str(f90_output_root),
            str(results_root),
        )
        for case in cases
    ]

    print(f"Running {len(cases)} case(s) with {args.jobs} worker(s)")
    print(f"Results root: {results_root}")
    print(f"JAX output: {jax_output_root}")
    print(f"Fortran output: {f90_output_root}")
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
            f"[{result.status:14}] {result.case:22} "
            f"jax={result.jax_elapsed_s:.1f}s f90={result.fortran_elapsed_s:.1f}s "
            f"total={result.elapsed_s:.1f}s{ts_info}"
        )

    if args.jobs == 1:
        for task in tasks:
            _emit(_worker(task))
    else:
        try:
            with mp.Pool(processes=args.jobs) as pool:
                for result in pool.imap_unordered(_worker, tasks):
                    _emit(result)
        except (PermissionError, OSError) as exc:
            print(f"WARNING: multiprocessing unavailable ({exc}); falling back to serial execution.")
            for task in tasks:
                _emit(_worker(task))

    results.sort(key=lambda r: r.case)
    summary = {
        "total_cases": len(results),
        "match": sum(r.status == "match" for r in results),
        "diff": sum(r.status == "diff" for r in results),
        "jax_failed": sum(r.status == "jax_failed" for r in results),
        "fortran_failed": sum(r.status == "fortran_failed" for r in results),
        "elapsed_s": time.time() - start,
        "cases": [asdict(r) for r in results],
    }

    summary_path = results_root / SUMMARY_FILENAME
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    final_bindiff_log = results_root / FINAL_BINDIFF_LOG_FILENAME
    final_diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-v", str(args.bindiff_verbose),
        "-t", str(args.bindiff_threshold),
        "-pt", str(args.bindiff_threshold),
        str(jax_output_root),
        str(f90_output_root),
    ]
    final_bindiff_rc = _run_and_log(final_diff_cmd, repo_root, final_bindiff_log)
    final_ts_list = _parse_earliest_timesteps(final_bindiff_log) if final_bindiff_rc != 0 else []
    final_avg_ts = sum(final_ts_list) / len(final_ts_list) if final_ts_list else -1.0

    print("\nSummary:")
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=2))
    if final_avg_ts >= 0:
        print(f"Average earliest diff timestep (across all cases): {final_avg_ts:.1f}")
    print(f"Detailed results: {summary_path}")
    print(f"Final bindiff log: {final_bindiff_log} (rc={final_bindiff_rc})")

    if summary["diff"] > 0 or summary["jax_failed"] > 0 or summary["fortran_failed"] > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
