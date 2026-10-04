#!/usr/bin/env python3
"""Run selected SCM cases with both drivers and diff outputs."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path


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

RESULTS_DIRNAME = Path("output") / "python_driver_test_results"
PYTHON_OUTPUT_DIRNAME = "python_output"
FORTRAN_OUTPUT_DIRNAME = "fortran_output"
SUMMARY_FILENAME = "case_compare_summary.json"
FINAL_BINDIFF_LOG_FILENAME = "final_bindiff.log"
HR_SPEC = "C8/0.2:0.8/4" # used to run multicol mode by varying C8

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tuner.system_defaults import default_max_workers as default_workers
from utilities.create_case_namelist import parse_forwarded_args
from run_scripts.run_scm_all import positive_int, split_values

@dataclass
class CaseResult:
    case: str
    status: str
    python_rc: int
    fortran_rc: int
    bindiff_rc: int
    python_elapsed_s: float
    fortran_elapsed_s: float
    elapsed_s: float
    case_dir: str
    note: str = ""
    avg_diff_timestep: float = -1.0


def _average_earliest_timestep(report_path: Path, case: str | None = None) -> float:
    """Read bindiff's summary without parsing its human-readable log."""
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        result = report["cases"][case] if case is not None else report
        value = result["average_earliest_timestep"]
        return float(value) if value is not None else -1.0
    except (OSError, ValueError, KeyError, TypeError):
        return -1.0


def _run_and_log(cmd: list[str], cwd: Path, log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(shlex.quote(part) for part in cmd) + "\n\n")
        proc = subprocess.run(cmd, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT)
        log.write(f"\n[exit_code] {proc.returncode}\n")
    return proc.returncode


def _tail(path: Path, n: int = 25) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-n:])


def _run_case(
    case: str,
    repo_root: Path,
    stats_file: Path,
    max_iters: int | None,
    bindiff_threshold: float,
    py_out: Path,
    f90_out: Path,
    results_root: Path,
    run_scm_args: tuple[str, ...] = (),
) -> CaseResult:
    run_scm = repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = repo_root / "run_scripts" / "run_bindiff_all.py"

    start = time.time()

    common_args = [
        str(run_scm),
        "-stats", str(stats_file),
        "-multicol", HR_SPEC,
    ]
    if max_iters is not None:
        common_args += ["-max_iters", str(max_iters)]

    common_args.extend(run_scm_args)
    py_cmd = [sys.executable, *common_args, "-python", "-output_dir", str(py_out), case]
    f90_cmd = [sys.executable, *common_args, "-output_dir", str(f90_out), case]

    py_log = results_root / f"{case}_run_python.log"
    f90_log = results_root / f"{case}_run_fortran.log"
    diff_log = results_root / f"{case}_bindiff.log"

    py_start = time.time()
    py_rc = _run_and_log(py_cmd, repo_root, py_log)
    py_elapsed = time.time() - py_start
    if py_rc != 0:
        return CaseResult(
            case=case,
            status="python_failed",
            python_rc=py_rc,
            fortran_rc=-1,
            bindiff_rc=-1,
            python_elapsed_s=py_elapsed,
            fortran_elapsed_s=0.0,
            elapsed_s=time.time() - start,
            case_dir=str(results_root),
            note=_tail(py_log),
        )

    f90_start = time.time()
    f90_rc = _run_and_log(f90_cmd, repo_root, f90_log)
    f90_elapsed = time.time() - f90_start
    if f90_rc != 0:
        return CaseResult(
            case=case,
            status="fortran_failed",
            python_rc=py_rc,
            fortran_rc=f90_rc,
            bindiff_rc=-1,
            python_elapsed_s=py_elapsed,
            fortran_elapsed_s=f90_elapsed,
            elapsed_s=time.time() - start,
            case_dir=str(results_root),
            note=_tail(f90_log),
        )

    # The JSON report supplies the timing summary; -case avoids prefix
    # collisions in the shared output directories.
    diff_report = results_root / f"{case}_bindiff.json"
    diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-verbose", "2",
        "-case", case,
        "-threshold", str(bindiff_threshold),
        "-percent_threshold", str(bindiff_threshold),
        "-result_json", str(diff_report),
        str(py_out),
        str(f90_out),
    ]
    diff_report.unlink(missing_ok=True)
    diff_rc = _run_and_log(diff_cmd, repo_root, diff_log)
    avg_ts = _average_earliest_timestep(diff_report, case) if diff_rc != 0 else -1.0

    status = "match" if diff_rc == 0 else "diff"
    return CaseResult(
        case=case,
        status=status,
        python_rc=py_rc,
        fortran_rc=f90_rc,
        bindiff_rc=diff_rc,
        python_elapsed_s=py_elapsed,
        fortran_elapsed_s=f90_elapsed,
        elapsed_s=time.time() - start,
        case_dir=str(results_root),
        note=f"Bindiff report: {diff_report}",
        avg_diff_timestep=avg_ts,
    )


def _worker(task: tuple) -> CaseResult:
    (
        case,
        repo_root,
        stats_file,
        max_iters,
        bindiff_threshold,
        py_output_root,
        f90_output_root,
        results_root,
        run_scm_args,
    ) = task

    return _run_case(
        case=case,
        repo_root=Path(repo_root),
        stats_file=Path(stats_file),
        max_iters=max_iters,
        bindiff_threshold=bindiff_threshold,
        py_out=Path(py_output_root),
        f90_out=Path(f90_output_root),
        results_root=Path(results_root),
        run_scm_args=run_scm_args,
    )


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run SCM cases with both Fortran and Python drivers, then compare "
            "outputs using run_bindiff_all.py."
        ),
        add_help=False, allow_abbrev=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        '-workers', dest='jobs',
        type=positive_int,
        default=default_workers(),
        help="Maximum concurrent case pairs (default: half the available logical CPUs).",
        metavar='N',
    )
    parser.add_argument(
        '-stats', dest='stats',
        default="input/stats/all_stats.in",
        help="Stats file passed to run_scm.py.",
    )
    parser.add_argument(
        '-max_iters', dest='max_iters',
        type=int,
        default=None,
        help="Override max_iters for all cases (default: use per-case values from DEFAULT_CASES).",
    )
    parser.add_argument(
        '-bindiff_verbose', dest='bindiff_verbose',
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="Verbosity level for the final combined run_bindiff_all.py run.",
    )
    parser.add_argument(
        '-bindiff_threshold', dest='bindiff_threshold',
        type=float,
        default=1.0e-7,
        help="Difference threshold passed to run_bindiff_all.py via -threshold.",
    )
    parser.add_argument(
        '-keep_existing', dest='keep_existing',
        action="store_true",
        help="Do not delete existing output dirs before rerun.",
    )
    parser.add_argument(
        '-cases', dest='cases',
        nargs="+",
        default=None,
        help="Case names to run (default is the curated no-micro/no-BUGS/no-sponge/no-SILHS set).",
    )
    args, forwarded = parse_forwarded_args(parser, argv)
    reserved = {"-jax", "-python", "-exe", "-driver_test", "-gdb",
                "-output_dir", "-multicol", "-override", "-install_dir"}
    for argument in forwarded:
        if argument.split("=", 1)[0] in reserved:
            parser.error(f"{argument} is controlled by the comparison harness")
    if args.max_iters is not None and args.max_iters < 1:
        parser.error("-max_iters must be positive")
    args.cases = split_values(args.cases) if args.cases is not None else None
    if args.cases == []:
        parser.error("-cases must contain at least one case name")
    args.run_scm_args = tuple(forwarded)
    return args


def main() -> int:
    args = parse_args()

    # Build case -> max_iters mapping.
    cases = list(args.cases) if args.cases else list(DEFAULT_CASES.keys())
    case_iters = {
        case: args.max_iters if args.max_iters is not None
              else DEFAULT_CASES.get(case)
        for case in cases
    }

    repo_root = Path(__file__).resolve().parents[1]
    results_root = (repo_root / RESULTS_DIRNAME).resolve()
    py_output_root = results_root / PYTHON_OUTPUT_DIRNAME
    f90_output_root = results_root / FORTRAN_OUTPUT_DIRNAME
    stats_file = (repo_root / args.stats).resolve()

    if not stats_file.exists():
        print(f"ERROR: stats file not found: {stats_file}")
        return 2

    run_scm = repo_root / "run_scripts" / "run_scm.py"
    run_bindiff = repo_root / "run_scripts" / "run_bindiff_all.py"
    if not run_scm.exists() or not run_bindiff.exists():
        print("ERROR: could not locate run_scripts/run_scm.py or run_bindiff_all.py")
        return 2

    if results_root.exists() and not args.keep_existing:
        shutil.rmtree(results_root)
    py_output_root.mkdir(parents=True, exist_ok=True)
    f90_output_root.mkdir(parents=True, exist_ok=True)

    tasks = [
        (
            case,
            str(repo_root),
            str(stats_file),
            case_iters[case],
            args.bindiff_threshold,
            str(py_output_root),
            str(f90_output_root),
            str(results_root),
            args.run_scm_args,
        )
        for case in cases
    ]

    print(f"Running {len(cases)} case(s) with {args.jobs} worker(s)")
    print(f"Results root: {results_root}")
    print(f"Python output: {py_output_root}")
    print(f"Fortran output: {f90_output_root}")
    print(f"Stats file: {stats_file}")

    start = time.time()
    results: list[CaseResult] = []

    def _emit(result: CaseResult) -> None:
        results.append(result)
        ts_info = ""
        if result.avg_diff_timestep >= 0:
            ts_info = f" avg_diff_ts={result.avg_diff_timestep:.1f}"
        print(
            f"[{result.status:14}] {result.case:22} "
            f"py={result.python_elapsed_s:.1f}s f90={result.fortran_elapsed_s:.1f}s "
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
        "python_failed": sum(r.status == "python_failed" for r in results),
        "fortran_failed": sum(r.status == "fortran_failed" for r in results),
        "elapsed_s": time.time() - start,
        "cases": [asdict(r) for r in results],
    }

    summary_path = results_root / SUMMARY_FILENAME
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    final_bindiff_log = results_root / FINAL_BINDIFF_LOG_FILENAME
    final_bindiff_report = results_root / "final_bindiff.json"
    final_diff_cmd = [
        sys.executable,
        str(run_bindiff),
        "-verbose", str(args.bindiff_verbose),
        "-threshold", str(args.bindiff_threshold),
        "-percent_threshold", str(args.bindiff_threshold),
        "-result_json", str(final_bindiff_report),
        str(py_output_root),
        str(f90_output_root),
    ]
    final_bindiff_report.unlink(missing_ok=True)
    final_bindiff_rc = _run_and_log(final_diff_cmd, repo_root, final_bindiff_log)
    final_avg_ts = _average_earliest_timestep(final_bindiff_report) if final_bindiff_rc != 0 else -1.0

    print("\nSummary:")
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=2))
    if final_avg_ts >= 0:
        print(f"Average earliest diff timestep (across all cases): {final_avg_ts:.1f}")
    print(f"Detailed results: {summary_path}")
    print(f"Final bindiff log: {final_bindiff_log} (rc={final_bindiff_rc})")

    # non-zero when there are diffs or failures
    if summary["diff"] > 0 or summary["python_failed"] > 0 or summary["fortran_failed"] > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()
    raise SystemExit(main())
