#!/usr/bin/env python3
"""Explore SCM stability at several timesteps with either standalone runner.

This is the existing native sweep around run_scripts/run_scm.py, which runs
src/clubb_standalone.F90 and src/clubb_driver.F90 (run_clubb). Passing -jax
selects their JAX counterparts within the same workflow.

For each case, use the native debug level 0, disable statistics output and set
equal main/radiation timesteps. Report each successful run and stop that case
at its first failed timestep. Continue with other cases and retain the native
exploratory exit status; this test reports stability limits rather than a
convergence threshold. Optional case/timestep lists, a step limit and output
root support bounded backend checks without changing the native defaults.
Four case workers run expensive cases first; each case's timesteps remain in
order. Complete case reports are printed as workers finish.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import sys
import threading
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run_scripts.run_scm import extract_jax_options
from utilities.output_paths import resolve_output_dir


CLUBB_ROOT = Path(__file__).resolve().parents[1]
RUN_SCM = CLUBB_ROOT / "run_scripts" / "run_scm.py"
DEFAULT_CASES = (
    "arm", "arm_97", "astex_a209", "atex", "bomex", "cgils_s6", "cgils_s11", "cgils_s12",
    "clex9_nov02", "clex9_oct14", "dycoms2_rf01", "dycoms2_rf01_fixed_sst",
    "dycoms2_rf02_do", "dycoms2_rf02_ds", "dycoms2_rf02_nd", "dycoms2_rf02_so",
    "fire", "gabls2", "gabls3", "gabls3_night", "jun25_altocu", "lba", "mc3e", "mpace_a",
    "mpace_b", "mpace_b_silhs", "nov11_altocu", "rico", "rico_silhs", "twp_ice", "wangara",
)

# Scheduling hints from the full CPU sweep in Jenkins build #22, in seconds
# across all five timesteps. Unmeasured cases retain their input order.
CASE_COST_SECONDS = {
    "rico": 440, "lba": 410, "nov11_altocu": 304, "clex9_oct14": 299,
    "dycoms2_rf02_do": 225, "dycoms2_rf02_ds": 219, "rico_silhs": 113,
    "gabls2": 106,
}


def parse_args() -> argparse.Namespace:
    """Retain the native sweep defaults and allow a focused JAX run."""
    parser = argparse.ArgumentParser(description=__doc__, add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        "-cases", default=",".join(DEFAULT_CASES), help="Comma-separated cases",
    )
    parser.add_argument(
        "-timesteps", default="600,1200,1800,2400,3000",
        help="Comma-separated dt_main/dt_rad values in seconds",
    )
    parser.add_argument("-max_iters", type=int, help="Optional iteration limit for each run")
    parser.add_argument(
        "-workers", type=int, default=4,
        help="Concurrent case workers; each case's timesteps remain ordered. Default: 4.",
    )
    parser.add_argument(
        "-output_root", dest="out_root", type=Path,
        help="Output root for separate case/timestep directories; bare names go under output/",
    )
    parser.add_argument(
        "-jax", action="store_true",
        help="Use JAX; an attached -jax=VALUE is forwarded unchanged to the launcher",
    )
    argv, jax_options, jax_occurrences = extract_jax_options(sys.argv[1:])
    args = parser.parse_args(argv)
    if jax_occurrences > 1:
        parser.error("-jax may be specified only once")
    args.jax_options = jax_options
    args.cases = [case.strip() for case in args.cases.split(",") if case.strip()]
    if not args.cases:
        parser.error("-cases must include at least one case")
    if len(args.cases) != len(set(args.cases)):
        parser.error("-cases must not repeat a case")
    try:
        args.timesteps = [int(value.strip()) for value in args.timesteps.split(",")]
    except ValueError:
        parser.error("-timesteps must contain comma-separated integers")
    if not args.timesteps or any(dt <= 0 for dt in args.timesteps):
        parser.error("-timesteps must contain positive integers")
    if args.max_iters is not None and args.max_iters <= 0:
        parser.error("-max_iters must be positive")
    if args.workers <= 0:
        parser.error("-workers must be positive")
    if args.out_root is not None:
        try:
            args.out_root = resolve_output_dir(args.out_root)
        except ValueError as exc:
            parser.error(str(exc))
    return args


def run_case(case: str, args: argparse.Namespace, stop: threading.Event) -> str:
    """Run one ordered sweep and return its complete, non-interleaved report."""
    lines = [f"---------------- Running {case} ----------------"]
    for dt in args.timesteps:
        if stop.is_set():
            break
        cmd = [
            sys.executable, str(RUN_SCM),
            "-debug", "0", "-tout", "0",
            "-dt_main", str(dt), "-dt_rad", str(dt),
        ]
        if args.out_root is not None:
            cmd.extend(["-output_dir", str(args.out_root / case / f"dt{dt}")])
        if args.jax:
            cmd.append("-jax" if args.jax_options is None else f"-jax={args.jax_options}")
        if args.max_iters is not None:
            cmd.extend(["-max_iters", str(args.max_iters)])
        cmd.append(case)
        result = subprocess.run(
            cmd,
            cwd=CLUBB_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
        if result.returncode != 0:
            lines.append(result.stdout.rstrip())
            lines.append(f"---- FAIL @ dt = {dt}")
            break
        lines.append(f"--- PASS @ dt = {dt}")
    return "\n".join(lines)


def main() -> int:
    """Run expensive cases first and report each completed case as one block."""
    args = parse_args()
    cases = sorted(args.cases, key=lambda case: -CASE_COST_SECONDS.get(case, 0))
    workers = min(args.workers, len(cases))
    print(f"Running {len(cases)} case sweeps with {workers} workers", flush=True)
    stop = threading.Event()
    # Threads only orchestrate isolated model subprocesses; they do not import JAX.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_case, case, args, stop) for case in cases]
        try:
            for future in as_completed(futures):
                print(future.result(), flush=True)
        except BaseException:
            stop.set()
            for future in futures:
                future.cancel()
            raise
    # As in the native script, an unstable timestep identifies the case's
    # limit; it does not make this exploratory test a failed Jenkins build.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
