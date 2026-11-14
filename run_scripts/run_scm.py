#!/usr/bin/env python3
"""Prepare and run CLUBB cases, groups or named lists through the shared case runner.

Use this command to select native Fortran, the Python SCM or JAX with the same
model options. The namelist utility owns case preparation; clubb_case_runner
owns runtime selection, process execution and output capture.

Usage:
    ./run_scripts/run_scm.py bomex
    ./run_scripts/run_scm.py -jax=cpu -max_iters 20 -output_dir jax_bomex bomex
    ./run_scripts/run_scm.py -cases priority_cases -workers 4
    ./run_scripts/run_scm.py -cases clex_group
    run_case() prepares and runs one case; main() expands selections and returns
    a combined failure status after all scheduled runs complete.

Understanding:
    Output paths follow utilities/output_paths.py. The case writes stats and
    timing output according to its settings; console output is also captured
    in <case>_log. Single runs bypass worker pools and preparation subprocesses.
    Scripts that need several full runs in one process can use the runner
    directly. The JAX launcher owns its environment and device modifiers.
"""

from pathlib import Path
import subprocess
import sys
import time

CLUBB_ROOT = str(Path(__file__).resolve().parents[1])
if CLUBB_ROOT not in sys.path:
    sys.path.insert(0, CLUBB_ROOT)

from utilities.clubb_case_runner import clubb_case_runner, parse_case_arguments  # noqa: E402
from input.scripts.create_case_namelist import create_case_namelist_from_args  # noqa: E402
from utilities.output_paths import resolve_output_dir  # noqa: E402


def run_case(args, output_dir, *, echo=True):
    """Prepare and execute one concrete case directly; return its exit status."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        namelist = create_case_namelist_from_args(args, output_dir)
        with clubb_case_runner(args, echo=echo) as runner:
            result = runner.run_and_wait(namelist, log_file=output_dir / f"{args.case_name}_log")
        if result["state"] != "completed":
            print(f"ERROR -- {args.case_name}: {result['error']}", file=sys.stderr)
            return 1
        expect_stats = ((args.stats or "").strip().lower() != "none") and args.tout != 0
        if expect_stats and not (output_dir / f"{args.case_name}_stats.nc").is_file():
            print(f"WARNING: stats output not found for {args.case_name}")
        return 0
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR -- {args.case_name}: {exc}", file=sys.stderr)
        return 1


def main():
    """Resolve selections before scheduling; bypass the pool for a single run."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from copy import copy
    from input.scripts.input_catalog import load_catalog, resolve_selection, catalog_index
    from input.scripts.create_case_namelist import case_run_name

    args = parse_case_arguments(sys.argv[1:], multiple_cases=True)
    try:
        cases = resolve_selection("cases", args.case_selections)
        output_dir = resolve_output_dir(args.out_dir)
        flag_catalog = load_catalog("flags")
        _, groups, _ = catalog_index("flags", flag_catalog)
        sweep = args.flags in {*groups, *flag_catalog.get("lists", {})}
        flags = resolve_selection("flags", args.flags) if sweep else [args.flags]
        tasks = []
        for flag in flags:
            for case in cases:
                concrete = copy(args)
                concrete.case_selector = case
                concrete.case_name = case_run_name(case)
                concrete.flags = flag
                destination = output_dir / flag.replace("/", "__") if sweep else output_dir
                tasks.append((concrete, destination))
        # Reject aliases that would race on the same generated input or model output.
        destinations = [(str(path), item.case_name) for item, path in tasks]
        if len(set(destinations)) != len(destinations):
            raise ValueError("Selected cases have colliding output names")
        if len(tasks) == 1:
            return run_case(*tasks[0])
        workers = min(args.workers, len(tasks))
        print(f"Running {len(tasks)} cases with up to {workers} workers")
        status = 0
        if workers == 1:
            results = ((item.case_name, run_case(item, path, echo=args.show_output))
                       for item, path in tasks)
            for name, code in results:
                print(f"{'ERROR' if code else 'COMPLETE'} -- {name}")
                status = max(status, code)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(run_case, item, path, echo=args.show_output): item.case_name
                           for item, path in tasks}
                for future in as_completed(futures):
                    name = futures[future]
                    try:
                        code = future.result()
                    except Exception as exc:
                        print(f"ERROR -- {name}: {exc}", file=sys.stderr)
                        code = 1
                    print(f"{'ERROR' if code else 'COMPLETE'} -- {name}")
                    status = max(status, code)
        return status
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    if any(token == '-python' or token.split('=', 1)[0] == '-multicol' for token in sys.argv[1:]):
        from utilities.setup_python_venv import ensure_python_venv
        ensure_python_venv()
    start_time = time.perf_counter()
    try:
        exit_code = main()
    finally:
        elapsed_time = time.perf_counter() - start_time
        print('-' * 50)
        print(f'run_scm.py total runtime: {elapsed_time:.2f} s')
        print('-' * 50)
    sys.exit(exit_code)
