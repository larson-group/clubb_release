#!/usr/bin/env python3
"""Prepare and run one CLUBB case through the shared case runner.

Use this command to select native Fortran, the Python SCM or JAX with the same
model options. The namelist utility owns case preparation; clubb_case_runner
owns runtime selection, process execution and output capture.

Usage:
    ./run_scripts/run_scm.py bomex
    ./run_scripts/run_scm.py -jax=cpu -max_iters 20 -output_dir jax_bomex bomex
    main() parses command-line arguments, prepares an aggregate namelist,
    waits for the case, and returns zero on success or one on a run failure.

Understanding:
    Output paths follow utilities/output_paths.py. The case writes stats and
    timing output according to its settings; console output is also captured
    in <case>_log. This command runs one case and closes the runner afterward.
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
from utilities.create_case_namelist import create_case_namelist_from_args  # noqa: E402
from utilities.output_paths import resolve_output_dir  # noqa: E402


def main():
    """Read CLI options, generate one aggregate namelist and run it to completion.

    Return 0 on model success or 1 on preparation/runtime failure. Write the
    generated input and console log under the resolved output directory; the
    selected model writes its configured output. Argument errors exit via argparse.
    """
    args = parse_case_arguments(sys.argv[1:])
    try:
        output_dir = resolve_output_dir(args.out_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {output_dir}")
        # Prepare complete inputs before the runner starts the selected runtime.
        namelist = create_case_namelist_from_args(args, output_dir)
        print(f"=================== Running {args.case_name} ===================")
        # One full case here; scripts needing reuse keep this context open longer.
        with clubb_case_runner(args, echo=True) as runner:
            result = runner.run_and_wait(namelist, log_file=output_dir / f"{args.case_name}_log")
        if result['state'] != 'completed':
            print(result['error'], file=sys.stderr)
            return 1
        # Missing optional output is reported separately from model failure.
        expect_stats = ((args.stats or "").strip().lower() != "none") and args.tout != 0
        stats_file = output_dir / f"{args.case_name}_stats.nc"
        if expect_stats and not stats_file.is_file():
            print(f"WARNING: stats output not found: {stats_file}")
        return 0
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    if '-python' in sys.argv[1:]:
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
