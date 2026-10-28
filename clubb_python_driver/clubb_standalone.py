#!/usr/bin/env python3
"""Run complete Python CLUBB cases from aggregate namelists.

Mirrors src/clubb_standalone.F90: each request initializes, advances and cleans
up its model. -worker keeps the process alive for further filenames without
retaining initialized model allocations.

Usage:
    Pass a generated namelist path to run once, or -worker to read one path per
    stdin line until EOF. --quiet suppresses per-timestep progress.
    - main(): select one-shot or worker execution from the command line.
    - worker_main(l_stdout=True): read paths and acknowledge completed cases.
    - run_case(namelist_path, l_stdout=True): run the full lifecycle.

Understanding:
    Generate aggregate inputs with input/scripts/create_case_namelist.py; a model
    fragment alone does not include parameter, physics and stats settings.
    Workers print CLUBB_WORKER_READY once and CLUBB_WORKER_DONE after each
    successful case's output closes. Exceptions exit without acknowledging
    success; utilities/clubb_case_runner.py captures the diagnostic and status.
    Python exceptions replace the native standalone's error-container checks.
    The compiled Python API is required; the existing Python bootstrap prepares
    dependencies when invoked as a command. Select the installed API through the
    case runner or set PYTHONPATH before invoking this module directly.
    This Python SCM does not write timing summaries.
"""
import sys
import time
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

from clubb_python_driver.clubb_case_initalization import (
    clean_up_clubb,
    init_clubb_case,
)
from clubb_python_driver.advance_clubb_to_end import advance_clubb_to_end


def main():
    """Parse frontend options and select one complete case or the filename loop."""
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print("Usage: python -m clubb_python_driver.clubb_standalone <namelist_path> [--quiet]")
        print("  namelist_path: path to the complete runner namelist")
        print("  -worker: read full-case namelist paths from stdin until EOF")
        print("  --quiet: suppress per-timestep output")
        sys.exit(0 if '--help' in sys.argv else 1)

    l_stdout = '--quiet' not in sys.argv
    if sys.argv[1] == '-worker':
        worker_main(l_stdout)
    else:
        run_case(sys.argv[1], l_stdout)


# ------------------------------------------------------------------------------
def worker_main(l_stdout=True):
    """Read namelist paths from stdin until EOF and flush readiness/completion.

    Mirrors worker_main in src/clubb_standalone.F90. l_stdout controls progress
    for every request. A failed run raises before its completion marker.
    """
    print('CLUBB_WORKER_READY', flush=True)
    for line in sys.stdin:
        run_case(line.rstrip('\r\n'), l_stdout)
        print('CLUBB_WORKER_DONE', flush=True)


# ------------------------------------------------------------------------------
def run_case(namelist_path, l_stdout=True):
    """Run namelist_path's full lifecycle and configured output; return None.

    Mirrors run_case's run_clubb call in src/clubb_standalone.F90, using explicit
    Python state. l_stdout controls per-timestep progress. Cleanup runs after
    successful advancement or an advancement exception.
    """
    t0 = time.time()
    # A worker request gets the same fresh setup as an independent invocation.
    state = init_clubb_case(namelist_path)
    try:
        advance_clubb_to_end(state, l_stdout=l_stdout)
    finally:
        # Release this case's resources before the worker accepts another path.
        clean_up_clubb(state)
    elapsed = time.time() - t0

    print(f"Completed {state['ifinal']} timesteps in {elapsed:.1f}s")


if __name__ == '__main__':
    main()
