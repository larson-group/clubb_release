#!/usr/bin/env python3
"""Run complete JAX CLUBB cases from aggregate namelists.

Mirrors src/clubb_standalone.F90: each request initializes, advances and cleans
up its model. -worker keeps the process alive for further filenames without
retaining initialized model allocations.

Usage:
    Pass a generated namelist path to run once, or -worker to read one path per
    stdin line until EOF. --quiet suppresses per-timestep progress.
    - main(): select one-shot or worker execution from the command line.
    - worker_main(l_stdout=True): read paths and acknowledge completed cases.
    - run_case(namelist_path, l_stdout=True, wait_for_completion=False): run the
      full lifecycle, optionally waiting for submitted JAX work.

Understanding:
    Generate aggregate inputs with utilities/create_case_namelist.py; a model
    fragment alone does not include parameter, physics and stats settings.
    Workers print CLUBB_WORKER_READY once and CLUBB_WORKER_DONE after each
    successful case's output closes. Exceptions exit without acknowledging
    success; utilities/clubb_case_runner.py captures the diagnostic and status.
    Python exceptions replace the native standalone's error-container checks.
    JAX's run_jax.py launcher owns environment/device selection. Worker requests
    wait for outstanding device computations and callbacks before cleanup and
    acknowledgment; the one-shot path preserves its existing asynchronous
    behavior. Compiled programs can survive cleanup of the per-case state.
"""
import os
import sys
import time

import jax
from clubb_jax.src.clubb_case_initalization import (
    clean_up_clubb,
    init_clubb_case,
    set_case_initial_conditions,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end


def main():
    """Check CLI/device selection and run one case or the persistent filename loop."""
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print("Usage: python -m clubb_jax.src.clubb_standalone <namelist_path> [--quiet]")
        print("  namelist_path: path to *_model.in file")
        print("  -worker: read full-case namelist paths from stdin until EOF")
        print("  --quiet: suppress per-timestep output")
        sys.exit(0 if '--help' in sys.argv else 1)

    if os.environ.get('CLUBB_JAX_PROFILE') == 'gpu':
        devices = jax.devices()
        accelerator = os.environ.get('CLUBB_JAX_ACCELERATOR', 'cuda13')
        expected_backend = 'metal' if accelerator == 'metal' else 'gpu'
        if jax.default_backend().lower() != expected_backend:
            raise RuntimeError(
                f'The GPU profile did not initialize the {expected_backend} backend.'
            )
        print('==> JAX GPU backend initialized; visible devices:', flush=True)
        for device in devices:
            print(f'    JAX GPU {device.id}: {device.device_kind}', flush=True)

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
        run_case(line.rstrip('\r\n'), l_stdout, wait_for_completion=True)
        print('CLUBB_WORKER_DONE', flush=True)


# ------------------------------------------------------------------------------
def run_case(namelist_path, l_stdout=True, *, wait_for_completion=False):
    """Run namelist_path's full lifecycle and configured output; return None.

    Mirrors run_case's run_clubb call in src/clubb_standalone.F90, using explicit
    Python state. l_stdout controls per-timestep progress. Cleanup runs after
    successful advancement or an advancement exception.
    Runtime batches reset initial conditions before each parameter slice.
    wait_for_completion fences device work/callbacks for worker acknowledgment.
    """
    t0 = time.time()
    # A worker request gets the same fresh setup as an independent invocation.
    state = init_clubb_case(namelist_path)
    try:
        num_batches = state['total_param_sets'] // state['ngrdcol']
        for batch_num in range(1, num_batches + 1):
            set_case_initial_conditions(state, batch_num=batch_num)
            advance_clubb_to_end(state, l_stdout=l_stdout)
        if wait_for_completion:
            # A worker acknowledgment must include outstanding device work.
            jax.block_until_ready(state)
            jax.effects_barrier()
    finally:
        # Release this case's resources before the worker accepts another path.
        clean_up_clubb(state)
    elapsed = time.time() - t0

    print(f"Completed {state['ifinal']} timesteps in {elapsed:.1f}s")


if __name__ == '__main__':
    main()
