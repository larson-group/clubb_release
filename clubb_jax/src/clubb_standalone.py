#!/usr/bin/env python3
"""Python CLUBB standalone — replaces clubb_standalone.F90.

Usage:
    python -m clubb_jax.src.clubb_standalone input/case_setups/bomex_model.in
    python -m clubb_jax.src.clubb_standalone input/case_setups/bomex_model.in --quiet
"""
import os
import sys
import time

import jax
from clubb_jax.src.clubb_case_initalization import (
    clean_up_clubb,
    init_clubb_case,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print("Usage: python -m clubb_jax.src.clubb_standalone <namelist_path> [--quiet]")
        print("  namelist_path: path to *_model.in file")
        print("  --quiet: suppress per-timestep output")
        sys.exit(0 if '--help' in sys.argv else 1)

    namelist_path = sys.argv[1]
    l_stdout = '--quiet' not in sys.argv

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

    t0 = time.time()
    state = init_clubb_case(namelist_path)
    advance_clubb_to_end(state, l_stdout=l_stdout)
    clean_up_clubb(state)
    elapsed = time.time() - t0

    print(f"Completed {state['ifinal']} timesteps in {elapsed:.1f}s")


if __name__ == '__main__':
    main()
