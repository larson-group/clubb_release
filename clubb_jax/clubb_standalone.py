#!/usr/bin/env python3
"""Python CLUBB standalone — replaces clubb_standalone.F90.

Usage:
    python -m clubb_jax.clubb_standalone input/case_setups/bomex_model.in
    python -m clubb_jax.clubb_standalone input/case_setups/bomex_model.in --quiet
"""
import sys
import time
from clubb_jax.clubb_case_initalization import (
    clean_up_clubb,
    init_clubb_case,
)
from clubb_jax.advance_clubb_to_end import advance_clubb_to_end


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print("Usage: python -m clubb_jax.clubb_standalone <namelist_path> [--quiet]")
        print("  namelist_path: path to *_model.in file")
        print("  --quiet: suppress per-timestep output")
        sys.exit(0 if sys.argv[1:] and sys.argv[1] in ('-h', '--help') else 1)

    namelist_path = sys.argv[1]
    l_stdout = '--quiet' not in sys.argv

    t0 = time.time()
    state = init_clubb_case(namelist_path)
    advance_clubb_to_end(state, l_stdout=l_stdout)
    clean_up_clubb(state)
    elapsed = time.time() - t0

    print(f"Completed {state['ifinal']} timesteps in {elapsed:.1f}s")


if __name__ == '__main__':
    main()
