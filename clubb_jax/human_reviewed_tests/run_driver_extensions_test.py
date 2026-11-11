#!/usr/bin/env python3
"""Check JAX driver behavior that the clubb_driver bindiff runs cannot see.

Run from the repository root with:
    python3 clubb_jax/human_reviewed_tests/run_driver_extensions_test.py

The script initializes one short BOMEX case: eight 60-second timesteps,
statistics output disabled, and four parameter sets that differ only in C8
(0.2, 0.4, 0.6, 0.8). The model runs two columns at a time, so the four sets
form two runtime batches: batch 1 holds C8 = 0.2, 0.4 and batch 2 holds
C8 = 0.6, 0.8. state['clubb_params_all'] stores all four sets and
state['clubb_params'] the two sets of the active batch.

Two independent checks run on that case:

1. Split-window equivalence (check_split_window_equivalence). Running steps
   1-8 in one call must give bit-identical final fields to running steps 1-4
   and then 5-8 in two calls. The JAX loss driver (clubb_loss_driver.py,
   following src/clubb_loss_driver.F90) advances the model window by window
   this way, so a split must not change the result.
2. Batch selection (check_batch_selection). Resetting with batch_num=2 must
   load parameter rows 3-4 into the active parameters, and the run must then
   finish and give different results from batch 1. Both runs compared by the
   clubb_driver bindiff use the same batch selection code, so they would
   agree even if it picked the wrong rows; this check compares against the
   stored parameter table directly.

Both checks compare the final thlm, rtm, wp2, rtp2 and thlp2. Exact equality
is expected because both paths run the same code on the same device. The
native reinitialization and repeated-advance sequence of
src/clubb_driver_test.F90 is mirrored separately by
clubb_jax/src/clubb_driver_test.py and checked by the clubb_driver bindiff.
"""
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == '__main__':
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

import numpy as np
from input.scripts.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import (
    init_clubb_case, set_case_initial_conditions, clean_up_clubb,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end

# Final model fields compared by both checks.
FIELDS = ('thlm', 'rtm', 'wp2', 'rtp2', 'thlp2')


def copy_fields(state):
    """Return host copies of the compared FIELDS from the current model state."""
    return {name: np.asarray(state[name]).copy() for name in FIELDS}


def check_split_window_equivalence(state):
    """Check that steps 1-8 in one call equal steps 1-4 then 5-8 in two calls.

    Both runs use batch 1 (C8 = 0.2, 0.4). Returns the batch 1 fields from the
    uninterrupted run, which check_batch_selection uses as its reference.
    """
    # Reference run: reset to the initial conditions of batch 1 and advance
    # through all eight timesteps in one call.
    set_case_initial_conditions(state, batch_num=1)
    advance_clubb_to_end(state, l_stdout=False)
    batch1_fields = copy_fields(state)

    # Windowed run: reset to the same initial conditions and parameters, then
    # advance through the same eight timesteps as two windows, the way the loss
    # driver does. The second call continues from the state the first left.
    set_case_initial_conditions(state, batch_num=1)
    advance_clubb_to_end(state, False, itime_start=1, itime_end=4)
    advance_clubb_to_end(state, False, itime_start=5, itime_end=8)
    assert not state['err_info'].is_fatal()

    # Same parameters and steps, so the split must not change any bit.
    for name in FIELDS:
        np.testing.assert_array_equal(state[name], batch1_fields[name], err_msg=name)
    return batch1_fields


def check_batch_selection(state, batch1_fields):
    """Check that batch 2 loads parameter rows 3-4 and changes the results.

    batch1_fields holds the final fields of the uninterrupted batch 1 run.
    """
    # Reset to the initial conditions and select batch 2 (C8 = 0.6, 0.8). The
    # active parameters must be rows 3-4 (clubb_params_all[2:4]) of the stored
    # table of all four sets.
    set_case_initial_conditions(state, batch_num=2)
    np.testing.assert_array_equal(state['clubb_params'], state['clubb_params_all'][2:4])

    # Run batch 2 through all eight timesteps. Different C8 values must change
    # at least one field; identical results would mean batch 2 had no effect.
    advance_clubb_to_end(state, l_stdout=False)
    assert not state['err_info'].is_fatal()
    assert any(np.any(np.asarray(state[name]) != batch1_fields[name]) for name in FIELDS)


def main():
    """Initialize the BOMEX case once, run both checks, then clean up."""
    with tempfile.TemporaryDirectory(prefix='clubb-lifecycle-') as folder:
        # 4 C8 values run 2 columns at a time: 2 batches of 2 parameter sets.
        path = create_case_namelist_file(
            'bomex', Path(folder), stats='none', debug='-1',
            multicol='C8/0.2:0.8/4', batch_size=2, max_iters=8,
            dt_main=60, dt_rad=60,
        )
        state = init_clubb_case(str(path))
        try:
            assert state['total_param_sets'] == 4 and state['ngrdcol'] == 2
            batch1_fields = check_split_window_equivalence(state)
            check_batch_selection(state, batch1_fields)
        finally:
            clean_up_clubb(state)
    print('JAX extensions: window equivalence and runtime batches pass')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
