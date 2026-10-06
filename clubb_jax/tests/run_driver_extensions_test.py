#!/usr/bin/env python3
"""Additional JAX window and parameter-batch checks beyond the native driver test.

Run from the repository root with:
    python3 clubb_jax/tests/run_driver_extensions_test.py

This JAX-only extension uses eight 60-second BOMEX timesteps with statistics
output disabled, four distinct C8 parameter sets and runtime width two. It
requires exact equality of five final state fields between uninterrupted and
successive 1:4/5:8 windows. Selecting the second batch must activate the saved
parameter rows, complete without a fatal error and change a checked field.

These assertions exercise src/clubb_driver.F90's reset/batch behavior and
src/clubb_loss_driver.F90's successive-window calls. The native reinitialization
and double-timestep sequence is mirrored separately by
clubb_jax/src/clubb_driver_test.py; its comparison uses saved NetCDF statistics.
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
from utilities.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import (
    init_clubb_case, set_case_initial_conditions, clean_up_clubb,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end


def main():
    with tempfile.TemporaryDirectory(prefix='clubb-lifecycle-') as folder:
        path = create_case_namelist_file(
            'bomex', Path(folder), stats='none', debug='-1',
            multicol='C8/0.2:0.8/4', batch_size=2, max_iters=8,
            dt_main=60, dt_rad=60,
        )
        state = init_clubb_case(str(path))
        fields = ('thlm', 'rtm', 'wp2', 'rtp2', 'thlp2')
        try:
            assert state['total_param_sets'] == 4 and state['ngrdcol'] == 2
            set_case_initial_conditions(state, batch_num=1)
            advance_clubb_to_end(state, l_stdout=False)
            expected = {name: np.asarray(state[name]).copy() for name in fields}
            set_case_initial_conditions(state, batch_num=1)
            advance_clubb_to_end(state, False, itime_start=1, itime_end=4)
            advance_clubb_to_end(state, False, itime_start=5, itime_end=8)
            assert not state['err_info'].is_fatal()
            for name in fields:
                np.testing.assert_array_equal(state[name], expected[name], err_msg=name)
            set_case_initial_conditions(state, batch_num=2)
            np.testing.assert_array_equal(
                state['clubb_params'], state['clubb_params_all'][2:4],
            )
            advance_clubb_to_end(state, l_stdout=False)
            assert not state['err_info'].is_fatal()
            assert any(np.any(np.asarray(state[name]) != expected[name]) for name in fields)
        finally:
            clean_up_clubb(state)
    print('JAX extensions: window equivalence and runtime batches pass')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
