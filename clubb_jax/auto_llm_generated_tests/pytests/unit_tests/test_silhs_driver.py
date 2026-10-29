"""SILHS subtimestep statistics keep the previous window and count it once."""

import jax
import jax.numpy as jnp
import numpy as np
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats


def test_subtimestep_statistics_preserve_previous_window_and_count_once():
    previous = JaxStats.empty(l_sample=True, names=("rrm_auto",), ncol=2, max_nlev=3)
    previous = previous.update("rrm_auto", jnp.full((2, 3), 10.0))

    def run(st):
        for value in (1.0, 2.0, 6.0):
            st = st.update("rrm_auto", jnp.full((2, 3), value))
        return st.average_subtimesteps(previous, 3)

    stats = jax.jit(run)(previous)
    np.testing.assert_array_equal(stats.buffers[0], 13.0)
    np.testing.assert_array_equal(stats.nsamples[0], 2)
