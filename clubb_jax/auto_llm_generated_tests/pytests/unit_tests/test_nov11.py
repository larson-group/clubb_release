from types import SimpleNamespace

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.Benchmark_cases.nov11 import nov11_altocu_rtm_adjust


def test_nov11_adjustment_is_jittable_and_preserves_dtype():
    gr = SimpleNamespace(
        zt=jnp.array([[2800.0, 3001.0], [3000.0, 3101.0]], dtype=jnp.float32),
        zm=jnp.array([[0.0], [100.0]], dtype=jnp.float32),
    )
    rtm = jnp.ones((2, 2), dtype=jnp.float32)
    adjust = jax.jit(
        lambda values, time: nov11_altocu_rtm_adjust(2, gr, time, 0.0, 300.0, values)
    )

    before = adjust(rtm, jnp.array(3599.0, dtype=jnp.float32))
    during = adjust(rtm, jnp.array(3600.0, dtype=jnp.float32))

    assert before.dtype == rtm.dtype
    assert during.dtype == rtm.dtype
    np.testing.assert_array_equal(before, np.ones((2, 2), dtype=np.float32))
    np.testing.assert_allclose(
        during,
        np.array([[1.0, 0.89], [1.0, 0.89]], dtype=np.float32),
        rtol=0.0,
        atol=0.0,
    )
