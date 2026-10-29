"""Source fatal diagnostics retain JIT-compatible, per-column status."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.SILHS.est_kessler_microphys_module import calc_estimate
from clubb_jax.src.SILHS.silhs_importance_sample_module import (
    importance_sampling_driver,
    limit_category_weights,
)


@pytest.mark.parametrize("debug_level", [-1, 0, 2])
def test_impossible_weight_transfer_reports_fatal_status(monkeypatch, debug_level):
    monkeypatch.setattr(error_code, "_debug_level", debug_level)
    real = jnp.full(8, 0.125)
    prescribed = jnp.zeros(8)
    _, error = jax.jit(limit_category_weights)(real, prescribed)
    assert bool(error)


def test_weight_limiter_error_reaches_importance_driver(monkeypatch):
    from clubb_jax.src.SILHS import silhs_importance_sample_module as sampler

    # Inject the local transfer result, so this tests status ownership rather
    # than relying on an invalid probability allocation reaching that branch.
    monkeypatch.setattr(error_code, "_debug_level", -1)
    monkeypatch.setattr(
        sampler, "limit_category_weights", lambda real, prescribed: (prescribed, jnp.array(True))
    )
    samples = jnp.array([0.125, 0.375, 0.625, 0.875])
    run = jax.jit(lambda x: importance_sampling_driver(
        4,                   # In
        0.2, 0.6,            # In
        0.3,                  # In
        0.4, 0.7,            # In
        1, True,             # In
        True, False, False,  # In
        x, x, x,             # InOut
        jax.random.PRNGKey(2),  # In
    ))
    assert bool(run(samples)[-1])


@pytest.mark.parametrize("fraction", ["mixture", "cloud1", "cloud2"])
@pytest.mark.parametrize("invalid", [-0.125, 1.125])
def test_kessler_invalid_fraction_is_local_to_column(fraction, invalid):
    fractions = {
        "mixture": jnp.full((2, 1), 0.5),
        "cloud1": jnp.full((2, 1), 0.25),
        "cloud2": jnp.full((2, 1), 0.75),
    }
    fractions[fraction] = fractions[fraction].at[0, 0].set(invalid)
    rc = jnp.array([0.0, 0.0004, 0.0008, 0.0006])[None, :, None]
    rc = jnp.broadcast_to(rc, (2, 4, 1))
    component = jnp.broadcast_to(jnp.array([1, 2, 1, 2])[None, :, None], rc.shape)
    run = jax.jit(lambda mixture, cloud1, cloud2: calc_estimate(
        4, mixture, cloud1, cloud2, rc,  # In
        component, jnp.ones_like(rc),    # In
        False,                           # In
        1.0e-3, 0.2e-3,                 # In
    ))
    estimate, error = run(fractions["mixture"], fractions["cloud1"], fractions["cloud2"])
    np.testing.assert_array_equal(error[:, 0], [True, False])
    np.testing.assert_allclose(estimate[:, 0], 3.0e-7, rtol=0.0, atol=1.0e-22)


def test_kessler_empty_sample_reports_source_stop():
    empty = jnp.empty((2, 0, 1))
    _, error = calc_estimate(
        0, jnp.full((2, 1), 0.5),     # In
        jnp.ones((2, 1)), jnp.ones((2, 1)), empty,  # In
        empty.astype(jnp.int32), empty,  # In
        False,                           # In
        1.0, 0.0,                       # In
    )
    assert np.all(np.asarray(error))
