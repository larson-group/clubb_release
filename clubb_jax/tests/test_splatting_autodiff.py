"""Splatting gradients at clipped frequencies and a zero tuning coefficient."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.advance_helper_module import wp23_term_splat_lhs
from clubb_jax.src.CLUBB_core.grid_class import setup_grid, zm2zt, zm2zt2zm

jax.config.update("jax_enable_x64", True)


def _weighted_loss(outputs):
    return sum(
        jnp.sum(jnp.linspace(1.0, 2.0, y.shape[1])[None, :] * (y + y**2))
        for y in outputs[:2]
    )


@pytest.mark.parametrize("coefficient", [0.0, 0.25])
@pytest.mark.parametrize("frequency_squared", [-4e-4, 0.0])
def test_clipped_frequency_has_zero_forward_and_reverse_derivatives(
    coefficient, frequency_squared,
):
    # Spacing > 60 m makes each averaging window contain only its own level,
    # with no below-ground samples. All frequencies stay on the clipped branch.
    gr = setup_grid(ngrdcol=1, deltaz=100.0, zm_init=0.0, zm_top=800.0, grid_type=1)
    n2 = jnp.full((1, gr.nzm), frequency_squared)
    density = jnp.linspace(0.7, 1.2, gr.nzm)[None, :]

    def fn(x):
        return wp23_term_splat_lhs(
            gr.nzm, gr.nzt, 1, gr, jnp.array([coefficient]),
            x, jnp.full_like(x, 100.0), density,
        )[:2]

    direction = jnp.linspace(-0.5, 1.0, gr.nzm)[None, :]
    outputs, tangents = jax.jit(lambda x, v: jax.jvp(fn, (x,), (v,)))(n2, direction)
    gradient = jax.jit(jax.grad(lambda x: _weighted_loss(fn(x))))(n2)
    for value in (*outputs, *tangents, gradient):
        np.testing.assert_array_equal(value, jnp.zeros_like(value))


@pytest.mark.parametrize("coefficient", [0.0, 0.25])
def test_mixed_columns_preserve_values_and_coefficient_sensitivity(coefficient):
    gr = setup_grid(ngrdcol=2, deltaz=40.0, zm_init=0.0, zm_top=1200.0, grid_type=1)
    profile = jnp.where(jnp.arange(gr.nzm) < gr.nzm // 2, -3e-4, 8e-4)
    n2 = jnp.stack((profile, profile[::-1]))
    density = jnp.broadcast_to(jnp.linspace(1.2, 0.7, gr.nzm), n2.shape)
    coefficients = jnp.array([coefficient, 0.5])

    def fn(x, rho, c):
        return wp23_term_splat_lhs(
            gr.nzm, gr.nzt, 2, gr, c, x, jnp.full_like(x, 100.0), rho,
        )

    outputs = fn(n2, density, coefficients)
    # Preserve the original square-root formula and interpolation for the
    # forward answer, including positive and negative averaged frequencies.
    averaged = outputs[2]
    assert jnp.any(averaged < 0) and jnp.any(averaged > 0)
    assert jnp.min(jnp.abs(averaged)) > 1e-6
    clipped = jnp.sqrt(jnp.maximum(0.0, averaged))
    expected_wp2 = coefficients[:, None] * zm2zt2zm(gr.nzm, gr.nzt, 2, gr, clipped)
    expected_wp3 = 1.5 * coefficients[:, None] * zm2zt(gr.nzm, gr.nzt, 2, gr, clipped)
    np.testing.assert_allclose(outputs[0], expected_wp2, rtol=2e-15, atol=1e-16)
    np.testing.assert_allclose(outputs[1], expected_wp3, rtol=2e-15, atol=1e-16)

    def loss(x, rho, c):
        return _weighted_loss(fn(x, rho, c))

    direction = (
        1e-4 * jnp.cos(jnp.arange(n2.size)).reshape(n2.shape),
        0.1 * jnp.sin(jnp.arange(n2.size)).reshape(n2.shape),
        jnp.array([0.7, -0.3]),
    )
    inputs = (n2, density, coefficients)
    _, tangent = jax.jit(lambda x, v: jax.jvp(loss, x, v))(inputs, direction)
    gradients = jax.jit(jax.grad(loss, argnums=(0, 1, 2)))(*inputs)
    for gradient in gradients:
        assert np.isfinite(gradient).all()
    eps = 1e-4
    plus = tuple(x + eps * v for x, v in zip(inputs, direction))
    minus = tuple(x - eps * v for x, v in zip(inputs, direction))
    finite_difference = (loss(*plus) - loss(*minus)) / (2 * eps)
    np.testing.assert_allclose(tangent, finite_difference, rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(
        sum(jnp.sum(g * v) for g, v in zip(gradients, direction)),
        finite_difference, rtol=1e-7, atol=1e-10,
    )
    if coefficient == 0.0:
        np.testing.assert_array_equal(outputs[0][0], 0.0)
        np.testing.assert_array_equal(outputs[1][0], 0.0)
        np.testing.assert_array_equal(gradients[0][0], 0.0)
        np.testing.assert_array_equal(gradients[1][0], 0.0)
        # A shortcut returning zero when C=0 would lose the derivative for
        # turning splatting on. The expected slope is the unit-C response.
        unit = fn(n2, density, jnp.ones(2))
        expected_slope = sum(
            jnp.sum(jnp.linspace(1.0, 2.0, y.shape[1]) * y[0]) for y in unit[:2]
        )
        assert expected_slope > 0
        np.testing.assert_allclose(gradients[2][0], expected_slope, rtol=1e-13)
