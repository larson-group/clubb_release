"""Forward/reverse derivatives of conservative fills, including inactive windows."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.fill_holes import (
    fill_holes_global,
    fill_holes_sliding_window,
    fill_holes_vertical,
)

jax.config.update("jax_enable_x64", True)


@pytest.mark.parametrize("fill", [fill_holes_global, fill_holes_sliding_window])
@pytest.mark.parametrize("threshold", [0.0, 0.125])
def test_at_threshold_is_identity_in_both_ad_modes(fill, threshold):
    # A whole column at the floor used to evaluate a discarded 0/0. Include
    # an active column in the same call to exercise per-column masking.
    field = jnp.full((2, 20), threshold)
    field = field.at[1, :].add(0.5).at[1, 8].set(threshold - 0.25)
    rho = jnp.ones_like(field)
    dz = jnp.ones_like(field)
    tangent = jnp.linspace(0.2, 1.0, 20)[None, :]
    tangent = jnp.concatenate((tangent, jnp.zeros_like(tangent)))
    weights = jnp.arange(1.0, 21.0)[None, :]

    def fn(x):
        return fill(20, 2, threshold, 0, 19, dz, rho, x)

    result, jvp = jax.jit(lambda x, v: jax.jvp(fn, (x,), (v,)))(field, tangent)
    grad = jax.jit(jax.grad(lambda x: jnp.sum(fn(x) * weights)))(field)
    np.testing.assert_array_equal(result[0], field[0])
    np.testing.assert_allclose(jvp, tangent, rtol=0, atol=0)
    assert np.isfinite(grad).all()
    np.testing.assert_allclose(grad[0], weights[0], rtol=0, atol=0)


@pytest.mark.parametrize("fill_type", [1, 2])
@pytest.mark.parametrize("grid_dir", [1, -1])
@pytest.mark.parametrize("has_hole", [False, True])
def test_zero_windows_preserve_mass_and_its_derivatives(fill_type, grid_dir, has_hole):
    # Sliding windows consisting only of zeros occur even when the rest of
    # the column is nonzero. A distant hole also exercises the active fill.
    field = jnp.array([[0.0] * 8 + [0.5] * 12])
    if has_hole:
        field = field.at[0, 8].set(-0.25)
    rho = jnp.linspace(0.7, 1.2, 20)[None, :]
    dz = jnp.linspace(20.0, 80.0, 20)[None, :]
    if grid_dir == -1:
        field, rho, dz = (jnp.flip(x, axis=1) for x in (field, rho, dz))
        dz = -dz
    lower, upper = (0, 19) if grid_dir == 1 else (19, 0)
    weights = rho * dz
    tangent = jnp.linspace(0.2, 1.0, 20)[None, :]

    def fn(x):
        return fill_holes_vertical(
            20, 1, 0.0, lower, upper, dz, rho, grid_dir, fill_type, x,
        )

    result, jvp = jax.jit(lambda x, v: jax.jvp(fn, (x,), (v,)))(field, tangent)
    grad = jax.jit(jax.grad(lambda x: jnp.sum(fn(x) * weights)))(field)
    assert np.isfinite(jvp).all()
    np.testing.assert_allclose(grad, weights, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(jnp.sum(result * weights), jnp.sum(field * weights), rtol=1e-13)
    np.testing.assert_allclose(jnp.sum(jvp * weights), jnp.sum(tangent * weights), rtol=1e-13)
    if not has_hole:
        np.testing.assert_array_equal(result, field)
        np.testing.assert_array_equal(jvp, tangent)
    else:
        assert np.min(result) >= -1e-14


@pytest.mark.parametrize("fill_type", [1, 2])
@pytest.mark.parametrize("mean_sign", [1.0, -1.0])
def test_active_fill_derivatives_match_finite_differences(fill_type, mean_sign):
    # Exercise both average-above-floor and insufficient-mass branches.
    field = mean_sign * jnp.array([[0.5] * 8 + [-0.3] + [0.6] * 11])
    rho = jnp.linspace(0.8, 1.2, 20)[None, :]
    dz = jnp.linspace(30.0, 70.0, 20)[None, :]
    weights = jnp.linspace(1.0, 3.0, 20)[None, :]
    tangent = jnp.linspace(0.1, 0.9, 20)[None, :]

    def fn(x):
        return fill_holes_vertical(20, 1, 0.0, 0, 19, dz, rho, 1, fill_type, x)

    def loss(x):
        # A spatially weighted nonlinear loss detects incorrect redistribution
        # derivatives that a conserved-mass objective alone would miss.
        return jnp.sum(weights * fn(x) ** 2)

    _, jvp = jax.jit(lambda x, v: jax.jvp(loss, (x,), (v,)))(field, tangent)
    grad = jax.jit(jax.grad(loss))(field)
    step = 1e-6
    fd = (loss(field + step * tangent) - loss(field - step * tangent)) / (2 * step)
    assert np.isfinite(grad).all()
    np.testing.assert_allclose(jvp, fd, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(jnp.sum(grad * tangent), fd, rtol=1e-7, atol=1e-9)


@pytest.mark.parametrize("fill", [fill_holes_global, fill_holes_sliding_window])
def test_known_conservative_forward_answer(fill):
    # One full five-level window: clipping adds one unit, so the four donors
    # each give back 0.25. Boundaries outside the fill range stay untouched.
    field = jnp.array([[9.0, -1.0, 2.0, 2.0, 2.0, 2.0, 8.0]])
    result = fill(7, 1, 0.0, 1, 5, jnp.ones_like(field), jnp.ones_like(field), field)
    np.testing.assert_allclose(result, [[9.0, 0.0, 1.75, 1.75, 1.75, 1.75, 8.0]], rtol=1e-15)


def test_horizontal_tke_fill_inactive_division_and_active_derivatives():
    from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_wp2_from_horz_tke

    # No donor energy, sufficient donors, and insufficient donors, respectively.
    x = jnp.array([[.2, .1, .1], [.02, .05, .7], [-.4, .2, .15]])
    def fn(x):
        wp, up, vp = fill_holes_wp2_from_horz_tke(
            3, 1, .1, 0, 2, x[:, 0][None, :], x[:, 1][None, :], x[:, 2][None, :])
        return jnp.stack((wp[0], up[0], vp[0]), axis=1)
    weights = jnp.arange(1., 10.).reshape(3, 3)
    loss = lambda x: jnp.sum(weights * fn(x)**2)
    grad = jax.jit(jax.grad(loss))(x)
    assert np.isfinite(grad).all()
    # Stay away from the exact donor/receiver clipping boundaries for FD.
    direction = jnp.array([[.03, 0., 0.], [0., .01, -.03], [.02, -.01, .02]])
    _, tangent = jax.jvp(loss, (x,), (direction,))
    np.testing.assert_allclose(tangent, jnp.vdot(grad, direction), atol=1e-12)
    h = 1e-5
    np.testing.assert_allclose(tangent, (loss(x+h*direction)-loss(x-h*direction))/(2*h), rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(jnp.sum(fn(x), axis=1), jnp.sum(x, axis=1), atol=1e-15)
