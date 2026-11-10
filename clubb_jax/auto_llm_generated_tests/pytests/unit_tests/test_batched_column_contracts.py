"""Check column independence, ordered reductions and interpolation exits from the Fortran refactor."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.advance_helper_module import vertical_avg, vertical_integral
from clubb_jax.src.CLUBB_core.interpolation import zlinterp


def test_vertical_reductions_keep_distinct_columns():
    """Different density weights produce independent column integrals and averages."""
    rho = jnp.array([[1.0, 2.0, 1.0], [2.0, 1.0, 3.0]])
    field = jnp.array([[2.0, 4.0, 6.0], [10.0, 20.0, 30.0]])
    dz = jnp.ones((2, 3))
    integral = vertical_integral(2, 3, rho, field, dz)
    average = vertical_avg(2, 3, rho, field, dz)
    assert integral.shape == average.shape == (2,)
    np.testing.assert_array_equal(integral, [16.0, 130.0])
    np.testing.assert_allclose(average, [4.0, 130.0 / 6.0], rtol=1e-15)
    with pytest.raises(ValueError, match="complete column/level"):
        vertical_integral(1, 3, rho[0], field[0], dz[0])


def test_vertical_integral_masks_inactive_nan_columns_and_gradients():
    """An inactive column contributes zero even when its input fields contain NaNs."""
    rho = jnp.array([[1.0, 2.0], [jnp.nan, jnp.nan]])
    field = jnp.array([[3.0, 4.0], [jnp.nan, jnp.nan]])
    dz = jnp.array([[1.0, 1.0], [jnp.nan, jnp.nan]])
    active = jnp.array([True, False])
    reduce = jax.jit(lambda value: vertical_integral(2, 2, rho, value, dz, active))
    np.testing.assert_array_equal(reduce(field), [11.0, 0.0])
    gradient = jax.grad(lambda value: jnp.sum(reduce(value)))(field)
    np.testing.assert_array_equal(gradient, [[1.0, 2.0], [0.0, 0.0]])


def test_vertical_integral_keeps_input_level_order():
    """Cancellation follows the Fortran left-to-right level accumulation."""
    field = jnp.array([[1.0e16, 1.0, -1.0e16], [1.0, -1.0, 3.0]], dtype=jnp.float64)
    weights = jnp.ones_like(field)
    np.testing.assert_array_equal(vertical_integral(2, 3, weights, field, weights), [0.0, 3.0])


def test_zlinterp_keeps_column_output_limits_and_first_above_source_exit():
    """A later in-range target stays zero after an above-source exit in that column."""
    source = jnp.array([[0.0, 1.0, 2.0], [10.0, 20.0, 30.0]])
    values = jnp.array([[0.0, 2.0, 4.0], [1.0, 3.0, 5.0]])
    target = jnp.array([[0.5, 3.0, 1.5, 0.0], [15.0, 25.0, 30.0, 10.0]])
    result = jax.jit(lambda: zlinterp(2, 4, 3, target, source, values, jnp.array([4, 2])))()
    assert result.shape == (2, 4)
    np.testing.assert_array_equal(result, [[1.0, 0.0, 0.0, 0.0], [2.0, 4.0, 0.0, 0.0]])


def test_zlinterp_inactive_column_does_not_read_nan_values():
    """An output-limited column has zero values and sensitivities despite NaN source values."""
    source = jnp.array([[0.0, 1.0], [0.0, 1.0]])
    values = jnp.array([[2.0, 4.0], [jnp.nan, jnp.nan]])
    target = jnp.array([[0.5], [0.5]])
    interpolate = lambda value: zlinterp(2, 1, 2, target, source, value, jnp.array([1, 0]))
    np.testing.assert_array_equal(interpolate(values), [[3.0], [0.0]])
    np.testing.assert_array_equal(jax.grad(lambda value: jnp.sum(interpolate(value)))(values), [[0.5, 0.5], [0.0, 0.0]])
