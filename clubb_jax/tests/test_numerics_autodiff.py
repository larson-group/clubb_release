"""Clipped-root values and the explicit zero-slope boundary convention."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped
from clubb_jax.src.CLUBB_core import (
    pdf_utilities, pdf_closure_module, adg1_adg2_3d_luhar_pdf,
)


@pytest.mark.parametrize('dtype', [jnp.float32, jnp.float64])
def test_clipped_root_values_and_both_derivatives(dtype):
    x = jnp.array([-1., 0., 1e-12, .25, 1., 4.], dtype=dtype)
    expected = jnp.array([0., 0., 5e5, 1., .5, .25], dtype=dtype)
    value, tangent = jax.jit(lambda x: jax.jvp(sqrt_clipped, (x,), (jnp.ones_like(x),)))(x)
    gradient = jax.jit(jax.grad(lambda x: jnp.sum(sqrt_clipped(x))))(x)
    np.testing.assert_array_equal(value, jnp.sqrt(jnp.maximum(x, 0)))
    np.testing.assert_allclose(tangent, expected, rtol=1e-6)
    np.testing.assert_allclose(gradient, expected, rtol=1e-6)
    assert jnp.isnan(sqrt_clipped(jnp.array(jnp.nan, dtype=dtype)))
    # A disabled contribution must still have sensitivity to its coefficient.
    assert jax.grad(lambda c: c * sqrt_clipped(jnp.array(4., dtype=dtype)))(0.) == 2.


def test_pdf_helpers_share_boundary_convention():
    for module in (pdf_utilities, pdf_closure_module, adg1_adg2_3d_luhar_pdf):
        assert jax.jit(jax.grad(module.sqrt_clipped))(0.) == 0.
        assert jax.jit(jax.grad(module.sqrt_clipped))(-1.) == 0.
