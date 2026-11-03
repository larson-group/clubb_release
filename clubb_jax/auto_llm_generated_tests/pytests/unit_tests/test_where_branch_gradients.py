#!/usr/bin/env python3
"""Check that unselected jnp.where branches do not turn gradients into NaN.

jnp.where evaluates both of its branches, and autodiff multiplies the zero weight of
the unselected branch by that branch's derivative. If the branch evaluates sqrt(0) or
1/0, the derivative is infinite and the gradient becomes NaN even though the returned
value is correct. These tests differentiate the routines at those inputs and also
check that the values on the selected branch still follow the source formula.
"""

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import _reciprocal_condition_number_estimate
from clubb_jax.src.CLUBB_core.pdf_utilities import corr_LL2NN, corr_NL2NN, corr_NN2LL


def test_corr_NL2NN_gradient_is_finite_when_sigma_is_zero():
    """A zero y variance selects corr_x_y, whose gradient in that variance must be finite."""
    gradient = jax.grad(lambda y_var: corr_NL2NN(0.5, 0.0, y_var))(0.0)
    assert np.isfinite(gradient)


def test_corr_LL2NN_gradient_is_finite_when_sigma_is_zero():
    """Zero variances select corr_x_y, whose gradient in those variances must be finite."""
    gradient = jax.grad(lambda x_var: corr_LL2NN(0.5, 0.0, 0.0, x_var, 0.0))(0.0)
    assert np.isfinite(gradient)


def test_corr_NN2LL_gradient_is_finite_when_sigma_is_zero():
    """Zero variances select corr_x_y_n, whose gradient in those variances must be finite."""
    gradient = jax.grad(lambda x_var: corr_NN2LL(0.5, 0.0, 0.0, x_var, 0.0))(0.0)
    assert np.isfinite(gradient)


def test_corr_routines_keep_source_values_when_sigma_is_positive():
    """Positive variances use the source lognormal correlation formulas."""
    corr, sigma_x, sigma_y, x_var, y_var = 0.5, 0.3, 0.4, 0.2, 0.1
    assert np.isclose(corr_NL2NN(corr, sigma_y, y_var), corr * np.sqrt(y_var) / sigma_y, rtol=1.0e-14)
    assert np.isclose(
        corr_LL2NN(corr, sigma_x, sigma_y, x_var, y_var),
        np.log(1.0 + corr * np.sqrt(x_var * y_var)) / (sigma_x * sigma_y),
        rtol=1.0e-14,
    )
    assert np.isclose(
        corr_NN2LL(corr, sigma_x, sigma_y, x_var, y_var),
        (np.exp(sigma_x * sigma_y * corr) - 1.0) / np.sqrt(x_var * y_var),
        rtol=1.0e-14,
    )


def test_condition_estimate_gradient_is_finite_for_a_zero_matrix():
    """A zero 1x1 matrix has no defined condition number, so the result is zero with a finite gradient."""
    def estimate(scale):
        return _reciprocal_condition_number_estimate(scale * jnp.ones((1, 1)))

    assert estimate(0.0) == 0.0
    assert np.isfinite(jax.grad(estimate)(0.0))
