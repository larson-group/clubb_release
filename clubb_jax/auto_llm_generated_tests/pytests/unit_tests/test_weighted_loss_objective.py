"""Check shared objective gradients using orthogonal synthetic profile patterns.

The two patterns have known unit variance and zero cross covariance. Their
amplitudes and a constant offset give independent analytic loss/derivative
oracles without running CLUBB or importing its Fortran bindings.
"""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped
from utilities.loss_metrics import case_rmse, profile_comparison, root_mean_square


@pytest.mark.parametrize('weights', [(1., 1., 1.), (2., .25, 3.), (.5, 4., .2)])
@pytest.mark.parametrize('parameters', [(.7, .4, .3), (0., 0., .3)])
def test_weighted_profile_gradient_matches_orthogonal_pattern_oracle(weights, parameters):
    """Check smooth profiles and the explicit zero-norm convention at flat models."""
    truth = jnp.array([-1., 1., -1., 1.])
    independent = jnp.array([-1., -1., 1., 1.])

    def objective(parameters):
        """Construct a profile from two known orthogonal patterns and a mean offset."""
        alpha, beta, offset = parameters
        model = alpha * truth + beta * independent + offset
        return profile_comparison(model, truth, weights, array_module=jnp,
                                  sqrt=sqrt_clipped)['metrics']['rmse']

    alpha, beta, offset = parameters
    amplitude = np.hypot(alpha, beta)
    wm, ws, wc = weights
    expected_squared = wm * offset**2 + ws * (amplitude - 1.)**2 + 2 * wc * (amplitude - alpha)
    expected = np.sqrt(expected_squared)
    amplitude_gradient = np.array([alpha, beta]) / amplitude if amplitude else np.zeros(2)
    gradient = np.r_[2 * ws * (amplitude - 1.) * amplitude_gradient
                     + 2 * wc * (amplitude_gradient - [1., 0.]), 2 * wm * offset] / (2 * expected)
    value, actual = jax.jit(jax.value_and_grad(objective))(jnp.asarray(parameters))
    assert actual.shape == (3,)
    np.testing.assert_allclose(value, expected, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(actual, gradient, rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize('weights', [(1., 1., 1.), (2., .25, 3.), (.5, 4., .2)])
def test_flat_benchmark_retains_variability_error_and_its_gradient(weights):
    """Undefined exported shape diagnostics must not suppress actual model variance."""
    model = jnp.array([-1., 1., -1., 1.]) + .3
    truth = jnp.zeros(4)

    def objective(model):
        """Compare to a constant benchmark using the established scale-one convention."""
        return profile_comparison(model, truth, weights, array_module=jnp,
                                  sqrt=sqrt_clipped)['metrics']['rmse']

    wm, ws, _ = weights
    expected = np.sqrt(wm * .3**2 + ws)
    value, gradient = jax.value_and_grad(objective)(model)
    expected_gradient = (wm * .3 + ws * (np.asarray(model) - .3)) / (4 * expected)
    assert gradient.shape == model.shape
    np.testing.assert_allclose(value, expected, rtol=1e-12)
    np.testing.assert_allclose(gradient, expected_gradient, rtol=1e-12, atol=1e-13)
    metrics = profile_comparison(model, truth, weights, array_module=jnp, sqrt=sqrt_clipped)['metrics']
    assert np.isnan([metrics['correlation'], metrics['std_ratio']]).all()


def test_exact_match_has_finite_zero_loss_gradient():
    """Use the existing zero derivative convention at an exact RMSE match."""
    truth = jnp.array([-1., 1., 2., 0.])
    objective = lambda model: profile_comparison(model, truth, array_module=jnp,
                                                 sqrt=sqrt_clipped)['metrics']['rmse']
    value, gradient = jax.value_and_grad(objective)(truth)
    np.testing.assert_array_equal(value, 0.)
    np.testing.assert_array_equal(gradient, jnp.zeros_like(truth))


def test_tiny_pattern_error_remains_visible_when_reported_correlation_rounds_to_one():
    """Use profile differences to retain a tiny orthogonal shape perturbation."""
    truth = jnp.array([-1., 1., -1., 1.])
    pattern = jnp.array([-1., -1., 1., 1.])

    def objective(amplitude):
        """Measure only the shape contribution of an orthogonal perturbation."""
        return profile_comparison(truth + amplitude * pattern, truth, (1., .1, 1.),
                                  array_module=jnp, sqrt=sqrt_clipped)['metrics']['rmse']

    amplitude = 1e-10
    represented_amplitude = np.mean((np.asarray(truth + amplitude * pattern) - np.asarray(truth))
                                    * np.asarray(pattern))
    value, gradient = jax.value_and_grad(objective)(amplitude)
    np.testing.assert_allclose(value, represented_amplitude, rtol=1e-8, atol=1e-20)
    np.testing.assert_allclose(gradient, 1., rtol=1e-8)


def test_case_and_case_collection_rms_gradients_match_weighted_norm():
    """Reduce windows/fields and then cases while preserving independent candidates."""
    values = jnp.array([[[1., 2.], [3., 4.]], [[2., 3.], [4., 5.]]])
    field_weights = [1., 3.]
    weights = np.broadcast_to(np.asarray(field_weights)[None, :, None], values.shape)
    expected_case = np.sqrt(np.sum(weights * np.asarray(values)**2, axis=(0, 1)) / 8)
    objective = lambda errors: jnp.sum(case_rmse(errors, field_weights,
                                                 array_module=jnp, sqrt=sqrt_clipped))
    actual_case = case_rmse(values, field_weights, array_module=jnp, sqrt=sqrt_clipped)
    gradient = jax.grad(objective)(values)
    assert actual_case.shape == (2,) and gradient.shape == values.shape
    np.testing.assert_allclose(actual_case, expected_case, rtol=1e-12)
    np.testing.assert_allclose(gradient, weights * np.asarray(values) / (8 * expected_case), rtol=1e-12)
    np.testing.assert_allclose(root_mean_square(jnp.stack([actual_case, actual_case]), axis=0,
                                               array_module=jnp, sqrt=sqrt_clipped), expected_case, rtol=1e-12)


def test_rms_keeps_large_finite_penalties_finite_and_excludes_zero_weight_items():
    """Scale before squaring and avoid contamination by explicitly excluded values."""
    values = jnp.array([1e30, 1e30], dtype=jnp.float32)
    assert bool(jnp.isfinite(root_mean_square(values, array_module=jnp, sqrt=sqrt_clipped)))
    np.testing.assert_allclose(root_mean_square(jnp.array([2., jnp.nan]), jnp.array([1., 0.]),
                                                array_module=jnp, sqrt=sqrt_clipped), 2.)
