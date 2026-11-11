"""Check differentiation of the existing time-averaged stats loss."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.clubb_loss_driver import (
    calculate_field_loss,
    invalid_loss_penalty,
    loss_field_type,
    loss_request_type,
)


def test_sampled_profile_loss_has_gradient_through_time_average():
    """Compare the two-sample stats loss and its gradient with analytic values."""
    initial = JaxStats.empty(l_sample=True, names=('wp2',), ncol=1, max_nlev=2)
    field = loss_field_type(
        clubb_var_name='wp2', benchmark_var_name='wp2', stats_var_id=initial.name_to_slot['wp2'],
        k_min=0, k_max=1, truth_profile=jnp.array([[1., 2.]]),
    )
    request = loss_request_type(fields=[field])

    def loss(c8):
        stats = initial.update('wp2', c8 * jnp.array([[1., 2.]]))
        stats = stats.begin_timestep(l_sample=True, l_reset_accumulators=False)
        stats = stats.update('wp2', c8 * jnp.array([[3., 4.]]))
        return calculate_field_loss(request, stats, 0)[0][0, 0]

    value, gradient = jax.jit(jax.value_and_grad(loss))(2.)
    # The benchmark [1, 2] has population standard deviation 0.5.
    np.testing.assert_allclose(value, 2 * np.sqrt(12.5), rtol=1e-13)
    np.testing.assert_allclose(gradient, 18. / np.sqrt(12.5), rtol=1e-13)


def test_exact_profile_match_has_finite_zero_gradient():
    """Use zero as the gradient convention at the nondifferentiable exact match."""
    initial = JaxStats.empty(l_sample=True, names=('wp2',), ncol=1, max_nlev=2)
    field = loss_field_type(
        clubb_var_name='wp2', benchmark_var_name='wp2', stats_var_id=initial.name_to_slot['wp2'],
        k_min=0, k_max=1, truth_profile=jnp.array([[1., 2.]]),
    )
    request = loss_request_type(fields=[field])

    def loss(scale):
        stats = initial.update('wp2', scale * jnp.array([[1., 2.]]))
        return calculate_field_loss(request, stats, 0)[0][0, 0]

    value, gradient = jax.jit(jax.value_and_grad(loss))(1.)
    np.testing.assert_array_equal(value, 0.)
    np.testing.assert_array_equal(gradient, 0.)


def test_failed_column_does_not_change_healthy_rmse_gradient():
    """A nonfinite neighbor is penalized while the healthy column keeps its gradient."""
    initial = JaxStats.empty(l_sample=True, names=('wp2',), ncol=2, max_nlev=2)
    field = loss_field_type(
        clubb_var_name='wp2', benchmark_var_name='wp2', stats_var_id=initial.name_to_slot['wp2'],
        k_min=0, k_max=1, truth_profile=jnp.array([[1., 2.]]),
    )
    request = loss_request_type(fields=[field])

    def scores(profiles):
        return calculate_field_loss(request, initial.update('wp2', profiles), 0)[0][0]

    profiles = jnp.array([[1.2, 2.4], [jnp.nan, 2.]])
    values = scores(profiles)
    gradient = jax.grad(lambda x: jnp.sum(scores(x)))(profiles)
    np.testing.assert_allclose(values[0], 2 * np.sqrt(.1), rtol=1e-13)
    np.testing.assert_array_equal(values[1], invalid_loss_penalty)
    np.testing.assert_allclose(gradient[0], [2. / np.sqrt(10.), 4. / np.sqrt(10.)], rtol=1e-13)


@pytest.mark.parametrize('target', ['loss', 'metrics'])
@pytest.mark.parametrize('selected', [['C8'], ['C7', 'C8'], None])
@pytest.mark.parametrize('mode', ['auto', 'forward', 'reverse'])
@pytest.mark.parametrize('field_count', [1, 2])
def test_loss_driver_batches_parameter_directions_and_retains_primal_stats(monkeypatch, selected, mode, field_count, target):
    """Check analytic derivatives, parameter order and unbatched stats for two columns.

    Substitute a polynomial scoring kernel for model advancement so the public
    driver batching contract has an independent, inexpensive numerical oracle.
    """
    from types import SimpleNamespace
    from clubb_jax.src import clubb_loss_driver as driver
    from utilities.loss_metrics import LOSS_METRIC_NAMES, case_rmse

    params = jnp.zeros((2, len(driver.PARAM_NAMES)))
    params = params.at[:, driver.PNAME_IDX['C8']].set(jnp.array([.4, .9]))
    params = params.at[:, driver.PNAME_IDX['C7']].set(jnp.array([.2, .7]))
    request = SimpleNamespace(state={'clubb_params': params}, fields=['wp2', 'thlm'][:field_count])
    monkeypatch.setattr(driver, 'active_request', request)
    calls = []
    weights = jnp.array([[1., 2.], [3., 4.]])[:, :field_count]

    def score(candidates, *, checkpoint=False):
        """Score two independent polynomial columns and retain their primal profiles."""
        calls.append(checkpoint)
        x, y = candidates[:, driver.PNAME_IDX['C8']], candidates[:, driver.PNAME_IDX['C7']]
        profile = jnp.stack([x, y], axis=1)
        stats = JaxStats.empty(l_sample=True, names=('wp2',), ncol=2, max_nlev=2).update('wp2', profile)
        request.state['_loss_stats_snapshots'] = [(2, (stats,))]
        metrics = tuple((i + 1) * weights[:, :, None] * (x * x + 3 * x * y + y)[None, None, :]
                     for i in range(len(LOSS_METRIC_NAMES)))
        return (*metrics, case_rmse(metrics[0], array_module=jnp))

    monkeypatch.setattr(driver, 'clubb_get_loss_for_params', score)
    progress = []
    result = driver.evaluate_loss(params, grad=True, gradient_parameters=selected,
                                  gradient_mode=mode, gradient_target=target,
                                  progress=lambda *counts: progress.append(counts))
    parameters = driver.PARAM_NAMES if selected is None else selected
    assert result['gradient_parameters'] == list(parameters)
    assert result['gradient_mode'] == (('reverse' if len(parameters) > (1 if target == 'loss' else field_count) else 'forward') if mode == 'auto' else mode)
    assert calls == [result['gradient_mode'] == 'reverse']
    assert progress == [(0, len(parameters)), (len(parameters), len(parameters))]
    x, y = np.array([.4, .9]), np.array([.2, .7])
    expected = np.stack([2 * x + 3 * y if name == 'C8' else 3 * x + 1 if name == 'C7'
                         else np.zeros(2) for name in parameters], axis=-1)
    for i, name in enumerate(LOSS_METRIC_NAMES):
        values = result['metrics'][name]
        assert values.shape == (2, field_count, 2)
        np.testing.assert_allclose(values, (i + 1) * np.asarray(weights)[:, :, None]
                                   * (x * x + 3 * x * y + y)[None, None, :], rtol=1e-13)
        if target == 'metrics':
            derivative = result['gradient'][name]
            assert derivative.shape == (2, field_count, 2, len(parameters))
            np.testing.assert_allclose(derivative, (i + 1) * np.asarray(weights)[:, :, None, None]
                                   * expected[None, None, :, :], rtol=1e-13)
    if target == 'loss':
        assert set(result['gradient']) == {'loss'}
        factor = np.sqrt(np.mean(np.asarray(weights) ** 2))
        np.testing.assert_allclose(result['gradient']['loss'], factor * expected, rtol=1e-13)
        assert result['gradient']['loss'].shape == (2, len(parameters))
    stats = request.state['_jax_stats']
    assert stats is request.state['_loss_stats_snapshots'][-1][1][-1]
    assert stats.buffers[0].shape == (1, 2, 2)
    np.testing.assert_allclose(stats.buffers[0][0], np.stack([x, y], axis=1), rtol=1e-13)


def test_unknown_gradient_mode_is_rejected_before_advancement(monkeypatch):
    """An invalid mode fails before invoking model advancement, even without gradients."""
    from clubb_jax.src import clubb_loss_driver as driver

    def advance(_):
        """Fail the test if validation accidentally starts a model evaluation."""
        raise AssertionError('Model advancement must not run')

    monkeypatch.setattr(driver, 'clubb_get_loss_for_params', advance)
    for grad in (False, True):
        with pytest.raises(ValueError, match='gradient_mode'):
            driver.evaluate_loss(jnp.zeros((1, len(driver.PARAM_NAMES))), grad=grad, gradient_mode='other')


@pytest.mark.parametrize('benchmark', [[0., 1., 2.], [1., 1., 1.]])
def test_flat_profile_taylor_pullbacks_stay_finite_for_all_metric_seeds(benchmark):
    """A flat model's inactive correlation root cannot contaminate other pullbacks.

    Compare the full Taylor Jacobian in both modes at the established zero-
    variance boundary convention, for variable and constant benchmarks.
    """
    from clubb_jax.src.clubb_loss_driver import calculate_profile_loss

    truth = jnp.array(benchmark)
    model = jnp.zeros(3)
    score = lambda profile: calculate_profile_loss(profile, truth, (1., 1., 1.))
    forward, reverse = jax.jacfwd(score)(model), jax.jacrev(score)(model)
    for actual, expected in zip(reverse, forward):
        assert np.isfinite(actual).all()
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
