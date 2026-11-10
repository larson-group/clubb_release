"""Undefined Taylor diagnostics must not remove usable loss values or gradients."""
from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src import clubb_loss_driver as driver
from clubb_jax.src.CLUBB_core.stats_netcdf import JaxStats
from utilities.loss_metrics import case_rmse


@pytest.mark.parametrize('mode', ['forward', 'reverse'])
@pytest.mark.parametrize('target', ['loss', 'metrics'])
@pytest.mark.parametrize('flat_reference', [False, True])
def test_flat_diagnostics_preserve_case_and_metric_pullbacks(monkeypatch, mode, target, flat_reference):
    """Use analytic three-level profiles through the public two-column AD adapter.

    A flat reference has varying model amplitudes and positive loss gradients.
    A flat model against a variable reference retains a negative loss gradient,
    despite undefined correlation and a nondifferentiable std-ratio boundary.
    """
    amplitudes = [.25, .5] if flat_reference else [0., .5]
    truth = jnp.ones(3) if flat_reference else jnp.array([1., 2., 3.])
    offset = 1. if flat_reference else 2.
    pattern = jnp.array([-1., 0., 1.])
    params = jnp.zeros((2, len(driver.PARAM_NAMES))).at[:, driver.PNAME_IDX['C8']].set(jnp.array(amplitudes))
    request = SimpleNamespace(state={'clubb_params': params}, fields=['wp2'])
    monkeypatch.setattr(driver, 'active_request', request)

    def score(candidates, *, checkpoint=False):
        """Score actual profile math; substitute only the costly model trajectory."""
        model = offset + candidates[:, driver.PNAME_IDX['C8'], None] * pattern
        stats = JaxStats.empty(l_sample=True, names=('wp2',), ncol=2, max_nlev=3).update('wp2', model)
        request.state['_loss_stats_snapshots'] = [(2, (stats,))]
        metrics = tuple(values[None, None, :] for values in driver.calculate_profile_loss(model, truth, (1., 1., 1.)))
        return (*metrics, case_rmse(metrics[0], array_module=jnp, sqrt=driver.sqrt_clipped))

    monkeypatch.setattr(driver, 'clubb_get_loss_for_params', score)
    result = driver.evaluate_loss(params, grad=True, gradient_parameters=['C8'], gradient_mode=mode, gradient_target=target)
    expected = np.full(2, np.sqrt(2./3.) if flat_reference else -1.)
    np.testing.assert_allclose(result['loss'], np.asarray(amplitudes)*np.sqrt(2./3.) if flat_reference else 1.-np.asarray(amplitudes))
    if target == 'loss':
        assert result['gradient']['loss'].shape == (2, 1)
        np.testing.assert_allclose(result['gradient']['loss'][:, 0], expected, rtol=1e-13)
    else:
        np.testing.assert_allclose(result['gradient']['rmse'][0, 0, :, 0], expected, rtol=1e-13)
        np.testing.assert_array_equal(result['gradient']['bias_norm'], 0.)
        correlation, std = result['gradient']['correlation'], result['gradient']['std_ratio']
        if flat_reference:
            assert np.isnan(correlation).all() and np.isnan(std).all()
        else:
            assert np.isnan(correlation[0, 0, 0, 0]) and np.isnan(std[0, 0, 0, 0])
            np.testing.assert_allclose([correlation[0, 0, 1, 0], std[0, 0, 1, 0]], [0., 1.], atol=1e-13)
    assert np.isfinite(result['metrics']['rmse']).all()
    assert np.isnan(result['metrics']['correlation'][0, 0, 0])
    if flat_reference:
        assert np.isnan(result['metrics']['std_ratio']).all()
    else:
        assert result['metrics']['std_ratio'][0, 0, 0] == 0.


@pytest.mark.parametrize('mode', ['forward', 'reverse'])
def test_rounded_constant_model_keeps_the_mean_bias_gradient(mode):
    """Flat detection must not replace derivatives of the profile mean."""
    import jax
    truth = jnp.full(87, .1)
    score = lambda x: driver.calculate_profile_loss(jnp.full(87, x), truth, (1., 1., 1.))[0]
    gradient = (jax.jacfwd if mode == 'forward' else jax.jacrev)(score)(.2)
    np.testing.assert_allclose(score(.2), .1, rtol=1e-13)
    np.testing.assert_allclose(gradient, 1., rtol=1e-13)
    metrics = driver.calculate_profile_loss(jnp.full(87, .2), truth, (1., 1., 1.))
    assert np.isnan(metrics[1]) and np.isnan(metrics[2])
