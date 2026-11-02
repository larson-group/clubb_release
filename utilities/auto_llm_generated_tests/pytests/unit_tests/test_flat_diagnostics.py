"""Check undefined shape diagnostics without weakening loss validation or JSON."""
import json

import numpy as np
import pytest

from utilities.loss_metrics import calculate_profile_loss_metrics, loss_metric_values_valid, loss_result_for_json


@pytest.mark.parametrize('model,truth,expected_std', [([2., 2., 2.], [1., 2., 3.], 0.), ([1., 2., 3.], [2., 2., 2.], np.nan), ([2., 2., 2.], [2., 2., 2.], np.nan)])
def test_flat_profiles_export_undefined_without_changing_rmse(model, truth, expected_std):
    """A flat profile leaves a finite direct RMSE and undefined shape information."""
    result = calculate_profile_loss_metrics(model, truth)
    scale = np.std(truth) or 1.
    np.testing.assert_allclose(result['rmse'], np.sqrt(np.mean((np.array(model)-truth)**2))/scale)
    assert np.isnan(result['correlation'])
    np.testing.assert_equal(result['std_ratio'], expected_std)
    assert all(loss_metric_values_valid(name, value) for name, value in result.items())


def test_nonfinite_diagnostics_are_null_in_strict_json_and_nan_in_numeric_data():
    """Sanitizing scientific JSON leaves the caller's numeric arrays intact."""
    numeric = {'metrics': {'correlation': np.array([np.nan, 1.]), 'rmse': np.array([0., 2.])}, 'gradient': {'std_ratio': np.float64(np.nan)}}
    converted = json.loads(json.dumps(loss_result_for_json(numeric), allow_nan=False))
    assert converted['metrics']['correlation'] == [None, 1.]
    assert converted['metrics']['rmse'] == [0., 2.]
    assert converted['gradient']['std_ratio'] is None
    assert np.isnan(numeric['metrics']['correlation'][0])
    for name in ('correlation', 'std_ratio', 'rmse', 'bias_norm'):
        assert not loss_metric_values_valid(name, np.inf)
        assert bool(loss_metric_values_valid(name, np.nan)) == (name in ('correlation', 'std_ratio'))


@pytest.mark.parametrize('count', [3, 87, 100])
def test_rounded_constant_profiles_use_the_existing_unit_normalization(count):
    """A nonzero constant has zero variance even if its arithmetic mean rounds."""
    result = calculate_profile_loss_metrics(np.full(count, .2), np.full(count, .1))
    np.testing.assert_allclose(result['rmse'], .1, rtol=1e-14)
    np.testing.assert_allclose(result['bias_norm'], .1, rtol=1e-14)
    assert np.isnan(result['correlation']) and np.isnan(result['std_ratio'])


def test_small_real_variation_is_not_classified_as_a_flat_profile():
    """The flat check uses exact equality, with no arbitrary variance threshold."""
    profile = np.array([.1, np.nextafter(.1, 1.), np.nextafter(.1, 0.)])
    result = calculate_profile_loss_metrics(profile, profile)
    np.testing.assert_allclose([result['rmse'], result['correlation'], result['std_ratio']], [0., 1., 1.], atol=1e-15)
