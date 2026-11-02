"""Check the standard-deviation-scaled profile RMSE against independent values."""
import numpy as np

from utilities.loss_metrics import calculate_rmse


def test_rmse_uses_mean_and_benchmark_stddev():
    """Duplicating levels leaves the mean error unchanged."""
    model = np.array([1., 2.])
    benchmark = np.array([4., 6.])
    expected = np.sqrt(25. / 2.)
    np.testing.assert_allclose(calculate_rmse(model, benchmark), expected)
    np.testing.assert_allclose(
        calculate_rmse(np.tile(model, 2), np.tile(benchmark, 2)), expected,
    )


def test_rmse_uses_unit_fallback_for_flat_benchmark():
    """A flat reference uses the established unit normalization."""
    np.testing.assert_allclose(calculate_rmse([1., 3.], [2., 2.]), 1.)


def test_contribution_weights_are_positive_and_not_rescaled():
    """Reject zero/nonfinite contributions and preserve positive magnitudes."""
    import pytest
    from utilities.loss_metrics import normalize_loss_weights
    assert normalize_loss_weights({'mean': 2., 'std': .3, 'corr': 1.7}) == {'mean': 2., 'std': .3, 'corr': 1.7}
    assert normalize_loss_weights() == {'mean': 1., 'std': 1., 'corr': 1.}
    for bad in (0., -1., float('nan'), float('inf'), None):
        with pytest.raises(ValueError, match='positive'):
            normalize_loss_weights({'corr': bad})
    with pytest.raises(ValueError, match='Unknown'):
        normalize_loss_weights({'variance': 1.})
