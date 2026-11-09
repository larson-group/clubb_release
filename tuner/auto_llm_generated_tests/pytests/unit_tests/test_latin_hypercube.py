"""Check stratification, linked ranges and current checkpoint continuation of LHS."""

from collections import deque
import pickle

import numpy as np
import pytest

from tuner.algorithm_registry import parse_strategy, prepare_strategy
from tuner.tuning_strategy import build_tuning_strategy


def _strategy(seed=17):
    """Build the registered sampler with one linked coordinate and one fixed value."""
    request = {'strategy': parse_strategy('latin_hypercube:17'), 'backend': 'fortran',
               'batch_size': 4, 'max_workers': 2, 'cases': ['synthetic'],
               'parameter_ranges': [{'name': 'A', 'targets': ['A', 'B'], 'min': 2., 'max': 6.},
                                    {'name': 'C', 'min': 8., 'max': 8.}]}
    prepare_strategy(request)
    assert request['total_samples'] == 17 and not request.get('grad', False)
    return build_tuning_strategy(request['strategy'], param_names=['A', 'B', 'C', 'untuned'],
                                 default_params_row=[3., 3., 8., 99.],
                                 parameter_ranges=request['parameter_ranges'], seed=seed)


@pytest.mark.parametrize('capacity', [1, 4, 20])
def test_lhs_design_survives_batching_feedback_and_resume(capacity):
    """Every stratum appears once even when feedback and checkpoint timing change."""
    strategy = _strategy()
    samples = []
    while not strategy.is_exhausted():
        queue = deque()
        strategy.fill(queue, capacity)
        samples.extend(queue)
        strategy.tell([dict(sample, total_loss=42.) for sample in reversed(queue)])
        strategy = pickle.loads(pickle.dumps(strategy))
    assert len(samples) == strategy.estimated_sample_count() == 17
    values = np.array([sample['param_row'] for sample in samples])
    assert values.shape == (17, 4)
    np.testing.assert_array_equal(values[:, 0], values[:, 1])
    np.testing.assert_array_equal(values[:, 2:], np.tile([8., 99.], (17, 1)))
    np.testing.assert_array_equal(np.sort(np.floor((values[:, 0]-2.)/4.*17).astype(int)), np.arange(17))
    reference, changed = deque(), deque()
    _strategy().fill(reference, 17)
    _strategy(seed=18).fill(changed, 17)
    assert samples == list(reference)
    assert samples != list(changed)
    assert [sample['sample_id'] for sample in samples] == list(range(17))
