"""Check range-derived Resolve grids through request preparation and construction.

The shared unit wrapper and the tuner's first Jenkins unit stage own these
checks. Unequal ranges, fixed coordinates and explicit overrides must give the
same budget and actual points with either prepared or empty algorithm options.
No model runtime is needed.
"""

from collections import deque
import itertools

import pytest

from tuner.algorithm_registry import parse_strategy, prepare_strategy
from tuner.tuning_strategy import build_tuning_strategy


@pytest.mark.parametrize('bounds,option,expected_spacing,axes', [
    ([(0., 2.), (10., 16.), (7., 7.)], None, 1.,
     [[0., 1., 2.], [10., 11., 12., 13., 14., 15., 16.], [7.]]),
    ([(0., 2.), (10., 16.)], 2., 2., [[0., 2.], [10., 12., 14., 16.]]),
    ([(7., 7.), (3., 3.)], None, 1., [[7.], [3.]]),
])
def test_resolve_default_spacing_matches_previews_and_physical_grid(bounds, option, expected_spacing, axes):
    """Check the smallest active width, explicit overrides and fixed-only grids."""
    names = [f'x{i}' for i in range(len(bounds))]
    ranges = [dict(name=name, min=low, max=high) for name, (low, high) in zip(names, bounds)]
    config = parse_strategy('resolve' if option is None else f'resolve:{option}')
    request = {'strategy': config, 'parameter_ranges': ranges}
    prepare_strategy(request)
    assert request['strategy']['options']['spacing'] == expected_spacing
    for strategy_config in (config, request['strategy']):
        strategy = build_tuning_strategy(strategy_config, param_names=names,
                                         default_params_row=[low for low, _ in bounds],
                                         parameter_ranges=ranges)
        assert strategy.grid_values == axes
        pending = deque()
        strategy.fill(pending, request['total_samples'] + 1)
        assert len(pending) == strategy.estimated_sample_count() == request['total_samples']
        assert strategy.is_exhausted()
        assert [sample['param_row'] for sample in pending] == [list(point) for point in itertools.product(*axes)]
