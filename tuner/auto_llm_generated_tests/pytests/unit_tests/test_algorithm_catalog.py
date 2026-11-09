"""Check JSON discovery and current strategy checkpoint continuation.

These shared Python unit checks need no model runtime. The catalog extension
and typed Dash rendering are covered by the existing analytic-Adam and Dash
suites; this module checks metadata-only loading and current pickle continuation.
"""

from collections import deque
import json
from pathlib import Path
import pickle
import subprocess
import sys

import pytest

from tuner.algorithm_registry import ALGORITHMS, get_algorithm, _read_catalog, prepare_strategy
from tuner.tuning_scheduler import TuningScheduler
from tuner.tuning_strategy import build_tuning_strategy


def test_catalog_listing_requires_no_algorithm_or_scientific_imports():
    """List real algorithm options in a fresh interpreter without site packages."""
    root = Path(__file__).resolve().parents[4]
    script = (
        "import json,sys; from tuner.algorithm_registry import list_algorithms; "
        "print(json.dumps(list_algorithms())); "
        "assert 'numpy' not in sys.modules and 'jax' not in sys.modules; "
        "assert not any(name.startswith('tuner.algorithms.') for name in sys.modules)"
    )
    result = subprocess.run([sys.executable, '-S', '-c', script], cwd=root,
                            check=True, capture_output=True, text=True)
    metadata = json.loads(result.stdout)
    assert {entry['name'] for entry in metadata} == set(ALGORITHMS)
    assert any(entry['requires_gradients'] for entry in metadata)
    assert {entry['category'] for entry in metadata} == {'optimizer', 'sampler'}
    assert all(entry['options'] for entry in metadata)


@pytest.mark.parametrize('change,message', [
    ({'options': {'rate': {'type': 'unknown'}}}, 'Unsupported option type'),
    ({'options': {'rate': {'type': 'float', 'default': -1., 'minimum': 0.}}}, '>='),
    ({'requires_gradients': 'false'}, 'must be a boolean'),
    ({'cli_order': ['undeclared']}, 'undeclared'),
    ({'category': 'search'}, 'category must be optimizer or sampler'),
    ({'category': None}, 'category must be optimizer or sampler'),
])
def test_invalid_catalog_entries_fail_before_building_controls(monkeypatch, change, message):
    """Reject invalid categories, types, defaults and positional metadata at lookup."""
    entry = {'implementation': 'not_imported:Strategy', 'category': 'sampler',
             'options': {'rate': {'type': 'float', 'default': .1}}, **change}
    monkeypatch.setitem(ALGORITHMS, 'invalid_catalog_entry', entry)
    get_algorithm.cache_clear()
    try:
        with pytest.raises(ValueError, match=message):
            get_algorithm('invalid_catalog_entry')
    finally:
        get_algorithm.cache_clear()


def test_catalog_root_must_be_an_algorithm_mapping(tmp_path):
    """A JSON array must not silently become a different discovery contract."""
    path = tmp_path / 'algorithms.json'
    path.write_text('[]')
    with pytest.raises(ValueError, match='JSON object'):
        _read_catalog(path)


@pytest.mark.parametrize('name,options', [
    ('random', {'max_samples': 5}), ('simann', {'max_iters': 3}),
    ('adam', {'max_updates': 2}), ('adam_gradient', {'max_updates': 2}),
])
def test_scheduler_resumes_current_strategy_and_chain_state(tmp_path, name, options):
    """Restore current checkpoints through the scheduler and compare every proposal.

    Include pending evaluations and reversed feedback so adaptive chain state,
    sample IDs and RNG continuation must survive. Grid instances are excluded:
    the baseline itertools.product state is not pickleable on Python 3.14.
    """
    request = {'strategy': {'name': name, 'options': options}, 'backend': 'jax',
               'max_workers': 2, 'batch_size': 4, 'cases': ['case'],
               'parameter_ranges': [{'name': 'A', 'min': 0., 'max': 2.}]}
    prepare_strategy(request)
    strategy = build_tuning_strategy(request['strategy'], param_names=['A'],
                                     default_params_row=[1.5], parameter_ranges=request['parameter_ranges'], seed=7)
    pending = deque()
    strategy.fill(pending, 2)
    path = tmp_path / 'resume.pkl'
    path.write_bytes(pickle.dumps({'version': 1, 'strategy': strategy, 'pending_samples': list(pending)}))
    # Exercise only the real checkpoint reader; no workers or SCM runs are needed.
    scheduler = object.__new__(TuningScheduler)
    scheduler.checkpoint_path, scheduler.cases, scheduler.baselines = path, ['case'], {}
    scheduler._restore_resume_checkpoint()
    restored = scheduler.strategy
    assert type(restored) is type(strategy)
    resumed_pending = scheduler.pending_samples
    for _ in range(100):
        assert list(pending) == list(resumed_pending)
        completed = [{'sample_id': sample['sample_id'], 'total_loss': (sample['param_row'][0] - .4)**2,
                      'total_gradient': [2. * (sample['param_row'][0] - .4)], 'gradient_parameters': ['A']}
                     for sample in reversed(pending)]
        strategy.tell(completed)
        restored.tell(completed)
        pending.clear()
        resumed_pending.clear()
        assert strategy.is_exhausted() == restored.is_exhausted()
        if strategy.is_exhausted():
            break
        strategy.fill(pending, 32)
        restored.fill(resumed_pending, 32)
    else:
        pytest.fail('Bounded strategies did not exhaust after restoring their checkpoint')
