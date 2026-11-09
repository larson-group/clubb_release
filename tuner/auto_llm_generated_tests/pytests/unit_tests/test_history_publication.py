"""Check count/time publication boundaries without waiting for a real model."""

import numpy as np
import pytest

from tuner import sample_history
from tuner.algorithm_registry import get_algorithm, normalize_options, parse_strategy, prepare_strategy
from tuner.sample_history import SampleHistoryWriter, sample_history_paths


def _writer(path):
    """Create one-field history using the real atomic NPZ writer."""
    return SampleHistoryWriter(job_dir=path, param_names=['A'], case_names=['case'],
                               field_names=['field'], metric_names=['rmse'],
                               case_configs=[], case_window_counts={'case': 1})


def _entry(index):
    """Provide a distinct completed candidate for duplicate/loss checks."""
    return {'sample_id': index, 'batch_id': index, 'all_params': {'A': float(index)},
            'case_loss': {'case': float(index)},
            'field_metrics': {'case': {'field': {'subwindows': [{'rmse': float(index)}]}}}}


def _ids(path):
    """Read all published rows in file order without retaining open NPZ handles."""
    result = []
    for file in sample_history_paths(path):
        with np.load(file) as data:
            result.extend(data['sample_id'].tolist())
    return result


def test_first_set_is_published_immediately_and_starts_the_timer(tmp_path, monkeypatch):
    """An arbitrary spin-up duration must not cause an empty file or early timer."""
    now = [0.]
    monkeypatch.setattr(sample_history.time, 'monotonic', lambda: now[0])
    writer = _writer(tmp_path)
    now[0] = 1000.
    writer.flush_if_due()
    writer.append([])
    assert not sample_history_paths(tmp_path)
    assert writer.last_flush_monotonic is None
    writer.append([_entry(0), _entry(1)])
    assert _ids(tmp_path) == [0, 1]
    assert writer.last_flush_monotonic == 1000.
    now[0] = 1019.999
    writer.append([_entry(2)])
    assert _ids(tmp_path) == [0, 1]
    now[0] = 1020.
    # No new append: the scheduler's timer tick must expose the pending row.
    writer.flush_if_due()
    assert _ids(tmp_path) == [0, 1, 2]
    now[0] = 2000.
    writer.flush_if_due()
    assert len(sample_history_paths(tmp_path)) == 2


def test_five_hundred_samples_publish_before_time_limit_and_reset_clock(tmp_path, monkeypatch):
    """Count-based publication includes exactly 500 rows after the first set."""
    now = [0.]
    monkeypatch.setattr(sample_history.time, 'monotonic', lambda: now[0])
    writer = _writer(tmp_path)
    writer.append([_entry(0)])
    now[0] = 1.
    writer.append([_entry(i) for i in range(1, 500)])
    assert _ids(tmp_path) == [0]
    now[0] = 2.
    writer.append([_entry(500)])
    assert _ids(tmp_path) == list(range(501))
    writer.append([_entry(501)])
    now[0] = 21.999
    writer.flush_if_due()
    assert _ids(tmp_path) == list(range(501))
    now[0] = 22.
    writer.flush_if_due()
    assert _ids(tmp_path) == list(range(502))
    writer.close()
    assert len(sample_history_paths(tmp_path)) == 3


def test_closing_and_resuming_preserve_each_pending_sample_once(tmp_path, monkeypatch):
    """Final flush and a restarted writer must neither overwrite nor duplicate rows."""
    monkeypatch.setattr(sample_history.time, 'monotonic', lambda: 0.)
    writer = _writer(tmp_path)
    writer.append([_entry(0)])
    writer.append([_entry(1)])
    writer.close()
    writer.close()
    resumed = _writer(tmp_path)
    resumed.append([_entry(2), _entry(3)])
    resumed.close()
    assert _ids(tmp_path) == [0, 1, 2, 3]
    assert len(sample_history_paths(tmp_path)) == 3


def test_catalog_defaults_and_gradient_chain_option_follow_execution():
    """Default budgets agree across methods and no public chain override remains."""
    assert normalize_options('random', {})['max_samples'] == 200
    assert normalize_options('adam', {})['max_updates'] == 30
    assert normalize_options('adam_gradient', {})['max_updates'] == 30
    spec = get_algorithm('adam_gradient')
    assert spec.label == 'Adam (grad)'
    assert 'chain_count' not in {option.name for option in spec.options}
    with pytest.raises(ValueError, match='derived scheduling'):
        parse_strategy('adam_gradient:chain_count=1')
    request = {'backend': 'jax', 'strategy': {'name': 'adam_gradient', 'options': {'chain_count': 1}},
               'cases': ['case'], 'max_workers': 3, 'batch_size': 2, 'parameter_ranges': []}
    prepare_strategy(request)
    assert request['strategy']['options']['chain_count'] == 6
    assert request['total_samples'] == 186
