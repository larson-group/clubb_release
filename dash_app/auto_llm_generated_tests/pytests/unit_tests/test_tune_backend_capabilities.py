"""Exercise backend incompatibility warnings and launch gating without model jobs."""

import json

from dash import Dash
import pytest

from dash_app.tune_tab.algorithm_controls import algorithm_backend_error, option_values
from dash_app.tune_tab.callbacks_display import register_display_callbacks
from tuner.algorithm_registry import ALGORITHMS, get_algorithm


@pytest.fixture
def render_state():
    """Invoke the registered UI callback with a complete valid draft and inspect outputs by owner."""
    app = Dash(__name__)
    register_display_callbacks(app)
    entry = next(item for item in app.callback_map.values()
                 if item['callback'].__wrapped__.__name__ == 'render_tuning_state')
    assert any(item['id'] == 'compile-run-implementation' for item in entry['inputs'])

    def render(algorithm='adam_gradient', backend='fortran', **changes):
        """Allow one setting to change while preserving valid cases, ranges and options."""
        arguments = dict(
            status={}, top_results=[], best_results=[], loss_runs={}, active_job={},
            workspace_selection={'mode': 'draft'}, strategy_mode=algorithm,
            loss_weights={}, aggregation_mode='quantile_weighted',
            _aggregation_weight_1=.1, _aggregation_weight_2=.4, _aggregation_weight_3=.4, _aggregation_weight_4=.1,
            selected_config='default', _tunable_configs=['default'],
            strategy_option_values=option_values(algorithm), tunable_names=['C8'],
            _batch_size=4, _max_workers=2, _scm_override='', selected_case_names=['bomex'],
            time_start_values=[0], time_end_values=[600], average_time_values=[300],
            altitude_min_values=[20], altitude_max_values=[2940],
            member_ids=[{'row': 0, 'member': 0}], member_values=['C8'], range_row_order=[0],
            min_values=[.2], max_values=[.8], run_implementation=backend, config_button_ids=[])
        arguments.update(changes)
        values = entry['callback'].__wrapped__(**arguments)
        assert len(values) == len(entry['output'])
        return {(output.component_id['type'] if isinstance(output.component_id, dict) else output.component_id,
                 output.component_property): value for output, value in zip(entry['output'], values)}

    return render


@pytest.mark.parametrize('algorithm', ['adam_gradient', 'adam', 'random'])
@pytest.mark.parametrize('backend', ['fortran', 'jax'])
def test_draft_start_and_warning_follow_selected_backend(render_state, algorithm, backend):
    """Changing only the backend must block gradient methods and preserve value-only methods."""
    result = render_state(algorithm, backend)
    incompatible = algorithm == 'adam_gradient' and backend == 'fortran'
    assert result['tune-start-button', 'disabled'] is incompatible
    index = list(ALGORITHMS).index(algorithm)
    title = result['tune-algorithm-button', 'title'][index]
    summary = result['tune-algorithm-summary', 'children']
    if incompatible:
        assert 'Select JAX' in title and 'Select JAX' in summary
        assert result['tune-algorithm-button', 'style'][index]['backgroundColor'] == '#991b1b'
    else:
        assert 'Select JAX' not in title and 'Select JAX' not in summary
    # Incompatible methods remain selectable so their explanation/options can be read.
    assert result['tune-algorithm-button', 'disabled'][index] is False


def test_new_catalog_gradient_method_gets_the_same_warning(monkeypatch, render_state):
    """A new JSON registration inherits backend gating without a method-specific UI branch."""
    monkeypatch.setitem(ALGORITHMS, 'another_grad', dict(ALGORITHMS['adam_gradient'], label='Another gradient method'))
    get_algorithm.cache_clear()
    try:
        assert 'Another gradient method' in algorithm_backend_error('another_grad', 'fortran')
        assert algorithm_backend_error('another_grad', 'jax') == ''
        result = render_state('another_grad', 'fortran')
        assert result['tune-start-button', 'disabled'] is True
        assert 'Select JAX' in result['tune-algorithm-summary', 'children']
    finally:
        get_algorithm.cache_clear()


@pytest.mark.parametrize('saved_backend,blocked', [('jax', False), ('fortran', True)])
def test_continue_uses_saved_backend_instead_of_shared_chooser(tmp_path, render_state, saved_backend, blocked):
    """Readonly continuation validates its immutable runtime even when another tab changes the chooser."""
    (tmp_path / 'request.json').write_text(json.dumps({'backend': saved_backend}))
    result = render_state(
        backend='fortran' if saved_backend == 'jax' else 'jax',
        status={'state': 'stopped', 'job_dir': str(tmp_path)},
        workspace_selection={'mode': 'readonly'}, selected_case_names=[])
    assert result['tune-start-button', 'disabled'] is blocked
    assert ('Select JAX' in result['tune-algorithm-summary', 'children']) is blocked
