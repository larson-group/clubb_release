"""Registry-driven Tune controls, hydration and broker requests without model runs."""

from dash import Dash
import pytest

from tuner.algorithm_registry import ALGORITHMS, _read_catalog, get_algorithm
from dash_app.services import TuneRequest
from dash_app.tune_tab.algorithm_controls import (
    algorithm_summary, build_algorithm_controls, option_ids, option_values, options_from_values,
)
from dash_app.tune_tab.callbacks_display import register_display_callbacks
from dash_app.tune_tab.callbacks_runs import agent_request_to_tune_controls
from dash_app.shared import actions


def _components(value):
    """Walk the actual layout tree to inspect registered input types and values."""
    if isinstance(value, (tuple, list)):
        for child in value:
            yield from _components(child)
    elif hasattr(value, 'to_plotly_json'):
        yield value
        yield from _components(getattr(value, 'children', None))


def test_a_new_registration_renders_typed_controls_and_round_trips(tmp_path, monkeypatch):
    """An extra list entry supplies numeric, string and boolean editors automatically."""
    import json
    catalog = tmp_path / "algorithms.json"
    catalog.write_text(json.dumps({"synthetic_ui": {
        "implementation": "unused_dash_algorithm:Strategy",
        "category": "sampler",
        "options": {"count": {"type": "int", "default": 3},
                    "step": {"type": "float", "default": .1},
                    "text": {"type": "string", "default": "hello"},
                    "enabled": {"type": "bool", "default": True}}}}))
    # Rendering reads metadata without importing the implementation.
    monkeypatch.setitem(ALGORITHMS, 'synthetic_ui', _read_catalog(catalog)['synthetic_ui'])
    get_algorithm.cache_clear()
    try:
        saved = {'count': 7, 'step': .03, 'text': 'saved', 'enabled': False}
        nodes = list(_components(build_algorithm_controls('synthetic_ui', saved)))
        inputs = {component.id['name']: component for component in nodes
                  if isinstance(getattr(component, 'id', None), dict)
                  and component.id.get('type') == 'tune-algorithm-option'
                  and component.id['algorithm'] == 'synthetic_ui'}
        assert inputs['count'].type == 'number' and inputs['count'].step == 1
        assert inputs['step'].type == 'number' and inputs['step'].step == 'any'
        assert inputs['text'].type == 'text'
        assert inputs['enabled'].value is False
        labels = {component.htmlFor for component in nodes if hasattr(component, 'htmlFor')}
        assert {json.dumps(control.id, sort_keys=True, separators=(',', ':'))
                for control in inputs.values()} <= labels
        assert options_from_values('synthetic_ui', option_values('synthetic_ui', saved)) == saved
        assert options_from_values('synthetic_ui', option_values('synthetic_ui', {'text': ''}))['text'] == ''
    finally:
        get_algorithm.cache_clear()


def test_saved_gradient_algorithm_hydrates_all_its_options():
    """Reloading a real canonical request preserves AD requirements and tuning values."""
    request = {'backend': 'jax', 'cases': ['bomex'], 'selected_fields': ['wp2'],
               'strategy': {'name': 'adam_gradient', 'options': {'max_updates': 4, 'beta1': .8}},
               'parameter_ranges': [{'name': 'C8', 'min': .2, 'max': .8}]}
    controls = agent_request_to_tune_controls(request, {'bomex': {'clubb_fields': ['wp2']}})
    assert controls['strategy'] == 'adam_gradient'
    options = options_from_values('adam_gradient', controls['strategy_options'])
    assert options['max_updates'] == 4 and options['beta1'] == .8
    assert controls['runtime']['implementation'] == 'jax'
    typed = TuneRequest(request_id='gradient-ui-test', backend='jax', cases=['bomex'],
                        strategy='adam_gradient', strategy_options={'max_updates': 4},
                        parameter_ranges=[{'name': 'C8', 'min': .2, 'max': .8}])
    assert typed.grad is True and typed.strategy_options['max_updates'] == 4
    with pytest.raises(ValueError, match='jax'):
        TuneRequest(request_id='gradient-ui-test', cases=['bomex'], strategy='adam_gradient',
                    parameter_ranges=[{'name': 'C8', 'min': .2, 'max': .8}])


def test_native_request_forwards_every_declared_option_without_method_adapter(monkeypatch):
    """The broker entry point preserves generic options including less common Adam settings."""
    captured = {}
    monkeypatch.setattr(actions, 'launch_tuning', lambda *args, **kw: captured.update(kw) or {})
    options = {'max_updates': 4, 'chain_count': 3, 'beta1': .8, 'epsilon': 1.e-6}
    actions.launch_tuning_request({'backend': 'jax', 'strategy': {'name': 'adam_gradient', 'options': options}})
    assert captured['strategy_options'] == options


def test_gradient_chain_count_has_no_control_and_previews_capacity_budget():
    """Chains are derived internally; saved explicit counts do not create a control."""
    values = option_values('adam_gradient', {'max_updates': 2})
    assert 'chain_count' not in options_from_values('adam_gradient', values)
    ranges = [{'name': 'C8', 'min': .2, 'max': .8}]
    summary = algorithm_summary('adam_gradient', values, 4, 3, ['bomex', 'atex'], ranges, backend='jax')
    assert summary == '36 candidate evaluations · requires JAX gradients'
    explicit = option_values('adam_gradient', {'max_updates': 2, 'chain_count': 1})
    assert 'chain_count' not in options_from_values('adam_gradient', explicit)
    assert algorithm_summary('adam_gradient', explicit, 4, 3, ['bomex', 'atex'], ranges, backend='jax') == (
        '36 candidate evaluations · requires JAX gradients')


def test_blank_resolve_spacing_previews_the_range_derived_default():
    """A blank grid control derives its budget from the selected physical ranges."""
    values = option_values('resolve')
    assert options_from_values('resolve', values) == {}
    ranges = [{'name': 'A', 'min': 0., 'max': 2.}, {'name': 'B', 'min': 10., 'max': 16.}]
    assert algorithm_summary('resolve', values, 1, 1, ['bomex'], ranges, backend='fortran') == (
        '21 candidate evaluations')


def test_registered_controls_lock_with_the_rest_of_a_saved_workspace():
    """A readonly revision locks every generated input and shows its selected panel."""
    app = Dash(__name__)
    register_display_callbacks(app)
    entry = next(item for item in app.callback_map.values()
                 if item['callback'].__wrapped__.__name__ == 'render_tuning_state')
    callback = entry['callback'].__wrapped__
    result = callback(
        status={}, top_results=[], best_results=[], loss_runs={}, active_job={},
        workspace_selection={'mode': 'readonly'}, strategy_mode='adam_gradient',
        loss_weights={}, aggregation_mode='quantile_weighted',
        _aggregation_weight_1=.1, _aggregation_weight_2=.4, _aggregation_weight_3=.4, _aggregation_weight_4=.1,
        selected_config='default', _tunable_configs=['default'],
        strategy_option_values=option_values('adam_gradient', {'max_updates': 2}),
        tunable_names=['C8'], _batch_size=1, _max_workers=1, _scm_override='',
        selected_case_names=['bomex'], time_start_values=[0], time_end_values=[600],
        average_time_values=[300], altitude_min_values=[20], altitude_max_values=[2940],
        member_ids=[{'row': 0, 'member': 0}], member_values=['C8'], range_row_order=[0],
        min_values=[.2], max_values=[.8], run_implementation='jax', config_button_ids=[])
    assert len(result) == len(entry['output'])
    by_id = {output.component_id['type']: value for output, value in zip(entry['output'], result)
             if isinstance(output.component_id, dict) and output.component_property == 'disabled'}
    assert by_id['tune-algorithm-option'] == [True] * len(option_ids())
    assert by_id['tune-algorithm-button'] == [True] * len(ALGORITHMS)


def test_generic_options_preserve_the_existing_agent_random_budget_limit():
    """Canonical agent options cannot bypass the public API's bounded random budget."""
    assert actions._normalize_tune_strategy('random', [], strategy_options={})['options']['max_samples'] == 200
    with pytest.raises(ValueError):
        actions._normalize_tune_strategy('random', [], strategy_options={'max_samples': 1001})
    assert actions._normalize_tune_strategy('random', [], strategy_options={'max_samples': 20})['options']['max_samples'] == 20
