"""Exercise Tune's setting cards and algorithm-derived gradients without a server."""

from types import SimpleNamespace

from dash import no_update
import pytest

from dash_app.tune_tab import callbacks_runs, callbacks_settings
from dash_app.tune_tab.layout import build_top_controls
from dash_app.tune_tab.tab import build_initial_tune_state
from dash_app.tune_tab.callbacks_display import _load_sample_history, _sample_history_signature
from tuner import sample_history
from utilities.loss_metrics import LOSS_METRIC_NAMES


class _Callbacks:
    """Collect the functions registered by the real settings callback owner."""

    def __init__(self):
        """Keep callbacks by name for direct invocation with a synthetic click."""
        self.functions = {}

    def callback(self, *args, **kwargs):
        """Retain the function while accepting Dash's normal callback declaration."""
        def register(function):
            """Store one actual callback without starting a browser or server."""
            self.functions[function.__name__] = function
            return function
        return register


def _walk(component):
    """Walk the actual component tree, including list-valued children."""
    if isinstance(component, (list, tuple)):
        for child in component:
            yield from _walk(child)
    elif hasattr(component, 'to_plotly_json'):
        yield component
        yield from _walk(getattr(component, 'children', None))


@pytest.mark.parametrize('section', ['algorithms', 'gradient', 'loss', 'aggregation'])
def test_settings_help_opens_and_closes_shared_notecards(monkeypatch, section):
    """Each question-mark opens its dialog, and mounting a zero-click Close cannot dismiss it."""
    callbacks = _Callbacks()
    callbacks_settings.register_settings_callbacks(callbacks)
    show = callbacks.functions['show_settings_help']
    context = SimpleNamespace(triggered_id={'type': 'tune-settings-help-open', 'section': section},
                              triggered=[{'value': 1}])
    monkeypatch.setattr(callbacks_settings, 'callback_context', context)
    card = show([1], [])
    dialog = next(node for node in _walk(card) if getattr(node, 'role', None) == 'dialog')
    assert dialog.to_plotly_json()['props']['aria-modal'] == 'true'
    close_id = {'type': 'tune-settings-help-close', 'section': section}
    assert any(getattr(node, 'id', None) == close_id for node in _walk(card))
    context.triggered_id = close_id
    context.triggered = [{'value': 0}]
    assert show([1], [0]) is no_update
    context.triggered = [{'value': 1}]
    assert show([1], [1]) == ''


def test_settings_rows_have_associated_labels_and_no_optional_gradient_checkbox():
    """All settings help buttons are reachable and form controls retain proper labels."""
    data = build_initial_tune_state()
    nodes = list(_walk(build_top_controls(data)))
    assert not any(getattr(node, 'id', None) == 'tune-grad' for node in nodes)
    help_ids = [node.id for node in nodes if isinstance(getattr(node, 'id', None), dict)
                and node.id.get('type') == 'tune-settings-help-open']
    assert {item['section'] for item in help_ids} == {'algorithms', 'gradient', 'loss', 'aggregation'}
    for row in nodes:
        if str(getattr(row, 'className', '')).startswith('tune-setting-row'):
            label, control = row.children
            assert label.htmlFor == control.id
    mode = next(node for node in nodes if getattr(node, 'id', None) == 'tune-gradient-mode')
    assert mode.disabled is True
    data['strategy_mode'] = 'adam_gradient'
    mode = next(node for node in _walk(build_top_controls(data))
                if getattr(node, 'id', None) == 'tune-gradient-mode')
    assert mode.disabled is False


@pytest.mark.parametrize('algorithm,gradients', [('random', False), ('adam', False), ('adam_gradient', True)])
def test_native_requests_derive_gradients_from_algorithm(monkeypatch, algorithm, gradients):
    """Native Tune cannot request optional gradients for a loss-value algorithm."""
    monkeypatch.setattr(callbacks_runs, 'load_tunable_names', lambda config: ['C8'])
    request = callbacks_runs.build_request_payload(
        case_names=['bomex'], time_start_values=[0], time_end_values=[600], average_time_values=[300],
        altitude_min_values=[0], altitude_max_values=[3000], selected_fields=['wp2'],
        param_values=['C8'], min_values=[.2], max_values=[.8], batch_size=4, max_workers=2,
        strategy_mode=algorithm, loss_weights={'mean': 1., 'std': 1., 'corr': 1.},
        aggregation_mode='quantile_weighted', strategy_options={}, selected_config='default',
        scm_override='', implementation='jax', gradient_mode='reverse')
    assert request['grad'] is gradients
    if gradients:
        assert request['strategy']['options']['chain_count'] == 8


def test_landscape_observes_partial_history_after_timer_tick(tmp_path, monkeypatch):
    """Slow evaluations must not hide pending samples from Dash's existing history reader."""
    now = [1000.]
    monkeypatch.setattr(sample_history.time, 'monotonic', lambda: now[0])
    writer = sample_history.SampleHistoryWriter(
        job_dir=tmp_path, param_names=['C8'], case_names=['bomex'], field_names=['wp2'],
        metric_names=list(LOSS_METRIC_NAMES), case_configs=[], case_window_counts={'bomex': 1})

    def entry(index):
        """Provide labeled metrics and a distinct parameter value for one candidate."""
        return {'sample_id': index, 'batch_id': index, 'all_params': {'C8': .2 + index * .1},
                'case_loss': {'bomex': 1. / (index + 1)},
                'field_metrics': {'bomex': {'wp2': {'subwindows': [
                    {metric: 1. / (index + 1) for metric in LOSS_METRIC_NAMES}]}}}}

    assert _load_sample_history(tmp_path) is None
    writer.append([entry(0)])
    first_signature = _sample_history_signature(tmp_path)
    assert first_signature
    assert _load_sample_history(tmp_path)['sample_id'].tolist() == [0]
    now[0] = 1010.
    writer.append([entry(1)])
    assert _sample_history_signature(tmp_path) == first_signature
    now[0] = 1020.
    writer.flush_if_due()
    assert _sample_history_signature(tmp_path) != first_signature
    history = _load_sample_history(tmp_path)
    assert history['sample_id'].tolist() == [0, 1]
    assert history['all_params'].shape == (2, 1)
    assert history['loss_metrics'].shape == (2, 1, len(LOSS_METRIC_NAMES))


def test_algorithm_help_reads_current_guide_and_resolves_repository_links(tmp_path, monkeypatch):
    """Reopening help reads edited documentation and keeps its source links usable."""
    from dash import dcc
    from dash_app.tune_tab import layout

    guide = tmp_path / 'tuner' / 'algorithms' / 'README.md'
    guide.parent.mkdir(parents=True)
    guide.write_text('# Methods\n\n[Random](random.py) and [tuner](../README.md#quick-start).', encoding='utf-8')
    monkeypatch.setattr(layout, 'REPO_ROOT', tmp_path)
    monkeypatch.setattr(layout, 'source_provenance', lambda root: {'commit': 'abc123'})
    card = layout.settings_help_notecard('algorithms')
    markdown = next(node for node in _walk(card) if isinstance(node, dcc.Markdown))
    assert markdown.children == (
        '# Methods\n\n[Random](https://github.com/larson-group/clubb/blob/abc123/tuner/algorithms/random.py) '
        'and [tuner](https://github.com/larson-group/clubb/blob/abc123/tuner/README.md#quick-start).')
    assert markdown.link_target == '_blank'
    guide.write_text('# Updated methods\n\nNew shared explanation.', encoding='utf-8')
    card = layout.settings_help_notecard('algorithms')
    markdown = next(node for node in _walk(card) if isinstance(node, dcc.Markdown))
    assert markdown.children == guide.read_text(encoding='utf-8')
