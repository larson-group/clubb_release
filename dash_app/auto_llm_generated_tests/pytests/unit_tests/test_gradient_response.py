"""Check the gradient page's request, plotting and broker-task contracts."""
from concurrent.futures import Future
import json
from types import SimpleNamespace

from dash import Dash
import pytest

from dash_app.misc_tab.gradient_response import callbacks, multi_callbacks
from dash_app.plot_tab.tasks import PlotTasks
from dash_app.services.gradient_response import validate_request


def test_response_figure_anchors_tangent_at_baseline():
    """The plotted tangent uses the map baseline value and derivative."""
    figure = callbacks.response_figure('rmse', 'RMSE', {
        'values': [0., 1., 2.], 'baseline': 0., 'parameter': 'C8',
        'metrics': {'rmse': {'values': [2., 5., 10.], 'baseline': 5., 'gradient': 3.}},
    })
    assert list(figure.data[0].y) == [2., 5., 10.]
    assert list(figure.data[1].y) == [5., 8., 11.]
    assert list(figure.data[2].x) == [0.]
    assert list(figure.data[2].y) == [5.]


@pytest.mark.parametrize('minimum,maximum', [(1., 1.), (2., 1.), (float('nan'), 1.), (-1., 1.)])
def test_invalid_range_is_rejected_before_a_model_run(minimum, maximum):
    """Nonfinite, reversed and out-of-bounds ranges cannot launch a worker."""
    with pytest.raises(ValueError):
        validate_request('C8', minimum, maximum, 'wp2', 'bomex', 0, 600)


def test_page_submits_and_displays_only_the_current_request(monkeypatch):
    """Map inputs submit a request, and changes clear the previous result."""
    app = Dash(__name__)
    multi_callbacks.register_callbacks(app)
    update = next(entry['callback'].__wrapped__ for entry in app.callback_map.values()
                  if entry['callback'].__wrapped__.__name__ == 'update_map')
    context = SimpleNamespace(triggered_id='gradient-multi-go')
    monkeypatch.setattr(multi_callbacks, 'ctx', context)
    calls = []
    monkeypatch.setattr(multi_callbacks, 'perform_action', lambda name, payload, **kwargs: calls.append((name, payload)))
    min_ids = [{'type': 'gradient-range-min', 'index': p} for p in ('C8', 'C7')]
    max_ids = [{'type': 'gradient-range-max', 'index': p} for p in ('C8', 'C7')]
    result, status, disabled, polling_disabled, request = update(1, 0, 'bomex', 0, 600, ['wp2', 'thlm'], ['C8', 'C7'], ['0.3', '0.2'], ['0.7', '0.8'], 'auto', 1., 1., 1., min_ids, max_ids, 'browser1', None)
    assert calls[0][0] == 'plot_task_submit'
    assert calls[0][1]['args'] == ['bomex', 0, 600, ['wp2', 'thlm'], [{'name': 'C8', 'min': .3, 'max': .7}, {'name': 'C7', 'min': .2, 'max': .8}]]
    assert result is None and disabled and not polling_disabled and 'Queued' in status
    context.triggered_id = 'gradient-multi-poll'
    changed = update(1, 1, 'bomex', 0, 600, ['wp2', 'thlm'], ['C8', 'C7'], ['0.3', '0.2'],
                     ['0.7', '0.8'], 'reverse', 1., 1., 1., min_ids, max_ids, 'browser1', request)
    assert changed[0] is None and 'Inputs changed' in changed[1] and changed[3]
    assert len(calls) == 1
    context.triggered_id = 'gradient-multi-fields'
    result, _, disabled, polling_disabled, _ = update(1, 1, 'bomex', 0, 600, ['thlm'], ['C8', 'C7'], ['0.3', '0.2'], ['0.7', '0.8'], 'auto', 1., 1., 1., min_ids, max_ids, 'browser1', request)
    assert result is None and not disabled and polling_disabled


def test_gradient_work_has_its_own_bounded_broker_lane():
    """A gradient run does not block the existing figure lane or lose its result."""
    executors = []
    class Executor:
        def __init__(self):
            self.futures = []
            executors.append(self)
        def submit(self, function, request):
            future = Future()
            self.futures.append(future)
            return future
        def shutdown(self, **kwargs):
            pass
    queue = PlotTasks(Executor)
    try:
        queue.submit({'scope': 'browser:profile', 'revision': 1, 'task': 'profile', 'args': []})
        queue.submit({'scope': 'explorer:gradient', 'revision': 1, 'task': 'gradient_matrix', 'args': ['bomex', 0, 600, ['wp2'], [{'name': 'C8', 'min': .3, 'max': .7}]]})
        assert len(executors) == 2
        executors[1].futures[0].set_result(json.dumps({'gradient': 3.}))
        assert queue.poll('explorer:gradient', 1)['result'] == {'gradient': 3.}
        assert queue.poll('browser:profile', 1)['state'] == 'running'
    finally:
        queue.close()


def test_sample_tangents_use_original_points_and_omit_unavailable_slopes():
    """Short segments preserve the sampled curve and cached midpoint tangent."""
    result = {'values': [0., 1., 2.], 'baseline': 1., 'parameter': 'C8',
              'metrics': {'rmse': {'values': [2., 5., 10.], 'baseline': 5., 'gradient': 3.,
                                           'sample_gradients': [1., 2., None]}}}
    figure = callbacks.response_figure('rmse', 'RMSE', result)
    assert list(figure.data[0].y) == [2., 5., 10.]
    assert list(figure.data[1].y) == [2., 5., 8.]
    segments = figure.data[3]
    assert segments.name == 'Sample gradients'
    assert list(segments.x) == [0., .35, None, .65, 1.35, None]
    assert list(segments.y) == pytest.approx([2., 2.35, None, 4.3, 5.7, None])
    assert segments.line.color != figure.data[1].line.color


def test_all_sample_request_is_bound_to_the_completed_sweep(monkeypatch):
    """The button requires a completed sweep and an older sweep cannot finish into it."""
    app = Dash(__name__)
    multi_callbacks.register_callbacks(app)
    functions = {entry['callback'].__wrapped__.__name__: entry['callback'].__wrapped__
                 for entry in app.callback_map.values()}
    mapping = {'loss_weights': {'mean': 1., 'std': 1., 'corr': 1.}, 'requested_gradient_mode': 'reverse', 'map_revision': 2, 'case': 'bomex', 'time_range_seconds': [0, 600],
               'fields': ['wp2'], 'parameter_ranges': [{'name': 'C8', 'min': .3, 'max': .7}]}
    detail = {'map_revision': 2, 'sweep_revision': 3, 'parameter': 'C8', 'field': 'wp2'}
    controls = functions['show_point_controls']
    assert controls(mapping, None, True) == (True, {'display': 'none'})
    assert controls(mapping, detail, True) == (False, {'display': 'flex'})
    assert controls(mapping, detail, False)[0]
    context = SimpleNamespace(triggered_id='gradient-multi-points-go')
    monkeypatch.setattr(multi_callbacks, 'ctx', context)
    calls = []
    monkeypatch.setattr(multi_callbacks, 'perform_action', lambda name, payload, **kw: calls.append((name, payload)))
    update = functions['update_point_gradients']
    data, status, polling_disabled, request = update(1, 0, mapping, detail, 'browser', None)
    assert data is None and not polling_disabled and '16 samples' in status
    assert calls[0][1]['task'] == 'gradient_matrix_response_gradients'
    assert calls[0][1]['gradient_mode'] == 'reverse'
    assert calls[0][1]['args'] == ['C8', 'wp2', 'bomex', 0, 600, ['wp2'], mapping['parameter_ranges']]
    context.triggered_id = 'gradient-multi-points-poll'
    result = update(1, 1, mapping, dict(detail, sweep_revision=4), 'browser', request)
    assert result[0] is None and result[2]
    assert len(calls) == 1


def test_multicolumn_derivative_request_uses_all_samples_in_one_loss_call(monkeypatch):
    """Request one parameter direction for all 16 columns and preserve their derivative order."""
    import numpy as np
    from dash_app.services import gradient_response, gradient_matrix
    from utilities.loss_metrics import LOSS_METRIC_NAMES

    monkeypatch.setattr(gradient_matrix, 'validate_matrix_request', lambda *args: args)
    calls = []

    def evaluate(*args, **kwargs):
        """Return distinct per-column derivatives from the simulated loss driver."""
        calls.append((args, kwargs))
        return {'loss_weights': {'mean': 1., 'std': 1., 'corr': 1.}, 'field_names': ['wp2'], 'metrics': {name: np.ones((1, 1, 16)) for name in LOSS_METRIC_NAMES},
                'gradient': {name: np.arange(16.).reshape(1, 1, 16, 1) for name in LOSS_METRIC_NAMES},
                'gradient_mode': 'reverse', 'timings': {'total': 1.}, 'runtime': {'reused': True}}

    monkeypatch.setattr(gradient_response, 'evaluate_experiment', evaluate)
    result = gradient_response.calculate_response('C8', 'wp2', 'bomex', 0, 600, ['wp2'],
                                                  [{'name': 'C8', 'min': .3, 'max': .7}], grad=True, gradient_mode='reverse')
    assert len(calls) == 1 and calls[0][0][4] == 16
    assert calls[0][1]['grad'] and calls[0][1]['gradient_parameters'] == ['C8']
    assert calls[0][1]['gradient_mode'] == 'reverse'
    assert result['gradient_mode'] == 'reverse'
    assert not calls[0][1]['profiles']
    assert result['values'] == np.linspace(.3, .7, 16).tolist()
    assert result['metrics']['rmse']['gradients'] == list(range(16))
    assert 'profiles' not in result


@pytest.mark.parametrize('values', [[2., 5., 10.], [0., 0., 0.], [-1.e-9, -1.e-9, -1.e-9]])
def test_response_range_uses_sampled_values_not_tangents(values):
    """Extreme map and sample slopes do not stretch the sampled response axis."""
    metric = {'values': values, 'baseline': values[1], 'gradient': 1.e6,
              'sample_gradients': [1.e6, -1.e6, 1.e6]}
    result = {'values': [0., 1., 2.], 'baseline': 1., 'parameter': 'C8',
              'metrics': {'rmse': metric}}
    figure = callbacks.response_figure('rmse', 'RMSE', result)
    low, high = figure.layout.yaxis.range
    assert low < min(values) <= max(values) < high
    assert high - low < max(1., max(map(abs, values))) * 3
    metric['gradient'] = None
    metric.pop('sample_gradients')
    assert tuple(callbacks.response_figure('rmse', 'RMSE', result).layout.yaxis.range) == (low, high)
