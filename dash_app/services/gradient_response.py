"""Evaluate and format Dash gradient explorer experiments in retained workers.

execute_gradient_request() dispatches map, 16-column response and sample-gradient
requests from the Dash broker and returns their plot data as JSON.
calculate_response() formats four metrics and optional model/benchmark profiles
for the selected field; derivative requests return one slope per sample.

All model work goes through evaluate_experiment(), which retains a shared loss
runner for a fixed case, averaging interval, field set, batch width and loss
weights. close_session() releases it when the configuration changes, an
initialization/evaluation fails or the worker exits. The loss runner owns model
initialization, statistics and benchmark data. Imports here stay independent of
Dash UI packages and defer model dependencies until evaluation.
"""
from __future__ import annotations

import atexit
import json
import math
import os
import tempfile
import time

SAMPLE_COUNT = 16


def validate_request(parameter, minimum, maximum, field, case, start, end):
    """Return the validated parameter range, field, case and averaging seconds."""
    from tuner.loss_backend import available_loss_fields
    from tuner.parameter_ranges import tunable_parameter_defaults
    from tuner.case_defaults import available_case_names
    from utilities.clubb_settings_validation import parameter_hard_bounds, parse_finite_number

    if not isinstance(parameter, str) or parameter not in tunable_parameter_defaults():
        raise ValueError('Select a valid tuning parameter.')
    if not isinstance(field, str) or field not in available_loss_fields():
        raise ValueError('Select a valid loss field.')
    if not isinstance(case, str) or case not in available_case_names():
        raise ValueError('Select a case with an LES benchmark.')
    try:
        start_seconds, end_seconds = float(start), float(end)
    except (TypeError, ValueError) as exc:
        raise ValueError('Enter averaging start and end times in seconds.') from exc
    if (not math.isfinite(start_seconds) or not math.isfinite(end_seconds)
            or not start_seconds.is_integer() or not end_seconds.is_integer()
            or start_seconds < 0 or start_seconds >= end_seconds):
        raise ValueError('Averaging times must be whole seconds with 0 ≤ start < end.')
    low, high = parse_finite_number(minimum), parse_finite_number(maximum)
    if low is None or high is None or low >= high:
        raise ValueError(f'{parameter} needs finite numbers with minimum less than maximum.')
    bounds = parameter_hard_bounds().get(parameter, {})
    if (bounds.get('min') is not None and low < bounds['min']) or (
        bounds.get('max') is not None and high > bounds['max']
    ):
        raise ValueError(f'The {parameter} range exceeds its model bounds.')
    return str(parameter), low, high, str(field), case, int(start_seconds), int(end_seconds)


_session = None
_configuration = None
_directory = None


def close_session():
    """Release the loss session retained by the hosting process."""
    global _session, _directory, _configuration
    if _session is not None:
        _session.close()
        _session = None
    if _directory is not None:
        _directory.cleanup()
        _directory = None
    _configuration = None


atexit.register(close_session)


def evaluate_experiment(case, start, end, fields, batch_size, parameters, *,
                        grad=False, gradient_parameters=None, gradient_mode='auto', profiles=False, loss_weights=None, progress_file=None):
    """Evaluate a request directly in the persistent JAX worker's loss session."""
    global _session, _directory, _configuration
    from tuner.case_defaults import read_case_defaults
    from tuner.status import atomic_write_json
    from utilities.clubb_loss_runner import clubb_loss_runner, resolve_gradient_mode, validate_gradient_mode

    from utilities.loss_metrics import normalize_loss_weights
    loss_weights = normalize_loss_weights(loss_weights)
    validate_gradient_mode(gradient_mode)
    started = time.monotonic()
    configuration = (case, start, end, tuple(fields), batch_size, tuple(loss_weights.items()))
    reused = _session is not None and configuration == _configuration
    if progress_file:
        atomic_write_json(progress_file, {'step': 'runtime', 'step_started_at': time.time(), 'durations': {}})
    if not reused:
        close_session()
        _directory = tempfile.TemporaryDirectory(prefix='clubb-explorer-session-')
        _session = clubb_loss_runner('jax', case, fields, _directory.name, batch_size=batch_size,
            case_defaults=read_case_defaults(case, {'time_average_range': [start, end], 'num_time_windows': 1}),
            loss_weights=loss_weights)
        try:
            _session.init()
        except BaseException:
            close_session()
            raise
        _configuration = configuration
    runtime = time.monotonic() - started
    params = _session.default_params.copy()
    for name, values in parameters.items():
        params[:, _session.param_names.index(name)] = values
    evaluation_started = time.monotonic()
    mode = resolve_gradient_mode(gradient_mode, len(gradient_parameters or _session.param_names), len(fields))

    def progress(completed, total):
        """Publish completed derivative directions through the broker's progress file."""
        if progress_file:
            atomic_write_json(progress_file, {'step': 'gradient', 'completed': completed,
                'gradient_mode': mode,
                'total': total, 'step_started_at': time.time() - (time.monotonic() - evaluation_started),
                'durations': {'runtime': runtime}})

    if progress_file and not grad:
        atomic_write_json(progress_file, {'step': 'samples', 'step_started_at': time.time(),
                                          'durations': {'runtime': runtime}})
    try:
        result = _session.evaluate(params, grad=grad, gradient_parameters=gradient_parameters,
                                   gradient_mode=gradient_mode, gradient_target='metrics',
                                   profiles=profiles, progress=progress)
    except BaseException:
        close_session()
        raise
    result['field_names'] = _session.field_names
    result['timings'] = {'runtime': runtime, 'total': time.monotonic() - started}
    result['runtime'] = {'worker_pid': os.getpid(), 'reused': reused}
    return result


def calculate_response(parameter, field, case, start, end, fields, parameter_ranges, *, grad=False, gradient_mode='auto', loss_weights=None, progress_file=None):
    """Evaluate 16 samples, optionally returning the swept-parameter derivative at each.

    grad=True evaluates all columns together through the loss driver's gradient
    interface. Derivatives have the same sample order as values. Forward-only
    calls include profiles; derivative calls augment an already displayed sweep.
    """
    import numpy as np
    from utilities.loss_metrics import loss_metric_values_valid, loss_result_for_json
    from .gradient_matrix import validate_matrix_request

    case, start, end, fields, parameter_ranges = validate_matrix_request(case, start, end, fields, parameter_ranges)
    if field not in fields or parameter not in [item['name'] for item in parameter_ranges]:
        raise ValueError('Select a field and parameter from the current gradient map.')
    selected = next(item for item in parameter_ranges if item['name'] == parameter)
    values = np.linspace(selected['min'], selected['max'], SAMPLE_COUNT)
    updates = {item['name']: item['min'] / 2 + item['max'] / 2 for item in parameter_ranges}
    updates[parameter] = values.tolist()
    data = evaluate_experiment(case, start, end, fields, SAMPLE_COUNT, updates, grad=grad,
                               gradient_parameters=[parameter] if grad else None,
                               gradient_mode=gradient_mode,
                               profiles=not grad, loss_weights=loss_weights, progress_file=progress_file)
    fi = data['field_names'].index(field)
    metrics = {name: {'values': np.asarray(array)[0, fi].tolist()} for name, array in data['metrics'].items()}
    if any(not loss_metric_values_valid(name, item['values']).all() for name, item in metrics.items()) or max(metrics['rmse']['values']) >= 1.e30:
        raise RuntimeError('A parameter sample failed or produced invalid metrics')
    from tuner.case_defaults import read_case_defaults
    result = {'case': case, 'field': field, 'parameter': parameter, 'time_range_seconds': [start, end],
        'altitude_range_m': read_case_defaults(case)['altitude_comparison_range'], 'values': values.tolist(),
        'loss_weights': data['loss_weights'], 'metrics': metrics, 'timings': data['timings'], 'runtime': data['runtime']}
    if grad:
        result['gradient_mode'] = data['gradient_mode']
        for name in metrics:
            metrics[name]['gradients'] = [float(value) if np.isfinite(value) else None
                                          for value in data['gradient'][name][0, fi, :, 0]]
    else:
        profile = data['profiles'][field]
        result['profiles'] = {'heights_m': profile['heights_m'].tolist(), 'units': profile['units'],
                              'benchmark': profile['benchmark'][0].tolist(), 'model': profile['model'][0].tolist()}
    return loss_result_for_json(result)


def execute_gradient_request(request):
    """Execute a broker-validated task type and return its plot data as JSON."""
    if request['task'] == 'gradient_matrix':
        from .gradient_matrix import calculate_matrix
        result = calculate_matrix(*request['args'], gradient_mode=request.get('gradient_mode', 'auto'),
                                  loss_weights=request.get('loss_weights'), progress_file=request['_progress_file'])
    elif request['task'] in {'gradient_matrix_response', 'gradient_matrix_response_gradients'}:
        result = calculate_response(*request['args'], grad=request['task'] == 'gradient_matrix_response_gradients',
                                    gradient_mode=request.get('gradient_mode', 'auto'),
                                    loss_weights=request.get('loss_weights'), progress_file=request['_progress_file'])
    else:
        raise ValueError('Unknown gradient task')
    return json.dumps(result, allow_nan=False)
