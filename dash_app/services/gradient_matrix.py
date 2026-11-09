"""Prepare field diagnostics and gradients for the Dash gradient explorer.

validate_matrix_request() checks field and parameter selections without loading
the model. calculate_matrix() evaluates midpoint parameters in a persistent
loss session and returns four metric values per field. It requests either
field-by-parameter metric gradients or one combined RMSE-loss gradient per
parameter. The loss driver owns differentiation; this module formats plot data
and timings for the browser.
"""
from __future__ import annotations

from .gradient_response import evaluate_experiment, validate_request, validate_parameter_range


def validate_matrix_request(case, start, end, fields, parameter_ranges):
    """Validate the comparison settings and each selected parameter range."""
    from tuner.loss_backend import available_loss_fields
    from tuner.parameter_ranges import default_parameter_ranges
    from utilities.clubb_settings_validation import parameter_hard_bounds

    ranges = default_parameter_ranges()
    if not isinstance(parameter_ranges, list) or not parameter_ranges:
        raise ValueError('Select one or more parameters and enter their ranges.')
    if not isinstance(fields, list) or not fields or len(set(fields)) != len(fields):
        raise ValueError('Select one or more distinct fields.')
    available_fields = set(available_loss_fields())
    if any(not isinstance(f, str) or f not in available_fields for f in fields):
        raise ValueError('Select valid loss fields.')
    normalized = []
    seen = set()
    hard_bounds = parameter_hard_bounds()
    for item in parameter_ranges:
        if not isinstance(item, dict):
            raise ValueError('Enter a range for every selected parameter.')
        name = item.get('name')
        if not isinstance(name, str) or name not in ranges or name in seen:
            raise ValueError('Select distinct valid tuning parameters.')
        seen.add(name)
        if not normalized:
            name, low, high, _, case, start, end = validate_request(
                name, item.get('min'), item.get('max'), fields[0], case, start, end)
        else:
            low, high = validate_parameter_range(name, item.get('min'), item.get('max'), hard_bounds)
        normalized.append({'name': name, 'min': low, 'max': high})
    return case, start, end, fields, normalized


def calculate_matrix(case, start, end, fields, parameter_ranges, *, gradient_mode='auto', gradient_target='metrics',
                     loss_weights=None, progress_file=None):
    """Return field metrics and either their gradients or the combined-loss gradient."""
    import numpy as np
    from utilities.loss_metrics import loss_metric_values_valid, loss_result_for_json
    from utilities.clubb_loss_runner import validate_gradient_target
    validate_gradient_target(gradient_target)
    case, start, end, fields, parameter_ranges = validate_matrix_request(case, start, end, fields, parameter_ranges)
    parameters = [item['name'] for item in parameter_ranges]
    baseline = [item['min'] / 2 + item['max'] / 2 for item in parameter_ranges]
    data = evaluate_experiment(case, start, end, fields, 1, dict(zip(parameters, baseline)),
                              grad=True, gradient_parameters=parameters, gradient_mode=gradient_mode,
                              gradient_target=gradient_target, loss_weights=loss_weights, progress_file=progress_file)
    metrics = {}
    for name, array in data['metrics'].items():
        values = np.asarray(array, dtype=float)[0, :, 0]
        if not loss_metric_values_valid(name, values).all() or (name == 'rmse' and np.any(values >= 1.e30)):
            raise RuntimeError('The baseline model run failed or produced invalid metrics')
        metrics[name] = {'values': loss_result_for_json(values)}
        if gradient_target == 'metrics':
            metrics[name]['gradients'] = [[float(v) if np.isfinite(v) else None for v in row]
                                          for row in data['gradient'][name][0, :, 0]]
    from tuner.case_defaults import read_case_defaults
    return {'case': case, 'time_range_seconds': [start, end], 'fields': fields,
        'parameters': parameters, 'parameter_ranges': parameter_ranges, 'baseline': baseline,
        'altitude_range_m': read_case_defaults(case)['altitude_comparison_range'],
        'loss_weights': data['loss_weights'], 'metrics': metrics, 'timings': data['timings'], 'runtime': data['runtime'],
        'requested_gradient_mode': gradient_mode, 'gradient_mode': data['gradient_mode'],
        'gradient_target': gradient_target,
        'loss': float(np.asarray(data['loss'])[0]) if gradient_target == 'loss' else None,
        'loss_gradients': [float(v) if np.isfinite(v) else None for v in np.asarray(data['gradient']['loss'])[0]]
            if gradient_target == 'loss' else None}
