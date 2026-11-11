"""Connect the gradient maps and selected-cell response to Dash."""
import math
import time

from dash import ALL, Input, Output, State, ctx, no_update
import plotly.graph_objects as go

from dash_app.misc_tab.gradient_response.callbacks import METRICS, profile_figure, progress_text, response_figure
from dash_app.services.gradient_matrix import validate_matrix_request
from dash_app.services.gradient_response import validate_request
from dash_app.shared.broker_client import perform_action
from dash_app.shared.plot_theme import figure_theme
from tuner.case_defaults import read_case_defaults
from utilities.loss_metrics import normalize_loss_weights
from utilities.clubb_settings_validation import parse_finite_number
from utilities.clubb_loss_runner import validate_gradient_mode
from .multi_layout import parameter_range_rows


def matrix_figure(data, metric, theme="dark"):
    """Display gradients scaled by 10% of each entered parameter interval."""
    figure = go.Figure()
    label = next(label for name, label, _ in METRICS if name == metric)
    figure.update_layout(**figure_theme(theme), margin=dict(l=105, r=35, t=55, b=100),
                         title=dict(text=label, x=0, y=0.99, yanchor='top'),
                         xaxis_title='Parameter', yaxis_title='Field')
    if not data or metric not in data['metrics']:
        figure.add_annotation(text='Calculate a map to compare gradients.', showarrow=False,
                              x=.5, y=.5, xref='paper', yref='paper')
        return figure
    ranges = {item['name']: item for item in data['parameter_ranges']}
    steps = [(ranges[p]['max'] - ranges[p]['min']) * 0.1 for p in data['parameters']]
    derivatives = data['metrics'][metric]['gradients']
    values = data['metrics'][metric]['values']
    matrix = [[None if v is None else v * steps[j] for j, v in enumerate(row)] for row in derivatives]
    finite = [abs(v) for row in matrix for v in row if v is not None and math.isfinite(v)]
    limit = max(finite, default=0) or 1
    custom = [[[('undefined' if values[i] is None else 'unavailable') if derivatives[i][j] is None else f'{derivatives[i][j]:.5g}',
                steps[j], data['baseline'][j], 'undefined' if values[i] is None else f'{values[i]:.5g}']
               for j in range(len(data['parameters']))] for i in range(len(data['fields']))]
    figure.add_trace(go.Heatmap(
        x=data['parameters'], y=data['fields'], z=matrix, customdata=custom,
        colorscale='RdBu_r', zmin=-limit, zmax=limit, zmid=0,
        colorbar=dict(title='Δ metric'), hoverongaps=False,
        hovertemplate=('Field: %{y}<br>Parameter: %{x}<br>Baseline: %{customdata[2]:.5g}'
                       '<br>Metric: %{customdata[3]}<br>Raw gradient: %{customdata[0]}'
                       '<br>Positive step: %{customdata[1]:.5g}<br>Predicted change: %{z:.5g}<extra></extra>'),
    ))
    unavailable = [(i, j) for i, row in enumerate(derivatives) for j, slope in enumerate(row) if slope is None]
    if unavailable:
        labels = [('undefined' if values[i] is None else
                   'gradient<br>undefined' if metric == 'std_ratio' and values[i] == 0.0 else
                   'gradient<br>unavailable') for i, j in unavailable]
        figure.add_scatter(x=[data['parameters'][j] for i, j in unavailable],
                           y=[data['fields'][i] for i, j in unavailable],
                           mode='text', text=labels, showlegend=False,
                           textfont=dict(color='#9ca3af', size=12),
                           customdata=[custom[i][j] for i, j in unavailable],
                           hovertemplate=('Field: %{y}<br>Parameter: %{x}<br>Metric: %{customdata[3]}'
                                          '<br>Gradient: %{customdata[0]}<extra></extra>'))
    figure.update_yaxes(autorange='reversed')
    return figure


def matrix_progress(progress, finished=False):
    """Format the active map stage and elapsed time."""
    progress = progress or {}
    if finished:
        return (f"Complete in {(progress.get('timings') or {}).get('total', 0):.1f} s"
                f" · {progress['gradient_mode']} mode. Click a cell to explore it.")
    step = progress.get('step')
    label = {'runtime': 'Preparing JAX', 'preparing': 'Preparing benchmark and case',
             'gradient': 'Calculating gradients', 'finalizing': 'Finishing', 'complete': 'Writing results'}.get(step, 'Starting')
    if step == 'gradient':
        label += f" · {progress.get('gradient_mode', 'auto')} mode · {progress.get('total', 0)} parameters"
    started = progress.get('step_started_at')
    elapsed = max(0., time.time() - started) if isinstance(started, (float, int)) else 0.
    return f'{label} · {elapsed:.1f} s'


def current_sweep(mapping, detail):
    """Return whether the completed sweep belongs to the currently displayed map."""
    return bool(mapping and detail and detail.get('map_revision') == mapping['map_revision']
                and 'sweep_revision' in detail)


def register_callbacks(app):
    """Register independent map, selection and detail callbacks."""
    @app.callback(Output('gradient-multi-start', 'value'), Output('gradient-multi-end', 'value'),
                  Input('gradient-multi-case', 'value'), prevent_initial_call=True)
    def case_time_range(case):
        """Set the selected case's usual comparison window in seconds."""
        return tuple(read_case_defaults(case)['time_average_range']) if case else (no_update, no_update)

    @app.callback(Output('gradient-multi-ranges', 'children'), Input('gradient-multi-parameters', 'value'),
                  State({'type': 'gradient-range-min', 'index': ALL}, 'id'),
                  State({'type': 'gradient-range-min', 'index': ALL}, 'value'),
                  State({'type': 'gradient-range-max', 'index': ALL}, 'id'),
                  State({'type': 'gradient-range-max', 'index': ALL}, 'value'))
    def show_parameter_ranges(parameters, min_ids, min_values, max_ids, max_values):
        """Keep edited endpoints when parameters are added or removed."""
        existing = {}
        for ids, values, key in ((min_ids, min_values, 'min'), (max_ids, max_values, 'max')):
            for item, value in zip(ids or [], values or []):
                existing.setdefault(item['index'], {})[key] = value
        return parameter_range_rows(parameters or [], existing)

    @app.callback(Output('gradient-multi-range-error', 'children'),
                  Input('gradient-multi-parameters', 'value'),
                  Input({'type': 'gradient-range-min', 'index': ALL}, 'value'),
                  Input({'type': 'gradient-range-max', 'index': ALL}, 'value'),
                  State({'type': 'gradient-range-min', 'index': ALL}, 'id'),
                  State({'type': 'gradient-range-max', 'index': ALL}, 'id'))
    def show_range_error(parameters, min_values, max_values, min_ids, max_ids):
        """Show invalid parameter endpoints before a map is submitted."""
        mins = {item['index']: value for item, value in zip(min_ids or [], min_values or [])}
        maxs = {item['index']: value for item, value in zip(max_ids or [], max_values or [])}
        for name in parameters or []:
            if name not in mins or name not in maxs:
                continue
            try:
                validate_request(name, mins[name], maxs[name], 'wp2', 'bomex', 0, 600)
            except ValueError as exc:
                return str(exc)
        return ''

    @app.callback(Output('gradient-multi-result', 'data'), Output('gradient-multi-status', 'children'),
                  Output('gradient-multi-go', 'disabled'), Output('gradient-multi-poll', 'disabled'),
                  Output('gradient-multi-request', 'data'),
                  Input('gradient-multi-go', 'n_clicks'), Input('gradient-multi-poll', 'n_intervals'),
                  Input('gradient-multi-case', 'value'), Input('gradient-multi-start', 'value'),
                  Input('gradient-multi-end', 'value'), Input('gradient-multi-fields', 'value'),
                  Input('gradient-multi-parameters', 'value'),
                  Input({'type': 'gradient-range-min', 'index': ALL}, 'value'),
                  Input({'type': 'gradient-range-max', 'index': ALL}, 'value'),
                  Input('gradient-multi-mode', 'value'),
                  *[Input(f'gradient-multi-weight-{name}', 'value') for name in ('mean', 'std', 'corr')],
                  State({'type': 'gradient-range-min', 'index': ALL}, 'id'),
                  State({'type': 'gradient-range-max', 'index': ALL}, 'id'),
                  State('gradient-multi-scope', 'data'), State('gradient-multi-request', 'data'))
    def update_map(clicks, _ticks, case, start, end, fields, parameters, min_values, max_values, gradient_mode, mean_weight, std_weight, corr_weight, min_ids, max_ids, scope, previous):
        """Submit and poll the map, clearing results when its inputs change."""
        previous = previous or {}
        trigger = ctx.triggered_id
        min_by_name = {item['index']: value for item, value in zip(min_ids or [], min_values or [])}
        max_by_name = {item['index']: value for item, value in zip(max_ids or [], max_values or [])}
        parameter_ranges = [{'name': name, 'min': min_by_name.get(name), 'max': max_by_name.get(name)}
                            for name in parameters or []]
        if trigger == 'gradient-multi-go' and clicks:
            try:
                args = list(validate_matrix_request(case, start, end, fields, parameter_ranges))
                loss_weights = normalize_loss_weights({'mean': mean_weight, 'std': std_weight, 'corr': corr_weight})
                validate_gradient_mode(gradient_mode)
                request = {'scope': f'{scope}:map', 'revision': int(previous.get('revision', 0)) + 1,
                           'selection': 0, 'task': 'gradient_matrix', 'args': args, 'gradient_mode': gradient_mode,
                           'loss_weights': loss_weights}
                perform_action('plot_task_submit', request)
                return None, 'Queued. The first run may take a few minutes.', True, False, request
            except (ValueError, RuntimeError) as exc:
                return None, str(exc), False, True, previous
        if trigger == 'gradient-multi-poll' and previous.get('scope'):
            try:
                current_ranges = [{'name': item['name'], 'min': parse_finite_number(item['min']),
                                   'max': parse_finite_number(item['max'])} for item in parameter_ranges]
                current_args = [case, parse_finite_number(start), parse_finite_number(end), fields, current_ranges]
                if (current_args != previous['args'] or gradient_mode != previous['gradient_mode']
                        or normalize_loss_weights({'mean': mean_weight, 'std': std_weight, 'corr': corr_weight}) != previous['loss_weights']):
                    return None, 'Inputs changed; calculate a new map.', False, True, previous
                result = perform_action('plot_task_poll', {'scope': previous['scope'], 'revision': previous['revision']}, ensure_running=False)
            except (ValueError, RuntimeError) as exc:
                return no_update, str(exc), False, True, previous
            state = result.get('state')
            if state == 'queued':
                return no_update, 'Queued behind another multi-field calculation.', True, False, previous
            if state == 'running':
                return no_update, matrix_progress(result.get('progress')), True, False, previous
            if state == 'ready':
                data = result['result']
                data['map_revision'] = previous['revision']
                return data, matrix_progress(data, finished=True), False, True, previous
            return None, str(result.get('error') or 'The request is no longer available.'), False, True, previous
        return None, '', False, True, previous

    @app.callback(*[Output(f'gradient-multi-map-{element}', 'figure') for _, _, element in METRICS],
                  *[Output(f'gradient-multi-map-{element}-container', 'style') for _, _, element in METRICS],
                  Input('gradient-multi-result', 'data'), Input('theme-store', 'data'))
    def show_map(data, theme):
        """Show all four metric maps for the current baseline."""
        height = max(388, 150 + 32 * len(data['fields'])) if data else 388
        return (*[matrix_figure(data, metric, theme) for metric, _, _ in METRICS],
                *[{'height': f'{height}px'} for _ in METRICS])

    @app.callback(Output('gradient-multi-detail-result', 'data'), Output('gradient-multi-detail-status', 'children'),
                  Output('gradient-multi-detail-poll', 'disabled'), Output('gradient-multi-detail-request', 'data'),
                  *[Input(f'gradient-multi-map-{element}', 'clickData') for _, _, element in METRICS],
                  Input('gradient-multi-result', 'data'),
                  Input('gradient-multi-detail-poll', 'n_intervals'),
                  State('gradient-multi-scope', 'data'), State('gradient-multi-detail-request', 'data'),
                  prevent_initial_call=True)
    def update_detail(click_rmse, click_correlation, click_bias, click_std, data, _ticks, scope, previous):
        """Run and poll the independent response sweep for the clicked map cell."""
        previous = previous or {}
        trigger = ctx.triggered_id
        if trigger == 'gradient-multi-result':
            return None, '', True, previous
        clicks = dict(zip((f'gradient-multi-map-{element}' for _, _, element in METRICS),
                          (click_rmse, click_correlation, click_bias, click_std)))
        click = clicks.get(trigger)
        if click and data and click.get('points'):
            point = click['points'][0]
            parameter, field = point.get('x'), point.get('y')
            if parameter not in data['parameters'] or field not in data['fields']:
                return None, 'Select a valid map cell.', True, previous
            baseline = data['baseline'][data['parameters'].index(parameter)]
            note = f'{field} × {parameter}; sweep centered on baseline {baseline:.5g}.'
            try:
                args = [parameter, field, data['case'], *data['time_range_seconds'],
                        data['fields'], data['parameter_ranges']]
                request = {'scope': f'{scope}:detail', 'revision': int(previous.get('revision', 0)) + 1,
                           'selection': 0, 'task': 'gradient_matrix_response', 'args': args,
                           'note': note, 'map_revision': data['map_revision'], 'loss_weights': data['loss_weights']}
                perform_action('plot_task_submit', request)
                return None, note + ' Queued.', False, request
            except (ValueError, RuntimeError) as exc:
                return None, str(exc), True, previous
        if (trigger == 'gradient-multi-detail-poll' and previous.get('scope') and data
                and previous.get('map_revision') == data.get('map_revision')):
            try:
                result = perform_action('plot_task_poll', {'scope': previous['scope'], 'revision': previous['revision']}, ensure_running=False)
            except RuntimeError as exc:
                return no_update, str(exc), True, previous
            state = result.get('state')
            if state == 'queued':
                return no_update, previous['note'] + ' Queued.', False, previous
            if state == 'running':
                return no_update, previous['note'] + '\n' + progress_text(result.get('progress')), False, previous
            if state == 'ready':
                detail = result['result']
                detail['map_revision'] = previous['map_revision']
                detail['sweep_revision'] = previous['revision']
                return detail, previous['note'] + '\n' + progress_text({'timings': detail['timings']}, finished=True), True, previous
            return None, str(result.get('error') or 'The response is no longer available.'), True, previous
        return None, '', True, previous

    @app.callback(Output('gradient-multi-points-go', 'disabled'),
                  Output('gradient-multi-points-controls', 'style'),
                  Input('gradient-multi-result', 'data'), Input('gradient-multi-detail-result', 'data'),
                  Input('gradient-multi-points-poll', 'disabled'))
    def show_point_controls(mapping, detail, polling_disabled):
        """Offer the all-sample calculation only after the current sweep completes."""
        available = current_sweep(mapping, detail)
        return not available or not polling_disabled, {'display': 'flex' if available else 'none'}

    @app.callback(Output('gradient-multi-points-result', 'data'),
                  Output('gradient-multi-points-status', 'children'),
                  Output('gradient-multi-points-poll', 'disabled'),
                  Output('gradient-multi-points-request', 'data'),
                  Input('gradient-multi-points-go', 'n_clicks'),
                  Input('gradient-multi-points-poll', 'n_intervals'),
                  Input('gradient-multi-result', 'data'), Input('gradient-multi-detail-result', 'data'),
                  State('gradient-multi-scope', 'data'), State('gradient-multi-points-request', 'data'),
                  prevent_initial_call=True)
    def update_point_gradients(clicks, _ticks, mapping, detail, scope, previous):
        """Request a multicolumn derivative pass and discard results from older sweeps."""
        previous = previous or {}
        trigger = ctx.triggered_id
        if trigger in {'gradient-multi-result', 'gradient-multi-detail-result'} or not current_sweep(mapping, detail):
            return None, '', True, previous
        if trigger == 'gradient-multi-points-go' and clicks:
            request = {'scope': f'{scope}:points', 'revision': int(previous.get('revision', 0)) + 1,
                       'selection': 0, 'task': 'gradient_matrix_response_gradients',
                       'gradient_mode': mapping['requested_gradient_mode'],
                       'loss_weights': mapping['loss_weights'],
                       'args': [detail['parameter'], detail['field'], mapping['case'], *mapping['time_range_seconds'],
                                mapping['fields'], mapping['parameter_ranges']],
                       'map_revision': detail['map_revision'], 'sweep_revision': detail['sweep_revision']}
            try:
                perform_action('plot_task_submit', request)
                return None, 'Queued gradients for all 16 samples.', False, request
            except (ValueError, RuntimeError) as exc:
                return None, str(exc), True, previous
        if (trigger == 'gradient-multi-points-poll' and previous.get('scope')
                and previous.get('map_revision') == detail['map_revision']
                and previous.get('sweep_revision') == detail['sweep_revision']):
            try:
                result = perform_action('plot_task_poll', {'scope': previous['scope'], 'revision': previous['revision']},
                                        ensure_running=False)
            except RuntimeError as exc:
                return no_update, str(exc), True, previous
            state = result.get('state')
            if state == 'queued':
                return no_update, 'Queued behind another gradient calculation.', False, previous
            if state == 'running':
                return no_update, progress_text(result.get('progress')), False, previous
            if state == 'ready':
                points = result['result']
                points.update(map_revision=previous['map_revision'], sweep_revision=previous['sweep_revision'])
                message = progress_text({'timings': points['timings']}, finished=True)
                message += f" · {points['gradient_mode']} mode"
                message += '\nPurple segments show the gradient at each sample.'
                if any(slope is None for metric, _, _ in METRICS for slope in points['metrics'][metric]['gradients']):
                    message += ' Unavailable slopes are omitted.'
                return points, message, True, previous
            return None, str(result.get('error') or 'The sample gradients are no longer available.'), True, previous
        return None, '', True, previous

    @app.callback(*[Output(f'gradient-multi-curve-{element}', 'figure') for _, _, element in METRICS],
                  Output('gradient-multi-profiles', 'figure'),
                  Input('gradient-multi-detail-result', 'data'), Input('gradient-multi-result', 'data'),
                  Input('gradient-multi-points-result', 'data'), Input('theme-store', 'data'))
    def show_detail(data, mapping, points, theme):
        """Overlay cached midpoint and optional all-sample tangents on the existing sweep."""
        if data and mapping and data.get('map_revision') == mapping.get('map_revision'):
            field_index = mapping['fields'].index(data['field'])
            param_index = mapping['parameters'].index(data['parameter'])
            result = dict(data, baseline=mapping['baseline'][param_index], metrics={
                metric: dict(data['metrics'][metric], baseline=mapping['metrics'][metric]['values'][field_index],
                             gradient=mapping['metrics'][metric]['gradients'][field_index][param_index])
                for metric, _, _ in METRICS
            })
            if (points and current_sweep(mapping, data) and points.get('map_revision') == data['map_revision']
                    and points.get('sweep_revision') == data['sweep_revision']
                    and points['parameter'] == data['parameter'] and points['field'] == data['field']
                    and points['values'] == data['values']):
                for metric, _, _ in METRICS:
                    result['metrics'][metric]['sample_gradients'] = points['metrics'][metric]['gradients']
        else:
            result = None
        return (*[response_figure(metric, label, result, theme) for metric, label, _ in METRICS], profile_figure(result, theme))
