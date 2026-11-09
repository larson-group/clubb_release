"""Build the field-by-parameter gradient explorer."""
import uuid
import sys

from dash import dcc, html
from dash_app.shared.components import styled_dropdown
from tuner.parameter_ranges import default_parameter_ranges
from tuner.case_defaults import read_case_defaults


def parameter_range_rows(parameters, existing=None):
    """Render editable parameter endpoints, preserving already-entered values."""
    defaults = default_parameter_ranges()
    existing = existing or {}
    for bounds in defaults.values():
        if bounds['min'] == sys.float_info.min and bounds['default'] > 0:
            bounds['min'] = max(bounds['min'], bounds['default'] / 4)
    return [html.Div([
        html.Div(parameter, className='gradient-range-name'),
        html.Div([
            html.Div([html.Label('Min'), dcc.Input(id={'type': 'gradient-range-min', 'index': parameter},
                                                   type='text', value=existing.get(parameter, {}).get('min', repr(defaults[parameter]['min'])), debounce=True)]),
            html.Div([html.Label('Max'), dcc.Input(id={'type': 'gradient-range-max', 'index': parameter},
                                                   type='text', value=existing.get(parameter, {}).get('max', repr(defaults[parameter]['max'])), debounce=True)]),
        ], className='gradient-response-range'),
    ], className='gradient-range-row') for parameter in parameters if parameter in defaults]


def build_multi_layout(params, fields, cases):
    """Return selection controls, a gradient map and its response detail plots."""
    default_case = 'bomex'
    default_parameters = ['C8', 'Lscale_tke_multiplier', 'beta']
    default_start, default_end = read_case_defaults(default_case)['time_average_range']
    return html.Section([
        dcc.Store(id='gradient-multi-scope', data=uuid.uuid4().hex),
        dcc.Store(id='gradient-multi-request'),
        dcc.Store(id='gradient-multi-result'),
        dcc.Store(id='gradient-multi-detail-request'),
        dcc.Store(id='gradient-multi-detail-result'),
        dcc.Store(id='gradient-multi-points-request'),
        dcc.Store(id='gradient-multi-points-result'),
        dcc.Interval(id='gradient-multi-poll', interval=1500, disabled=True),
        dcc.Interval(id='gradient-multi-detail-poll', interval=1500, disabled=True),
        dcc.Interval(id='gradient-multi-points-poll', interval=1500, disabled=True),
        html.Div([
            html.Div('JAX EXPLORER', className='gradient-response-eyebrow'),
            html.H2('Loss and gradient'),
            html.P('Compare field metrics and gradients across selected parameters. Select a field to explore its response.'),
        ], className='gradient-multi-header'),
        html.Div([
            html.Div([
                html.Div('EXPERIMENT', className='gradient-response-eyebrow'),
                html.Label('Case', htmlFor='gradient-multi-case'),
                styled_dropdown(id='gradient-multi-case', options=[{'label': c, 'value': c} for c in cases], value=default_case, clearable=False),
                html.Label('Parameters', htmlFor='gradient-multi-parameters'),
                styled_dropdown(id='gradient-multi-parameters', options=[{'label': p, 'value': p} for p in params],
                                value=default_parameters, multi=True),
                html.Div(parameter_range_rows(default_parameters), id='gradient-multi-ranges'),
                html.Div(id='gradient-multi-range-error', role='status', className='gradient-range-error'),
                html.Label('Fields', htmlFor='gradient-multi-fields'),
                styled_dropdown(id='gradient-multi-fields', options=[{'label': f, 'value': f} for f in fields],
                                value=['wp2', 'cloud_frac'], multi=True),
                html.Label('Taylor RMSE weights'),
                html.Div([
                    html.Div([html.Label(label), dcc.Input(
                        id=f'gradient-multi-weight-{name}', type='number', value=1,
                        min=0, step='any', debounce=True,
                    )]) for name, label in (('mean', 'Mean'), ('std', 'Std'), ('corr', 'Pattern'))
                ], className='gradient-response-range'),
                html.P('Strictly positive. Unit weights recover normalized RMSE. The other three diagnostics keep their definitions.',
                       className='gradient-response-help'),
                html.P('Undefined shape diagnostics are marked explicitly; the RMSE remains available for flat profiles.',
                       className='gradient-response-help'),
                html.Label('Gradient mode' , htmlFor='gradient-multi-mode'),
                styled_dropdown(id='gradient-multi-mode', options=[
                    {'label': 'Auto', 'value': 'auto'},
                    {'label': 'Forward', 'value': 'forward'},
                    {'label': 'Backward (reverse)', 'value': 'reverse'},
                ], value='auto', clearable=False),
                html.Label('Gradient target', htmlFor='gradient-multi-target'),
                styled_dropdown(id='gradient-multi-target', options=[
                    {'label': 'Each field and metric', 'value': 'metrics'},
                    {'label': 'Combined RMSE loss', 'value': 'loss'},
                ], value='metrics', clearable=False),
                html.P('Combined RMSE loss differentiates one scalar across all selected fields. All four field metrics remain visible.',
                       className='gradient-response-help'),
                html.P('Auto uses backward differentiation when parameters outnumber fields, or when more than one parameter is selected for combined loss.',
                       className='gradient-response-help'),
                html.Div([
                    html.Div([html.Label('Average start (s)', htmlFor='gradient-multi-start'),
                              dcc.Input(id='gradient-multi-start', type='number', value=default_start, min=0, step=1, debounce=True)]),
                    html.Div([html.Label('Average end (s)', htmlFor='gradient-multi-end'),
                              dcc.Input(id='gradient-multi-end', type='number', value=default_end, min=0, step=1, debounce=True)]),
                ], className='gradient-response-range'),
                html.Button('Calculate maps', id='gradient-multi-go', n_clicks=0, className='gradient-response-go clubb-activity-stripes'),
                html.P('Each parameter is set to the midpoint of its range. Gradient colors show the predicted change for a positive step of 10% of that range.', className='gradient-response-help'),
                html.P('The model runs from the case’s initial time through the selected averaging end.', className='gradient-response-help'),
                html.Div(id='gradient-multi-status', role='status', className='gradient-response-status'),
            ], className='gradient-response-controls'),
            html.Div([
                html.P('Select a field in a metric chart to calculate its response.', className='gradient-response-help'),
                html.Div([dcc.Graph(id='gradient-multi-combined', config={'displaylogo': False, 'responsive': True})],
                         id='gradient-multi-combined-container', className='gradient-response-chart gradient-multi-map-card'),
                html.Div([
                    html.Div([dcc.Graph(id=f'gradient-multi-map-{metric}', config={'displaylogo': False, 'responsive': True})],
                             id=f'gradient-multi-map-{metric}-container', className='gradient-response-chart gradient-multi-map-card')
                    for metric in ('rmse', 'correlation', 'bias-norm', 'std-ratio')
                ], className='gradient-response-metric-grid'),
                html.Div(id='gradient-multi-error', role='alert', className='gradient-range-error'),
                html.Div(id='gradient-multi-detail-status', role='status', className='gradient-response-status gradient-multi-detail-status'),
                html.Div([
                    html.Button('Calculate all maps', id='gradient-multi-points-go', n_clicks=0, disabled=True,
                                title='Calculate the swept-parameter gradient at all 16 sample values.',
                                className='gradient-response-points-go clubb-activity-stripes'),
                    html.Div(id='gradient-multi-points-status', role='status', className='gradient-response-status'),
                ], id='gradient-multi-points-controls', className='gradient-multi-points-controls', style={'display': 'none'}),
                html.Div([
                    html.Div([dcc.Graph(id=f'gradient-multi-curve-{metric}', config={'displaylogo': False, 'responsive': True})],
                             className='gradient-response-chart')
                    for metric in ('rmse', 'correlation', 'bias-norm', 'std-ratio')
                ], className='gradient-response-metric-grid'),
                html.Div([dcc.Graph(id='gradient-multi-profiles', config={'displaylogo': False, 'responsive': True})],
                         className='gradient-response-chart gradient-response-profile'),
            ], className='gradient-response-results'),
        ], className='gradient-response-grid'),
    ])
