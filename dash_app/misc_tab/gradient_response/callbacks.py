"""Connect gradient, profile and timing views to the broker-owned run."""
from __future__ import annotations

import time
import math

import plotly.graph_objects as go
from plotly.colors import sample_colorscale

from dash_app.shared.plot_theme import figure_theme


METRICS = (
    ('rmse', 'RMSE', 'rmse'),
    ('correlation', 'Correlation', 'correlation'),
    ('bias_norm', 'Normalized bias', 'bias-norm'),
    ('std_ratio', 'Standard-deviation ratio', 'std-ratio'),
)
STEPS = (
    ('runtime', 'Preparing JAX'),
    ('preparing', 'Preparing benchmark and case'),
    ('samples', 'Evaluating 16 parameter values together'),
    ('gradient', 'Calculating gradients for all 16 samples'),
    ('finalizing', 'Finishing'),
)


def response_figure(metric, label, result=None, theme="dark"):
    """Plot valid metric samples and tangents, identifying undefined diagnostics."""
    figure = go.Figure()
    figure.update_layout(
        **figure_theme(theme), margin=dict(l=60, r=20, t=90, b=55),
        title=dict(text=label, x=0, y=0.99, yanchor='top'), xaxis_title='Parameter value', yaxis_title=label,
        legend=dict(orientation='h', y=1.02, yanchor='bottom', x=0),
    )
    if result is None:
        figure.add_annotation(text='Select a map cell to calculate.', showarrow=False, x=.5, y=.5,
                              xref='paper', yref='paper')
        return figure
    values = result['values']
    selected = result['metrics'][metric]
    baseline, value, gradient = result['baseline'], selected['baseline'], selected['gradient']
    finite_values = [y for y in selected['values'] if y is not None and math.isfinite(y)]
    figure.add_scatter(x=values, y=selected['values'], name='CLUBB', mode='lines+markers', connectgaps=False,
                       line=dict(color='#28b7ac' if theme == 'dark' else '#087f78', width=3), marker=dict(size=7))
    baseline_defined = value is not None and math.isfinite(value)
    if baseline_defined and gradient is not None and math.isfinite(gradient):
        tangent = [value + gradient * (x - baseline) for x in values]
        figure.add_scatter(x=values, y=tangent, name='Map gradient', mode='lines',
                           line=dict(color='#df704f', width=2, dash='dash'))
        title = f'{label}  ·  slope {gradient:.5g}'
    elif not baseline_defined:
        title = f'{label}  ·  undefined at midpoint'
    elif result.get('gradient_target') == 'loss':
        title = f'{label}  ·  value only'
    elif metric == 'std_ratio' and value == 0.0:
        title = f'{label}  ·  gradient undefined at flat model'
    else:
        title = f'{label}  ·  gradient unavailable'
    if baseline_defined:
        figure.add_scatter(x=[baseline], y=[value], name='Map baseline', mode='markers',
                           marker=dict(color='#df704f', size=9, symbol='diamond'))
    slopes = selected.get('sample_gradients')
    if slopes is not None:
        half_width = min(b - a for a, b in zip(values, values[1:])) * .35
        xs, ys, hover = [], [], []
        for x, y, slope in zip(values, selected['values'], slopes):
            if y is None or not math.isfinite(y) or slope is None or not math.isfinite(slope):
                continue
            left, right = max(values[0], x - half_width), min(values[-1], x + half_width)
            xs.extend([left, right, None])
            ys.extend([y + slope * (left - x), y + slope * (right - x), None])
            hover.extend([[x, y, slope], [x, y, slope], None])
        figure.add_scatter(x=xs, y=ys, customdata=hover, name='Sample gradients', mode='lines',
                           connectgaps=False, line=dict(color='#c69afa' if theme == 'dark' else '#7550b3', width=3),
                           hovertemplate=('Parameter: %{customdata[0]:.5g}<br>Metric: %{customdata[1]:.5g}'
                                          '<br>Gradient: %{customdata[2]:.5g}<extra></extra>'))
        figure.update_layout(margin=dict(t=110))
    missing = len(selected['values']) - len(finite_values)
    if missing:
        reason = ('Correlation requires variation in both profiles.' if metric == 'correlation'
                  else 'Std ratio requires variation in the benchmark.')
        text = ('Undefined: ' + reason if not finite_values
                else f'{missing} undefined samples shown as gaps. ' + reason)
        figure.add_annotation(text=text, showarrow=False, x=.5, y=.5 if not finite_values else .02,
                              xref='paper', yref='paper')
    if finite_values:
        low, high = min(finite_values), max(finite_values)
        padding = 0.05 * (high - low if high != low else abs(low) or 1.0)
        figure.update_yaxes(range=[low - padding, high + padding])
    figure.update_layout(xaxis_title=result['parameter'], title=dict(text=title, x=0))
    return figure


def profile_figure(result=None, theme="dark"):
    """Plot the time-averaged model profiles with the benchmark profile."""
    figure = go.Figure()
    figure.update_layout(
        **figure_theme(theme), title=dict(text='Time-averaged profiles', x=0),
        margin=dict(l=70, r=25, t=55, b=60),
        xaxis_title='Field value', yaxis_title='Height (m)',
    )
    if result is None:
        figure.add_annotation(text='Profiles will appear after the run.', showarrow=False, x=.5, y=.5,
                              xref='paper', yref='paper')
        return figure
    profiles = result['profiles']
    colors = sample_colorscale('Viridis', [0.25 + 0.7 * i / (len(result['values']) - 1)
                                           for i in range(len(result['values']))])
    for value, profile, color in zip(result['values'], profiles['model'], colors):
        figure.add_scatter(x=profile, y=profiles['heights_m'], mode='lines',
                           name=f"{result['parameter']} = {value:.5g}",
                           line=dict(color=color, width=1.5))
    figure.add_scatter(x=profiles['benchmark'], y=profiles['heights_m'], name='LES', mode='lines',
                       line=dict(color='#edf3fa' if theme == 'dark' else '#172d3e', width=3, dash='dash'))
    figure.update_layout(xaxis_title=f"{result['field']} ({profiles['units']})")
    return figure


def progress_text(progress, *, finished=False):
    """Format the active step and elapsed time plus completed stage durations."""
    progress = progress or {}
    durations = progress.get('durations') or progress.get('timings') or {}
    if finished:
        heading = f"Complete in {durations.get('total', sum(durations.values())):.1f} s."
    else:
        step = progress.get('step')
        label = 'Writing results' if step == 'complete' else dict(STEPS).get(step, 'Starting')
        if step == 'gradient' and progress.get('gradient_mode'):
            label += f" · {progress['gradient_mode']} mode"
        started = progress.get('step_started_at')
        elapsed = max(0., time.time() - started) if isinstance(started, (int, float)) else 0.
        heading = f'{label} · {elapsed:.1f} s'
    completed = [f'{label}: {durations[key]:.1f} s' for key, label in STEPS if key in durations]
    return heading + ('\n' + '\n'.join(completed) if completed else '')
