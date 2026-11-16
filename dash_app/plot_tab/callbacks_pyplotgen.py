"""Plot-tab callbacks for broker-owned PyPlotGen exports."""

from __future__ import annotations

import time
import subprocess

from dash import ALL, ClientsideFunction, Input, Output, State, callback_context, dcc, html, no_update

from dash_app.shared.broker_client import perform_action
from dash_app.shared.notecard import notecard
from dash_app.persistence import WORKSPACE_TOKEN
from .async_callbacks import task_callback
from .pyplotgen_runtime import pyplotgen_help, read_pyplotgen_progress


ACTIVE_STATES = {"queued", "submitting", "running", "stopping"}


def _current_job(broker_snapshot, action):
    broker_job = dict((broker_snapshot or {}).get("pyplotgen") or {})
    action_job = dict((action or {}).get("job") or {})
    if (
        (action or {}).get("kind") == "stopping"
        and action_job.get("run_id") == broker_job.get("run_id")
        and str(broker_job.get("state") or "") in ACTIVE_STATES
    ):
        return action_job
    if float(action_job.get("started_at") or 0) > float(broker_job.get("started_at") or 0):
        return action_job
    return broker_job or action_job


def register_pyplotgen_callbacks(app):
    """Wire help-derived export controls to the existing broker-owned export lifecycle."""
    @task_callback(app, 'pyplotgen_help',
        Output('plots-pyplotgen-help', 'data'),
        Input('plots-pyplotgen-options', 'id'),
    )
    def load_pyplotgen_help(_mounted):
        """Load the CLI's help in a worker so tab mounting does not wait for imports."""
        try:
            return pyplotgen_help()
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            return {'error': f'Could not read PyPlotGen -h: {exc}', 'options': [], 'text': ''}

    @app.callback(
        Output('plots-pyplotgen-options', 'options'),
        Output('plots-pyplotgen-help-status', 'children'),
        Input('plots-pyplotgen-help', 'data'),
    )
    def update_pyplotgen_options(help_data):
        """Populate searchable option names and hover descriptions from the CLI help."""
        data = help_data or {}
        return ([{'label': option['flag'], 'value': option['flag'], 'title': option['description']}
                 for option in data.get('options', [])], data.get('error') or '')

    @app.callback(
        Output('plots-pyplotgen-option-values', 'children'),
        Input('plots-pyplotgen-options', 'value'),
        Input('plots-pyplotgen-help', 'data'),
        State({'type': 'plots-pyplotgen-option-value', 'flag': ALL}, 'value'),
        State({'type': 'plots-pyplotgen-option-value', 'flag': ALL}, 'id'),
    )
    def update_pyplotgen_value_fields(selected, help_data, values, ids):
        """Show value fields only for selected options whose help lists a value hint."""
        previous = {item['flag']: value for item, value in zip(ids or [], values or [])}
        catalog = {option['flag']: option for option in (help_data or {}).get('options', [])}
        fields = []
        for flag in selected or []:
            metavar = catalog.get(flag, {}).get('metavar')
            if not metavar:
                continue
            fields.append(html.Label([
                html.Span(flag + ' ' + metavar),
                dcc.Input(id={'type': 'plots-pyplotgen-option-value', 'flag': flag},
                          type='text', value=previous.get(flag, ''), placeholder=metavar,
                          persistence=WORKSPACE_TOKEN, persistence_type='local'),
            ], title=catalog[flag]['description']))
        return fields

    @app.callback(
        Output('plots-pyplotgen-help-modal', 'children'),
        Input('plots-pyplotgen-help-open', 'n_clicks'),
        Input({'type': 'plots-pyplotgen-help-close', 'index': ALL}, 'n_clicks'),
        Input('plots-pyplotgen-help', 'data'),
        State('plots-pyplotgen-help-modal', 'children'),
        prevent_initial_call=True,
    )
    def show_pyplotgen_help(open_clicks, close_clicks, help_data, current):
        """Show PyPlotGen's verbatim -h text and a short explanation of gallery export."""
        trigger = callback_context.triggered_id
        if isinstance(trigger, dict) and trigger.get('type') == 'plots-pyplotgen-help-close':
            return '' if any(value for value in close_clicks or []) else no_update
        if trigger == 'plots-pyplotgen-help-open' and not open_clicks:
            return no_update
        if trigger != 'plots-pyplotgen-help-open' and not current:
            return no_update
        data = help_data or {}
        return notecard('PyPlotGen', [
            html.P('PyPlotGen creates a static gallery of CLUBB plots from the active output folders. '
                   'The completed gallery opens in a new tab.'),
            html.P('Search for options, press Enter to select one, and repeat to add more. '
                   'Fill in any value fields below the selector; quote paths containing spaces. '
                   'Dash supplies -clubb and -output_dir from the active outputs and export folder.'),
            html.Pre(data.get('error') or data.get('text') or 'Loading PyPlotGen -h…',
                     className='plots-pyplotgen-help-text'),
        ], {'type': 'plots-pyplotgen-help-close', 'index': 0}, subtitle='pyplotgen/pyplotgen.py -h', size='large')

    @app.callback(
        Output("plots-pyplotgen-action", "data"),
        Input("plots-pyplotgen-run", "n_clicks"),
        State("plots-output-dirs", "data"),
        State('plots-pyplotgen-options', 'value'),
        State({'type': 'plots-pyplotgen-option-value', 'flag': ALL}, 'value'),
        State({'type': 'plots-pyplotgen-option-value', 'flag': ALL}, 'id'),
        State("dashboard-broker-jobs", "data"),
        State("plots-pyplotgen-action", "data"),
        prevent_initial_call=True,
    )
    def launch_pyplotgen(clicks, output_dirs, options, values, ids, broker_snapshot, current_action):
        """Forward selected flags/values, or cancel an active export through its broker owner."""
        if not clicks:
            return no_update
        try:
            current_job = _current_job(broker_snapshot, current_action)
            if str(current_job.get("state") or "") in {"queued", "submitting", "running"}:
                response = perform_action("stop_pyplotgen_request", {}, internal=True)
                return {"kind": "stopping", "at": time.time(), "job": dict(response.get("job") or {})}
            response = perform_action(
                "launch_pyplotgen_request",
                {"output_dirs": list(output_dirs or []), 'options': list(options or []),
                 'option_values': {item['flag']: value for item, value in zip(ids or [], values or [])}},
                internal=True,
            )
            return {"kind": "started", "at": time.time(), "job": dict(response.get("job") or {})}
        except (OSError, RuntimeError, ValueError) as exc:
            return {"kind": "error", "at": time.time(), "message": str(exc)}

    @app.callback(
        Output("plots-pyplotgen-progress-interval", "disabled"),
        Input("dashboard-broker-jobs", "data"),
        Input("plots-pyplotgen-action", "data"),
    )
    def toggle_pyplotgen_progress(broker_snapshot, action):
        return str(_current_job(broker_snapshot, action).get("state") or "") not in ACTIVE_STATES

    @app.callback(
        Output("plots-pyplotgen-run", "disabled"),
        Output("plots-pyplotgen-run", "className"),
        Output("plots-pyplotgen-run-label", "children"),
        Output("plots-pyplotgen-status", "children"),
        Input("plots-pyplotgen-progress-interval", "n_intervals"),
        Input("dashboard-broker-jobs", "data"),
        Input("plots-pyplotgen-action", "data"),
    )
    def render_pyplotgen_status(_ticks, broker_snapshot, action):
        action = dict(action or {})
        job = _current_job(broker_snapshot, action)
        state = str(job.get("state") or "")
        if action.get("kind") == "error":
            return False, "plots-pyplotgen-button plots-pyplotgen-button--error", "Generate PyPlotGen", html.Span(
                str(action.get("message") or "PyPlotGen could not be started"),
                className="plots-pyplotgen-error",
            )
        if state in ACTIVE_STATES:
            stopping = state == "stopping"
            label = "Stopping PyPlotGen…" if stopping else "Cancel PyPlotGen"
            progress = read_pyplotgen_progress(str(job.get("log_path") or ""))
            current, total = progress or (0, 0)
            percentage = min(100.0, 100.0 * current / total) if total else 0.0
            progress_bar = html.Div(
                [
                    html.Div(className="plots-pyplotgen-progress-fill", style={"width": f"{percentage:.1f}%"}),
                    html.Span(f"{current} / {total}" if total else "0 / ?", className="plots-pyplotgen-progress-text"),
                ],
                className="plots-pyplotgen-progress",
            )
            return stopping, "plots-pyplotgen-button plots-pyplotgen-button--active", label, progress_bar
        if state == "finished":
            html_path = str(job.get("html_path") or "")
            return (
                False,
                "plots-pyplotgen-button plots-pyplotgen-button--finished",
                "Generate PyPlotGen",
                html.A(
                    html_path,
                    id="plots-pyplotgen-html-path",
                    href=str(job.get("html_url") or ""),
                    target="_blank",
                    className="plots-pyplotgen-result-link",
                ),
            )
        if state in {"stopped", "cancelled"}:
            return False, "plots-pyplotgen-button", "Generate PyPlotGen", "PyPlotGen cancelled"
        if state == "error":
            detail = str(job.get("error") or f"PyPlotGen exited with status {job.get('returncode')}")
            return False, "plots-pyplotgen-button plots-pyplotgen-button--error", "Generate PyPlotGen", html.Span(
                detail, className="plots-pyplotgen-error"
            )
        return False, "plots-pyplotgen-button", "Generate PyPlotGen", ""

    app.clientside_callback(
        ClientsideFunction(namespace="plotsPyplotgen", function_name="reserveWindow"),
        Output("plots-pyplotgen-popup", "data"),
        Input("plots-pyplotgen-run", "n_clicks"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        ClientsideFunction(namespace="plotsPyplotgen", function_name="openCompleted"),
        Output("plots-pyplotgen-opened-run", "data"),
        Input("dashboard-broker-jobs", "data"),
        Input("plots-pyplotgen-action", "data"),
        State("plots-pyplotgen-popup", "data"),
        State("plots-pyplotgen-opened-run", "data"),
    )
