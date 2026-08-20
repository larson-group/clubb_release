"""Plot-tab callbacks for broker-owned PyPlotGen exports."""

from __future__ import annotations

import time

from dash import ClientsideFunction, Input, Output, State, html, no_update

from dash_app.shared.broker_client import perform_action
from .pyplotgen_runtime import read_pyplotgen_progress


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
    @app.callback(
        Output("plots-pyplotgen-action", "data"),
        Input("plots-pyplotgen-run", "n_clicks"),
        State("plots-loaded-output-dirs", "data"),
        State("dashboard-broker-jobs", "data"),
        State("plots-pyplotgen-action", "data"),
        prevent_initial_call=True,
    )
    def launch_pyplotgen(clicks, output_dirs, broker_snapshot, current_action):
        if not clicks:
            return no_update
        try:
            current_job = _current_job(broker_snapshot, current_action)
            if str(current_job.get("state") or "") in {"queued", "submitting", "running"}:
                response = perform_action("stop_pyplotgen_request", {}, internal=True)
                return {"kind": "stopping", "at": time.time(), "job": dict(response.get("job") or {})}
            response = perform_action(
                "launch_pyplotgen_request",
                {"output_dirs": list(output_dirs or [])},
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
