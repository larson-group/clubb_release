"""Hidden Dash polling for transient MCP browser handoffs and job recovery."""

from __future__ import annotations

from dash import Input, Output, State, dcc, html, no_update

from dash_app.shared.activity import acknowledge_ui_request, claim_ui_request, read_activity


def dashboard_handoff():
    """Return non-visible state used by semantic MCP browser handoffs."""

    return html.Div(
        [
            dcc.Store(id="dashboard-request", data=None),
            dcc.Store(id="dashboard-broker-jobs", data={}),
            dcc.Store(id="dashboard-handoff-last-action", data=0),
            dcc.Interval(id="dashboard-handoff-interval", interval=1500, n_intervals=0),
        ],
        style={"display": "none"},
    )


def register_dashboard_handoff_callbacks(app):
    """Poll durable broker state without maintaining agent presence or chat."""

    @app.callback(
        Output("dashboard-request", "data"),
        Output("dashboard-broker-jobs", "data"),
        Output("dashboard-tabs", "value", allow_duplicate=True),
        Output("dashboard-handoff-last-action", "data"),
        Input("dashboard-handoff-interval", "n_intervals"),
        State("dashboard-handoff-last-action", "data"),
        State("dashboard-broker-jobs", "data"),
        prevent_initial_call=True,
    )
    def refresh_dashboard_handoff(_ticks, last_action_id, current_jobs):
        request = claim_ui_request(last_action_id)
        snapshot = read_activity()
        jobs = dict(snapshot.get("jobs") or {})
        request_id = int((request or {}).get("id") or 0)
        try:
            seen_id = int(last_action_id or 0)
        except (TypeError, ValueError):
            seen_id = 0
        jobs_changed = jobs != dict(current_jobs or {})
        next_request = request if request and request_id > seen_id else no_update
        next_tab = request.get("tab", "tutorial") if request and request_id > seen_id else no_update
        next_action_id = request_id if request and request_id > seen_id else no_update
        if request and request_id > seen_id:
            is_plot_state_request = (
                request.get("tab") == "plots"
                and request.get("operation") in {"set_view", "add_budget", "remove"}
            )
            if not is_plot_state_request:
                acknowledge_ui_request(request_id)
        return (
            next_request,
            jobs if jobs_changed else no_update,
            next_tab,
            next_action_id,
        )
