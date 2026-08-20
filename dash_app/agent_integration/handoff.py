"""Hidden Dash polling for transient MCP browser handoffs and job recovery."""

from __future__ import annotations

from dash import ClientsideFunction, Input, Output, State, dcc, html, no_update

from dash_app.shared.activity import acknowledge_ui_request, broker_jobs, claim_ui_request, count_active_jobs


ACTIVE_STATES = {"queued", "submitting", "running", "stopping"}


def dashboard_tab_labels(jobs):
    """Describe broker work in the tab where its controls and output live."""
    jobs = dict(jobs or {})

    compile_job = dict(jobs.get("compile") or {})
    compile_state = str(compile_job.get("state") or "")
    compile_kind = str((compile_job.get("job") or {}).get("kind") or compile_job.get("kind") or "")
    compile_label = "Compile"
    if compile_state in ACTIVE_STATES:
        action = "rebuilding" if compile_kind == "rebuild" else "compiling"
        if compile_state == "stopping":
            action = "stopping"
        compile_label = f"Compile · {action}"

    run_summary = dict(jobs.get("run_summary") or {})
    run_parts = []
    for state in ("running", "queued", "stopping"):
        count = int(run_summary.get(state) or 0)
        if count:
            run_parts.append(f"{count} {state}")
    run_label = "Run" + (" · " + " · ".join(run_parts) if run_parts else "")

    profile = dict(jobs.get("profile") or {})
    profile_state = str(profile.get("state") or "")
    profile_label = "Profile"
    if profile_state in ACTIVE_STATES:
        profile_label = "Profile · stopping" if profile_state == "stopping" else "Profile · timing"

    tune = dict(jobs.get("tune") or {})
    tune_state = str(tune.get("state") or "")
    loss_count = sum(
        str((record or {}).get("state") or "") in ACTIVE_STATES
        for record in (jobs.get("loss_runs") or {}).values()
    )
    tune_parts = []
    if tune_state in ACTIVE_STATES:
        tune_parts.append("stopping" if tune_state == "stopping" else "tuning")
    if loss_count:
        tune_parts.append(f"{loss_count} result run{'s' if loss_count != 1 else ''}")
    tune_label = "Tune" + (" · " + " · ".join(tune_parts) if tune_parts else "")
    pyplotgen = dict(jobs.get("pyplotgen") or {})
    pyplotgen_state = str(pyplotgen.get("state") or "")
    plots_label = "Plots"
    if pyplotgen_state in ACTIVE_STATES:
        plots_label = "Plots · stopping" if pyplotgen_state == "stopping" else "Plots · generating"
    return compile_label, run_label, profile_label, tune_label, plots_label


def dashboard_handoff(base_title: str = "CLUBB Dash"):
    """Return non-visible state used by semantic MCP browser handoffs."""

    return html.Div(
        [
            dcc.Store(id="dashboard-request", data=None),
            dcc.Store(id="dashboard-broker-jobs", data={}),
            dcc.Store(id="dashboard-active-process-count", data=0),
            dcc.Store(id="dashboard-base-title", data=str(base_title)),
            dcc.Store(id="dashboard-title-signal", data=None),
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
        Output("dashboard-active-process-count", "data"),
        Output("dashboard-tabs", "value", allow_duplicate=True),
        Output("dashboard-handoff-last-action", "data"),
        Input("dashboard-handoff-interval", "n_intervals"),
        State("dashboard-handoff-last-action", "data"),
        State("dashboard-broker-jobs", "data"),
        State("dashboard-active-process-count", "data"),
        prevent_initial_call=True,
    )
    def refresh_dashboard_handoff(_ticks, last_action_id, current_jobs, current_active_count):
        request = claim_ui_request(last_action_id)
        jobs = broker_jobs()
        request_id = int((request or {}).get("id") or 0)
        try:
            seen_id = int(last_action_id or 0)
        except (TypeError, ValueError):
            seen_id = 0
        jobs_changed = jobs != dict(current_jobs or {})
        active_count = count_active_jobs(jobs)
        count_changed = active_count != int(current_active_count or 0)
        next_request = request if request and request_id > seen_id else no_update
        next_tab = (
            request.get("tab", "tutorial")
            if request and request_id > seen_id and not request.get("preserve_tab")
            else no_update
        )
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
            active_count if count_changed else no_update,
            next_tab,
            next_action_id,
        )

    app.clientside_callback(
        ClientsideFunction(namespace="dashboardTitle", function_name="sync"),
        Output("dashboard-title-signal", "data"),
        Input("dashboard-active-process-count", "data"),
        State("dashboard-base-title", "data"),
    )

    @app.callback(
        Output("dashboard-tab-compile", "label"),
        Output("dashboard-tab-run", "label"),
        Output("dashboard-tab-profile", "label"),
        Output("dashboard-tab-tune", "label"),
        Output("dashboard-tab-plots", "label"),
        Input("dashboard-broker-jobs", "data"),
    )
    def update_tab_activity_labels(jobs):
        return dashboard_tab_labels(jobs)
