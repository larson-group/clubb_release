"""Callbacks for run-tab case and stats selection state."""

import json

from dash import ALL, ClientsideFunction, Input, Output, State

from .layout import stats_button_style


def register_selection_callbacks(app, case_groups):
    """Register callbacks that manage selected cases and stats-file choices."""

    @app.callback(
        Output({"type": "run-stats-button", "name": ALL}, "style"),
        Output("run-stats-custom-open", "style"),
        Input("run-selected-stats-file", "data"),
        State({"type": "run-stats-button", "name": ALL}, "id"),
    )
    def update_stats_button_styles(selected_stats, ids):
        """Highlight the currently selected stats-file button."""
        if not ids:
            return [], stats_button_style(False)
        return [
            stats_button_style(button_id.get("name") == selected_stats)
            for button_id in ids
        ], stats_button_style(selected_stats not in {"all", "standard", "core"})

    app.clientside_callback(
        """function(clicks, selected) {
            const ctx = dash_clientside.callback_context;
            const event = (ctx.triggered || [])[0] || {};
            const id = ctx.triggered_id;
            if (!(Number(event.value) > 0) || !id || !id.name) return dash_clientside.no_update;
            const next = Array.from(selected || []);
            const index = next.indexOf(id.name);
            if (index < 0) next.push(id.name); else next.splice(index, 1);
            return next;
        }""",
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-case-button", "name": ALL}, "n_clicks"),
        State("run-selected-cases", "data"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        """function(clicks) {
            const ctx = dash_clientside.callback_context;
            const event = (ctx.triggered || [])[0] || {};
            const id = ctx.triggered_id;
            if (!(Number(event.value) > 0) || !id) return dash_clientside.no_update;
            return GROUPS[id.name] || dash_clientside.no_update;
        }""".replace("GROUPS", json.dumps(case_groups)),
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-group-button", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        "function(clicks) { return clicks > 0 ? [] : dash_clientside.no_update; }",
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input("run-deselect", "n_clicks"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace",
            function_name="clearRunTransientState",
        ),
        Output("run-clear-persistence-signal", "data"),
        Input("run-clear", "n_clicks"),
        prevent_initial_call=True,
    )
