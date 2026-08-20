"""Callbacks for run-tab case and stats selection state."""

from dash import ALL, ClientsideFunction, Input, Output, State, callback_context, no_update

from .layout import stats_button_style


def register_selection_callbacks(app, case_groups):
    """Register callbacks that manage selected cases and stats-file choices."""

    @app.callback(
        Output({"type": "run-stats-button", "name": ALL}, "style"),
        Input("run-selected-stats-file", "data"),
        State({"type": "run-stats-button", "name": ALL}, "id"),
    )
    def update_stats_button_styles(selected_stats, ids):
        """Highlight the currently selected stats-file button."""
        if not ids:
            return []
        return [
            stats_button_style(button_id.get("name") == selected_stats)
            for button_id in ids
        ]

    @app.callback(
        Output("run-selected-stats-file", "data", allow_duplicate=True),
        Input({"type": "run-stats-button", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def select_stats_file(_n_clicks):
        """Persist the selected stats-file name from the clicked button."""
        trigger_id = callback_context.triggered_id
        if isinstance(trigger_id, dict):
            return trigger_id.get("name", no_update)
        return no_update

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-case-button", "name": ALL}, "n_clicks"),
        State("run-selected-cases", "data"),
        prevent_initial_call=True,
    )
    def select_case(_n_clicks, selected_cases):
        """Toggle one case in the next-submission selection."""
        trigger_id = callback_context.triggered_id
        case_name = trigger_id.get("name") if isinstance(trigger_id, dict) else None
        if not case_name:
            return no_update
        selected = list(selected_cases or [])
        if case_name in selected:
            selected.remove(case_name)
        else:
            selected.append(case_name)
        return selected

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-group-button", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def select_case_group(_n_clicks):
        """Replace the current selection with one predefined case group."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict):
            return no_update
        group_cases = case_groups.get(trigger_id.get("name") or "", [])
        return list(group_cases) if group_cases else no_update

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input("run-deselect", "n_clicks"),
        prevent_initial_call=True,
    )
    def deselect_all(_deselect_clicks):
        """Clear the next-submission case selection."""
        return []

    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace",
            function_name="clearRunTransientState",
        ),
        Output("run-clear-persistence-signal", "data"),
        Input("run-clear", "n_clicks"),
        prevent_initial_call=True,
    )
