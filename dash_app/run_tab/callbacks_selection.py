"""Callbacks for run-tab case and stats selection state."""

from dash import ALL, Input, Output, State, callback_context, no_update

from .layout import case_button_style, stats_button_style
from .runtime import snapshot_status_lists


def register_selection_callbacks(app, case_groups):
    """Register callbacks that manage selected cases and stats-file choices."""

    @app.callback(
        Output({"type": "run-case-button", "name": ALL}, "style"),
        Input("run-selected-cases", "data"),
        Input("run-completed-cases", "data"),
        Input("run-failed-cases", "data"),
        Input("run-running-cases", "data"),
        Input("run-queued-cases", "data"),
        State({"type": "run-case-button", "name": ALL}, "id"),
    )
    def update_case_button_styles(selected_cases, completed_cases, failed_cases, running_cases, queued_cases, ids):
        """Color case buttons based on selected, running, completed, and failed state."""
        if not ids:
            return []
        completed_global, failed_global = snapshot_status_lists()
        selected = set(selected_cases or [])
        completed = set(completed_cases or []) | set(completed_global)
        failed = set(failed_cases or []) | set(failed_global)
        running = set((running_cases or {}).keys())
        queued = {item.get("case") for item in (queued_cases or []) if item.get("case")}
        styles = []
        for case_id in ids:
            name = case_id.get("name")
            if name in failed:
                color = "#dc2626"
            elif name in completed:
                color = "#16a34a"
            elif name in running or name in queued:
                color = "linear-gradient(135deg, #dc2626 0%, #16a34a 100%)"
            else:
                color = "#2563eb"
            styles.append(case_button_style(color, name in selected))
        return styles

    @app.callback(
        Output({"type": "run-stats-button", "name": ALL}, "style"),
        Input("run-selected-stats-file", "data"),
        State({"type": "run-stats-button", "name": ALL}, "id"),
    )
    def update_stats_button_styles(selected_stats, ids):
        """Highlight the currently selected stats-file button."""
        if not ids:
            return []
        return [stats_button_style(btn_id.get("name") == selected_stats) for btn_id in ids]

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
        State("run-running-cases", "data"),
        prevent_initial_call=True,
    )
    def select_case(_n_clicks, selected_cases, running_cases):
        """Toggle one case in the selection list unless that case is already running."""
        trigger_id = callback_context.triggered_id
        case_name = trigger_id.get("name") if isinstance(trigger_id, dict) else None
        if not case_name or case_name in (running_cases or {}):
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
        group_name = trigger_id.get("name")
        group_cases = case_groups.get(group_name or "", [])
        return list(group_cases) if group_cases else no_update

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input("run-deselect", "n_clicks"),
        prevent_initial_call=True,
    )
    def deselect_all(_n_clicks):
        """Clear all currently selected cases."""
        return []
