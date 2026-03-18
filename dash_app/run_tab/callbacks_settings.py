"""Callbacks for run-tab flag/parameter editing and multicol UI state."""

from dash import (
    ALL,
    Input,
    MATCH,
    Output,
    Patch,
    State,
    ClientsideFunction,
    callback_context,
    no_update,
)

from .layout import build_multicol_row, field_style
from .state import RUN_FINALIZED, RUN_LOCK, RUN_STATUS


def blank_multicol_row(row_id):
    """Return one empty multicol hypergrid row record."""
    return {"id": row_id, "param": "", "min": "", "max": "", "npoints": ""}


def invalidate_cached_run_results(completed_cases, failed_cases):
    """Clear cached run status and completed/failed stores after a settings edit."""
    with RUN_LOCK:
        if RUN_STATUS or RUN_FINALIZED:
            RUN_STATUS.clear()
            RUN_FINALIZED.clear()

    completed_update = [] if completed_cases else no_update
    failed_update = [] if failed_cases else no_update
    return completed_update, failed_update


def register_settings_callbacks(app):
    """Register callbacks that track dirty state and changed-field styling."""

    app.clientside_callback(
        ClientsideFunction(
            namespace="runTab",
            function_name="syncParamRowClass",
        ),
        Output({"type": "run-param-container", "file": MATCH, "name": MATCH}, "className"),
        Input({"type": "run-param", "file": MATCH, "name": MATCH}, "value"),
        Input({"type": "run-param", "file": MATCH, "name": MATCH}, "disabled"),
        State({"type": "run-param", "file": MATCH, "name": MATCH}, "id"),
        State("run-defaults-by-key", "data"),
        prevent_initial_call=True,
    )

    @app.callback(
        Output("run-multicol-rows", "children"),
        Output("run-multicol-next-id", "data"),
        Output("run-multicol-row-order", "data"),
        Input("run-multicol-add", "n_clicks"),
        Input({"type": "run-hr-remove", "index": ALL}, "n_clicks"),
        State("run-multicol-next-id", "data"),
        State("run-multicol-row-order", "data"),
        State("run-tunable-names", "data"),
        prevent_initial_call=True,
    )
    def sync_multicol_rows(_add_clicks, _remove_clicks, next_id, row_order, tunable_names):
        """Add or remove multicol rows without mirroring row values through the server."""
        current_order = list(row_order or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "run-multicol-add":
            row_id = int(next_id or 0)
            patch = Patch()
            patch.append(build_multicol_row(blank_multicol_row(row_id), tunable_names or []))
            return patch, row_id + 1, current_order + [row_id]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-hr-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order:
                return no_update, no_update, no_update
            patch = Patch()
            del patch[current_order.index(remove_id)]
            current_order.remove(remove_id)
            return patch, no_update, current_order

        return no_update, no_update, no_update

    @app.callback(
        Output({"type": "run-flag-container", "name": ALL}, "style"),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-flag", "name": ALL}, "value"),
        State("run-defaults", "data"),
        State("run-flag-names", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def sync_flag_styles_and_invalidate(flag_values, defaults_data, flag_names_data, completed_cases, failed_cases):
        """Update only the affected flag row and invalidate stale run results."""
        flag_names = flag_names_data or []
        trigger_id = callback_context.triggered_id
        completed_update, failed_update = invalidate_cached_run_results(
            completed_cases, failed_cases
        )

        if not flag_names:
            return [], completed_update, failed_update

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-flag":
            triggered_name = trigger_id.get("name")
            if triggered_name in flag_names:
                flag_styles = [no_update] * len(flag_names)
                index = flag_names.index(triggered_name)
                current = bool((flag_values or [])[index]) if index < len(flag_values or []) else False
                changed = current != defaults_data["flags"].get(triggered_name)
                flag_styles[index] = field_style(changed)
                return flag_styles, completed_update, failed_update

        flag_styles = []
        for name, values in zip(flag_names, flag_values or []):
            current = bool(values)
            flag_styles.append(field_style(current != defaults_data["flags"].get(name)))
        return flag_styles, completed_update, failed_update

    @app.callback(
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        Input("run-multicol-add", "n_clicks"),
        Input({"type": "run-hr-remove", "index": ALL}, "n_clicks"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        Input({"type": "run-hr-min", "index": ALL}, "value"),
        Input({"type": "run-hr-max", "index": ALL}, "value"),
        Input({"type": "run-hr-npoints", "index": ALL}, "value"),
        Input("run-opt-max-iters", "value"),
        Input("run-opt-debug", "value"),
        Input("run-opt-dt-main", "value"),
        Input("run-opt-dt-rad", "value"),
        Input("run-opt-tout", "value"),
        Input("run-opt-out-dir", "value"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def mark_nonboolean_settings_dirty(
        _param_values,
        _multicol_add_clicks,
        _multicol_remove_clicks,
        _multicol_params,
        _multicol_mins,
        _multicol_maxes,
        _multicol_npoints,
        _opt_max_iters,
        _opt_debug,
        _opt_dt_main,
        _opt_dt_rad,
        _opt_tout,
        _opt_out_dir,
        completed_cases,
        failed_cases,
    ):
        """Invalidate stale run results on any non-boolean settings change."""
        if callback_context.triggered_id is None:
            return no_update, no_update
        return invalidate_cached_run_results(completed_cases, failed_cases)

    app.clientside_callback(
        ClientsideFunction(namespace="runTab", function_name="syncMulticolDisabled"),
        Output({"type": "run-param", "file": "tunable", "name": ALL}, "disabled"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "id"),
        State({"type": "run-param", "file": "tunable", "name": ALL}, "disabled"),
        prevent_initial_call=True,
    )
