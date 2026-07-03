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

from tunable_configs import tunable_config_names

from .config_state import build_tunable_config_state
from .layout import build_multicol_row, build_right_pane, field_style
from .state import RUN_FINALIZED, RUN_LOCK, RUN_STATUS


def blank_multicol_row(row_id):
    """Return one empty multicol hypergrid row record."""
    return {"id": row_id, "param": "", "min": "", "max": "", "npoints": "4"}


def default_multicol_row(row_id, tunable_names, default_ranges=None, selected_params=None):
    """Return a new multicol row initialized to the first useful tunable parameter."""
    selected = {str(value or "").strip() for value in (selected_params or [])}
    selected.discard("")
    available = [name for name in (tunable_names or []) if name]
    param = next((name for name in available if name not in selected), available[0] if available else "")
    row = blank_multicol_row(row_id)
    row["param"] = param
    derived = (default_ranges or {}).get(param, {})
    row["min"] = derived.get("min", "")
    row["max"] = derived.get("max", "")
    return row


def _fallback_rows_by_id(rows):
    """Return saved multicol rows indexed by row id."""
    return {row.get("id"): row for row in (rows or []) if isinstance(row, dict)}


def multicol_rows_from_values(
    row_order,
    param_values,
    min_values,
    max_values,
    npoint_values,
    tunable_names,
    default_ranges=None,
    fallback_rows=None,
):
    """Preserve live multicol rows across config changes."""
    available = set(tunable_names or [])
    fallback_by_id = _fallback_rows_by_id(fallback_rows)
    row_ids = list(row_order or [])
    if not row_ids and fallback_by_id:
        row_ids = [row.get("id") for row in (fallback_rows or []) if isinstance(row, dict)]
    rows = []
    for index, row_id in enumerate(row_ids):
        fallback = fallback_by_id.get(row_id, {})
        param = fallback.get("param", "") if index >= len(param_values or []) else param_values[index]
        param = str(param).strip() if param is not None else ""
        min_value = fallback.get("min", "") if index >= len(min_values or []) or min_values[index] is None else str(min_values[index])
        max_value = fallback.get("max", "") if index >= len(max_values or []) or max_values[index] is None else str(max_values[index])
        if param and param in available:
            derived = (default_ranges or {}).get(param, {})
            if not min_value.strip() and derived.get("min"):
                min_value = derived["min"]
            if not max_value.strip() and derived.get("max"):
                max_value = derived["max"]
        rows.append(
            {
                "id": row_id,
                "param": param,
                "min": min_value,
                "max": max_value,
                "npoints": fallback.get("npoints", "4") if index >= len(npoint_values or []) or npoint_values[index] is None else str(npoint_values[index]),
            }
        )
    if rows:
        return rows
    return [blank_multicol_row(0)]


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

    @app.callback(
        Output("run-selected-config", "data"),
        Output("run-defaults", "data"),
        Output("run-defaults-by-key", "data"),
        Output("run-flag-names", "data"),
        Output("run-param-meta", "data"),
        Output("run-tunable-names", "data"),
        Output("run-tunable-default-ranges", "data"),
        Output("run-right-pane", "children"),
        Output("run-multicol-next-id", "data", allow_duplicate=True),
        Output("run-multicol-row-order", "data", allow_duplicate=True),
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-config-button", "name": ALL}, "n_clicks"),
        State("run-tunable-configs", "data"),
        State("run-selected-config", "data"),
        State("run-multicol-next-id", "data"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State({"type": "run-hr-npoints", "index": ALL}, "value"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        prevent_initial_call=True,
    )
    def select_tunable_config(
        _clicks,
        configs,
        current_config,
        next_id,
        row_order,
        saved_multicol_rows,
        multicol_param_values,
        multicol_min_values,
        multicol_max_values,
        multicol_npoint_values,
        completed_cases,
        failed_cases,
    ):
        """Switch the run tab to another tunable config and rebuild config-dependent controls."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "run-config-button":
            return (no_update,) * 13

        config_name = str(trigger_id.get("name", "")).strip()
        if not config_name or config_name not in tunable_config_names(configs):
            return (no_update,) * 13
        if config_name == current_config:
            return (no_update,) * 13

        config_state = build_tunable_config_state(config_name, configs)
        multicol_rows = multicol_rows_from_values(
            row_order,
            multicol_param_values,
            multicol_min_values,
            multicol_max_values,
            multicol_npoint_values,
            config_state["tunable_names"],
            config_state["tunable_default_ranges"],
            saved_multicol_rows,
        )
        multicol_row_order = [row["id"] for row in multicol_rows]
        multicol_next_id = max([int(row_id) for row_id in multicol_row_order], default=-1) + 1
        multicol_next_id = max(int(next_id or 0), multicol_next_id)
        right_pane_data = {"tunable_configs": configs or [], "multicol_rows": multicol_rows, **config_state}
        completed_update, failed_update = invalidate_cached_run_results(completed_cases, failed_cases)

        return (
            config_state["selected_config"],
            config_state["defaults"],
            config_state["defaults_by_key"],
            config_state["flag_names"],
            config_state["param_meta"],
            config_state["tunable_names"],
            config_state["tunable_default_ranges"],
            build_right_pane(right_pane_data),
            multicol_next_id,
            multicol_row_order,
            multicol_rows,
            completed_update,
            failed_update,
        )

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
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Input("run-multicol-add", "n_clicks"),
        Input({"type": "run-hr-remove", "index": ALL}, "n_clicks"),
        State("run-multicol-next-id", "data"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        State("run-tunable-names", "data"),
        State("run-tunable-default-ranges", "data"),
        State({"type": "run-hr-param", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def sync_multicol_rows(add_clicks, remove_clicks, next_id, row_order, rows_state, tunable_names, default_ranges, param_values):
        """Add or remove multicol rows without mirroring row values through the server."""
        current_order = list(row_order or [])
        current_rows = list(rows_state or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "run-multicol-add":
            if not add_clicks:
                return no_update, no_update, no_update, no_update
            row_id = int(next_id or 0)
            patch = Patch()
            row = default_multicol_row(row_id, tunable_names or [], default_ranges or {}, param_values or [])
            patch.append(build_multicol_row(row, tunable_names or []))
            return patch, row_id + 1, current_order + [row_id], current_rows + [row]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "run-hr-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order:
                return no_update, no_update, no_update, no_update
            remove_pos = current_order.index(remove_id)
            if remove_pos >= len(remove_clicks or []) or not (remove_clicks or [])[remove_pos]:
                return no_update, no_update, no_update, no_update
            patch = Patch()
            del patch[remove_pos]
            current_order.remove(remove_id)
            updated_rows = [row for row in current_rows if row.get("id") != remove_id]
            return patch, no_update, current_order, updated_rows

        return no_update, no_update, no_update, no_update

    @app.callback(
        Output("run-multicol-rows-state", "data", allow_duplicate=True),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        Input({"type": "run-hr-min", "index": ALL}, "value"),
        Input({"type": "run-hr-max", "index": ALL}, "value"),
        Input({"type": "run-hr-npoints", "index": ALL}, "value"),
        State("run-multicol-row-order", "data"),
        State("run-multicol-rows-state", "data"),
        prevent_initial_call=True,
    )
    def persist_multicol_rows(param_values, min_values, max_values, npoint_values, row_order, saved_multicol_rows):
        """Persist live multicol values outside the right-pane subtree."""
        if callback_context.triggered_id is None:
            return no_update
        return multicol_rows_from_values(
            row_order,
            param_values,
            min_values,
            max_values,
            npoint_values,
            [],
            {},
            saved_multicol_rows,
        )

    @app.callback(
        Output({"type": "run-hr-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "run-hr-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
        State({"type": "run-hr-min", "index": ALL}, "value"),
        State({"type": "run-hr-max", "index": ALL}, "value"),
        State("run-multicol-row-order", "data"),
        State("run-tunable-default-ranges", "data"),
        prevent_initial_call=True,
    )
    def autofill_multicol_bounds(param_values, min_values, max_values, row_order, default_ranges):
        """Fill multicol min/max from the selected parameter's default value."""
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "run-hr-param":
            return no_update, no_update

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return no_update, no_update
        row_pos = current_order.index(row_id)

        values = list(param_values or [])
        if row_pos >= len(values):
            return no_update, no_update
        param_name = values[row_pos]
        derived_range = (default_ranges or {}).get(param_name)
        if not derived_range:
            return no_update, no_update

        updated_min = list(min_values or [])
        updated_max = list(max_values or [])
        while len(updated_min) < len(values):
            updated_min.append("")
        while len(updated_max) < len(values):
            updated_max.append("")

        updated_min[row_pos] = derived_range.get("min", "")
        updated_max[row_pos] = derived_range.get("max", "")
        return updated_min, updated_max

    @app.callback(
        Output("run-batch-size", "disabled"),
        Input({"type": "run-hr-param", "index": ALL}, "value"),
    )
    def sync_batch_size_disabled(multicol_param_values):
        """Enable batch size only when a multicol parameter is selected."""
        return not any(str(value or "").strip() for value in (multicol_param_values or []))

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
        Input("run-batch-size", "value"),
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
        _batch_size,
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
