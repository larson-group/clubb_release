"""Callbacks for tuning-tab selection and parameter-range editing."""

from __future__ import annotations

from dash import ALL, Input, Output, Patch, State, callback_context, no_update

from .discovery import available_fields_for_cases
from .layout import build_case_config_row, build_param_range_row
from .runtime import empty_status_payload

from tuner.taylor_metrics import AGGREGATION_MODE_NAMES, LOSS_MODE_NAMES


def blank_param_range_row(row_id):
    """Return one empty tuning-range row record."""
    return {"id": row_id, "param": "", "min": "", "max": ""}


def blank_case_config_row(row_id):
    """Return one empty case-setup row record."""
    return {
        "id": row_id,
        "case_name": "",
        "time_start": "",
        "time_end": "",
        "average_time_seconds": "",
    }


def case_config_row_from_defaults(row_id, case_name, case_data):
    """Return a case-setup row initialized from discovered defaults."""
    defaults = (case_data or {}).get(case_name, {})
    time_range = defaults.get("time_average_range", ["", ""])
    num_windows = int(defaults.get("num_time_windows", 1) or 1)
    average_time_seconds = ""
    if len(time_range) > 1 and num_windows > 0:
        average_time_seconds = int((int(time_range[1]) - int(time_range[0])) / num_windows)
    return {
        "id": row_id,
        "case_name": case_name,
        "time_start": time_range[0] if len(time_range) > 0 else "",
        "time_end": time_range[1] if len(time_range) > 1 else "",
        "average_time_seconds": average_time_seconds,
    }


def selected_case_names(case_values):
    """Return nonblank case names from case row dropdown values."""
    return [
        value.strip()
        for value in (case_values or [])
        if isinstance(value, str) and value.strip()
    ]


def case_options_by_row(case_values, case_names):
    """Return case dropdown options that exclude cases selected in other rows."""
    selected_by_row = selected_case_names(case_values)
    selected_counts = {}
    for name in selected_by_row:
        selected_counts[name] = selected_counts.get(name, 0) + 1

    options_by_row = []
    for current_value in case_values or []:
        current_name = current_value.strip() if isinstance(current_value, str) else ""
        row_options = []
        for name in case_names or []:
            if selected_counts.get(name, 0) == 0 or name == current_name:
                row_options.append({"label": name, "value": name})
        options_by_row.append(row_options)
    return options_by_row


def average_time_class_names(time_start_values, time_end_values, average_time_values):
    """Return CSS classes that mark invalid average-time entries."""
    classes = []
    for raw_start, raw_end, raw_average in zip(
        time_start_values or [],
        time_end_values or [],
        average_time_values or [],
    ):
        class_name = "tune-average-time-input"
        try:
            start_value = float(raw_start)
            end_value = float(raw_end)
            average_value = float(raw_average)
        except (TypeError, ValueError):
            classes.append(class_name)
            continue
        if int(start_value) != start_value or int(end_value) != end_value or int(average_value) != average_value:
            classes.append(class_name + " invalid")
            continue
        start = int(start_value)
        end = int(end_value)
        average = int(average_value)
        if end <= start or average < 1 or (end - start) % average != 0:
            class_name += " invalid"
        classes.append(class_name)
    return classes


def parameter_options_by_row(param_values, tunable_names):
    """Return dropdown options that exclude parameters selected in other rows."""
    selected_by_row = [
        value.strip()
        for value in (param_values or [])
        if isinstance(value, str) and value.strip()
    ]
    selected_counts = {}
    for name in selected_by_row:
        selected_counts[name] = selected_counts.get(name, 0) + 1

    options_by_row = []
    for current_value in param_values or []:
        current_name = current_value.strip() if isinstance(current_value, str) else ""
        row_options = []
        for name in tunable_names or []:
            if selected_counts.get(name, 0) == 0 or name == current_name:
                row_options.append({"label": name, "value": name})
        options_by_row.append(row_options)
    return options_by_row


def register_settings_callbacks(app):
    """Register callbacks that manage case defaults and range-row state."""

    @app.callback(
        Output("tune-strategy-mode", "data"),
        Input("tune-mode-random", "n_clicks"),
        Input("tune-mode-resolve", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_strategy_mode(_random_clicks, _resolve_clicks, active_job):
        """Persist the selected tuning strategy mode."""
        if active_job:
            return no_update
        trigger_id = callback_context.triggered_id
        if trigger_id == "tune-mode-random":
            return "random"
        if trigger_id == "tune-mode-resolve":
            return "resolve"
        return no_update

    @app.callback(
        Output("tune-loss-mode", "data"),
        Input("tune-loss-mode-scaled-rmse", "n_clicks"),
        Input("tune-loss-mode-centered-rmse-bias", "n_clicks"),
        Input("tune-loss-mode-taylor-components", "n_clicks"),
        Input("tune-loss-mode-taylor-components-squared", "n_clicks"),
        Input("tune-loss-mode-shape-first", "n_clicks"),
        Input("tune-loss-mode-bias-light-taylor", "n_clicks"),
        Input("tune-loss-mode-decomposed-taylor", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_loss_mode(
        _scaled_rmse_clicks,
        _centered_clicks,
        _components_clicks,
        _squared_clicks,
        _shape_first_clicks,
        _bias_light_clicks,
        _decomposed_clicks,
        active_job,
    ):
        """Persist the selected Taylor loss policy mode."""
        if active_job:
            return no_update
        mode_by_button = {
            "tune-loss-mode-scaled-rmse": "scaled_rmse",
            "tune-loss-mode-centered-rmse-bias": "centered_rmse_bias",
            "tune-loss-mode-taylor-components": "taylor_components",
            "tune-loss-mode-taylor-components-squared": "taylor_components_squared",
            "tune-loss-mode-shape-first": "shape_first",
            "tune-loss-mode-bias-light-taylor": "bias_light_taylor",
            "tune-loss-mode-decomposed-taylor": "decomposed_taylor",
        }
        selected = mode_by_button.get(callback_context.triggered_id)
        if selected in LOSS_MODE_NAMES:
            return selected
        return no_update

    @app.callback(
        Output("tune-aggregation-mode", "data"),
        Input("tune-aggregation-mean-max", "n_clicks"),
        Input("tune-aggregation-mean-worst-quantile", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_aggregation_mode(_mean_max_clicks, _worst_quantile_clicks, active_job):
        """Persist the selected loss aggregation policy mode."""
        if active_job:
            return no_update
        mode_by_button = {
            "tune-aggregation-mean-max": "mean_max",
            "tune-aggregation-mean-worst-quantile": "mean_worst_quantile",
        }
        selected = mode_by_button.get(callback_context.triggered_id)
        if selected in AGGREGATION_MODE_NAMES:
            return selected
        return no_update

    @app.callback(
        Output("tune-case-rows", "children"),
        Output("tune-case-next-id", "data"),
        Output("tune-case-row-order", "data"),
        Input("tune-case-add", "n_clicks"),
        Input({"type": "tune-case-remove", "index": ALL}, "n_clicks"),
        State("tune-case-next-id", "data"),
        State("tune-case-row-order", "data"),
        State({"type": "tune-case-name", "index": ALL}, "value"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def sync_case_rows(_add_clicks, _remove_clicks, next_id, row_order, case_values, case_data, active_job):
        """Add or remove case setup rows without mirroring row values through the server."""
        if active_job:
            return no_update, no_update, no_update

        current_order = list(row_order or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "tune-case-add":
            row_id = int(next_id or 0)
            case_names = sorted((case_data or {}).keys())
            row_options = case_options_by_row(list(case_values or []) + [""], case_names)[-1]
            patch = Patch()
            patch.append(build_case_config_row(blank_case_config_row(row_id), row_options))
            return patch, row_id + 1, current_order + [row_id]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "tune-case-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order:
                return no_update, no_update, no_update
            patch = Patch()
            del patch[current_order.index(remove_id)]
            current_order.remove(remove_id)
            return patch, no_update, current_order

        return no_update, no_update, no_update

    @app.callback(
        Output({"type": "tune-case-name", "index": ALL}, "options"),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
    )
    def sync_case_dropdown_options(case_values, case_data, active_job):
        """Remove already-selected cases from the other case dropdown rows."""
        if active_job:
            return no_update
        return case_options_by_row(case_values or [], sorted((case_data or {}).keys()))

    @app.callback(
        Output({"type": "tune-case-time-start", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-time-end", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-average-time", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State("tune-case-row-order", "data"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def autofill_case_defaults(
        case_values,
        time_start_values,
        time_end_values,
        average_time_values,
        row_order,
        case_data,
        active_job,
    ):
        """Fill editable case settings from case defaults when a case is selected."""
        if active_job:
            return no_update, no_update, no_update

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-case-name":
            return no_update, no_update, no_update

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return no_update, no_update, no_update
        row_pos = current_order.index(row_id)
        case_name = (case_values or [None])[row_pos] if row_pos < len(case_values or []) else None
        if not case_name or case_name not in (case_data or {}):
            return no_update, no_update, no_update

        defaults = case_config_row_from_defaults(row_id, case_name, case_data)
        outputs = [
            list(time_start_values or []),
            list(time_end_values or []),
            list(average_time_values or []),
        ]
        while any(len(values) < len(case_values or []) for values in outputs):
            for values in outputs:
                while len(values) < len(case_values or []):
                    values.append("")

        outputs[0][row_pos] = defaults["time_start"]
        outputs[1][row_pos] = defaults["time_end"]
        outputs[2][row_pos] = defaults["average_time_seconds"]
        return tuple(outputs)

    @app.callback(
        Output({"type": "tune-case-average-time", "index": ALL}, "className"),
        Input({"type": "tune-case-time-start", "index": ALL}, "value"),
        Input({"type": "tune-case-time-end", "index": ALL}, "value"),
        Input({"type": "tune-case-average-time", "index": ALL}, "value"),
    )
    def mark_average_time_validity(time_start_values, time_end_values, average_time_values):
        """Highlight average-time values that do not divide the selected window."""
        return average_time_class_names(time_start_values, time_end_values, average_time_values)

    @app.callback(
        Output("tune-range-rows", "children"),
        Output("tune-range-next-id", "data"),
        Output("tune-range-row-order", "data"),
        Input("tune-range-add", "n_clicks"),
        Input({"type": "tune-range-remove", "index": ALL}, "n_clicks"),
        State("tune-range-next-id", "data"),
        State("tune-range-row-order", "data"),
        State("tune-tunable-names", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def sync_range_rows(_add_clicks, _remove_clicks, next_id, row_order, tunable_names, active_job):
        """Add or remove tuning-range rows without mirroring row values through the server."""
        if active_job:
            return no_update, no_update, no_update

        current_order = list(row_order or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "tune-range-add":
            row_id = int(next_id or 0)
            patch = Patch()
            patch.append(build_param_range_row(blank_param_range_row(row_id), tunable_names or []))
            return patch, row_id + 1, current_order + [row_id]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "tune-range-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order:
                return no_update, no_update, no_update
            patch = Patch()
            del patch[current_order.index(remove_id)]
            current_order.remove(remove_id)
            return patch, no_update, current_order

        return no_update, no_update, no_update

    @app.callback(
        Output({"type": "tune-range-param", "index": ALL}, "options"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        Input("tune-range-row-order", "data"),
        State("tune-tunable-names", "data"),
    )
    def sync_parameter_options(param_values, _row_order, tunable_names):
        """Hide parameters already selected in other tuning rows."""
        return parameter_options_by_row(param_values or [], tunable_names or [])

    @app.callback(
        Output({"type": "tune-range-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-range-row-order", "data"),
        State("tune-tunable-default-ranges", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def autofill_tune_range(_param_values, min_values, max_values, row_order, default_ranges, active_job):
        """Fill min/max from the selected parameter's default value."""
        if active_job:
            return no_update, no_update

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-range-param":
            return no_update, no_update

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return no_update, no_update
        row_pos = current_order.index(row_id)

        param_values = list(_param_values or [])
        if row_pos >= len(param_values):
            return no_update, no_update
        param_name = param_values[row_pos]
        derived_range = (default_ranges or {}).get(param_name)
        if not derived_range:
            return no_update, no_update

        updated_min = list(min_values or [])
        updated_max = list(max_values or [])
        while len(updated_min) < len(param_values):
            updated_min.append("")
        while len(updated_max) < len(param_values):
            updated_max.append("")

        updated_min[row_pos] = derived_range.get("min", "")
        updated_max[row_pos] = derived_range.get("max", "")
        return updated_min, updated_max

    @app.callback(
        Output("tune-field-selector", "options"),
        Output("tune-field-selector", "value"),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
    )
    def sync_case_fields(case_values, case_data, active_job):
        """Refresh field options and default field selection when the case changes."""
        if active_job:
            return no_update, no_update

        case_names = selected_case_names(case_values)
        fields = available_fields_for_cases(case_names, case_data or {})
        defaults = (case_data or {}).get(case_names[0], {}) if case_names else {}
        return (
            [{"label": field_name, "value": field_name} for field_name in fields],
            [field_name for field_name in defaults.get("default_clubb_fields", []) if field_name in fields],
        )

    @app.callback(
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-case-add", "n_clicks"),
        Input({"type": "tune-case-remove", "index": ALL}, "n_clicks"),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        Input({"type": "tune-case-time-start", "index": ALL}, "value"),
        Input({"type": "tune-case-time-end", "index": ALL}, "value"),
        Input({"type": "tune-case-average-time", "index": ALL}, "value"),
        Input("tune-field-selector", "value"),
        Input("tune-range-add", "n_clicks"),
        Input({"type": "tune-range-remove", "index": ALL}, "n_clicks"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        Input({"type": "tune-range-min", "index": ALL}, "value"),
        Input({"type": "tune-range-max", "index": ALL}, "value"),
        Input("tune-batch-size", "value"),
        Input("tune-max-workers", "value"),
        Input("tune-strategy-mode", "data"),
        Input("tune-loss-mode", "data"),
        Input("tune-aggregation-mode", "data"),
        Input("tune-random-max-samples", "value"),
        Input("tune-resolve-spacing", "value"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def invalidate_tuning_state(
        _case_add_clicks,
        _case_remove_clicks,
        _case_names,
        _case_time_start,
        _case_time_end,
        _case_average_time,
        _selected_fields,
        _add_clicks,
        _remove_clicks,
        _param_values,
        _min_values,
        _max_values,
        _batch_size,
        _max_workers,
        _strategy_mode,
        _loss_mode,
        _aggregation_mode,
        _random_max_samples,
        _resolve_spacing,
        active_job,
    ):
        """Clear stale tuning results when inactive; active jobs keep polling."""
        if active_job:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update
        return {}, empty_status_payload(), [], [], {}, True, ""
