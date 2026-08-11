"""Callbacks for tuning-tab selection and parameter-range editing."""

from __future__ import annotations

from dash import ALL, Input, Output, Patch, State, callback_context, no_update

from .discovery import (
    available_fields_for_cases,
    load_tunable_default_ranges,
    load_tunable_names,
    tunable_config_names,
)
from .layout import (
    DEFAULT_AVERAGE_TIME_SECONDS,
    build_case_config_row,
    build_config_buttons,
    build_param_range_row,
)
from .runtime import empty_status_payload
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.clubb_settings_validation import is_independently_tunable

from tuner.taylor_metrics import LOSS_MODE_NAMES, TIME_WINDOW_AGGREGATION_SCOPES
from tuner.request import evaluate_tune_settings


def blank_param_range_row(row_id):
    """Return one empty tuning-range row record."""
    return {"id": row_id, "param": "", "targets": [], "min": "", "max": ""}


def blank_case_config_row(row_id):
    """Return one empty case-setup row record."""
    return {
        "id": row_id,
        "case_name": "",
        "time_start": "",
        "time_end": "",
        "average_time_seconds": DEFAULT_AVERAGE_TIME_SECONDS,
        "altitude_min": "",
        "altitude_max": "",
    }


def case_config_row_from_defaults(row_id, case_name, case_data):
    """Return a case-setup row initialized from discovered defaults."""
    defaults = (case_data or {}).get(case_name, {})
    time_range = defaults.get("time_average_range", ["", ""])
    altitude_range = defaults.get("altitude_comparison_range", ["", ""])
    average_time_seconds = int(defaults.get("average_time_seconds") or DEFAULT_AVERAGE_TIME_SECONDS)
    return {
        "id": row_id,
        "case_name": case_name,
        "time_start": time_range[0] if len(time_range) > 0 else "",
        "time_end": time_range[1] if len(time_range) > 1 else "",
        "average_time_seconds": average_time_seconds,
        "altitude_min": altitude_range[0] if len(altitude_range) > 0 else "",
        "altitude_max": altitude_range[1] if len(altitude_range) > 1 else "",
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


def wildcard_no_update(values):
    """Return one Dash ``no_update`` value for each matched ALL component."""
    return [no_update for _ in (values or [])]


def removal_click_is_real(click_values, row_order, row_id):
    """Return whether a pattern-matched remove control was actually clicked.

    Dash can invoke an ``ALL`` input when a freshly hydrated row is inserted.
    Its new remove button has ``n_clicks=0`` and must not be mistaken for a
    request to delete the row that the workspace loader just created.
    """
    try:
        position = list(row_order or []).index(row_id)
        return int((click_values or [])[position] or 0) > 0
    except (IndexError, TypeError, ValueError):
        return False


def parameter_options_by_row(param_values, tunable_names):
    """Return dropdown options that exclude parameters selected in other rows."""
    canonical_values = [
        canonical_tunable_parameter_name(value, tunable_names)
        for value in (param_values or [])
    ]
    selected_by_row = [
        value
        for value in canonical_values
        if value
    ]
    selected_counts = {}
    for name in selected_by_row:
        selected_counts[name] = selected_counts.get(name, 0) + 1

    options_by_row = []
    for current_name in canonical_values:
        row_options = []
        for name in tunable_names or []:
            if selected_counts.get(name, 0) == 0 or name == current_name:
                row_options.append({"label": name, "value": name})
        options_by_row.append(row_options)
    return options_by_row


def sanitize_param_values_for_config(param_values, min_values, max_values, tunable_names, default_ranges):
    """Clear parameter rows that are not available in the selected config."""
    valid_names = set(tunable_names or [])
    updated_params = list(param_values or [])
    updated_min = list(min_values or [])
    updated_max = list(max_values or [])
    while len(updated_min) < len(updated_params):
        updated_min.append("")
    while len(updated_max) < len(updated_params):
        updated_max.append("")

    for idx, raw_name in enumerate(updated_params):
        name = canonical_tunable_parameter_name(raw_name, valid_names)
        if not name:
            continue
        if name not in valid_names:
            updated_params[idx] = None
            updated_min[idx] = ""
            updated_max[idx] = ""
            continue
        updated_params[idx] = name
        derived = (default_ranges or {}).get(name, {})
        if not str(updated_min[idx] or "").strip() and derived.get("min"):
            updated_min[idx] = derived["min"]
        if not str(updated_max[idx] or "").strip() and derived.get("max"):
            updated_max[idx] = derived["max"]

    return updated_params, updated_min, updated_max


def register_settings_callbacks(app):
    """Register callbacks that manage case defaults and range-row state."""

    @app.callback(
        Output("tune-config-buttons", "children"),
        Input("tune-tunable-configs", "data"),
        Input("tune-selected-config", "data"),
    )
    def render_config_buttons(configs, selected_config):
        """Refresh Tune's named-config choices after Run saves a config."""
        return build_config_buttons(configs or [], selected_config)

    @app.callback(
        Output("tune-strategy-mode", "data"),
        Input("tune-mode-random", "n_clicks"),
        Input("tune-mode-resolve", "n_clicks"),
        Input("tune-mode-simann", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_strategy_mode(_random_clicks, _resolve_clicks, _simann_clicks, active_job):
        """Persist the selected tuning strategy mode."""
        if active_job:
            return no_update
        trigger_id = callback_context.triggered_id
        if trigger_id == "tune-mode-random":
            return "random"
        if trigger_id == "tune-mode-resolve":
            return "resolve"
        if trigger_id == "tune-mode-simann":
            return "simann"
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
        Output("tune-time-window-aggregation-scope", "data"),
        Input("tune-aggregation-overall", "n_clicks"),
        Input("tune-aggregation-by-case", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_aggregation_scope(_overall_clicks, _by_case_clicks, active_job):
        """Persist whether quantile weighting pools all windows or each case."""
        if active_job:
            return no_update
        mode_by_button = {
            "tune-aggregation-overall": "overall",
            "tune-aggregation-by-case": "by_case",
        }
        selected = mode_by_button.get(callback_context.triggered_id)
        if selected in TIME_WINDOW_AGGREGATION_SCOPES:
            return selected
        return no_update

    @app.callback(
        Output("tune-selected-config", "data"),
        Output("tune-tunable-names", "data"),
        Output("tune-tunable-default-ranges", "data"),
        Output({"type": "tune-range-param", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-config-button", "name": ALL}, "n_clicks"),
        State("tune-tunable-configs", "data"),
        State("tune-selected-config", "data"),
        State({"type": "tune-range-param", "index": ALL}, "value"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_tunable_config(
        _clicks,
        configs,
        current_config,
        param_values,
        min_values,
        max_values,
        active_job,
    ):
        """Persist the selected tunable config and refresh parameter metadata."""
        if active_job:
            return (
                no_update,
                no_update,
                no_update,
                wildcard_no_update(param_values),
                wildcard_no_update(min_values),
                wildcard_no_update(max_values),
            )

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-config-button":
            return (
                no_update,
                no_update,
                no_update,
                wildcard_no_update(param_values),
                wildcard_no_update(min_values),
                wildcard_no_update(max_values),
            )

        config_name = str(trigger_id.get("name", "")).strip()
        if not config_name or config_name not in tunable_config_names(configs):
            return (
                no_update,
                no_update,
                no_update,
                wildcard_no_update(param_values),
                wildcard_no_update(min_values),
                wildcard_no_update(max_values),
            )
        tunable_names = load_tunable_names(config_name)
        default_ranges = load_tunable_default_ranges(config_name)
        updated_params, updated_min, updated_max = sanitize_param_values_for_config(
            param_values,
            min_values,
            max_values,
            tunable_names,
            default_ranges,
        )
        return config_name, tunable_names, default_ranges, updated_params, updated_min, updated_max

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
            if remove_id not in current_order or not removal_click_is_real(_remove_clicks, current_order, remove_id):
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
            return wildcard_no_update(case_values)
        return case_options_by_row(case_values or [], sorted((case_data or {}).keys()))

    @app.callback(
        Output({"type": "tune-case-time-start", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-time-end", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-average-time", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-altitude-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-case-altitude-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-min", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-max", "index": ALL}, "value"),
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
        altitude_min_values,
        altitude_max_values,
        row_order,
        case_data,
        active_job,
    ):
        """Fill editable case settings from case defaults when a case is selected."""
        if active_job:
            return tuple(wildcard_no_update(case_values) for _ in range(5))

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-case-name":
            return tuple(wildcard_no_update(case_values) for _ in range(5))

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return tuple(wildcard_no_update(case_values) for _ in range(5))
        row_pos = current_order.index(row_id)
        case_name = (case_values or [None])[row_pos] if row_pos < len(case_values or []) else None
        if not case_name or case_name not in (case_data or {}):
            return tuple(wildcard_no_update(case_values) for _ in range(5))

        # A workspace loader supplies a complete saved window in the same
        # render that creates its case dropdown.  Preserve those explicit
        # values if Dash reports the newly inserted dropdown as an input
        # change; defaults are only appropriate for a genuinely blank added
        # row.
        hydrated_values = [
            time_start_values,
            time_end_values,
            average_time_values,
            altitude_min_values,
            altitude_max_values,
        ]
        if all(
            row_pos < len(values or []) and (values or [])[row_pos] not in (None, "")
            for values in hydrated_values
        ):
            return tuple(wildcard_no_update(case_values) for _ in range(5))

        defaults = case_config_row_from_defaults(row_id, case_name, case_data)
        outputs = [
            list(time_start_values or []),
            list(time_end_values or []),
            list(average_time_values or []),
            list(altitude_min_values or []),
            list(altitude_max_values or []),
        ]
        while any(len(values) < len(case_values or []) for values in outputs):
            for values in outputs:
                while len(values) < len(case_values or []):
                    values.append("")

        outputs[0][row_pos] = defaults["time_start"]
        outputs[1][row_pos] = defaults["time_end"]
        outputs[2][row_pos] = defaults["average_time_seconds"]
        outputs[3][row_pos] = defaults["altitude_min"]
        outputs[4][row_pos] = defaults["altitude_max"]
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
            if remove_id not in current_order or not removal_click_is_real(_remove_clicks, current_order, remove_id):
                return no_update, no_update, no_update
            patch = Patch()
            del patch[current_order.index(remove_id)]
            current_order.remove(remove_id)
            return patch, no_update, current_order

        return no_update, no_update, no_update

    @app.callback(
        Output("tune-settings-resolution", "data"),
        Output("tune-settings-resolution-note", "children"),
        Output("tune-settings-resolution-note", "style"),
        Input("tune-selected-config", "data"),
        Input("tune-scm-override", "value"),
    )
    def resolve_tune_controls(selected_config, override):
        """Expose mode restrictions before a Tune worker is started."""
        try:
            resolution = evaluate_tune_settings(selected_config or "default", override or "")
        except (OSError, ValueError) as exc:
            return {}, str(exc), {"display": "block", "marginTop": "5px", "marginBottom": "5px", "color": "#fecaca"}
        issues = [issue for issue in resolution.get("issues", []) if issue.get("severity") == "error"]
        if not issues:
            return resolution, "", {"display": "none"}
        return (
            resolution,
            " ".join(str(issue.get("message") or "") for issue in issues),
            {"display": "block", "marginTop": "5px", "marginBottom": "5px", "color": "#fecaca"},
        )

    @app.callback(
        Output({"type": "tune-range-param", "index": ALL}, "options"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        Input("tune-range-row-order", "data"),
        Input("tune-tunable-names", "data"),
        Input("tune-settings-resolution", "data"),
    )
    def sync_parameter_options(param_values, _row_order, tunable_names, resolution):
        """Hide parameters already selected in other tuning rows."""
        states = dict((resolution or {}).get("parameter_states") or {})
        active_names = [
            name for name in (tunable_names or [])
            if is_independently_tunable(states.get(name))
        ]
        return parameter_options_by_row(param_values or [], active_names)

    @app.callback(
        Output({"type": "tune-range-targets", "index": ALL}, "data"),
        Output({"type": "tune-range-link-label", "index": ALL}, "children"),
        Output({"type": "tune-range-link-label", "index": ALL}, "style"),
        Input({"type": "tune-range-param", "index": ALL}, "value"),
        State({"type": "tune-range-targets", "index": ALL}, "data"),
        State("tune-range-row-order", "data"),
        prevent_initial_call=True,
    )
    def unlink_manual_parameter_change(param_values, target_values, row_order):
        """Keep linked preset targets until the row's selector is changed.

        A range remains a single logical coordinate when its displayed first
        target is untouched.  Selecting a different parameter is an explicit
        manual edit, so that row becomes an ordinary one-target range.
        """
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-range-param":
            return no_update, no_update, no_update
        row_id = trigger_id.get("index")
        order = list(row_order or [])
        if row_id not in order:
            return no_update, no_update, no_update
        row_pos = order.index(row_id)
        values = list(param_values or [])
        existing = [list(item or []) for item in (target_values or [])]
        if row_pos >= len(values):
            return no_update, no_update, no_update
        while len(existing) < len(values):
            existing.append([])
        param = str(values[row_pos] or "").strip()
        current = [str(item).strip() for item in existing[row_pos] if str(item).strip()]
        if not current or current[0] != param:
            existing[row_pos] = [param] if param else []
        labels = [" = ".join(targets) if len(targets) > 1 else "" for targets in existing]
        styles = [
            {"fontSize": "12px", "opacity": 0.8, "whiteSpace": "nowrap", "display": "inline-block" if label else "none"}
            for label in labels
        ]
        return existing, labels, styles

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
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-range-param":
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

        row_id = trigger_id.get("index")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return wildcard_no_update(min_values), wildcard_no_update(max_values)
        row_pos = current_order.index(row_id)

        param_values = list(_param_values or [])
        if row_pos >= len(param_values):
            return wildcard_no_update(min_values), wildcard_no_update(max_values)
        param_name = param_values[row_pos]
        derived_range = (default_ranges or {}).get(param_name)
        if not derived_range:
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

        # As with case windows, do not overwrite an explicit saved range
        # merely because the loader inserted the corresponding dropdown.
        if (
            row_pos < len(min_values or [])
            and row_pos < len(max_values or [])
            and (min_values or [])[row_pos] not in (None, "")
            and (max_values or [])[row_pos] not in (None, "")
        ):
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

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
        State("tune-workspace-selection", "data"),
    )
    def sync_case_fields(case_values, case_data, active_job, workspace_selection):
        """Refresh field options and default field selection when the case changes."""
        # Hydrated saved revisions already supply their exact selected fields.
        # Do not replace them with the case's generic defaults merely because
        # their saved case rows were rebuilt during load.
        selection = dict(workspace_selection or {})
        # A preset is rendered as a fresh ``mode: new`` draft before it owns a
        # durable workspace ID.  Its newly inserted case rows still emit this
        # callback, so treating only an existing ID as hydrated let the
        # generic first-case defaults overwrite its deliberate field list.
        # Keep that protection only while those exact preset cases remain;
        # manually changing cases in a new draft must still refresh fields.
        preset_cases = [str(name).strip() for name in selection.get("preset_case_names", []) if str(name).strip()]
        is_preset_hydration = preset_cases and selected_case_names(case_values) == preset_cases
        if active_job or selection.get("workspace_id") or is_preset_hydration:
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
        Input({"type": "tune-case-altitude-min", "index": ALL}, "value"),
        Input({"type": "tune-case-altitude-max", "index": ALL}, "value"),
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
        Input("tune-time-window-aggregation-scope", "data"),
        Input("tune-aggregation-weight-1", "value"),
        Input("tune-aggregation-weight-2", "value"),
        Input("tune-aggregation-weight-3", "value"),
        Input("tune-aggregation-weight-4", "value"),
        Input("tune-selected-config", "data"),
        Input("tune-random-max-samples", "value"),
        Input("tune-resolve-spacing", "value"),
        State("tune-active-job", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def invalidate_tuning_state(
        _case_add_clicks,
        _case_remove_clicks,
        _case_names,
        _case_time_start,
        _case_time_end,
        _case_average_time,
        _case_altitude_min,
        _case_altitude_max,
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
        _aggregation_scope,
        _aggregation_weight_1,
        _aggregation_weight_2,
        _aggregation_weight_3,
        _aggregation_weight_4,
        _selected_config,
        _random_max_samples,
        _resolve_spacing,
        active_job,
        workspace_selection,
    ):
        """Clear stale results only for an unsaved editable configuration.

        Loading a durable revision intentionally updates the same controls.
        Those hydration updates must preserve its retained results, including
        after a Dash refresh, rather than behaving like a user edit.
        """
        if active_job or dict(workspace_selection or {}).get("workspace_id"):
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update
        return {}, empty_status_payload(), [], [], {}, True, ""
