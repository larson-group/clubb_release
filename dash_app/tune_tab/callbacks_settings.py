"""Callbacks for tuning-tab selection and parameter-range editing."""

from __future__ import annotations

import json
from pathlib import Path

from dash import ALL, Input, Output, Patch, State, callback_context, html, no_update

from .discovery import (
    available_tunable_configs,
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
from .state import REPO_ROOT
from dash_app.run_tab.layout import build_run_config_buttons
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.clubb_settings_validation import is_independently_tunable
from utilities.save_tunable_config import (
    build_tuned_config_overrides,
    normalize_config_name,
    save_tunable_config,
)

from tuner.taylor_metrics import LOSS_MODE_NAMES, TIME_WINDOW_AGGREGATION_SCOPES
from tuner.request import evaluate_tune_settings


def blank_param_range_row(row_id):
    """Return one empty tuning-range row record."""
    return {"id": row_id, "param": "", "targets": [], "min": "", "max": ""}


def tuned_result_params(result):
    """Return the selected tuning coordinates from a status or results row."""
    return dict((result or {}).get("params") or (result or {}).get("selected_params") or {})


def tuned_result_options(results):
    """Build compact selector labels for at most the retained top 16 results."""
    options = []
    for index, result in enumerate(list(results or [])[:16]):
        rank = result.get("rank", index + 1)
        try:
            loss = f"{float(result.get('total_loss')):.6g}"
        except (TypeError, ValueError):
            loss = "--"
        params = tuned_result_params(result)
        preview_values = []
        for name, value in list(params.items())[:4]:
            try:
                rendered = f"{float(value):.6g}"
            except (TypeError, ValueError):
                rendered = str(value)
            preview_values.append(f"{name}={rendered}")
        preview = ", ".join(preview_values)
        if len(params) > 4:
            preview += ", …"
        label = f"Rank {rank} · loss {loss}"
        if preview:
            label += f" · {preview}"
        options.append({"label": label, "value": str(index)})
    return options


def tune_request_for_save(status, *, output_roots=None):
    """Load the immutable request belonging to the displayed Tune result."""
    job_dir = Path(str((status or {}).get("job_dir") or "")).resolve()
    roots = [
        Path(root).resolve()
        for root in (
            output_roots
            or [Path(REPO_ROOT) / "output" / "tuner", Path(REPO_ROOT) / "output_tuner"]
        )
    ]
    if not any(job_dir.is_relative_to(root) for root in roots):
        raise ValueError("The displayed Tune result has no valid saved job directory.")
    request_path = job_dir / "request.json"
    try:
        request = json.loads(request_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("The completed Tune request could not be read.") from exc
    if not isinstance(request, dict) or not request.get("config"):
        raise ValueError("The completed Tune request does not identify its source config.")
    return request


def blank_locked_param_range_row(row_id):
    """Return one empty two-member equality-constrained range."""
    return {"id": row_id, "targets": ["", ""], "min": "", "max": "", "locked": True}


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


def parameter_options_by_row(param_values, tunable_names, parameter_states=None):
    """Return unique options for every visible ordinary or grouped member."""
    parameter_states = dict(parameter_states or {})
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
            state = parameter_states.get(name)
            option = {"label": name, "value": name}
            if state and not is_independently_tunable(state):
                option["disabled"] = True
                option["title"] = str(
                    state.get("reason") or "Unavailable with the current CLUBB settings."
                )
            elif selected_counts.get(name, 0) and name != current_name:
                option["disabled"] = True
                option["title"] = "Already selected by another parameter row or locked group."
            row_options.append(option)
        options_by_row.append(row_options)
    return options_by_row


def parameter_targets_by_row(member_ids, member_values, row_order):
    """Return visible member values grouped in range-row order."""
    grouped = {row_id: [] for row_id in (row_order or [])}
    for component_id, value in zip(member_ids or [], member_values or []):
        if not isinstance(component_id, dict):
            continue
        row_id = component_id.get("row")
        if row_id not in grouped:
            continue
        grouped[row_id].append((int(component_id.get("member", 0)), str(value or "").strip()))
    return [
        [value for _member, value in sorted(grouped.get(row_id, []))]
        for row_id in (row_order or [])
    ]


def parameter_rows_from_controls(member_ids, member_values, row_order, min_values, max_values):
    """Capture the visible parameter editor without a mirrored target Store."""
    targets_by_row = parameter_targets_by_row(member_ids, member_values, row_order)
    rows = []
    for position, row_id in enumerate(row_order or []):
        targets = targets_by_row[position] if position < len(targets_by_row) else [""]
        rows.append(
            {
                "id": row_id,
                "targets": targets or [""],
                "min": (min_values or [""])[position] if position < len(min_values or []) else "",
                "max": (max_values or [""])[position] if position < len(max_values or []) else "",
                "locked": len(targets) > 1,
            }
        )
    return rows


def parameter_group_specs(member_ids, member_values, row_order):
    """Return nonblank logical coordinate labels and physical targets."""
    specs = []
    for targets in parameter_targets_by_row(member_ids, member_values, row_order):
        clean = [target for target in targets if target]
        if clean:
            specs.append({"label": " = ".join(clean), "targets": clean})
    return specs


def sanitize_param_values_for_config(
    param_values,
    min_values,
    max_values,
    tunable_names,
    default_ranges,
    member_ids=None,
    row_order=None,
):
    """Clear parameter rows that are not available in the selected config."""
    valid_names = set(tunable_names or [])
    updated_params = list(param_values or [])
    updated_min = list(min_values or [])
    updated_max = list(max_values or [])
    row_count = len(row_order or []) if member_ids is not None else len(updated_params)
    while len(updated_min) < row_count:
        updated_min.append("")
    while len(updated_max) < row_count:
        updated_max.append("")

    for idx, raw_name in enumerate(param_values or []):
        name = canonical_tunable_parameter_name(raw_name, valid_names)
        if not name:
            continue
        if name not in valid_names:
            updated_params[idx] = None
            continue
        updated_params[idx] = name

    if member_ids is None:
        primary_names = updated_params
    else:
        primary_names = [
            targets[0] if targets else ""
            for targets in parameter_targets_by_row(member_ids, updated_params, row_order)
        ]
    for idx, name in enumerate(primary_names):
        if not name:
            if idx < len(updated_min):
                updated_min[idx] = ""
                updated_max[idx] = ""
            continue
        derived = (default_ranges or {}).get(name, {})
        if not str(updated_min[idx] or "").strip() and derived.get("min"):
            updated_min[idx] = derived["min"]
        if not str(updated_max[idx] or "").strip() and derived.get("max"):
            updated_max[idx] = derived["max"]

    return updated_params, updated_min, updated_max


def register_settings_callbacks(app):
    """Register callbacks that manage case defaults and range-row state."""

    @app.callback(
        Output("tune-config-save-modal", "className"),
        Output("tune-config-save-result", "options"),
        Output("tune-config-save-result", "value"),
        Output("tune-config-save-name", "value"),
        Output("tune-config-save-feedback", "children"),
        Output("tune-config-save-overwrite", "data"),
        Output("tune-config-save-submit", "children"),
        Output("tune-config-save-submit", "className"),
        Output("tune-config-save-status", "children"),
        Output("tune-tunable-configs", "data", allow_duplicate=True),
        Output("run-tunable-configs", "data", allow_duplicate=True),
        Output("run-config-buttons", "children", allow_duplicate=True),
        Input("tune-config-save-open", "n_clicks"),
        Input("tune-config-save-cancel", "n_clicks"),
        Input("tune-config-save-submit", "n_clicks"),
        State("tune-config-save-result", "value"),
        State("tune-config-save-name", "value"),
        State("tune-config-save-overwrite", "data"),
        State("tune-top-results", "data"),
        State("tune-best-results", "data"),
        State("tune-status", "data"),
        State("run-selected-config", "data"),
        prevent_initial_call=True,
    )
    def save_tuned_config(
        open_clicks,
        _cancel_clicks,
        _submit_clicks,
        selected_result,
        name,
        overwrite_name,
        top_results,
        best_results,
        status,
        run_selected_config,
    ):
        """Save one retained Tune result on top of its exact run configuration."""
        hidden = "shared-notecard-overlay tune-config-save-modal--hidden"
        visible = "shared-notecard-overlay"
        default_submit_class = "tune-config-save-submit"
        results = list(top_results or best_results or [])[:16]
        trigger_id = callback_context.triggered_id
        if trigger_id == "tune-config-save-open" and open_clicks:
            options = tuned_result_options(results)
            return (
                visible,
                options,
                options[0]["value"] if options else None,
                "",
                "",
                None,
                "Save config",
                default_submit_class,
                no_update,
                no_update,
                no_update,
                no_update,
            )
        if trigger_id == "tune-config-save-cancel":
            return (
                hidden, no_update, None, "", "", None, "Save config",
                default_submit_class, no_update, no_update, no_update, no_update,
            )
        if trigger_id != "tune-config-save-submit":
            return (no_update,) * 12

        try:
            save_name = normalize_config_name(name)
            result_index = int(selected_result)
            if result_index < 0 or result_index >= len(results):
                raise ValueError("Choose one of the retained tuning results.")
        except (TypeError, ValueError) as exc:
            feedback = html.Div(
                str(exc),
                className="tune-config-save-feedback tune-config-save-feedback--error",
            )
            return (
                visible, no_update, no_update, no_update, feedback, None,
                "Save config", default_submit_class, no_update, no_update, no_update, no_update,
            )

        configs = available_tunable_configs()
        exists = save_name in tunable_config_names(configs)
        if exists and overwrite_name != save_name:
            feedback = html.Div(
                f"A config named {save_name} already exists. Confirm to replace it.",
                className="tune-config-save-feedback tune-config-save-feedback--warning",
            )
            return (
                visible,
                no_update,
                no_update,
                no_update,
                feedback,
                save_name,
                "Overwrite config",
                "tune-config-save-submit tune-config-save-submit--danger",
                no_update,
                no_update,
                no_update,
                no_update,
            )

        chosen = results[result_index]
        try:
            request = tune_request_for_save(status)
            source_config = request["config"]
            grouped_overrides = build_tuned_config_overrides(
                source_config,
                request.get("override", ""),
                tuned_result_params(chosen),
            )
            rank = chosen.get("rank", result_index + 1)
            saved = save_tunable_config(
                source_config,
                grouped_overrides,
                save_name,
                f"Saved from Tune result rank {rank} in {Path(status['job_dir']).name}.",
                force=exists,
            )
            from tuner.request import settings_schema_for_tune_config

            settings_schema_for_tune_config.cache_clear()
        except (KeyError, OSError, TypeError, ValueError) as exc:
            feedback = html.Div(
                f"Could not save config: {exc}",
                className="tune-config-save-feedback tune-config-save-feedback--error",
            )
            return (
                visible, no_update, no_update, no_update, feedback, None,
                "Save config", default_submit_class, no_update, no_update, no_update, no_update,
            )

        configs = available_tunable_configs()
        action = "Overwrote" if saved["overwritten"] else "Saved"
        message = f"{action} {saved['name']} from Tune rank {rank}."
        return (
            hidden,
            no_update,
            None,
            "",
            "",
            None,
            "Save config",
            default_submit_class,
            message,
            configs,
            configs,
            build_run_config_buttons(configs, run_selected_config).children,
        )

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
        Input("tune-mode-adam", "n_clicks"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_strategy_mode(_random_clicks, _resolve_clicks, _simann_clicks, _adam_clicks, active_job):
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
        if trigger_id == "tune-mode-adam":
            return "adam"
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
        Output({"type": "tune-range-member", "row": ALL, "member": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-config-button", "name": ALL}, "n_clicks"),
        State("tune-tunable-configs", "data"),
        State("tune-selected-config", "data"),
        State({"type": "tune-range-member", "row": ALL, "member": ALL}, "id"),
        State({"type": "tune-range-member", "row": ALL, "member": ALL}, "value"),
        State("tune-range-row-order", "data"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def select_tunable_config(
        _clicks,
        configs,
        current_config,
        member_ids,
        param_values,
        row_order,
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
            member_ids,
            row_order,
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
        Input("tune-range-add-locked", "n_clicks"),
        Input({"type": "tune-range-remove", "index": ALL}, "n_clicks"),
        Input({"type": "tune-range-member-add", "index": ALL}, "n_clicks"),
        Input({"type": "tune-range-member-remove", "row": ALL, "member": ALL}, "n_clicks"),
        State("tune-range-next-id", "data"),
        State("tune-range-row-order", "data"),
        State({"type": "tune-range-member", "row": ALL, "member": ALL}, "id"),
        State({"type": "tune-range-member", "row": ALL, "member": ALL}, "value"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-tunable-names", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def sync_range_rows(
        _add_clicks,
        _add_locked_clicks,
        _remove_clicks,
        _member_add_clicks,
        _member_remove_clicks,
        next_id,
        row_order,
        member_ids,
        member_values,
        min_values,
        max_values,
        tunable_names,
        active_job,
    ):
        """Edit visible ordinary rows and equality-constrained groups."""
        if active_job:
            return no_update, no_update, no_update

        current_order = list(row_order or [])
        trigger_id = callback_context.triggered_id

        if trigger_id == "tune-range-add":
            row_id = int(next_id or 0)
            patch = Patch()
            patch.append(build_param_range_row(blank_param_range_row(row_id), tunable_names or []))
            return patch, row_id + 1, current_order + [row_id]

        if trigger_id == "tune-range-add-locked":
            row_id = int(next_id or 0)
            patch = Patch()
            patch.append(build_param_range_row(blank_locked_param_range_row(row_id), tunable_names or []))
            return patch, row_id + 1, current_order + [row_id]

        if isinstance(trigger_id, dict) and trigger_id.get("type") == "tune-range-remove":
            remove_id = trigger_id.get("index")
            if remove_id not in current_order or not removal_click_is_real(_remove_clicks, current_order, remove_id):
                return no_update, no_update, no_update
            patch = Patch()
            del patch[current_order.index(remove_id)]
            current_order.remove(remove_id)
            return patch, no_update, current_order

        if not isinstance(trigger_id, dict):
            return no_update, no_update, no_update
        trigger_value = callback_context.triggered[0].get("value") if callback_context.triggered else 0
        if not trigger_value:
            return no_update, no_update, no_update

        rows = parameter_rows_from_controls(
            member_ids,
            member_values,
            current_order,
            min_values,
            max_values,
        )
        row_id = trigger_id.get("index", trigger_id.get("row"))
        row = next((item for item in rows if item["id"] == row_id), None)
        if row is None:
            return no_update, no_update, no_update

        if trigger_id.get("type") == "tune-range-member-add":
            row["targets"].append("")
            row["locked"] = True
        elif trigger_id.get("type") == "tune-range-member-remove":
            member = int(trigger_id.get("member", -1))
            if member < 0 or member >= len(row["targets"]):
                return no_update, no_update, no_update
            del row["targets"][member]
            if len(row["targets"]) <= 1:
                row["targets"] = row["targets"] or [""]
                row["locked"] = False
        else:
            return no_update, no_update, no_update

        return (
            [build_param_range_row(item, tunable_names or []) for item in rows],
            no_update,
            no_update,
        )

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
        Output({"type": "tune-range-member", "row": ALL, "member": ALL}, "options"),
        Input({"type": "tune-range-member", "row": ALL, "member": ALL}, "value"),
        Input("tune-range-row-order", "data"),
        Input("tune-tunable-names", "data"),
        Input("tune-settings-resolution", "data"),
    )
    def sync_parameter_options(param_values, _row_order, tunable_names, resolution):
        """Show the catalog while disabling parameters inactive in this mode."""
        states = dict((resolution or {}).get("parameter_states") or {})
        return parameter_options_by_row(param_values or [], tunable_names, states)

    @app.callback(
        Output({"type": "tune-range-min", "index": ALL}, "value", allow_duplicate=True),
        Output({"type": "tune-range-max", "index": ALL}, "value", allow_duplicate=True),
        Input({"type": "tune-range-member", "row": ALL, "member": ALL}, "value"),
        State({"type": "tune-range-member", "row": ALL, "member": ALL}, "id"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State("tune-range-row-order", "data"),
        State("tune-tunable-default-ranges", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def autofill_tune_range(_param_values, member_ids, min_values, max_values, row_order, default_ranges, active_job):
        """Fill min/max from the selected parameter's default value."""
        if active_job:
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict) or trigger_id.get("type") != "tune-range-member":
            return wildcard_no_update(min_values), wildcard_no_update(max_values)

        row_id = trigger_id.get("row")
        current_order = list(row_order or [])
        if row_id not in current_order:
            return wildcard_no_update(min_values), wildcard_no_update(max_values)
        row_pos = current_order.index(row_id)

        target_rows = parameter_targets_by_row(member_ids, _param_values, current_order)
        targets = target_rows[row_pos] if row_pos < len(target_rows) else []
        param_name = targets[0] if targets else None
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
        while len(updated_min) < len(current_order):
            updated_min.append("")
        while len(updated_max) < len(current_order):
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
        Input("tune-range-add-locked", "n_clicks"),
        Input({"type": "tune-range-remove", "index": ALL}, "n_clicks"),
        Input({"type": "tune-range-member-add", "index": ALL}, "n_clicks"),
        Input({"type": "tune-range-member-remove", "row": ALL, "member": ALL}, "n_clicks"),
        Input({"type": "tune-range-member", "row": ALL, "member": ALL}, "value"),
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
        Input("tune-adam-max-updates", "value"),
        Input("tune-adam-learning-rate-percent", "value"),
        Input("tune-adam-perturbation-percent", "value"),
        Input("tune-adam-spsa-pairs", "value"),
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
        _add_locked_clicks,
        _remove_clicks,
        _member_add_clicks,
        _member_remove_clicks,
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
        _adam_max_updates,
        _adam_learning_rate_percent,
        _adam_perturbation_percent,
        _adam_spsa_pairs,
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
