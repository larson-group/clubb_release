"""Callbacks for starting, stopping, and polling dashboard tuning runs."""

from __future__ import annotations

import math
import time

from dash import ALL, Input, Output, State, ctx, html, no_update

from .discovery import (
    available_fields_for_cases,
    load_tunable_default_ranges,
    load_tunable_names,
    tunable_config_names,
)
from .layout import build_case_config_row, build_param_range_row
from .runtime import (
    active_tuning_job,
    active_job_exited,
    clear_finished_job,
    empty_status_payload,
    poll_loss_runs,
    read_tuning_results,
    read_tuning_status,
    start_loss_run,
    stop_tuning_job,
)
from dash_app.shared.tunable_configs import canonical_tunable_parameter_name
from utilities.create_case_namelist import normalize_override_string
from utilities.clubb_settings_validation import is_independently_tunable

from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_LOSS_MODE,
    DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    LOSS_MODE_NAMES,
    TIME_WINDOW_AGGREGATION_SCOPES,
    normalize_aggregation_weights,
)
from tuner.presets import apply_preset, get_preset, list_presets
from tuner.request import evaluate_tune_settings


STARTUP_STATUS_GRACE_SECONDS = 15.0


def _case_list(case_names):
    if isinstance(case_names, str):
        case_names = [case_names]
    return [str(case_name).strip() for case_name in (case_names or []) if str(case_name).strip()]


def _integer_seconds(raw_value):
    value = float(raw_value)
    if int(value) != value:
        raise ValueError("time values must be integer seconds")
    return int(value)


def _numeric_value(raw_value, label):
    try:
        value = float(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    return value


def _initializing_status(job_data):
    """Return immediate browser state while a broker-owned worker starts.

    The scheduler writes its own ``initializing`` status shortly after launch.
    Supplying that state here avoids a transient stale ``stopped``/``idle``
    view and makes the client poller active without requiring a page refresh.
    """
    status = empty_status_payload()
    status.update(
        {
            "state": "initializing",
            "job_dir": str((job_data or {}).get("job_dir") or ""),
            "config": (job_data or {}).get("config", "default"),
        }
    )
    return status


def _case_configs_from_rows(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    altitude_min_values,
    altitude_max_values,
):
    configs = []
    for raw_name, raw_t0, raw_t1, raw_average, raw_altitude_min, raw_altitude_max in zip(
        case_names or [],
        time_start_values or [],
        time_end_values or [],
        average_time_values or [],
        altitude_min_values or [],
        altitude_max_values or [],
    ):
        case_name = "" if raw_name is None else str(raw_name).strip()
        if not case_name:
            continue
        time_start = _integer_seconds(raw_t0)
        time_end = _integer_seconds(raw_t1)
        average_time_seconds = _integer_seconds(raw_average)
        altitude_min = _numeric_value(raw_altitude_min, "altitude min")
        altitude_max = _numeric_value(raw_altitude_max, "altitude max")
        if time_end <= time_start:
            raise ValueError("time end must be after time start")
        if average_time_seconds < 1:
            raise ValueError("average_time_seconds must be >= 1")
        if (time_end - time_start) % average_time_seconds != 0:
            raise ValueError("average_time_seconds must divide the time range evenly")
        if altitude_max < altitude_min:
            raise ValueError("altitude max must be >= altitude min")
        num_time_windows = (time_end - time_start) // average_time_seconds
        configs.append(
            {
                "case_name": case_name,
                "time_average_range": [time_start, time_end],
                "average_time_seconds": average_time_seconds,
                "num_time_windows": int(num_time_windows),
                "altitude_comparison_range": [altitude_min, altitude_max],
            }
        )
    return configs


def format_tuning_error(status, active_job):
    """Render a tuning failure with enough context to diagnose without opening JSON."""
    error_message = str((status or {}).get("error_message") or "Unknown tuning error.")
    log_path = (active_job or {}).get("log_path") or (status or {}).get("log_path")
    visible_lines = error_message.splitlines()[:80]
    if len(error_message.splitlines()) > len(visible_lines):
        visible_lines.append("... truncated in UI; see worker.log for the full traceback ...")

    children = [
        html.Div("Tuning job failed.", style={"fontWeight": "700", "marginBottom": "6px"}),
        html.Pre(
            "\n".join(visible_lines),
            style={
                "whiteSpace": "pre-wrap",
                "maxHeight": "360px",
                "overflowY": "auto",
                "padding": "8px",
                "border": "1px solid #fca5a5",
                "borderRadius": "4px",
                "backgroundColor": "#fff1f2",
                "color": "#7f1d1d",
            },
        ),
    ]
    if log_path:
        children.append(
            html.Div(
                f"Worker log: {log_path}",
                style={"marginTop": "6px", "fontFamily": "monospace", "fontSize": "12px"},
            )
        )
    return html.Div(children)


def build_validation_message(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    altitude_min_values,
    altitude_max_values,
    selected_fields,
    param_values,
    min_values,
    max_values,
    batch_size,
    max_workers,
    strategy_mode,
    loss_mode,
    aggregation_mode,
    random_max_samples,
    resolve_spacing,
    simann_max_iters,
    simann_initial_temp,
    simann_final_temp,
    case_data,
    selected_config,
    tunable_configs,
    tunable_names,
    aggregation_scope=DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    aggregation_weights=DEFAULT_AGGREGATION_WEIGHTS,
    range_targets=None,
):
    """Validate the requested tuning configuration and return an error string if needed."""
    selected_config = str(selected_config or "").strip()
    if not selected_config:
        return "Select a tunable config."
    if selected_config not in tunable_config_names(tunable_configs):
        return f"Unknown tunable config: {selected_config}"
    cases = _case_list(case_names)
    if not cases:
        return "Select at least one case."
    if len(set(cases)) != len(cases):
        return "Each tuning case can only be added once."
    for case_name in cases:
        if not (case_data or {}).get(case_name, {}):
            return f"No tuning defaults found for case {case_name}."
    try:
        case_configs = _case_configs_from_rows(
            case_names,
            time_start_values,
            time_end_values,
            average_time_values,
            altitude_min_values,
            altitude_max_values,
        )
    except TypeError:
        return "Each case row needs numeric time, average-time, and altitude values."
    except ValueError as exc:
        message = str(exc)
        if "divide" in message:
            return "Average time must divide the selected time range evenly."
        if "altitude max" in message:
            return "Altitude max must be >= altitude min."
        if "altitude" in message:
            return "Each case row needs numeric altitude min/max values."
        return "Each case row needs integer time and average-time values."
    if len(case_configs) != len(cases):
        return "Each case row needs a case, time range, average time, and altitude range."
    for config in case_configs:
        case_name = config["case_name"]
        time_start, time_end = config["time_average_range"]
        average_time_seconds = config["average_time_seconds"]
        altitude_min, altitude_max = config["altitude_comparison_range"]
        if time_end <= time_start:
            return f"{case_name} requires time end > time start."
        if average_time_seconds < 1:
            return f"{case_name} requires average time >= 1 second."
        if (time_end - time_start) % average_time_seconds != 0:
            return f"{case_name} requires average time to divide the selected time range evenly."
        if altitude_max < altitude_min:
            return f"{case_name} requires altitude max >= altitude min."
    if not selected_fields:
        return "Select at least one CLUBB field to tune."
    for case_name in cases:
        case_fields = set(((case_data or {}).get(case_name, {}) or {}).get("clubb_fields", []))
        missing_fields = [field_name for field_name in selected_fields if field_name not in case_fields]
        if missing_fields:
            return f"Field(s) not available for {case_name}: {', '.join(missing_fields)}"
    if batch_size in (None, ""):
        return "Enter a batch size."
    try:
        batch_size_value = float(batch_size)
    except (TypeError, ValueError):
        return "Batch size must be an integer."
    if int(batch_size_value) != batch_size_value:
        return "Batch size must be an integer."
    batch_size = int(batch_size_value)
    if batch_size < 1:
        return "Batch size must be >= 1."
    if max_workers in (None, ""):
        return "Enter a max worker count."
    try:
        max_workers_value = float(max_workers)
    except (TypeError, ValueError):
        return "Max workers must be an integer."
    if int(max_workers_value) != max_workers_value:
        return "Max workers must be an integer."
    if int(max_workers_value) < 1:
        return "Max workers must be >= 1."
    if strategy_mode not in {"random", "resolve", "simann"}:
        return "Select a tuning mode."
    loss_mode = loss_mode or DEFAULT_LOSS_MODE
    aggregation_mode = aggregation_mode or DEFAULT_AGGREGATION_MODE
    if loss_mode not in LOSS_MODE_NAMES:
        return "Select a loss mode."
    if aggregation_mode not in AGGREGATION_MODE_NAMES:
        return "Select an aggregation mode."
    if aggregation_scope not in TIME_WINDOW_AGGREGATION_SCOPES:
        return "Select an aggregation scope."
    try:
        normalize_aggregation_weights(aggregation_weights)
    except ValueError as exc:
        return str(exc)
    if strategy_mode == "random":
        if random_max_samples in (None, ""):
            return "Enter a random max sample count."
        try:
            max_samples_value = float(random_max_samples)
        except (TypeError, ValueError):
            return "Random max samples must be an integer."
        if int(max_samples_value) != max_samples_value:
            return "Random max samples must be an integer."
        if int(max_samples_value) < 1:
            return "Random max samples must be >= 1."
    if strategy_mode == "resolve":
        if resolve_spacing in (None, ""):
            return "Enter a resolve spacing."
        try:
            spacing_value = float(resolve_spacing)
        except (TypeError, ValueError):
            return "Resolve spacing must be numeric."
        if spacing_value <= 0.0:
            return "Resolve spacing must be > 0."
    if strategy_mode == "simann":
        if simann_max_iters in (None, ""):
            return "Enter a SimAnn max iteration count."
        try:
            max_iters_value = float(simann_max_iters)
        except (TypeError, ValueError):
            return "SimAnn max iterations must be an integer."
        if int(max_iters_value) != max_iters_value:
            return "SimAnn max iterations must be an integer."
        if int(max_iters_value) < 1:
            return "SimAnn max iterations must be >= 1."
        for label, raw_value in (
            ("initial temperature", simann_initial_temp),
            ("final temperature", simann_final_temp),
        ):
            if raw_value in (None, ""):
                return f"Enter a SimAnn {label}."
            try:
                value = float(raw_value)
            except (TypeError, ValueError):
                return f"SimAnn {label} must be numeric."
            if value <= 0.0:
                return f"SimAnn {label} must be > 0."

    ranges = []
    seen_names = set()
    owned_targets = set()
    for index, (param_name, min_text, max_text) in enumerate(zip(param_values or [], min_values or [], max_values or [])):
        param_name = canonical_tunable_parameter_name(param_name, tunable_names)
        min_text = "" if min_text is None else str(min_text).strip()
        max_text = "" if max_text is None else str(max_text).strip()
        if not any((param_name, min_text, max_text)):
            continue
        if not all((param_name, min_text, max_text)):
            return "Each parameter row must include parameter, min, and max."
        if param_name in seen_names:
            return f"Duplicate tuning parameter: {param_name}"
        if param_name not in set(tunable_names or []):
            return f"{param_name} is not available in config {selected_config}."
        seen_names.add(param_name)
        raw_targets = (range_targets or [])[index] if index < len(range_targets or []) else [param_name]
        targets = [str(target).strip() for target in raw_targets or [] if str(target).strip()]
        if not targets or targets[0] != param_name:
            targets = [param_name]
        duplicate_targets = owned_targets.intersection(targets)
        if duplicate_targets:
            return "A physical tuning parameter may only belong to one range: " + ", ".join(sorted(duplicate_targets))
        owned_targets.update(targets)
        missing_targets = [target for target in targets if target not in set(tunable_names or [])]
        if missing_targets:
            return f"Linked parameter(s) not available in config {selected_config}: {', '.join(missing_targets)}"
        try:
            min_value = float(min_text.replace("D", "E").replace("d", "e"))
            max_value = float(max_text.replace("D", "E").replace("d", "e"))
        except ValueError:
            return f"Invalid numeric range for {param_name}."
        if min_value > max_value:
            return f"{param_name} requires min <= max."
        ranges.append({"name": param_name, "targets": targets, "min": min_value, "max": max_value})

    if not ranges:
        return "Add at least one parameter range to tune."
    return ""


def build_request_payload(
    case_names,
    time_start_values,
    time_end_values,
    average_time_values,
    altitude_min_values,
    altitude_max_values,
    selected_fields,
    param_values,
    min_values,
    max_values,
    batch_size,
    max_workers,
    strategy_mode,
    loss_mode,
    aggregation_mode,
    random_max_samples,
    resolve_spacing,
    simann_max_iters,
    simann_initial_temp,
    simann_final_temp,
    selected_config,
    scm_override,
    aggregation_scope=DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    aggregation_weights=DEFAULT_AGGREGATION_WEIGHTS,
    range_targets=None,
):
    """Build the worker request payload from the live UI values."""
    batch_size = int(float(batch_size))
    max_workers = int(float(max_workers))
    selected_config = str(selected_config or "default").strip() or "default"
    tunable_names = load_tunable_names(selected_config)
    ranges = []
    for index, (param_name, min_text, max_text) in enumerate(zip(param_values or [], min_values or [], max_values or [])):
        param_name = canonical_tunable_parameter_name(param_name, tunable_names)
        min_text = "" if min_text is None else str(min_text).strip()
        max_text = "" if max_text is None else str(max_text).strip()
        if not all((param_name, min_text, max_text)):
            continue
        raw_targets = (range_targets or [])[index] if index < len(range_targets or []) else [param_name]
        targets = [str(target).strip() for target in raw_targets or [] if str(target).strip()]
        if not targets or targets[0] != param_name:
            targets = [param_name]
        ranges.append(
            {
                "name": param_name,
                "targets": targets,
                "min": float(min_text.replace("D", "E").replace("d", "e")),
                "max": float(max_text.replace("D", "E").replace("d", "e")),
            }
        )
    strategy_options = {}
    if strategy_mode == "random":
        strategy_options["max_samples"] = int(float(random_max_samples))
    if strategy_mode == "resolve":
        strategy_options["spacing"] = float(resolve_spacing)
    if strategy_mode == "simann":
        strategy_options["max_iters"] = int(float(simann_max_iters))
        strategy_options["initial_temp"] = float(simann_initial_temp)
        strategy_options["max_final_temp"] = float(simann_final_temp)
        strategy_options["chain_count"] = max(1, max_workers * batch_size)

    case_configs = _case_configs_from_rows(
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        altitude_min_values,
        altitude_max_values,
    )

    request = {
        "config": selected_config,
        "override": normalize_override_string(scm_override),
        "cases": [config["case_name"] for config in case_configs],
        "case_configs": case_configs,
        "selected_fields": list(selected_fields or []),
        "batch_size": int(batch_size),
        "max_workers": int(max_workers),
        "loss_mode": loss_mode or DEFAULT_LOSS_MODE,
        "aggregation_mode": aggregation_mode or DEFAULT_AGGREGATION_MODE,
        "time_window_aggregation_mode": aggregation_mode or DEFAULT_AGGREGATION_MODE,
        "time_window_aggregation_scope": aggregation_scope or DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
        "aggregation_weights": normalize_aggregation_weights(aggregation_weights),
        "parameter_ranges": ranges,
        "strategy": {
            "name": strategy_mode,
            "options": strategy_options,
        },
    }
    return request


def agent_request_to_tune_controls(request, case_data):
    """Translate a canonical agent request into the native Tune form controls.

    The worker request is the source of truth.  Rendering it through the same
    row builders used by the manual UI keeps the displayed ranges and case
    windows editable once a running job has stopped.
    """
    request = dict(request or {})
    selected_config = str(request.get("config") or "default").strip() or "default"
    resolution = evaluate_tune_settings(selected_config, normalize_override_string(request.get("override")))
    tunable_names = [
        name for name in load_tunable_names(selected_config)
        if is_independently_tunable(resolution.get("parameter_states", {}).get(name))
    ]
    default_ranges = load_tunable_default_ranges(selected_config)
    raw_case_configs = list(request.get("case_configs") or [])
    case_rows = []
    for index, raw in enumerate(raw_case_configs):
        raw = dict(raw or {})
        case_name = str(raw.get("case_name") or "").strip()
        discovered_defaults = dict((case_data or {}).get(case_name) or {})
        time_range = list(raw.get("time_average_range") or discovered_defaults.get("time_average_range") or ["", ""])
        altitude_range = list(raw.get("altitude_comparison_range") or discovered_defaults.get("altitude_comparison_range") or ["", ""])
        start = time_range[0] if len(time_range) > 0 else ""
        end = time_range[1] if len(time_range) > 1 else ""
        average = raw.get("average_time_seconds", discovered_defaults.get("average_time_seconds"))
        if average in (None, ""):
            try:
                window_count = max(1, int(raw.get("num_time_windows") or 1))
                average = (int(end) - int(start)) // window_count
            except (TypeError, ValueError):
                average = ""
        case_rows.append(
            {
                "id": index,
                "case_name": case_name,
                "time_start": start,
                "time_end": end,
                "average_time_seconds": average,
                "altitude_min": altitude_range[0] if len(altitude_range) > 0 else "",
                "altitude_max": altitude_range[1] if len(altitude_range) > 1 else "",
            }
        )

    parameter_rows = []
    for index, raw in enumerate(list(request.get("parameter_ranges") or [])):
        raw = dict(raw or {})
        name = canonical_tunable_parameter_name(raw.get("name"), tunable_names)
        parameter_rows.append(
            {
                "id": index,
                "param": name,
                "targets": raw.get("targets") or [name],
                "min": str(raw.get("min") if raw.get("min") is not None else ""),
                "max": str(raw.get("max") if raw.get("max") is not None else ""),
            }
        )

    selected_cases = [row["case_name"] for row in case_rows if row["case_name"]]
    fields = [str(value).strip() for value in (request.get("selected_fields") or []) if str(value).strip()]
    field_names = available_fields_for_cases(selected_cases, case_data or {})
    strategy = dict(request.get("strategy") or {})
    strategy_name = str(strategy.get("name") or "random").strip().lower()
    options = dict(strategy.get("options") or {})
    return {
        "case_rows": case_rows,
        "parameter_rows": parameter_rows,
        "selected_config": selected_config,
        "tunable_names": tunable_names,
        "default_ranges": default_ranges,
        "field_options": [{"label": name, "value": name} for name in field_names],
        "fields": [name for name in fields if name in field_names],
        "strategy": strategy_name,
        "random_max_samples": options.get("max_samples", 100),
        "resolve_spacing": options.get("spacing", 0.1),
        "simann_max_iters": options.get("max_iters", 200),
        "simann_initial_temp": options.get("initial_temp", 1.0),
        "simann_final_temp": options.get("max_final_temp", 1.0e-12),
        "batch_size": request.get("batch_size", 1),
        "max_workers": request.get("max_workers", 1),
        "loss_mode": request.get("loss_mode") or DEFAULT_LOSS_MODE,
        "aggregation_mode": request.get("aggregation_mode") or DEFAULT_AGGREGATION_MODE,
        "aggregation_scope": request.get("time_window_aggregation_scope") or DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
        "aggregation_weights": request.get("aggregation_weights") or list(DEFAULT_AGGREGATION_WEIGHTS),
        "override": request.get("override") or "",
    }


def initial_tune_control_outputs():
    """Return the ordinary, editable Tune defaults for a new workspace."""
    # Imported lazily to avoid the module-registration cycle: tab.py imports
    # this callback module while building the initial Dash layout.
    from .tab import build_initial_tune_state

    initial = build_initial_tune_state()
    case_rows = list(initial["initial_case_rows"])
    param_rows = list(initial["initial_param_rows"])
    case_options = [{"label": name, "value": name} for name in initial["cases"]]
    case_children = [build_case_config_row(row, case_options) for row in case_rows]
    range_children = [build_param_range_row(row, initial["tunable_names"]) for row in param_rows]
    return (
        case_children,
        len(case_rows),
        [row["id"] for row in case_rows],
        range_children,
        len(param_rows),
        [row["id"] for row in param_rows],
        initial["selected_config"],
        initial["tunable_names"],
        initial["tunable_default_ranges"],
        [{"label": field, "value": field} for field in initial["field_options"]],
        initial["selected_fields"],
        initial["strategy_mode"],
        initial["loss_mode"],
        initial["aggregation_mode"],
        initial["time_window_aggregation_scope"],
        *initial["aggregation_weights"],
        initial["random_max_samples"],
        initial["resolve_spacing"],
        initial["simann_max_iters"],
        initial["simann_initial_temp"],
        initial["simann_final_temp"],
        initial["batch_size"],
        initial["max_workers"],
        initial["scm_override"],
    )


def register_run_callbacks(app):
    """Register tuning start/stop/poll callbacks."""

    @app.callback(
        Output("tune-case-rows", "children", allow_duplicate=True),
        Output("tune-case-next-id", "data", allow_duplicate=True),
        Output("tune-case-row-order", "data", allow_duplicate=True),
        Output("tune-range-rows", "children", allow_duplicate=True),
        Output("tune-range-next-id", "data", allow_duplicate=True),
        Output("tune-range-row-order", "data", allow_duplicate=True),
        Output("tune-selected-config", "data", allow_duplicate=True),
        Output("tune-tunable-names", "data", allow_duplicate=True),
        Output("tune-tunable-default-ranges", "data", allow_duplicate=True),
        Output("tune-field-selector", "options", allow_duplicate=True),
        Output("tune-field-selector", "value", allow_duplicate=True),
        Output("tune-strategy-mode", "data", allow_duplicate=True),
        Output("tune-loss-mode", "data", allow_duplicate=True),
        Output("tune-aggregation-mode", "data", allow_duplicate=True),
        Output("tune-time-window-aggregation-scope", "data", allow_duplicate=True),
        Output("tune-aggregation-weight-1", "value", allow_duplicate=True),
        Output("tune-aggregation-weight-2", "value", allow_duplicate=True),
        Output("tune-aggregation-weight-3", "value", allow_duplicate=True),
        Output("tune-aggregation-weight-4", "value", allow_duplicate=True),
        Output("tune-random-max-samples", "value", allow_duplicate=True),
        Output("tune-resolve-spacing", "value", allow_duplicate=True),
        Output("tune-simann-max-iters", "value", allow_duplicate=True),
        Output("tune-simann-initial-temp", "value", allow_duplicate=True),
        Output("tune-simann-final-temp", "value", allow_duplicate=True),
        Output("tune-batch-size", "value", allow_duplicate=True),
        Output("tune-max-workers", "value", allow_duplicate=True),
        Output("tune-scm-override", "value", allow_duplicate=True),
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-interval", "n_intervals", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-workspace-selection", "data"),
        State("tune-case-data", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def load_selected_workspace(workspace_selection, case_data, current_active_job):
        """Load only an explicitly selected Tune revision into the dashboard.

        Workspace metadata changes (for example a rename) deliberately do not
        reach this callback.  A selected revision is immutable, and merely
        refreshing the workspace browser must never rehydrate its controls or
        replace a user's in-progress configuration.  Persisted selections are
        reconciled against the browser list by the dedicated workspace
        callback, which changes ``tune-workspace-selection`` only when the
        selected lineage genuinely changes.
        """
        selection = dict(workspace_selection or {})
        workspace_id = str(selection.get("workspace_id") or "")
        revision_id = str(selection.get("revision_id") or "")
        if selection.get("mode") == "new":
            preset_name = str(selection.get("preset") or "").strip()
            if preset_name:
                try:
                    controls = agent_request_to_tune_controls(apply_preset({"preset": preset_name}), case_data or {})
                except (OSError, TypeError, ValueError) as exc:
                    return (no_update,) * 34 + (f"Could not apply Tune preset: {exc}",)
                case_options = [{"label": name, "value": name} for name in sorted((case_data or {}).keys())]
                case_children = [build_case_config_row(row, case_options) for row in controls["case_rows"]]
                range_children = [build_param_range_row(row, controls["tunable_names"]) for row in controls["parameter_rows"]]
                return (
                    case_children, len(controls["case_rows"]), [row["id"] for row in controls["case_rows"]],
                    range_children, len(controls["parameter_rows"]), [row["id"] for row in controls["parameter_rows"]],
                    controls["selected_config"], controls["tunable_names"], controls["default_ranges"],
                    controls["field_options"], controls["fields"], controls["strategy"], controls["loss_mode"],
                    controls["aggregation_mode"], controls["aggregation_scope"], *controls["aggregation_weights"],
                    controls["random_max_samples"], controls["resolve_spacing"],
                    controls["simann_max_iters"], controls["simann_initial_temp"], controls["simann_final_temp"],
                    # Presets define the scientific experiment, not local
                    # scheduling.  Preserve the user's batch/worker choices.
                    no_update, no_update, controls["override"], {}, empty_status_payload(),
                    [], [], {}, True, 0, f"Applied preset {preset_name}. Review this editable draft, then Start to save and run it.",
                )
            # A new workspace has no durable request yet.  Clear any prior
            # selected execution's retained plots/status, but leave the normal
            # editable control defaults in place for the user to configure.
            return (
                *initial_tune_control_outputs(),
                {},
                empty_status_payload(),
                [],
                [],
                {},
                True,
                0,
                "Configure this new workspace, then Start to create its original revision.",
            )
        if not workspace_id or not revision_id:
            return (no_update,) * 35
        try:
            from dash_app.shared.broker_client import perform_action

            loaded = perform_action(
                "load_tuning_workspace",
                {"workspace_id": workspace_id, "revision_id": revision_id},
                internal=True,
            )
        except Exception as exc:
            return (no_update,) * 34 + (f"Could not load saved Tune revision: {exc}",)
        request = dict(loaded.get("request") or {})
        job = dict(loaded.get("job") or {})
        if not job or not request:
            return (no_update,) * 34 + ("Saved Tune revision has no readable request.",)
        job.update({"workspace_id": workspace_id, "revision_id": revision_id})
        execution_state = str((loaded.get("execution") or {}).get("state") or "draft")
        message_prefix = (
            f"Loaded {workspace_id} / {revision_id} ({execution_state})."
            if selection.get("mode") != "draft"
            else f"Editing new revision {workspace_id} / {revision_id}."
        )
        job["broker_managed"] = True
        try:
            controls = agent_request_to_tune_controls(request, case_data or {})
        except (OSError, TypeError, ValueError) as exc:
            return (no_update,) * 34 + (f"Could not display agent tuning request: {exc}",)
        case_names = [row["case_name"] for row in controls["case_rows"]]
        case_options = [
            {"label": name, "value": name}
            for name in sorted((case_data or {}).keys())
        ]
        case_children = [build_case_config_row(row, case_options) for row in controls["case_rows"]]
        range_children = [build_param_range_row(row, controls["tunable_names"]) for row in controls["parameter_rows"]]
        status = read_tuning_status(job.get("status_path"))
        retained_results = read_tuning_results(job.get("results_path"))
        current_job = dict(current_active_job or {})
        same_live_job = (
            bool(current_job)
            and str(current_job.get("job_dir") or "") == str(job.get("job_dir") or "")
        )
        # Immediately after Start, the worker can still have its initial
        # ``draft`` status.  Do not erase the live dashboard job in that
        # narrow interval, or Stop/polling become unavailable before the
        # worker writes ``running``.
        display_job = (
            current_job
            if same_live_job
            else job
            if str(status.get("state") or "") in {"initializing", "running", "stopping"}
            else {}
        )
        if str((workspace_selection or {}).get("mode") or "") == "draft":
            message = message_prefix + " This revision is editable until it is started."
        else:
            message = message_prefix + " The controls show the exact worker request; create a new revision to change settings."
        return (
            case_children,
            len(case_names),
            list(range(len(case_names))),
            range_children,
            len(controls["parameter_rows"]),
            list(range(len(controls["parameter_rows"]))),
            controls["selected_config"],
            controls["tunable_names"],
            controls["default_ranges"],
            controls["field_options"],
            controls["fields"],
            controls["strategy"],
            controls["loss_mode"],
            controls["aggregation_mode"],
            controls["aggregation_scope"],
            *controls["aggregation_weights"],
            controls["random_max_samples"],
            controls["resolve_spacing"],
            controls["simann_max_iters"],
            controls["simann_initial_temp"],
            controls["simann_final_temp"],
            controls["batch_size"],
            controls["max_workers"],
            controls["override"],
            display_job,
            status,
            status.get("top_results", []),
            retained_results.get("best_results", []),
            retained_results.get("best_results_by_case", {}),
            False,
            0,
            message,
        )

    @app.callback(
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-interval", "n_intervals"),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Output("tune-right-tabs", "value", allow_duplicate=True),
        Input("tune-start-button", "n_clicks"),
        State({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-min", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-max", "index": ALL}, "value"),
        State("tune-field-selector", "value"),
        State({"type": "tune-range-param", "index": ALL}, "value"),
        State({"type": "tune-range-min", "index": ALL}, "value"),
        State({"type": "tune-range-max", "index": ALL}, "value"),
        State({"type": "tune-range-targets", "index": ALL}, "data"),
        State("tune-batch-size", "value"),
        State("tune-max-workers", "value"),
        State("tune-strategy-mode", "data"),
        State("tune-loss-mode", "data"),
        State("tune-aggregation-mode", "data"),
        State("tune-time-window-aggregation-scope", "data"),
        State("tune-aggregation-weight-1", "value"),
        State("tune-aggregation-weight-2", "value"),
        State("tune-aggregation-weight-3", "value"),
        State("tune-aggregation-weight-4", "value"),
        State("tune-random-max-samples", "value"),
        State("tune-resolve-spacing", "value"),
        State("tune-simann-max-iters", "value"),
        State("tune-simann-initial-temp", "value"),
        State("tune-simann-final-temp", "value"),
        State("tune-case-data", "data"),
        State("tune-selected-config", "data"),
        State("tune-tunable-configs", "data"),
        State("tune-tunable-names", "data"),
        State("tune-scm-override", "value"),
        State("tune-active-job", "data"),
        State("tune-workspace-selection", "data"),
        State("tune-status", "data"),
        prevent_initial_call=True,
    )
    def start_tuning(
        _n_clicks,
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        altitude_min_values,
        altitude_max_values,
        selected_fields,
        param_values,
        min_values,
        max_values,
        range_targets,
        batch_size,
        max_workers,
        strategy_mode,
        loss_mode,
        aggregation_mode,
        aggregation_scope,
        aggregation_weight_1,
        aggregation_weight_2,
        aggregation_weight_3,
        aggregation_weight_4,
        random_max_samples,
        resolve_spacing,
        simann_max_iters,
        simann_initial_temp,
        simann_final_temp,
        case_data,
        selected_config,
        tunable_configs,
        tunable_names,
        scm_override,
        active_job,
        workspace_selection,
        displayed_status,
    ):
        """Validate the tuning inputs, then launch the background worker."""
        active_status = read_tuning_status((active_job or {}).get("status_path")) if active_job else dict(displayed_status or {})
        if str(active_status.get("state") or "") in {"initializing", "running", "stopping"}:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, "A tuning job is already active.", no_update

        selection = dict(workspace_selection or {})
        workspace_id = str(selection.get("workspace_id") or "")
        revision_id = str(selection.get("revision_id") or "")
        selection_mode = str(selection.get("mode") or "")
        if selection_mode == "readonly" and workspace_id and revision_id:
            if str(active_status.get("state") or "") != "stopped":
                return no_update, no_update, no_update, no_update, no_update, no_update, no_update, "Started revisions are immutable; use Reset or New revision.", no_update
            try:
                from dash_app.shared.broker_client import perform_action

                result = perform_action(
                    "resume_tuning_revision",
                    {"workspace_id": workspace_id, "revision_id": revision_id},
                    internal=True,
                )
                resumed = dict(result.get("job") or {})
                if not resumed:
                    raise RuntimeError("dashboard broker did not return resumed tuning metadata")
                resumed["broker_managed"] = True
                resumed["browser_started_at_unix"] = time.time()
                return resumed, _initializing_status(resumed), [], [], {}, False, 1, "", "runtime"
            except Exception as exc:
                return no_update, no_update, no_update, no_update, no_update, no_update, no_update, str(exc), no_update

        validation_message = build_validation_message(
            case_names,
            time_start_values,
            time_end_values,
            average_time_values,
            altitude_min_values,
            altitude_max_values,
            selected_fields,
            param_values,
            min_values,
            max_values,
            batch_size,
            max_workers,
            strategy_mode,
            loss_mode,
            aggregation_mode,
            random_max_samples,
            resolve_spacing,
            simann_max_iters,
            simann_initial_temp,
            simann_final_temp,
            case_data,
            selected_config,
            tunable_configs,
            tunable_names,
            aggregation_scope=aggregation_scope,
            aggregation_weights=[aggregation_weight_1, aggregation_weight_2, aggregation_weight_3, aggregation_weight_4],
            range_targets=range_targets,
        )
        if validation_message:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, validation_message, no_update

        request = build_request_payload(
            case_names,
            time_start_values,
            time_end_values,
            average_time_values,
            altitude_min_values,
            altitude_max_values,
            selected_fields,
            param_values,
            min_values,
            max_values,
            batch_size,
            max_workers,
            strategy_mode,
            loss_mode,
            aggregation_mode,
            random_max_samples,
            resolve_spacing,
            simann_max_iters,
            simann_initial_temp,
            simann_final_temp,
            selected_config,
            scm_override,
            aggregation_scope=aggregation_scope,
            aggregation_weights=[aggregation_weight_1, aggregation_weight_2, aggregation_weight_3, aggregation_weight_4],
            range_targets=range_targets,
        )
        try:
            # The worker must be born in the durable broker, even when a user
            # pressed Tune's ordinary Start button. A Dash debug reload then
            # loses only this browser callback, never the worker/keepalive.
            from dash_app.shared.broker_client import perform_action

            if selection_mode == "draft" and workspace_id and revision_id:
                request["workspace_id"] = workspace_id
                request["revision_id"] = revision_id
            elif selection_mode == "new":
                # The name is UI metadata, deliberately kept out of the
                # scientific request that is written to request.json.
                request["workspace_display_name"] = str(selection.get("display_name") or "").strip()
                if selection.get("preset"):
                    request["preset"] = str(selection["preset"])
            result = perform_action("launch_tuning_request", {"request": request}, internal=True)
            active_job = dict(result.get("job") or {})
            if not active_job:
                raise RuntimeError("dashboard broker did not return tuning job metadata")
            active_job["broker_managed"] = True
            active_job["browser_started_at_unix"] = time.time()
            # A fresh workspace has no durable ID until the broker creates
            # it. Preserve the browser-only draft token so the workspace
            # lifecycle callback can seal only this exact new draft.
            if selection_mode == "new":
                active_job["draft_token"] = str(selection.get("draft_token") or "")
        except Exception as exc:
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, str(exc), no_update

        return active_job, _initializing_status(active_job), [], [], {}, False, 1, "", "runtime"

    @app.callback(
        Output("tune-workspace-selection", "data", allow_duplicate=True),
        Input({"type": "tune-preset", "name": ALL}, "n_clicks"),
        State("tune-active-job", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def apply_tune_preset(_clicks, active_job, selection):
        """Open a preset as a fresh, editable workspace draft."""
        if active_job:
            return no_update
        trigger = ctx.triggered_id
        if not isinstance(trigger, dict) or trigger.get("type") != "tune-preset":
            return no_update
        name = str(trigger.get("name") or "").strip()
        try:
            get_preset(name)
        except ValueError:
            return no_update
        return {
            "mode": "new",
            "preset": name,
            "preset_case_names": [item["case_name"] for item in get_preset(name)["case_configs"]],
            "display_name": f"{name}-{int(time.time())}",
            "draft_token": str(time.time_ns()),
        }

    @app.callback(
        Output({"type": "tune-preset", "name": ALL}, "disabled"),
        Input("tune-active-job", "data"),
        Input("tune-workspace-selection", "data"),
    )
    def disable_tune_presets(active_job, selection):
        selection = dict(selection or {})
        locked = bool(active_job) or str(selection.get("mode") or "") == "readonly"
        return [locked] * len(list_presets())

    @app.callback(
        Output("tune-interval", "disabled", allow_duplicate=True),
        Input("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def sync_tuning_poll_interval(active_job):
        """Keep polling coupled to the client-visible broker job handoff.

        A resumed revision is launched outside Dash and may update its status
        file after this callback returns.  Driving the interval from the
        active-job Store guarantees that the browser starts polling as soon as
        Continue publishes that job, rather than depending on timing among the
        Start callback's other multi-output updates.
        """
        return not bool(active_job)

    @app.callback(
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-stop-button", "n_clicks"),
        State("tune-active-job", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def stop_tuning(_n_clicks, active_job, workspace_selection):
        """Request a graceful stop and keep polling until the tuner exits."""
        job = dict(active_job or {})
        if not job:
            selection = dict(workspace_selection or {})
            workspace_id = str(selection.get("workspace_id") or "")
            revision_id = str(selection.get("revision_id") or "")
            if workspace_id and revision_id:
                try:
                    from dash_app.shared.broker_client import perform_action

                    loaded = perform_action(
                        "load_tuning_workspace",
                        {"workspace_id": workspace_id, "revision_id": revision_id},
                        internal=True,
                    )
                    job = dict(loaded.get("job") or {})
                except Exception as exc:
                    return no_update, no_update, no_update, no_update, f"Could not stop Tune revision: {exc}"
        if not job:
            return no_update, no_update, no_update, no_update, "No active Tune revision is loaded."
        status = read_tuning_status(job.get("status_path"))
        if str(status.get("state") or "") not in {"initializing", "running", "stopping"}:
            return no_update, no_update, no_update, no_update, "The selected Tune revision is no longer running."
        stop_tuning_job(job)
        return no_update, no_update, no_update, no_update, ""

    @app.callback(
        Output("tune-status", "data", allow_duplicate=True),
        Output("tune-top-results", "data", allow_duplicate=True),
        Output("tune-best-results", "data", allow_duplicate=True),
        Output("tune-best-results-by-case", "data", allow_duplicate=True),
        Output("tune-active-job", "data", allow_duplicate=True),
        Output("tune-interval", "disabled", allow_duplicate=True),
        Output("tune-runtime-error", "children", allow_duplicate=True),
        Input("tune-interval", "n_intervals"),
        State("tune-active-job", "data"),
        State("tune-status", "data"),
        prevent_initial_call=True,
    )
    def poll_tuning_status(_tick, active_job, displayed_status):
        """Poll the worker status file and retire finished jobs."""
        if not active_job:
            return no_update, no_update, no_update, no_update, no_update, True, no_update

        # A broker-managed job is leased by its detached broker.  Keeping the
        # old in-process heartbeat here would make two controllers update the
        # same control file after a Dash reload.
        if not active_job.get("broker_managed"):
            active_tuning_job(active_job).heartbeat()
        status = read_tuning_status(active_job.get("status_path"))
        top_results = status.get("top_results", [])
        retained_results = read_tuning_results(active_job.get("results_path"))
        best_results = retained_results.get("best_results", [])
        best_results_by_case = retained_results.get("best_results_by_case", {})
        state = status.get("state", "idle")
        browser_started_at = float(active_job.get("browser_started_at_unix", 0.0) or 0.0)
        startup_age = time.time() - browser_started_at if browser_started_at else float("inf")
        waiting_for_initial_status = (
            state == "stopped"
            and str((displayed_status or {}).get("state") or "") == "initializing"
            and startup_age < STARTUP_STATUS_GRACE_SECONDS
        )
        if waiting_for_initial_status:
            # ``resume`` starts a fresh process against a durable directory.
            # Its old stopped status is valid until the child imports and
            # writes initializing/running, so do not discard the live browser
            # handoff during that short window.
            return displayed_status, [], [], {}, active_job, False, ""
        validation_message = ""
        active_job_out = active_job
        interval_disabled = False

        exited = not active_job.get("broker_managed") and active_job_exited(active_job)
        if state in {"stopped", "error", "finished"} or exited:
            clear_finished_job(active_job)
            active_job_out = {}
            interval_disabled = True
            if state == "error":
                validation_message = format_tuning_error(status, active_job)

        return status, top_results, best_results, best_results_by_case, active_job_out, interval_disabled, validation_message

    @app.callback(
        Output("tune-loss-runs", "data", allow_duplicate=True),
        Output("tune-loss-run-interval", "disabled", allow_duplicate=True),
        Output("tune-loss-run-message", "children"),
        Input({"type": "tune-loss-run-button", "action": ALL}, "n_clicks_timestamp"),
        State("tune-top-results", "data"),
        State("tune-best-results", "data"),
        State({"type": "tune-case-name", "index": ALL}, "value"),
        State({"type": "tune-case-time-start", "index": ALL}, "value"),
        State({"type": "tune-case-time-end", "index": ALL}, "value"),
        State({"type": "tune-case-average-time", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-min", "index": ALL}, "value"),
        State({"type": "tune-case-altitude-max", "index": ALL}, "value"),
        State("tune-field-selector", "value"),
        State("tune-selected-config", "data"),
        State("tune-batch-size", "value"),
        State("tune-scm-override", "value"),
        State("tune-loss-runs", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def run_result_loss(
        _button_timestamps,
        top_results,
        best_results,
        case_names,
        time_start_values,
        time_end_values,
        average_time_values,
        altitude_min_values,
        altitude_max_values,
        selected_fields,
        selected_config,
        batch_size,
        scm_override,
        loss_runs,
        workspace_selection,
    ):
        """Run listed results through the window loss runner or normal SCM runner."""
        triggered_id = ctx.triggered_id
        if not isinstance(triggered_id, dict):
            return no_update, no_update, no_update
        trigger_value = ctx.triggered[0].get("value") if ctx.triggered else None
        if isinstance(trigger_value, list):
            values = trigger_value
        else:
            values = [trigger_value]
        try:
            click_timestamp = max(float(value) for value in values if value is not None)
        except (TypeError, ValueError):
            click_timestamp = -1.0
        if click_timestamp <= 0:
            return no_update, no_update, no_update

        action = triggered_id.get("action")
        if action not in {"window", "complete"}:
            return no_update, no_update, no_update
        rank_key = str(action)
        loss_runs = dict(loss_runs or {})
        existing_run = loss_runs.get(rank_key) or {}
        if float(existing_run.get("click_timestamp", -1.0) or -1.0) >= click_timestamp:
            return no_update, no_update, no_update
        if existing_run.get("state") == "running":
            return no_update, no_update, no_update

        scoreboard_results = list(top_results or best_results or [])
        param_sets = [result.get("params", {}) for result in scoreboard_results[:16]]
        if not param_sets:
            return no_update, no_update, "No result rows are available to run."
        try:
            if action == "window":
                case_configs = _case_configs_from_rows(
                    case_names,
                    time_start_values,
                    time_end_values,
                    average_time_values,
                    altitude_min_values,
                    altitude_max_values,
                )
                selected_case_names = [config["case_name"] for config in case_configs]
                run_case_configs = case_configs
                run_fields = list(selected_fields or [])
                message_prefix = "Started windowed loss run"
            else:
                selected_case_names = _case_list(case_names)
                run_case_configs = None
                run_fields = []
                message_prefix = "Started complete loss run"
            run_data = start_loss_run(
                selected_case_names,
                run_fields,
                param_sets,
                rank=action,
                case_configs=run_case_configs,
                run_mode=action,
                config=selected_config,
                batch_size=batch_size,
                override=scm_override,
                workspace_id=(workspace_selection or {}).get("workspace_id"),
                revision_id=(workspace_selection or {}).get("revision_id"),
                workspace_name=(workspace_selection or {}).get("display_name"),
            )
        except Exception as exc:
            return no_update, no_update, str(exc)

        run_data["click_timestamp"] = click_timestamp
        loss_runs[rank_key] = run_data
        return loss_runs, False, (
            f"{message_prefix} for {len(param_sets)} result rows "
            f"in {run_data['batch_count']} batch(es) of up to {run_data['batch_size']} "
            f"columns (pid {run_data['pid']}). Log: {run_data['log_path']}"
        )

    @app.callback(
        Output("tune-loss-runs", "data", allow_duplicate=True),
        Output("tune-loss-run-interval", "disabled", allow_duplicate=True),
        Input("tune-loss-run-interval", "n_intervals"),
        State("tune-loss-runs", "data"),
        prevent_initial_call=True,
    )
    def poll_result_loss_runs(_tick, loss_runs):
        """Poll ad-hoc loss runs and update row button states."""
        updated_runs, any_running = poll_loss_runs(loss_runs or {})
        return updated_runs, not any_running
