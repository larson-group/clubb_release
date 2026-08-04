import os

from dash import ALL, Input, Output, State, callback_context, html, no_update

from .benchmark_overlay import (
    clear_benchmark_caches,
    sanitize_enabled_sources,
)
from .layout import (
    benchmark_button,
    case_button,
    output_directory_options,
    selected_output_directory_chips,
)
from .plot_types.registry import PLOT_TYPES
from dash_app.services.profiles import (
    build_case_metadata as build_case_data,
    discover_output_directories,
    scan_case_outputs as scan_output_cases,
)
from dash_app.services.plots import apply_plot_request, resolve_benchmark_sources, toggle_benchmark_source
from .plot_types.shared import (
    clear_all_caches,
    duration_slider_marks,
    normalize_output_directory,
    ordered_case_names,
    snap_start_time_to_step,
    start_time_slider_marks,
    time_start_max_for_duration,
)
from .state import (
    clamp_float,
    clamp_height_range,
    average_length_bounds,
    default_average_length,
    empty_case_selection,
    initial_plot_state_for_case,
    remap_plot_types_for_case_mode,
)


def _is_positive_click(value):
    """Return true only for an actual Dash button click count."""
    return isinstance(value, (int, float)) and value > 0


def _is_same_case(current_case_data, case_name, output_dirs):
    """Return whether a reload keeps the same case and output directories."""
    current = current_case_data or {}
    return (
        current.get("name") == case_name
        and list(current.get("output_dirs") or []) == list(output_dirs or [])
    )


def _clear_plot_runtime_state(clear_shared=True, clear_benchmarks=True):
    if clear_shared:
        clear_all_caches()
    if clear_benchmarks:
        clear_benchmark_caches()
    for plot in PLOT_TYPES.values():
        clear_state = getattr(plot, "clear_render_state", None)
        if callable(clear_state):
            clear_state()


def register_case_callbacks(app):
    """Register callbacks that manage directories, cases, and case-driven resets."""
    @app.callback(
        Output("plots-extra-dir-control", "style"),
        Input("plots-show-extra-dir", "n_clicks"),
        State("plots-extra-dir-control", "style"),
        prevent_initial_call=True,
    )
    def show_extra_directory_input(clicks, current_style):
        """Reveal the optional raw-path escape hatch only when requested."""
        if not _is_positive_click(clicks):
            return no_update
        return {
            **(current_style or {}),
            "display": "flex",
            "marginTop": "8px",
            "alignItems": "center",
            "gap": "8px",
        }

    @app.callback(
        Output("plots-output-dir-picker", "options"),
        Output("plots-selected-output-dirs", "children"),
        Input("plots-output-dirs", "data"),
        Input("plots-refresh-cases", "n_clicks"),
    )
    def refresh_output_directory_picker(selected_dirs, _refresh_clicks):
        """Rescan addable folders and render compact selected-folder chips."""
        selected = [
            normalize_output_directory(str(path).strip())
            for path in selected_dirs or []
            if str(path or "").strip() and os.path.isdir(normalize_output_directory(str(path).strip()))
        ]
        records = discover_output_directories()
        return (
            output_directory_options(records, selected),
            selected_output_directory_chips(records, selected),
        )

    @app.callback(
        Output("plots-output-dirs", "data"),
        Output("plots-extra-dir-input", "value"),
        Output("plots-extra-dir-message", "children"),
        Input("plots-output-dir-picker", "value"),
        Input("plots-add-extra-dir", "n_clicks"),
        Input({"type": "plots-remove-output-dir", "path": ALL}, "n_clicks"),
        State("plots-output-dirs", "data"),
        State("plots-extra-dir-input", "value"),
        prevent_initial_call=True,
    )
    def build_output_dirs(selected_dir, _add_extra, _remove_clicks, current_output_dirs, extra_dir):
        """Add or remove one output folder while retaining a compact selection."""
        selected = []
        seen = set()
        for raw_entry in current_output_dirs or []:
            candidate = str(raw_entry or "").strip()
            if not candidate:
                continue
            normalized = normalize_output_directory(candidate)
            if not os.path.isdir(normalized):
                continue
            if normalized in seen:
                continue
            seen.add(normalized)
            selected.append(normalized)

        trigger = callback_context.triggered_id
        if trigger == "plots-output-dir-picker":
            candidate = str(selected_dir or "").strip()
            if candidate:
                normalized = normalize_output_directory(candidate)
                if os.path.isdir(normalized) and normalized not in seen:
                    selected.append(normalized)
        elif trigger == "plots-add-extra-dir":
            candidate = str(extra_dir or "").strip()
            if not candidate:
                return no_update, no_update, "Enter a directory path before adding it."
            normalized = normalize_output_directory(candidate)
            if not os.path.isdir(normalized):
                return no_update, no_update, f"Not found: {normalized}"
            if normalized not in seen:
                selected.append(normalized)
            else:
                return no_update, no_update, "That folder is already selected."
        elif isinstance(trigger, dict) and trigger.get("type") == "plots-remove-output-dir":
            removed = normalize_output_directory(str(trigger.get("path") or ""))
            selected = [path for path in selected if path != removed]
        else:
            return no_update, no_update, no_update
        if selected == list(current_output_dirs or []):
            return no_update, no_update, no_update
        return selected, "" if trigger == "plots-add-extra-dir" else no_update, ""

    @app.callback(
        Output("plots-case-button-container", "children"),
        Input("plots-output-dirs", "data"),
        Input("plots-case-data", "data"),
        Input("plots-refresh-cases", "n_clicks"),
    )
    def render_case_buttons(output_dirs, case_data, _refresh_clicks):
        """Render the case buttons for the active directory set and selection."""
        cases = scan_output_cases(output_dirs)
        selected_name = case_data.get("name") if case_data else None
        available_names = ordered_case_names(cases.keys())
        if not available_names:
            mode = "compare" if len(output_dirs or []) > 1 else "single"
            return [html.Div(f"No common cases found for {mode} mode.")]
        return [case_button(name, bool(cases.get(name)), selected=(name == selected_name)) for name in available_names]

    @app.callback(
        Output("plots-case-data", "data"),
        Output("plots-enabled-benchmarks", "data"),
        Output("plots-plot-order", "data"),
        Output("plots-plot-state", "data"),
        Output("plots-next-id", "data"),
        Output("plots-selected-column", "data"),
        Output("plots-column-mode", "value"),
        Output("plots-global-time-range", "min"),
        Output("plots-global-time-range", "max"),
        Output("plots-global-time-range", "value"),
        Output("plots-global-time-range", "step"),
        Output("plots-global-time-range", "marks"),
        Output("plots-global-time-point", "min"),
        Output("plots-global-time-point", "max"),
        Output("plots-global-time-point", "value"),
        Output("plots-global-time-point", "step"),
        Output("plots-global-time-point", "marks"),
        Output("plots-global-height-range", "min"),
        Output("plots-global-height-range", "max"),
        Output("plots-global-height-range", "value"),
        Output("plots-global-height-range", "marks"),
        Output("plots-global-height-range", "step"),
        Output("plots-time-override", "data", allow_duplicate=True),
        Input({"type": "plots-case-button", "name": ALL}, "n_clicks"),
        Input("plots-output-dirs", "data"),
        Input("plots-refresh-cases", "n_clicks"),
        Input("dashboard-request", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        State("plots-case-data", "data"),
        State("plots-enabled-benchmarks", "data"),
        State("plots-selected-column", "data"),
        State("plots-column-mode", "value"),
        State("plots-global-time-range", "value"),
        State("plots-global-time-point", "value"),
        State("plots-global-height-range", "value"),
        # Opening a large multi-column case performs substantial NetCDF I/O.
        # Run it in Diskcache's separate process: Flask remains safe and the
        # user can continue navigating while metadata is assembled.  The
        # callback stays data-only: process-local Plot caches must never be
        # cleared or marked from this worker.
        background=True,
        interval=200,
        prevent_initial_call=True,
    )
    def select_case(
        _clicks,
        output_dirs,
        _refresh_clicks,
        agent_request,
        plot_order,
        plot_state,
        next_id,
        current_case_data,
        current_enabled_benchmarks,
        current_column,
        current_column_mode,
        current_average_minutes,
        current_start_time,
        current_height_range,
    ):
        """Select a case and refresh global controls without resetting same-case reloads."""
        trigger = callback_context.triggered_id
        if trigger in ("plots-output-dirs", "plots-refresh-cases"):
            cases = scan_output_cases(output_dirs)
            available_names = ordered_case_names(cases.keys())
            if not available_names:
                return (*empty_case_selection(), None)
            preferred_name = (current_case_data or {}).get("name")
            case_name = preferred_name if preferred_name in cases else available_names[0]
        elif trigger == "dashboard-request":
            if (agent_request or {}).get("tab") != "plots" or (agent_request or {}).get("operation") not in {"set_view", "add_budget"}:
                return (no_update,) * 23
            requested_output_dir = str((agent_request or {}).get("output_dir") or "").strip()
            if requested_output_dir:
                output_dirs = [requested_output_dir]
            cases = scan_output_cases(output_dirs)
            requested_case = str((agent_request or {}).get("case") or "")
            if requested_case not in cases:
                return (no_update,) * 23
            case_name = requested_case
        elif isinstance(trigger, dict):
            case_name = trigger.get("name")
        else:
            return (no_update,) * 23
        files = scan_output_cases(output_dirs).get(case_name, [])
        if not case_name or not files:
            return (no_update,) * 23
        same_case = _is_same_case(current_case_data, case_name, output_dirs)
        case_data = build_case_data(case_name, files, output_dirs)
        case_data["preserve_plot_view"] = bool(same_case)
        updated_order = list(plot_order or [])
        updated_state = remap_plot_types_for_case_mode(plot_state, case_data)
        updated_next_id = int(next_id or 0)
        if trigger == "dashboard-request":
            request = dict(agent_request or {})
            if request.get("operation") == "add_budget" and not same_case:
                updated_order, updated_state, updated_next_id = initial_plot_state_for_case(case_data)
            transition = apply_plot_request(
                case_data,
                request,
                updated_order,
                updated_state,
                updated_next_id,
            )
            updated_order = transition.plot_order
            updated_state = transition.plot_state
            updated_next_id = transition.next_id
        if not updated_order and not updated_state and updated_next_id == 0:
            updated_order, updated_state, updated_next_id = initial_plot_state_for_case(case_data)
        slider_min, slider_max, slider_step = average_length_bounds(case_data)
        default_duration = default_average_length(slider_min, slider_max)
        requested_average = (agent_request or {}).get("average_minutes") if trigger == "dashboard-request" else None
        active_duration = clamp_float(
            requested_average if requested_average is not None else (current_average_minutes if same_case else default_duration),
            slider_min,
            slider_max,
            default_duration,
        )
        start_min = float(case_data.get("time_slider_start_min_seconds", 0))
        start_max = time_start_max_for_duration(case_data, active_duration)
        requested_start = None
        if trigger == "dashboard-request":
            requested_start = (agent_request or {}).get("time_start_seconds")
            if requested_start is None:
                requested_start = (agent_request or {}).get("time_seconds")
        active_start = snap_start_time_to_step(
            case_data,
            clamp_float(
                requested_start if requested_start is not None else (current_start_time if same_case else case_data.get("default_time_start_seconds", start_min)),
                start_min,
                start_max,
                start_min,
            ),
            active_duration,
        )
        height_min = float(case_data.get("height_slider_min", 0.0))
        height_max = float(case_data.get("height_slider_max", 1.0))
        default_height_range = case_data.get("default_height_range") or [height_min, height_max]
        active_height_range = clamp_height_range(
            current_height_range if same_case else default_height_range,
            height_min,
            height_max,
            default_height_range,
        )
        height_marks = {height_min: f"{height_min:g}", height_max: f"{height_max:g}"}
        max_column = max(int(case_data.get("columns_len") or 1) - 1, 0)
        active_column = int(clamp_float(current_column if same_case else 0, 0, max_column, 0))
        active_column_mode = current_column_mode if same_case and current_column_mode in {"single", "all"} else "single"
        if trigger == "dashboard-request" and "benchmark_sources" in (agent_request or {}):
            enabled_benchmarks = resolve_benchmark_sources(
                case_data,
                (agent_request or {}).get("benchmark_sources"),
                strict=True,
            )
        else:
            enabled_benchmarks = sanitize_enabled_sources(case_data, current_enabled_benchmarks)
        time_override = None
        preset = str((agent_request or {}).get("window_preset") or "") if trigger == "dashboard-request" else ""
        if preset in {"loss", "pyplotgen"}:
            exact_start = case_data.get(f"{preset}_time_start_seconds")
            exact_duration = case_data.get(f"{preset}_time_duration_minutes")
            try:
                exact_start = float(exact_start)
                exact_duration = float(exact_duration)
            except (TypeError, ValueError):
                exact_start = None
                exact_duration = None
            if exact_start is not None and exact_duration is not None:
                time_override = {
                    "mode": preset,
                    "start_seconds": exact_start,
                    "duration_minutes": exact_duration,
                    "slider_start_seconds": float(active_start),
                    "slider_duration_minutes": float(active_duration),
                }
        return (
            case_data,
            enabled_benchmarks,
            updated_order,
            updated_state,
            updated_next_id,
            active_column,
            active_column_mode,
            slider_min,
            slider_max,
            active_duration,
            slider_step,
            duration_slider_marks(slider_min, slider_max, active_duration),
            start_min,
            start_max,
            active_start,
            max(1.0e-6, float(active_duration)) * 60.0,
            start_time_slider_marks(case_data, active_start, active_duration),
            height_min,
            height_max,
            active_height_range,
            height_marks,
            float(case_data.get("height_step", 1.0)),
            time_override,
        )

    @app.callback(
        Output("plots-benchmark-button-container", "children"),
        Input("plots-case-data", "data"),
        Input("plots-enabled-benchmarks", "data"),
    )
    def sync_benchmark_controls(case_data, enabled_benchmarks):
        """Render benchmark toggle buttons in the header for the active case."""
        available = set((case_data or {}).get("benchmarks", {}).get("available_sources") or [])
        selected = set(sanitize_enabled_sources(case_data, enabled_benchmarks))
        return [
            benchmark_button("sam", "SAM LES", "sam" in available, selected=("sam" in selected)),
            benchmark_button("coamps", "COAMPS LES", "coamps" in available, selected=("coamps" in selected)),
        ]

    @app.callback(
        Output("plots-enabled-benchmarks", "data", allow_duplicate=True),
        Input({"type": "plots-benchmark-button", "source": ALL}, "n_clicks_timestamp"),
        State("plots-case-data", "data"),
        State("plots-enabled-benchmarks", "data"),
        prevent_initial_call=True,
    )
    def update_enabled_benchmarks(_click_timestamps, case_data, current_sources):
        trigger = callback_context.triggered_id
        if not isinstance(trigger, dict) or trigger.get("type") != "plots-benchmark-button":
            return no_update
        triggered = callback_context.triggered[0] if callback_context.triggered else None
        if not triggered:
            return no_update
        try:
            if int(triggered.get("value") or -1) < 0:
                return no_update
        except (TypeError, ValueError):
            return no_update
        source = trigger.get("source")
        sanitized = toggle_benchmark_source(case_data, current_sources, source)
        if sanitized == list(current_sources or []):
            return no_update
        _clear_plot_runtime_state(clear_shared=False, clear_benchmarks=True)
        return sanitized

    @app.callback(
        Output("plots-add-budget", "disabled"),
        Output("plots-add-profile", "disabled"),
        Output("plots-add-timeseries", "disabled"),
        Output("plots-add-timeheight", "disabled"),
        Output("plots-add-subcolumn", "disabled"),
        Output("plots-add-pdf-contour", "disabled"),
        Input("plots-case-data", "data"),
        Input("plots-output-dirs", "data"),
    )
    def set_add_button_enabled_state(case_data, _output_dirs):
        """Enable add buttons only for plot families supported by the current case."""
        if not case_data:
            return True, True, True, True, True, True
        return (
            case_data.get("compare_mode") or not bool(case_data.get("budget_groups")),
            not bool(case_data.get("profile_vars")),
            not bool(case_data.get("timeseries_vars")),
            not PLOT_TYPES["timeheight"].supports_case_data(case_data),
            case_data.get("compare_mode") or not bool(case_data.get("subcolumn_vars")),
            not PLOT_TYPES["pdf_contour"].supports_case_data(case_data),
        )
