import os

from dash import ALL, Input, Output, State, callback_context, html, no_update

from .benchmark_overlay import (
    clear_benchmark_caches,
    sanitize_enabled_sources,
)
from .layout import (
    active_output_items,
    available_output_buttons,
    benchmark_button,
    case_button,
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


def _normalize_output_dirs(output_dirs):
    """Normalize and de-duplicate selection state without pruning missing paths."""
    normalized = []
    seen = set()
    for raw_entry in output_dirs or []:
        candidate = str(raw_entry or "").strip()
        if not candidate:
            continue
        path = normalize_output_directory(candidate)
        if path not in seen:
            normalized.append(path)
            seen.add(path)
    return normalized


def _update_output_dirs(output_dirs, action, path):
    """Apply one explicit add/remove action while preserving every other path."""
    selected = _normalize_output_dirs(output_dirs)
    normalized = normalize_output_directory(str(path or ""))
    if action == "add" and normalized not in selected:
        return [*selected, normalized]
    if action == "remove":
        return [candidate for candidate in selected if candidate != normalized]
    return selected


def _catalog_tracking_paths(records, selected_dirs):
    """Keep active paths and previously added external paths in the catalog."""
    tracked = _normalize_output_dirs(selected_dirs)
    for record in records or []:
        if record.get("catalog_origin") != "external":
            continue
        path = str(record.get("path") or "")
        if path and path not in tracked:
            tracked.append(path)
    return tracked


def _triggered_click_is_positive():
    triggered = callback_context.triggered[0] if callback_context.triggered else None
    return bool(triggered) and _is_positive_click(triggered.get("value"))


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
        Output("plots-output-refresh-interval", "disabled"),
        Input("dashboard-tabs", "value"),
    )
    def enable_output_refresh(active_tab):
        """Poll only while Plot is the active top-level tab."""
        return active_tab != "plots"

    @app.callback(
        Output("plots-output-catalog", "data"),
        Input("dashboard-tabs", "value"),
        Input("plots-output-refresh-interval", "n_intervals"),
        State("plots-output-dirs", "data"),
        State("plots-output-catalog", "data"),
        prevent_initial_call=True,
    )
    def refresh_output_catalog(active_tab, _intervals, selected_dirs, current_catalog):
        """Refresh on Plot activation and each active ten-second interval."""
        if active_tab != "plots":
            return no_update
        tracked = _catalog_tracking_paths(current_catalog, selected_dirs)
        return discover_output_directories(selected_dirs=tracked)

    @app.callback(
        Output("plots-output-menu-expanded", "data"),
        Input("plots-output-menu-toggle", "n_clicks_timestamp"),
        State("plots-output-menu-expanded", "data"),
        prevent_initial_call=True,
    )
    def toggle_output_menu(_timestamp, expanded):
        if not _triggered_click_is_positive():
            return no_update
        return not bool(expanded)

    @app.callback(
        Output("plots-available-output-list", "children"),
        Output("plots-active-output-list", "children"),
        Output("plots-output-menu", "className"),
        Input("plots-output-catalog", "data"),
        Input("plots-output-dirs", "data"),
        Input("plots-output-menu-expanded", "data"),
    )
    def render_output_controls(records, selected_dirs, expanded):
        selected = _normalize_output_dirs(selected_dirs)
        known_paths = {str(record.get("path")) for record in (records or [])}
        if any(path not in known_paths for path in selected):
            records = discover_output_directories(selected_dirs=selected)
        menu_class = "plots-output-menu plots-output-menu--expanded" if expanded else "plots-output-menu"
        return (
            available_output_buttons(records, selected, expanded=bool(expanded)),
            active_output_items(records, selected),
            menu_class,
        )

    @app.callback(
        Output("plots-output-dirs", "data"),
        Output("plots-extra-dir-input", "value"),
        Output("plots-extra-dir-message", "children"),
        Output("plots-output-catalog", "data", allow_duplicate=True),
        Input({"type": "plots-add-output-dir", "path": ALL}, "n_clicks_timestamp"),
        Input("plots-add-extra-dir", "n_clicks_timestamp"),
        Input({"type": "plots-remove-output-dir", "path": ALL}, "n_clicks_timestamp"),
        State("plots-output-dirs", "data"),
        State("plots-extra-dir-input", "value"),
        State("plots-output-catalog", "data"),
        prevent_initial_call=True,
    )
    def build_output_dirs(_add_clicks, _add_extra, _remove_clicks, current_output_dirs, extra_dir, catalog):
        """Move one clicked output between the available menu and active tray."""
        if not _triggered_click_is_positive():
            return no_update, no_update, no_update, no_update
        selected = _normalize_output_dirs(current_output_dirs)
        seen = set(selected)

        trigger = callback_context.triggered_id
        if isinstance(trigger, dict) and trigger.get("type") == "plots-add-output-dir":
            normalized = normalize_output_directory(str(trigger.get("path") or ""))
            discovered = {str(record.get("path")) for record in (catalog or [])}
            if normalized not in discovered or normalized in seen:
                return no_update, no_update, no_update, no_update
            selected = _update_output_dirs(selected, "add", normalized)
        elif trigger == "plots-add-extra-dir":
            candidate = str(extra_dir or "").strip()
            if not candidate:
                return no_update, no_update, "Enter a directory path before adding it.", no_update
            normalized = normalize_output_directory(candidate)
            if not os.path.isdir(normalized):
                return no_update, no_update, f"Not found: {normalized}", no_update
            if normalized not in seen:
                selected = _update_output_dirs(selected, "add", normalized)
            else:
                return no_update, no_update, "That folder is already selected.", no_update
        elif isinstance(trigger, dict) and trigger.get("type") == "plots-remove-output-dir":
            removed = normalize_output_directory(str(trigger.get("path") or ""))
            selected = _update_output_dirs(selected, "remove", removed)
        else:
            return no_update, no_update, no_update, no_update
        if selected == list(current_output_dirs or []):
            return no_update, no_update, no_update, no_update
        tracked = _catalog_tracking_paths(catalog, selected)
        updated_catalog = discover_output_directories(selected_dirs=tracked)
        return selected, "" if trigger == "plots-add-extra-dir" else no_update, "", updated_catalog

    @app.callback(
        Output("plots-case-button-container", "children"),
        Input("plots-output-dirs", "data"),
        Input("plots-case-data", "data"),
    )
    def render_case_buttons(output_dirs, case_data):
        """Render the case buttons for the active directory set and selection."""
        cases = scan_output_cases(output_dirs)
        selected_name = case_data.get("name") if case_data else None
        available_names = ordered_case_names(cases.keys())
        if not available_names:
            return [html.Div("No cases found in the active outputs.")]
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
        if trigger == "plots-output-dirs":
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
            if not _triggered_click_is_positive():
                return (no_update,) * 23
            case_name = trigger.get("name")
        else:
            return (no_update,) * 23
        files = scan_output_cases(output_dirs).get(case_name, [])
        if not case_name or not files:
            return (no_update,) * 23
        case_data = build_case_data(case_name, files, output_dirs)
        same_case = _is_same_case(current_case_data, case_name, case_data.get("output_dirs"))
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
