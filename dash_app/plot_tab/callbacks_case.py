import os

from dash import ALL, Input, Output, State, callback_context, html, no_update

from .benchmark_overlay import clear_benchmark_caches, sanitize_enabled_sources
from .layout import benchmark_button, case_button, directory_entry
from .plot_types.registry import PLOT_TYPES
from .plot_types.shared import (
    build_case_data,
    clear_all_caches,
    normalize_output_directory,
    ordered_case_names,
    scan_output_cases,
)
from .state import (
    empty_case_selection,
    entry_list_or_default,
    live_dir_entries,
    remap_plot_types_for_case_mode,
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
        Output("plots-output-dir-entries", "data"),
        Input("plots-add-dir-row", "n_clicks"),
        Input({"type": "plots-remove-dir", "index": ALL}, "n_clicks"),
        State("plots-output-dir-entries", "data"),
        State({"type": "plots-dir-entry", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def mutate_dir_entries(_add_clicks, _remove_clicks, current_entries, live_values):
        """Add or remove directory rows while preserving the current typed values."""
        entries = live_dir_entries(current_entries, live_values)
        trigger = callback_context.triggered_id
        if trigger == "plots-add-dir-row":
            entries.append("")
            return entries
        if isinstance(trigger, dict) and trigger.get("type") == "plots-remove-dir":
            idx = int(trigger.get("index"))
            if 0 <= idx < len(entries):
                entries.pop(idx)
            if not entries:
                entries = [""]
            return entries
        return no_update

    @app.callback(
        Output("plots-dir-list", "children"),
        Input("plots-output-dir-entries", "data"),
    )
    def render_dir_list(dir_entries):
        """Render the current list of editable directory rows."""
        entries = entry_list_or_default(dir_entries)
        return [directory_entry(index, value) for index, value in enumerate(entries)]

    @app.callback(
        Output("plots-output-dirs", "data"),
        Input({"type": "plots-dir-entry", "index": ALL}, "value"),
        State("plots-output-dir-entries", "data"),
        State("plots-output-dirs", "data"),
    )
    def build_output_dirs(live_values, dir_entries, current_output_dirs):
        """Normalize the entered directories and keep only existing unique paths."""
        entries = live_dir_entries(dir_entries, live_values)
        valid_dirs = []
        seen = set()
        for raw_entry in entries:
            candidate = (raw_entry or "").strip()
            if not candidate:
                continue
            normalized = normalize_output_directory(candidate)
            if not os.path.isdir(normalized):
                continue
            if normalized in seen:
                continue
            seen.add(normalized)
            valid_dirs.append(normalized)
        if valid_dirs == list(current_output_dirs or []):
            return no_update
        return valid_dirs

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
        Output("plots-time-mode", "value"),
        Output("plots-global-time-range", "max"),
        Output("plots-global-time-range", "value"),
        Output("plots-global-time-range", "marks"),
        Output("plots-global-time-point", "max"),
        Output("plots-global-time-point", "value"),
        Output("plots-global-time-point", "marks"),
        Output("plots-global-height-range", "min"),
        Output("plots-global-height-range", "max"),
        Output("plots-global-height-range", "value"),
        Output("plots-global-height-range", "marks"),
        Output("plots-global-height-range", "step"),
        Input({"type": "plots-case-button", "name": ALL}, "n_clicks"),
        Input("plots-output-dirs", "data"),
        Input("plots-refresh-cases", "n_clicks"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        State("plots-case-data", "data"),
        State("plots-enabled-benchmarks", "data"),
        prevent_initial_call=True,
    )
    def select_case(_clicks, output_dirs, _refresh_clicks, plot_order, plot_state, next_id, current_case_data, current_enabled_benchmarks):
        """Select a case and reset the global controls for that case's dimensions."""
        trigger = callback_context.triggered_id
        if trigger in ("plots-output-dirs", "plots-refresh-cases"):
            cases = scan_output_cases(output_dirs)
            available_names = ordered_case_names(cases.keys())
            if not available_names:
                _clear_plot_runtime_state()
                return empty_case_selection()
            preferred_name = (current_case_data or {}).get("name")
            case_name = preferred_name if preferred_name in cases else available_names[0]
        elif isinstance(trigger, dict):
            case_name = trigger.get("name")
        else:
            return (no_update,) * 19
        files = scan_output_cases(output_dirs).get(case_name, [])
        if not case_name or not files:
            return (no_update,) * 19
        current_name = (current_case_data or {}).get("name")
        current_files = list((current_case_data or {}).get("files") or [])
        current_dirs = list((current_case_data or {}).get("output_dirs") or [])
        next_files = list(files)
        next_dirs = list(output_dirs or [])
        context_changed = (
            trigger == "plots-output-dirs"
            or current_name != case_name
            or current_files != next_files
            or current_dirs != next_dirs
        )
        if context_changed:
            _clear_plot_runtime_state()
        case_data = build_case_data(case_name, files, output_dirs)
        updated_order = list(plot_order or [])
        updated_state = remap_plot_types_for_case_mode(plot_state, case_data)
        updated_next_id = int(next_id or 0)
        if not updated_order and not updated_state and updated_next_id == 0 and case_data.get("profile_vars"):
            from .state import default_plot_state

            initial_id = 0
            updated_order = [initial_id]
            updated_state[str(initial_id)] = default_plot_state(case_data, initial_id, existing_state=updated_state)
            updated_next_id = 1
        slider_max = max(case_data.get("time_len") or 1, 1)
        marks = {1: "1", slider_max: str(slider_max)} if slider_max > 1 else {1: "1"}
        height_min = float(case_data.get("height_slider_min", 0.0))
        height_max = float(case_data.get("height_slider_max", 1.0))
        height_marks = {height_min: f"{height_min:g}", height_max: f"{height_max:g}"}
        enabled_benchmarks = sanitize_enabled_sources(case_data, current_enabled_benchmarks)
        return (
            case_data,
            enabled_benchmarks,
            updated_order,
            updated_state,
            updated_next_id,
            0,
            "single",
            "range",
            slider_max,
            case_data.get("default_time_range") or [1, slider_max],
            marks,
            slider_max,
            case_data.get("default_time_range", [1, slider_max])[0],
            marks,
            height_min,
            height_max,
            case_data.get("default_height_range") or [height_min, height_max],
            height_marks,
            float(case_data.get("height_step", 1.0)),
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
        """Treat benchmark-source selection as top-level plotting context."""
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
        current = list(sanitize_enabled_sources(case_data, current_sources))
        if source in current:
            current.remove(source)
        else:
            current.append(source)
        sanitized = sanitize_enabled_sources(case_data, current)
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
        Input("plots-case-data", "data"),
        Input("plots-output-dirs", "data"),
    )
    def set_add_button_enabled_state(case_data, _output_dirs):
        """Enable add buttons only for plot families supported by the current case."""
        if not case_data:
            return True, True, True, True, True
        return (
            case_data.get("compare_mode") or not bool(case_data.get("budget_groups")),
            not bool(case_data.get("profile_vars")),
            not bool(case_data.get("timeseries_vars")),
            case_data.get("compare_mode") or not bool(case_data.get("timeheight_vars")),
            case_data.get("compare_mode") or not bool(case_data.get("subcolumn_vars")),
        )
