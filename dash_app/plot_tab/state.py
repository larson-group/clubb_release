from dash import no_update

from .plot_types.profile_plot import PLOT as profile_plot
from .plot_types.registry import PLOT_TYPES
from .plot_types.shared import (
    OUTPUT_DIR,
    build_case_data,
    normalize_output_directory,
    ordered_case_names,
    scan_output_cases,
)

DEFAULT_OUTPUT_DIR = normalize_output_directory(OUTPUT_DIR)
DEFAULT_OUTPUT_ENTRY = "output/"
DEFAULT_PLAYBACK_INTERVAL_S = 1.0
PLAYBACK_INTERVAL_STEP_S = 0.1
MIN_PLAYBACK_INTERVAL_S = 0.1
MAX_PLAYBACK_INTERVAL_S = 5.0


def normalize_playback_interval(interval_s):
    """Clamp and snap playback speed to the supported 0.1 s interval grid."""
    snapped = round(float(interval_s) / PLAYBACK_INTERVAL_STEP_S) * PLAYBACK_INTERVAL_STEP_S
    return round(min(MAX_PLAYBACK_INTERVAL_S, max(MIN_PLAYBACK_INTERVAL_S, snapped)), 1)


def entry_list_or_default(dir_entries):
    """Return the stored directory entries or the default single output entry."""
    if dir_entries is None:
        return [DEFAULT_OUTPUT_ENTRY]
    entries = list(dir_entries)
    return entries if entries else [""]


def live_dir_entries(stored_entries, live_values):
    """Prefer the currently typed directory values when they match the stored rows."""
    entries = entry_list_or_default(stored_entries)
    if live_values is not None and len(live_values) == len(entries):
        return list(live_values)
    return entries


def empty_case_selection():
    """Return the callback payload for a UI state with no selectable cases."""
    return (
        None,
        [],
        no_update,
        no_update,
        no_update,
        0,
        "single",
        "range",
        1,
        [1, 1],
        {1: "1"},
        1,
        1,
        {1: "1"},
        0.0,
        1.0,
        [0.0, 1.0],
        {0.0: "0", 1.0: "1"},
        1.0,
    )


def format_column_values(values):
    """Summarize overplot parameter values as a compact range plus sample list."""
    formatted = [f"{value:g}" for value in values]
    if formatted:
        min_val = min(values)
        max_val = max(values)
        range_text = f"{min_val:g}-{max_val:g}"
    else:
        range_text = "[]"
    if len(formatted) > 6:
        shown = formatted[:2] + ["..."] + formatted[-2:]
    else:
        shown = formatted
    return f"{range_text} [{', '.join(shown)}]"


def default_plot_state(case_data, plot_id, plot_type=None, existing_state=None):
    """Build a default plot-state entry, preferring a variable not already in use."""
    module = PLOT_TYPES[plot_type or profile_plot.plot_type_id]
    state = module.make_default_state(case_data, plot_id)
    options = list(module.case_data_options(case_data))
    used_vars = {
        (entry or {}).get("var")
        for entry in (existing_state or {}).values()
        if (entry or {}).get("plot_type") == module.plot_type_id and (entry or {}).get("var")
    }
    current_var = state.get("var")
    if current_var in used_vars:
        for option in options:
            candidate = option.get("value")
            if candidate not in used_vars:
                state["var"] = candidate
                break
    state["size"] = state.get("size") or "normal"
    return state


def initialize_case_state(output_dirs=None):
    """Build the initial plots-tab state from the currently available output files."""
    directories = list(output_dirs or [DEFAULT_OUTPUT_DIR])
    cases = scan_output_cases(directories)
    ordered_names = ordered_case_names(cases.keys())
    if not ordered_names:
        return {
            "case_data": None,
            "plot_order": [],
            "plot_state": {},
            "next_id": 0,
            "enabled_benchmarks": [],
            "selected_column": 0,
            "column_mode": "single",
            "time_mode": "range",
            "time_slider_max": 1,
            "time_range": [1, 1],
            "time_marks": {1: "1"},
            "time_point_max": 1,
            "time_point": 1,
            "time_point_marks": {1: "1"},
            "height_min": 0.0,
            "height_max": 1.0,
            "height_range": [0.0, 1.0],
            "height_marks": {0.0: "0", 1.0: "1"},
            "height_step": 1.0,
        }
    case_name = ordered_names[0]
    case_data = build_case_data(case_name, cases[case_name], directories)
    plot_order = []
    plot_state = {}
    next_id = 0
    if case_data.get("profile_vars"):
        max_initial = min(3, len(case_data.get("profile_vars") or []))
        for plot_id in range(max_initial):
            plot_order.append(plot_id)
            plot_state[str(plot_id)] = default_plot_state(case_data, plot_id, existing_state=plot_state)
        next_id = max_initial
    slider_max = max(case_data.get("time_len") or 1, 1)
    marks = {1: "1", slider_max: str(slider_max)} if slider_max > 1 else {1: "1"}
    height_min = float(case_data.get("height_slider_min", 0.0))
    height_max = float(case_data.get("height_slider_max", 1.0))
    return {
        "case_data": case_data,
        "plot_order": plot_order,
        "plot_state": plot_state,
        "next_id": next_id,
        "enabled_benchmarks": [],
        "selected_column": 0,
        "column_mode": "single",
        "time_mode": "range",
        "time_slider_max": slider_max,
        "time_range": case_data.get("default_time_range") or [1, slider_max],
        "time_marks": marks,
        "time_point_max": slider_max,
        "time_point": case_data.get("default_time_range", [1, slider_max])[0],
        "time_point_marks": marks,
        "height_min": height_min,
        "height_max": height_max,
        "height_range": case_data.get("default_height_range") or [height_min, height_max],
        "height_marks": {height_min: f"{height_min:g}", height_max: f"{height_max:g}"},
        "height_step": float(case_data.get("height_step", 1.0)),
    }


def build_initial_plots_state(output_dirs=None):
    """Public wrapper for the initial plots-tab state construction."""
    return initialize_case_state(output_dirs=output_dirs)


def remap_plot_types_for_case_mode(plot_state, _case_data):
    """Return plot state unchanged after profile/compare consolidation."""
    return dict(plot_state or {})
