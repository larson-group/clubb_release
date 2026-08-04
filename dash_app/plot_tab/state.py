import os

from dash import no_update

from .plot_types.profile_plot import PLOT as profile_plot
from .plot_types.registry import PLOT_TYPES
from dash_app.services.profiles import build_case_metadata as build_case_data, scan_case_outputs as scan_output_cases
from .plot_types.shared import (
    OUTPUT_DIR,
    normalize_output_directory,
    ordered_case_names,
    duration_slider_marks,
    start_time_slider_marks,
    snap_start_time_to_step,
    time_start_max_for_duration,
)

DEFAULT_OUTPUT_DIR = normalize_output_directory(os.path.join(OUTPUT_DIR, "dash_default"))
DEFAULT_PLAYBACK_INTERVAL_S = 1.0
PLAYBACK_INTERVAL_STEP_S = 0.1
MIN_PLAYBACK_INTERVAL_S = 0.1
MAX_PLAYBACK_INTERVAL_S = 5.0
MAX_AVERAGE_LENGTH_MINUTES = 240.0
DEFAULT_PLOT_VAR_PRIORITY = (
    "cloud_frac",
    "wp2",
    "wp3",
    "thlm",
    "rcm",
    "rtm",
    "wpthlp",
    "wprtp",
    "rcm",
    "thlp2",
)
INITIAL_DEFAULT_PLOT_COUNT = 3
MAX_INITIAL_DEFAULT_PLOT_COUNT = 4
INITIAL_DEFAULT_PLOT_TYPES = (
    profile_plot.plot_type_id,
    "timeseries",
    "subcolumn",
    "timeheight",
)


def normalize_playback_interval(interval_s):
    """Clamp and snap playback speed to the supported 0.1 s interval grid."""
    snapped = round(float(interval_s) / PLAYBACK_INTERVAL_STEP_S) * PLAYBACK_INTERVAL_STEP_S
    return round(min(MAX_PLAYBACK_INTERVAL_S, max(MIN_PLAYBACK_INTERVAL_S, snapped)), 1)


def clamp_float(value, low, high, fallback):
    """Return a finite float clamped to [low, high], or fallback when invalid."""
    try:
        result = float(value)
    except (TypeError, ValueError):
        result = float(fallback)
    low = float(low)
    high = float(high)
    if high < low:
        low, high = high, low
    return max(low, min(result, high))


def default_average_length(slider_min, slider_max):
    """Return the default averaging length in minutes."""
    return clamp_float(slider_min, slider_min, slider_max, slider_min)


def average_length_bounds(case_data):
    """Return averaging slider bounds capped at the supported maximum."""
    slider_min = max(1.0e-6, float((case_data or {}).get("time_slider_duration_min_minutes") or 1))
    raw_max = max(slider_min, float((case_data or {}).get("time_slider_duration_max_minutes") or slider_min))
    slider_max = max(slider_min, min(raw_max, MAX_AVERAGE_LENGTH_MINUTES))
    slider_step = max(1.0e-6, float((case_data or {}).get("time_slider_duration_step_minutes") or slider_min))
    return slider_min, slider_max, slider_step


def clamp_height_range(value, height_min, height_max, fallback):
    """Return a two-value height range clamped to the slider bounds."""
    if isinstance(value, (list, tuple)) and len(value) == 2:
        raw_low, raw_high = value
    elif isinstance(fallback, (list, tuple)) and len(fallback) == 2:
        raw_low, raw_high = fallback
    else:
        raw_low, raw_high = height_min, height_max
    low = clamp_float(raw_low, height_min, height_max, height_min)
    high = clamp_float(raw_high, height_min, height_max, height_max)
    if high < low:
        low, high = high, low
    return [low, high]


def empty_case_selection():
    """Return the callback payload for a UI state with no selectable cases."""
    return (
        {},
        [],
        no_update,
        no_update,
        no_update,
        0,
        "single",
        1,
        1,
        1,
        1,
        {},
        0,
        1,
        1,
        1,
        None,
        0.0,
        1.0,
        [0.0, 1.0],
        {0.0: "0", 1.0: "1"},
        1.0,
    )


def format_column_values(values):
    """Summarize overplot parameter values as a compact range plus sample list."""
    formatted = [f"{value:.3g}" for value in values]
    if formatted:
        min_val = min(values)
        max_val = max(values)
        range_text = f"{min_val:.3g}-{max_val:.3g}"
    else:
        range_text = "[]"
    if len(formatted) > 16:
        shown = formatted[:4] + ["..."] + formatted[-4:]
    else:
        shown = formatted
    return f"{range_text} [{', '.join(shown)}]"


def _unique_values(values):
    """Return values in their original order, dropping repeats and nulls."""
    seen = set()
    unique = []
    for value in values:
        if value is None or value in seen:
            continue
        seen.add(value)
        unique.append(value)
    return unique


def _option_values(options):
    """Return dropdown option values in display order."""
    return [option.get("value") for option in options if option.get("value") is not None]


def _preferred_plot_var(module, options, plot_id, used_vars=None, current_var=None):
    """Choose a default variable, prioritizing common CLUBB fields before fallback defaults."""
    values = _option_values(options)
    value_set = set(values)
    if not values:
        return None
    candidates = _unique_values(
        list(DEFAULT_PLOT_VAR_PRIORITY)
        + list(getattr(module, "default_vars", []))
        + [current_var, module.default_value(options, plot_id)]
    )
    used_vars = set(used_vars or [])
    for candidate in candidates:
        if candidate in value_set and candidate not in used_vars:
            return candidate
    for candidate in values:
        if candidate not in used_vars:
            return candidate
    return current_var if current_var in value_set else values[min(max(int(plot_id), 0), len(values) - 1)]


def _plot_type_values(case_data, plot_type):
    module = PLOT_TYPES[plot_type]
    return set(_option_values(module.case_data_options(case_data)))


def _plot_state_for_var(case_data, plot_id, plot_type, var_name):
    module = PLOT_TYPES[plot_type]
    state = module.make_default_state(case_data, plot_id)
    state["var"] = var_name
    state["size"] = state.get("size") or "normal"
    return state


def _append_plot(plot_order, plot_state, case_data, plot_type, var_name):
    plot_id = len(plot_order)
    plot_order.append(plot_id)
    plot_state[str(plot_id)] = _plot_state_for_var(case_data, plot_id, plot_type, var_name)


def _append_timeheight_fillers(plot_order, plot_state, case_data):
    if len(plot_order) >= INITIAL_DEFAULT_PLOT_COUNT:
        return
    timeheight_type = "timeheight"
    if not PLOT_TYPES[timeheight_type].supports_case_data(case_data):
        return
    timeheight_values = _plot_type_values(case_data, timeheight_type)
    existing_timeheight_vars = {
        entry.get("var")
        for entry in plot_state.values()
        if entry.get("plot_type") == timeheight_type
    }
    candidate_vars = _unique_values(entry.get("var") for entry in plot_state.values())
    for var_name in candidate_vars:
        if len(plot_order) >= MAX_INITIAL_DEFAULT_PLOT_COUNT:
            return
        if var_name not in timeheight_values or var_name in existing_timeheight_vars:
            continue
        _append_plot(plot_order, plot_state, case_data, timeheight_type, var_name)
        existing_timeheight_vars.add(var_name)


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
    selected_var = _preferred_plot_var(module, options, plot_id, used_vars=used_vars, current_var=state.get("var"))
    if selected_var is not None:
        state["var"] = selected_var
    state["size"] = state.get("size") or "normal"
    return state


def initial_plot_state_for_case(case_data):
    """Build the initial plot set for a case using curated variables before broad fallbacks."""
    plot_order = []
    plot_state = {}
    if not case_data:
        return plot_order, plot_state, 0

    plot_type_values = {
        plot_type: _plot_type_values(case_data, plot_type)
        for plot_type in INITIAL_DEFAULT_PLOT_TYPES
        if PLOT_TYPES[plot_type].supports_case_data(case_data)
    }
    for var_name in _unique_values(DEFAULT_PLOT_VAR_PRIORITY):
        for plot_type in INITIAL_DEFAULT_PLOT_TYPES:
            if var_name not in plot_type_values.get(plot_type, set()):
                continue
            _append_plot(plot_order, plot_state, case_data, plot_type, var_name)
            break
        if len(plot_order) >= INITIAL_DEFAULT_PLOT_COUNT:
            break

    if not plot_order and case_data.get("profile_vars"):
        for plot_id in range(min(INITIAL_DEFAULT_PLOT_COUNT, len(case_data.get("profile_vars") or []))):
            plot_order.append(plot_id)
            plot_state[str(plot_id)] = default_plot_state(case_data, plot_id, existing_state=plot_state)

    _append_timeheight_fillers(plot_order, plot_state, case_data)
    return plot_order, plot_state, len(plot_order)


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
            "time_slider_min": 1,
            "time_slider_max": 1,
            "time_slider_step": 1,
            "time_range": 1,
            "time_marks": {1: "1m"},
            "time_point_min": 0,
            "time_point_max": 1,
            "time_point": 1,
            "time_point_step": 1,
            "time_point_marks": {0: "0s", 1: "1s"},
            "height_min": 0.0,
            "height_max": 1.0,
            "height_range": [0.0, 1.0],
            "height_marks": {0.0: "0", 1.0: "1"},
            "height_step": 1.0,
        }
    case_name = ordered_names[0]
    case_data = build_case_data(case_name, cases[case_name], directories)
    plot_order, plot_state, next_id = initial_plot_state_for_case(case_data)
    slider_min, slider_max, slider_step = average_length_bounds(case_data)
    default_duration = default_average_length(slider_min, slider_max)
    default_start = snap_start_time_to_step(
        case_data,
        clamp_float(
            case_data.get("default_time_start_seconds", 0),
            case_data.get("time_slider_start_min_seconds", 0),
            time_start_max_for_duration(case_data, default_duration),
            case_data.get("time_slider_start_min_seconds", 0),
        ),
        default_duration,
    )
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
        "time_slider_min": slider_min,
        "time_slider_max": slider_max,
        "time_slider_step": slider_step,
        "time_range": default_duration,
        "time_marks": duration_slider_marks(slider_min, slider_max, default_duration),
        "time_point_min": case_data.get("time_slider_start_min_seconds", 0),
        "time_point_max": time_start_max_for_duration(case_data, default_duration),
        "time_point": default_start,
        "time_point_step": max(1.0e-6, float(default_duration)) * 60.0,
        "time_point_marks": start_time_slider_marks(case_data, default_start, default_duration),
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
