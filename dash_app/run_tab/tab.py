"""Thin entrypoint for the run tab package."""

from dash import dcc

from .callbacks_console import register_console_callbacks
from .callbacks_runs import register_run_callbacks
from .callbacks_selection import register_selection_callbacks
from .callbacks_settings import register_settings_callbacks
from .discovery import load_available_cases, load_case_groups, load_stats_choices
from .layout import (
    build_case_buttons,
    build_flag_controls,
    build_layout,
    build_param_sections,
    build_stats_buttons,
    compute_width_hints,
)
from .namelist import is_bool_value, is_true, normalize_numeric_display, read_namelist_entries
from .state import DEFAULT_STATS_NAME, FLAGS_FILE, NO_STATS_NAME, SILHS_FILE, TUNABLE_FILE


def build_initial_run_state():
    """Collect static case, stats, and parameter metadata for the initial run-tab layout."""
    cases = load_available_cases()
    case_groups = load_case_groups(cases)
    stats_files = load_stats_choices()
    if DEFAULT_STATS_NAME in stats_files:
        default_stats_name = DEFAULT_STATS_NAME
    elif stats_files:
        default_stats_name = stats_files[0]
    else:
        default_stats_name = DEFAULT_STATS_NAME

    flag_entries = read_namelist_entries(FLAGS_FILE)
    flag_bools = [entry for entry in flag_entries if is_bool_value(entry["value"])]
    flag_params = [entry for entry in flag_entries if entry not in flag_bools]
    tunable_entries = read_namelist_entries(TUNABLE_FILE)
    silhs_entries = read_namelist_entries(SILHS_FILE)

    flag_names = [entry["name"] for entry in flag_bools]
    all_config_names = (
        flag_names
        + [entry["name"] for entry in flag_params]
        + [entry["name"] for entry in tunable_entries]
        + [entry["name"] for entry in silhs_entries]
    )
    label_width_px, _value_width_px, right_pane_width_px = compute_width_hints(all_config_names)

    param_entries = (
        [{"file": "flags", **entry} for entry in flag_params]
        + [{"file": "tunable", **entry} for entry in tunable_entries]
        + [{"file": "silhs", **entry} for entry in silhs_entries]
    )
    defaults = {
        "flags": {entry["name"]: is_true(entry["value"]) for entry in flag_bools},
        "params": {
            "flags": {entry["name"]: entry["value"] for entry in flag_params},
            "tunable": {entry["name"]: entry["value"] for entry in tunable_entries},
            "silhs": {entry["name"]: entry["value"] for entry in silhs_entries},
        },
    }
    defaults_by_key = {
        f"{entry['file']}:{entry['name']}": normalize_numeric_display(entry["value"])
        for entry in param_entries
    }
    flag_controls = build_flag_controls(flag_bools, is_true)

    return {
        "cases": cases,
        "case_groups": case_groups,
        "default_stats_name": default_stats_name,
        "defaults": defaults,
        "defaults_by_key": defaults_by_key,
        "flag_names": flag_names,
        "tunable_names": [entry["name"] for entry in tunable_entries],
        "param_meta": [{"file": entry["file"], "name": entry["name"]} for entry in param_entries],
        "right_pane_width_px": right_pane_width_px,
        "case_buttons": build_case_buttons(cases),
        "stats_buttons": build_stats_buttons(stats_files, default_stats_name, NO_STATS_NAME),
        "param_sections": build_param_sections(flag_params, flag_controls, tunable_entries, silhs_entries, label_width_px, normalize_numeric_display),
    }


def build_tab(app):
    """Build the run tab and register its callback groups."""
    initial_state = build_initial_run_state()

    # Wire case and stats selection first because the remaining callbacks depend on these stores.
    register_selection_callbacks(app, initial_state["case_groups"])

    # Register settings synchronization before run lifecycle so dirty-state invalidation is in place.
    register_settings_callbacks(app)

    # Register run launch, cancel, clear, and polling callbacks that mutate process state.
    register_run_callbacks(app)

    # Register console rendering last because it depends on the selection and run-state stores above.
    register_console_callbacks(app, initial_state["cases"])

    return dcc.Tab(label="Run", value="run", children=build_layout(initial_state))
