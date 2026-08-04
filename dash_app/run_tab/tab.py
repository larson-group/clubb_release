"""Thin entrypoint for the run tab package."""

from dash import dcc

from .callbacks_console import register_console_callbacks
from .callbacks_runs import register_run_callbacks
from .callbacks_selection import register_selection_callbacks
from .callbacks_settings import register_settings_callbacks
from .config_state import build_tunable_config_state
from .discovery import load_available_cases, load_case_groups, load_stats_choices
from .layout import build_case_buttons, build_layout, build_stats_buttons
from .state import DEFAULT_STATS_NAME, NO_STATS_NAME
from dash_app.shared.tunable_configs import available_tunable_configs, default_tunable_config_name


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

    tunable_configs = available_tunable_configs()
    config_state = build_tunable_config_state(default_tunable_config_name(tunable_configs), tunable_configs)

    return {
        "cases": cases,
        "case_groups": case_groups,
        "default_stats_name": default_stats_name,
        "tunable_configs": tunable_configs,
        "case_buttons": build_case_buttons(cases),
        "stats_buttons": build_stats_buttons(stats_files, default_stats_name, NO_STATS_NAME),
        **config_state,
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
