"""Thin entrypoint for the tuning tab package."""

from dash import dcc

from tuner.system_defaults import default_max_workers
from tuner.taylor_metrics import (
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_LOSS_MODE,
    DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
)

from .callbacks_display import register_display_callbacks
from .callbacks_runs import register_run_callbacks
from .callbacks_settings import register_settings_callbacks
from .callbacks_workspaces import register_workspace_callbacks
from .discovery import (
    available_fields_for_cases,
    available_tunable_configs,
    default_tunable_config_name,
    load_case_defaults,
    load_tunable_default_ranges,
    load_tunable_names,
)
from .layout import DEFAULT_AVERAGE_TIME_SECONDS, build_layout
from .runtime import empty_status_payload


def build_initial_tune_state():
    """Collect static metadata for the initial tuning-tab layout."""
    case_data = load_case_defaults()
    tunable_configs = available_tunable_configs()
    selected_config = default_tunable_config_name(tunable_configs)
    tunable_names = load_tunable_names(selected_config)
    tunable_default_ranges = load_tunable_default_ranges(selected_config)
    cases = sorted(case_data.keys())
    # A fresh Tune workspace should be immediately runnable without relying
    # on alphabetical case ordering.  Revisions intentionally bypass this
    # default and hydrate the exact request they were branched from.
    selected_cases = ["arm"] if "arm" in case_data else [cases[0]] if cases else []
    selected_case_defaults = case_data.get(selected_cases[0], {}) if selected_cases else {}
    field_options = available_fields_for_cases(selected_cases, case_data)
    max_name_len = max((len(name) for name in tunable_names), default=16)
    right_pane_width_px = max(360, min(760, int(180 + max_name_len * 7.5)))
    status = empty_status_payload()
    status_text = "state: idle | samples: 0 | best smart loss: --"
    initial_parameter_names = ("C4", "C8", "C11")
    initial_param_rows = [
        {
            "id": index,
            "param": name,
            "min": tunable_default_ranges.get(name, {}).get("min", ""),
            "max": tunable_default_ranges.get(name, {}).get("max", ""),
        }
        for index, name in enumerate(initial_parameter_names)
        if name in tunable_default_ranges
    ]
    initial_case_rows = []
    if selected_cases:
        time_range = selected_case_defaults.get("time_average_range", ["", ""])
        altitude_range = selected_case_defaults.get("altitude_comparison_range", ["", ""])
        average_time_seconds = int(
            selected_case_defaults.get("average_time_seconds") or DEFAULT_AVERAGE_TIME_SECONDS
        )
        initial_case_rows.append(
            {
                "id": 0,
                "case_name": selected_cases[0],
                "time_start": time_range[0],
                "time_end": time_range[1],
                "average_time_seconds": average_time_seconds,
                "altitude_min": altitude_range[0] if len(altitude_range) > 0 else "",
                "altitude_max": altitude_range[1] if len(altitude_range) > 1 else "",
            }
        )

    return {
        "cases": cases,
        "case_data": case_data,
        "tunable_configs": tunable_configs,
        "selected_config": selected_config,
        "tunable_names": tunable_names,
        "tunable_default_ranges": tunable_default_ranges,
        "selected_cases": selected_cases,
        "selected_fields": [
            field_name
            for field_name in selected_case_defaults.get("default_clubb_fields", [])
            if field_name in field_options
        ],
        "batch_size": 8,
        "max_workers": default_max_workers(),
        "strategy_mode": "random",
        "loss_mode": DEFAULT_LOSS_MODE,
        "aggregation_mode": DEFAULT_AGGREGATION_MODE,
        "aggregation_weights": list(DEFAULT_AGGREGATION_WEIGHTS),
        "time_window_aggregation_scope": DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
        "scm_override": "",
        "random_max_samples": 2000,
        "resolve_spacing": 0.1,
        "simann_max_iters": 200,
        "simann_initial_temp": 1.0,
        "simann_final_temp": 1.0e-12,
        "initial_case_rows": initial_case_rows,
        "initial_param_rows": initial_param_rows,
        "field_options": field_options,
        "right_pane_width_px": right_pane_width_px,
        "status": status,
        "status_text": status_text,
    }


def build_tab(app):
    """Build the tuning tab and register its callback groups."""
    initial_state = build_initial_tune_state()

    register_settings_callbacks(app)
    register_run_callbacks(app)
    register_display_callbacks(app)
    register_workspace_callbacks(app)

    return dcc.Tab(label="Tune", value="tune", children=build_layout(initial_state))
