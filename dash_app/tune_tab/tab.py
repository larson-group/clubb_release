"""Thin entrypoint for the tuning tab package."""

from dash import dcc

from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_MODE

from .callbacks_display import register_display_callbacks
from .callbacks_runs import register_run_callbacks
from .callbacks_settings import register_settings_callbacks
from .discovery import available_fields_for_cases, load_case_defaults, load_tunable_default_ranges, load_tunable_names
from .layout import build_layout
from .runtime import empty_status_payload


def build_initial_tune_state():
    """Collect static metadata for the initial tuning-tab layout."""
    case_data = load_case_defaults()
    tunable_names = load_tunable_names()
    tunable_default_ranges = load_tunable_default_ranges()
    cases = sorted(case_data.keys())
    selected_cases = [cases[0]] if cases else []
    selected_case_defaults = case_data.get(selected_cases[0], {}) if selected_cases else {}
    field_options = available_fields_for_cases(selected_cases, case_data)
    max_name_len = max((len(name) for name in tunable_names), default=16)
    right_pane_width_px = max(360, min(760, int(180 + max_name_len * 7.5)))
    status = empty_status_payload()
    status_text = "state: idle | samples: 0 | best smart loss: --"
    initial_param_rows = [
        {"id": 0, "param": "C4", "min": "1", "max": "4"},
        {"id": 1, "param": "C8", "min": "0.1", "max": "0.9"},
        {"id": 2, "param": "C11", "min": "0.1", "max": "0.9"},
    ]
    initial_case_rows = []
    if selected_cases:
        time_range = selected_case_defaults.get("time_average_range", ["", ""])
        num_windows = int(selected_case_defaults.get("num_time_windows", 1) or 1)
        average_time_seconds = ""
        if len(time_range) > 1 and num_windows > 0:
            average_time_seconds = int((int(time_range[1]) - int(time_range[0])) / num_windows)
        initial_case_rows.append(
            {
                "id": 0,
                "case_name": selected_cases[0],
                "time_start": time_range[0],
                "time_end": time_range[1],
                "average_time_seconds": average_time_seconds,
            }
        )

    return {
        "cases": cases,
        "case_data": case_data,
        "tunable_names": tunable_names,
        "tunable_default_ranges": tunable_default_ranges,
        "selected_cases": selected_cases,
        "selected_fields": [
            field_name
            for field_name in selected_case_defaults.get("default_clubb_fields", [])
            if field_name in field_options
        ],
        "batch_size": 8,
        "max_workers": 10,
        "strategy_mode": None,
        "loss_mode": DEFAULT_LOSS_MODE,
        "aggregation_mode": DEFAULT_AGGREGATION_MODE,
        "random_max_samples": 100,
        "resolve_spacing": 0.1,
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

    return dcc.Tab(label="Tune", value="tune", children=build_layout(initial_state))
