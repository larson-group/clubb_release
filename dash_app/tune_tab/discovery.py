"""Discovery helpers for the Dash tuning tab."""

from __future__ import annotations

import os
import re

from dash_app.run_tab.discovery import load_available_cases
from dash_app.run_tab.namelist import read_namelist_entries
from dash_app.shared.tunable_configs import (
    available_tunable_configs as _available_tunable_configs,
    default_tunable_config_name as _default_tunable_config_name,
    tunable_config_names as _tunable_config_names,
    tunable_params_file_for_config,
)

from .state import REPO_ROOT


# Keep the historical Tune-discovery import surface while the canonical
# implementations live in shared configuration discovery.
available_tunable_configs = _available_tunable_configs
default_tunable_config_name = _default_tunable_config_name
tunable_config_names = _tunable_config_names

from tuner.case_defaults import (  # noqa: E402
    CASE_DEFAULTS_PATH,
    DEFAULT_LOSS_FIELDS,
    load_case_defaults as load_tuner_case_defaults,
)
from tuner.parameter_ranges import default_parameter_ranges  # noqa: E402
from utilities.benchmark_converter import supported_fields  # noqa: E402


TUNING_STATS_FILE = os.path.join(REPO_ROOT, "input", "stats", "all_tuning_stats.in")


def read_tuning_stats_fields():
    """Return CLUBB fields present in the tuner stats namelist."""
    if not os.path.isfile(TUNING_STATS_FILE):
        return []
    text = open(TUNING_STATS_FILE, "r", encoding="utf-8").read()
    fields = []
    for entry_text in re.findall(r'(?im)^\s*entry\(\d+\)\s*=\s*"([^"]+)"\s*$', text):
        field_name = entry_text.split("|", 1)[0].strip()
        if field_name:
            fields.append(field_name)
    return fields


def available_tuning_fields():
    """Return normalized fields the tuner can request and CLUBB can output."""
    return sorted(set(supported_fields()) & set(read_tuning_stats_fields()))


def _dash_case_defaults(case_name, defaults, fields):
    return {
        "case_name": case_name,
        "case_defaults_path": str(CASE_DEFAULTS_PATH),
        "les_stats_file": defaults["les_stats_file"],
        "default_clubb_fields": [name for name in DEFAULT_LOSS_FIELDS if name in fields],
        "clubb_fields": fields,
        "altitude_comparison_range": list(defaults["altitude_comparison_range"]),
        "time_average_range": list(defaults["time_average_range"]),
        "num_time_windows": int(defaults.get("num_time_windows", 1)),
    }


def load_case_defaults():
    """Return all tunable cases and their fixed tuning metadata."""
    case_defaults = {}
    all_defaults = load_tuner_case_defaults()
    fields = available_tuning_fields()
    for case_name in load_available_cases():
        defaults = all_defaults.get(case_name)
        if defaults is not None:
            case_defaults[case_name] = _dash_case_defaults(case_name, defaults, fields)
    return case_defaults


def load_tunable_names(config_name=None):
    """Return tunable parameter names from the selected tunable namelist."""
    return [entry["name"] for entry in read_namelist_entries(tunable_params_file_for_config(config_name))]


def _format_tune_range_value(value):
    """Format an auto-filled tune range endpoint compactly."""
    return f"{float(value):.6g}"


def load_tunable_default_ranges(config_name=None):
    """Return Fortran-bound-aware, UI-formatted default Tune ranges.

    The shared resolver uses the compiled Fortran hard-bound table whenever an
    endpoint exists and uses the legacy quarter-to-four-times-default envelope
    only on open sides.
    """
    return {
        name: {
            "default": _format_tune_range_value(item["default"]),
            "min": _format_tune_range_value(item["min"]),
            "max": _format_tune_range_value(item["max"]),
        }
        for name, item in default_parameter_ranges(config_name).items()
    }


def available_fields_for_case(case_defaults):
    """Return CLUBB fields that the normalized benchmark converter can provide."""
    return sorted(set((case_defaults or {}).get("clubb_fields", [])))


def available_fields_for_cases(case_names, case_data):
    """Return fields that are available for every selected case."""
    selected_cases = [name for name in (case_names or []) if name in (case_data or {})]
    if not selected_cases:
        return []

    field_sets = [
        set(available_fields_for_case((case_data or {}).get(case_name, {})))
        for case_name in selected_cases
    ]
    if not field_sets:
        return []
    common_fields = set.intersection(*field_sets)
    return sorted(common_fields)
