"""Build run-tab tunable-config dependent UI metadata."""

from tunable_configs import (
    available_tunable_configs,
    default_tunable_config_name,
    tunable_config_file,
    tunable_config_names,
)

from .layout import build_flag_controls, build_param_sections, compute_width_hints
from .namelist import is_bool_value, is_true, normalize_numeric_display, read_namelist_entries


def _parse_tunable_default(value):
    """Parse a scalar Fortran namelist value as a float."""
    try:
        return float(str(value).strip().replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return None


def _format_range_value(value):
    """Format an auto-filled range endpoint compactly."""
    return f"{float(value):.6g}"


def _tunable_default_ranges(tunable_entries):
    """Return default min/max ranges derived from tunable parameter defaults."""
    ranges = {}
    for entry in tunable_entries:
        default_value = _parse_tunable_default(entry.get("value"))
        if default_value is None:
            continue
        low = default_value / 4.0
        high = default_value * 4.0
        if low > high:
            low, high = high, low
        ranges[entry["name"]] = {
            "default": _format_range_value(default_value),
            "min": _format_range_value(low),
            "max": _format_range_value(high),
        }
    return ranges


def _select_config_name(config_name=None, configs=None):
    """Return a valid config name for the available config list."""
    configs = list(configs or available_tunable_configs())
    selected = str(config_name or "").strip()
    if selected in tunable_config_names(configs):
        return selected
    return default_tunable_config_name(configs)


def build_tunable_config_state(config_name=None, configs=None):
    """Collect config-dependent run-tab defaults, controls, and parameter metadata."""
    selected_config = _select_config_name(config_name, configs)
    flag_entries = read_namelist_entries(tunable_config_file(selected_config, "configurable_model_flags.in"))
    flag_bools = [entry for entry in flag_entries if is_bool_value(entry["value"])]
    flag_params = [entry for entry in flag_entries if entry not in flag_bools]
    tunable_entries = read_namelist_entries(tunable_config_file(selected_config, "tunable_parameters.in"))
    silhs_entries = read_namelist_entries(tunable_config_file(selected_config, "silhs_parameters.in"))

    flag_names = [entry["name"] for entry in flag_bools]
    all_config_names = (
        flag_names
        + [entry["name"] for entry in flag_params]
        + [entry["name"] for entry in tunable_entries]
        + [entry["name"] for entry in silhs_entries]
    )
    label_width_px, value_width_px, right_pane_width_px = compute_width_hints(all_config_names)

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
        "selected_config": selected_config,
        "defaults": defaults,
        "defaults_by_key": defaults_by_key,
        "flag_names": flag_names,
        "tunable_names": [entry["name"] for entry in tunable_entries],
        "tunable_default_ranges": _tunable_default_ranges(tunable_entries),
        "param_meta": [{"file": entry["file"], "name": entry["name"]} for entry in param_entries],
        "right_pane_width_px": right_pane_width_px,
        "param_sections": build_param_sections(
            flag_params,
            flag_controls,
            tunable_entries,
            silhs_entries,
            label_width_px,
            value_width_px,
            normalize_numeric_display,
        ),
    }
