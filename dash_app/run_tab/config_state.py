"""Build run-tab tunable-config dependent UI metadata."""

from dash_app.shared.tunable_configs import (
    available_tunable_configs,
    default_tunable_config_name,
    tunable_config_file,
    tunable_config_names,
)
from utilities.clubb_settings_validation import (
    build_settings_schema,
    evaluate_settings,
    is_independently_tunable,
)

from .layout import build_flag_controls, build_param_sections, compute_width_hints
from .namelist import is_bool_value, is_true, normalize_numeric_display, read_namelist_entries


def _format_range_value(value):
    """Format an auto-filled range endpoint compactly."""
    return f"{float(value):.6g}"


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
    settings_schema = build_settings_schema(
        defaults["flags"],
        defaults["params"],
        [{"file": entry["file"], "name": entry["name"]} for entry in param_entries],
    )
    settings_resolution = evaluate_settings(settings_schema)
    schema_bounds = settings_schema["hard_bounds"]
    active_flag_relationships = [
        relation
        for relation in settings_resolution["flag_relationships"]
        if all(member in set(flag_names) for member in relation["members"])
    ]
    flag_controls = build_flag_controls(flag_bools, is_true, active_flag_relationships)
    active_linked_groups = [
        list(group)
        for group in settings_resolution["linked_parameter_groups"]
        if all(name in {entry["name"] for entry in tunable_entries} for name in group)
    ]
    tunable_control_entries = [
        {
            **entry,
            "disabled": not is_independently_tunable(
                settings_resolution["parameter_states"].get(entry["name"])
            ),
        }
        for entry in tunable_entries
    ]

    return {
        "selected_config": selected_config,
        "tunable_names": [entry["name"] for entry in tunable_entries],
        "tunable_default_ranges": {
            name: {key: _format_range_value(value) for key, value in item.items()}
            for name, item in settings_resolution["tunable_default_ranges"].items()
        },
        "linked_parameter_groups": active_linked_groups,
        "settings_resolution": settings_resolution,
        "settings_schema": settings_schema,
        "tunable_entries": tunable_entries,
        "param_meta": [
            {
                "file": entry["file"],
                "name": entry["name"],
                **(
                    {
                        "min": schema_bounds.get(entry["name"], {}).get("min"),
                        "max": schema_bounds.get(entry["name"], {}).get("max"),
                    }
                    if entry["file"] == "tunable"
                    else {}
                ),
            }
            for entry in param_entries
        ],
        "right_pane_width_px": right_pane_width_px,
        "param_sections": build_param_sections(
            flag_params,
            flag_controls,
            tunable_control_entries,
            active_linked_groups,
            silhs_entries,
            label_width_px,
            value_width_px,
            normalize_numeric_display,
        ),
    }
