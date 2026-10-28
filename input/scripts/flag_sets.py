"""Helpers for running CLUBB under the flag sets described by a JSON config file.

A flag config maps flag-set names to namelist overrides, for example::

    {"flag1": {"saturation_formula": 2, "l_diffuse_rtm_and_thlm": true}}

Shared by run_scripts/run_clubb_w_varying_flags.py and
tests/run_jax_vs_fortran_cases.py so both render overrides identically.
"""

import json

# Name of the injected flag set that runs with no overrides at all.
DEFAULT_FLAG_SET_NAME = "default"


def read_flag_settings(path):
    """Load the JSON mapping from flag-set name to override values."""
    if not str(path).endswith(".json"):
        raise ValueError(f"Flag config file must be a JSON file: {path}")

    with open(path) as f:
        return json.load(f)


def get_flag_sets(skip_default, flag_dict):
    """Normalize JSON flag sets and optionally inject the unmodified default run."""
    flag_sets = {}

    if not skip_default:
        flag_sets[DEFAULT_FLAG_SET_NAME] = None

    for flag_set_name, overrides in flag_dict.items():
        if flag_set_name == DEFAULT_FLAG_SET_NAME:
            raise ValueError(f"'{DEFAULT_FLAG_SET_NAME}' may not be used as a flag set name.")
        flag_sets[flag_set_name] = overrides

    return flag_sets


def format_override_value(value):
    """Render Python values into Fortran-friendly override strings."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    return str(value)


def build_override_arg(overrides):
    """Serialize one flag set into the comma-delimited -override format."""
    if not overrides:
        return None
    return ",".join(
        f"{key}={format_override_value(value)}"
        for key, value in overrides.items()
    )
