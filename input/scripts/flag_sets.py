"""Resolve flag sweeps from the shared catalog and render late run overrides.

read_flag_settings() expands a catalog group, list or preset. get_flag_sets()
adds the unmodified default run unless explicitly skipped. build_override_arg()
and format_override_value() preserve the override interface used by comparison
runners, where sweep settings apply after curated case settings and before user
edits. Collection member names remain local; collisions across collections use
qualified output labels.
"""
from input.scripts.input_catalog import catalog_presets, load_catalog, resolve_selection

DEFAULT_FLAG_SET_NAME = "default"
DEFAULT_FLAG_GROUP = "run_bindiff_w_flags_config_core_flags_group"


def read_flag_settings(selection):
    """Expand a catalog selector into independently applied, named sparse patches."""
    names = resolve_selection("flags", selection)
    patches = catalog_presets("flags", load_catalog("flags"))
    local_names = [name.rsplit("/", 1)[-1] for name in names]
    result = {
        (name.replace("/", "__") if local_names.count(local) > 1 else local):
        {key: value for category in ("model", "silhs")
         for key, value in patches[name].get(category, {}).items()}
        for name, local in zip(names, local_names)
    }
    if len(result) != len(names):
        raise ValueError("Selected flag members have colliding output labels")
    from input.scripts.flags_json_to_namelist import flag_settings
    for name in names:
        flag_settings(name)
    return result


def get_flag_sets(skip_default, flag_dict):
    """Optionally prepend an unmodified baseline to the selected independent runs."""
    if DEFAULT_FLAG_SET_NAME in flag_dict:
        raise ValueError("'default' may not be used as a flag set name")
    return {**({} if skip_default else {DEFAULT_FLAG_SET_NAME: None}), **flag_dict}


def format_override_value(value):
    """Render JSON booleans and numeric selectors for direct namelist overrides."""
    return (".true." if value else ".false.") if isinstance(value, bool) else str(value)


def build_override_arg(overrides):
    """Render one sparse sweep patch, or None for the unmodified baseline."""
    return ",".join(f"{key}={format_override_value(value)}" for key, value in overrides.items()) if overrides else None
