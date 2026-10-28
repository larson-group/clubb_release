#!/usr/bin/env python3
"""Resolve and serialize unified model/SILHS flag JSON.

flag_settings(selection, json_file=None, config=None) applies defaults, shared
config settings and an explicit preset/file, then checks boolean/integer types.
flags_to_namelist(settings) serializes final configurable flag groups for
create_case_namelist.py. flags_json_to_namelist() resolves and serializes a
selection in memory. Custom files contain the same model/silhs override objects
as preset entries. CLUBB owns valid combinations of physics options.
"""
from __future__ import annotations

from pathlib import Path
import re

from input.scripts.namelist import load_json
from input.scripts.namelist import namelists_to_text

DEFAULT_FLAGS_JSON = Path(__file__).resolve().parents[1] / "clubb_config_flags.jsonc"


def flag_settings(selection=None, json_file=None, *, config=None):
    """Resolve both flag categories, requiring booleans or integer selectors."""
    from input.scripts.create_case_namelist import resolve_input_settings
    settings = resolve_input_settings("flags", selection, config=config, json_file=json_file)
    if set(settings) != {"model", "silhs"}:
        raise ValueError("Flags JSON requires only model and silhs objects")
    for category, values in settings.items():
        if not isinstance(values, dict):
            raise ValueError(f"Flag {category} must be an object")
        for name, value in values.items():
            if not re.fullmatch(r"[A-Za-z_]\w*", name) or type(value) not in (bool, int):
                raise ValueError(f"Flag {category}.{name} requires a boolean or integer")
    return settings


def flags_to_namelist(settings):
    """Emit both model and SILHS configurable flag groups."""
    return namelists_to_text({"configurable_clubb_flags_nl": settings["model"],
                              "configurable_silhs_flags_nl": settings["silhs"]})


def flags_json_to_namelist(config=None, json_file=None):
    """Serialize a flag preset or snapshot; an omitted config uses catalog defaults."""
    return flags_to_namelist(flag_settings(config, json_file))
