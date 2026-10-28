#!/usr/bin/env python3
"""Resolve and serialize unified model/SILHS parameter JSON.

parameter_settings(selection, json_file=None, config=None) reads defaults,
shared-config settings and an explicit preset or override file, then validates
numeric types and column dimensions. parameters_to_namelist(settings) serializes
final model/SILHS groups for create_case_namelist.py. The convenience function
parameters_json_to_namelist() resolves and serializes a selection in memory.

Custom files have the same model/silhs override objects as preset entries.
Multicolumn snapshots may also supply multicol dimensions. Omitted keys inherit
the selected base; null values are rejected. Physical bounds remain owned by CLUBB.
"""
from __future__ import annotations

import math
from pathlib import Path
import re

from input.scripts.namelist import load_json
from input.scripts.namelist import namelists_to_text

DEFAULT_PARAMETERS_JSON = Path(__file__).resolve().parents[1] / "clubb_params.jsonc"
SETTING_NAME = re.compile(r"[A-Za-z_]\w*\Z")


def parameter_settings(selection=None, json_file=None, *, config=None):
    """Resolve numeric model/SILHS parameters and optional multicolumn dimensions."""
    from input.scripts.create_case_namelist import resolve_input_settings
    settings = resolve_input_settings("parameters", selection, config=config, json_file=json_file)
    if set(settings) - {"model", "silhs", "multicol"} or not {"model", "silhs"} <= settings.keys():
        raise ValueError("Parameters JSON requires model and silhs objects; only multicol is optional")
    for category in ("model", "silhs"):
        if not isinstance(settings[category], dict):
            raise ValueError(f"Parameter {category} must be an object")
        for name, value in settings[category].items():
            values = value if isinstance(value, list) else [value]
            if not SETTING_NAME.fullmatch(name) or not values or any(
                isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item)
                for item in values
            ):
                raise ValueError(f"Parameter {category}.{name} requires finite numbers")
    if "multicol" in settings:
        dims = settings["multicol"]
        if not isinstance(dims, dict) or set(dims) - {"ngrdcol", "batch_size"} or "ngrdcol" not in dims:
            raise ValueError("multicol requires ngrdcol and optional batch_size")
        if any(type(value) is not int or value < 1 for value in dims.values()):
            raise ValueError("multicol dimensions must be positive integers")
        if dims.get("batch_size", 1) > dims["ngrdcol"]:
            raise ValueError("batch_size cannot exceed ngrdcol")
        for value in settings["model"].values():
            if isinstance(value, list) and len(value) != dims["ngrdcol"]:
                raise ValueError("Parameter column arrays must match ngrdcol")
    elif any(isinstance(value, list) for value in settings["model"].values()):
        raise ValueError("Parameter column arrays require multicol dimensions")
    return settings


def parameters_to_namelist(settings):
    """Emit model and SILHS groups, preceded by multicol dimensions when present."""
    groups = {"multicol_def": settings["multicol"]} if "multicol" in settings else {}
    groups.update(clubb_params_nl=settings["model"], silhs_params_nl=settings["silhs"])
    return namelists_to_text(groups)


def parameters_json_to_namelist(config=None, json_file=None):
    """Serialize a parameter preset or snapshot; an omitted config uses defaults."""
    return parameters_to_namelist(parameter_settings(config, json_file))
