"""Shared tunable-parameter config discovery for Dash tabs."""

from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
TUNABLE_CONFIG_ROOT = REPO_ROOT / "input"
TUNABLE_PARAMETER_RENAMES: dict[str, str] = {}


def available_tunable_configs():
    """List the built-in baseline and parameter/flag preset names."""
    from input.scripts.create_case_namelist import DEFAULT_CONFIG, load_catalog
    names = {DEFAULT_CONFIG} | set(load_catalog("parameters")["presets"]) | set(load_catalog("flags")["presets"])
    return [{"label": name, "value": name, "path": str(TUNABLE_CONFIG_ROOT)}
            for name in sorted(names, key=lambda name: (name != "default", name))]


def config_settings(config_name, kind):
    """Return resolved JSON for a shared config and input kind."""
    from input.scripts.create_case_namelist import resolve_input_settings
    return resolve_input_settings(kind, config=config_name or "default")


def config_entries(config_name, section):
    """Return editable entries from parameter/flag JSON, preserving catalog order."""
    from input.scripts.namelist import setting_entries
    kind = "flags" if section == "flags" else "parameters"
    values = config_settings(config_name, kind)
    values = {**values["model"], **values["silhs"]} if section == "flags" else values["silhs" if section == "silhs" else "model"]
    return setting_entries(values)


def default_tunable_config_name(configs=None):
    """Return the default selected config name from a config option list."""
    configs = list(configs or available_tunable_configs())
    values = [str(config.get("value", "")).strip() for config in configs]
    if "default" in values:
        return "default"
    return values[0] if values else "default"


def tunable_config_names(configs=None):
    """Return selectable config names from a config option list."""
    return {
        str(config.get("value", "")).strip()
        for config in (configs or available_tunable_configs())
        if str(config.get("value", "")).strip()
    }


def canonical_tunable_parameter_name(name, available_names):
    """Resolve one documented historic name when its successor is available.

    This never invents a parameter.  The selected configuration remains the
    authority: an alias is accepted only when its declared successor appears
    in that configuration's current namelist.
    """
    candidate = str(name or "").strip()
    available = set(available_names or [])
    if candidate in available:
        return candidate
    successor = TUNABLE_PARAMETER_RENAMES.get(candidate)
    if successor in available:
        return successor
    return candidate
