"""Shared tunable-parameter config discovery for Dash tabs."""

from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
TUNABLE_CONFIG_ROOT = REPO_ROOT / "input" / "parameter_and_flag_configs"
REQUIRED_TUNABLE_CONFIG_FILES = (
    "tunable_parameters.in",
    "configurable_model_flags.in",
    "silhs_parameters.in",
)

TUNABLE_PARAMETER_RENAMES: dict[str, str] = {}

from utilities.create_case_namelist import resolve_tunable_config_dir  # noqa: E402


def _complete_tunable_config(path: Path) -> bool:
    """Return whether a tunable config directory has all files needed by run_scm."""
    return path.is_dir() and all((path / filename).is_file() for filename in REQUIRED_TUNABLE_CONFIG_FILES)


def available_tunable_configs():
    """Return complete named tunable configs available to Dash tabs."""
    configs = []
    if TUNABLE_CONFIG_ROOT.is_dir():
        for path in sorted(TUNABLE_CONFIG_ROOT.iterdir(), key=lambda item: item.name):
            if _complete_tunable_config(path):
                configs.append({"label": path.name, "value": path.name, "path": str(path.resolve())})
    configs.sort(key=lambda item: (item["value"] != "default", item["value"]))
    return configs


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


def tunable_config_dir(config_name=None):
    """Return the resolved config directory for a config name or directory."""
    return Path(resolve_tunable_config_dir(config_name or "default"))


def tunable_config_file(config_name=None, filename="tunable_parameters.in"):
    """Return a file path inside the resolved config directory."""
    return str(tunable_config_dir(config_name) / os.fspath(filename))


def tunable_params_file_for_config(config_name=None):
    """Return the tunable parameter file for a config name or directory."""
    return tunable_config_file(config_name, "tunable_parameters.in")


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
