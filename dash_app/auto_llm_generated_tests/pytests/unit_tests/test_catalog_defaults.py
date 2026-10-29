"""Provisional Dash baseline selection check, owned by clubb_dash."""

from dash_app.shared.tunable_configs import (
    available_tunable_configs,
    config_settings,
    default_tunable_config_name,
)
from input.scripts import create_case_namelist as inputs


def test_dash_offers_baseline_without_a_default_preset(monkeypatch):
    """Keep baseline selectable when only a non-default parameter preset exists."""
    catalogs = {
        "parameters": {"defaults": {"model": {"C8": 0.5}, "silhs": {}},
                       "presets": {"experiment": {"model": {"C8": 0.2}}}},
        "flags": {"defaults": {"model": {"enabled": False}, "silhs": {}}, "presets": {}},
        "cases": {"defaults": {}, "presets": {}},
    }
    monkeypatch.setattr(inputs, "load_catalog", lambda kind, json_file=None: catalogs[kind])
    configs = available_tunable_configs()
    assert "default" in {config["value"] for config in configs}
    assert default_tunable_config_name(configs) == "default"
    assert config_settings("default", "parameters") == catalogs["parameters"]["defaults"]
    assert config_settings("default", "flags") == catalogs["flags"]["defaults"]
