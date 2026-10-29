"""Provisional checks for selecting catalog defaults without empty presets."""

import json

import pytest

from input.scripts import create_case_namelist as inputs
from input.scripts.flags_json_to_namelist import flags_json_to_namelist, flags_to_namelist
from input.scripts.parameters_json_to_namelist import parameters_json_to_namelist, parameters_to_namelist
from input.scripts.save_tunable_config import save_tunable_config


@pytest.fixture
def catalogs(tmp_path, monkeypatch):
    """Provide three override catalogs with a parameter-only preset and no default preset."""
    values = {
        "cases": {"model_setting": {"dt_main": 30}},
        "parameters": {"model": {"C8": 0.5}, "silhs": {}},
        "flags": {"model": {"enabled": False}, "silhs": {}},
    }
    for kind, filename in inputs.CATALOG_FILES.items():
        presets = {"experiment": {"model": {"C8": 0.2}}} if kind == "parameters" else {}
        (tmp_path / filename).write_text(json.dumps({"defaults": values[kind], "presets": presets}))
    monkeypatch.setattr(inputs, "INPUT_ROOT", tmp_path)
    return tmp_path, values


def test_shared_default_needs_no_preset_and_unknown_names_still_fail(catalogs):
    """Default resolves every baseline; named configs apply only where defined."""
    _, defaults = catalogs
    assert "default" in inputs.available_configs()
    for kind, expected in defaults.items():
        assert inputs.resolve_input_settings(kind, config="default") == expected
        with pytest.raises(ValueError, match="Unknown config"):
            inputs.resolve_input_settings(kind, config="typo")
    assert inputs.resolve_input_settings("parameters", config="experiment")["model"]["C8"] == 0.2
    assert inputs.resolve_input_settings("flags", config="experiment") == defaults["flags"]
    with pytest.raises(ValueError, match="Unknown flags preset"):
        inputs.resolve_input_settings("flags", "experiment")


@pytest.mark.parametrize("kind,convert,serialize", [
    ("parameters", parameters_json_to_namelist, parameters_to_namelist),
    ("flags", flags_json_to_namelist, flags_to_namelist),
])
def test_converter_omitted_selection_uses_catalog_defaults(catalogs, kind, convert, serialize):
    """Convenience converters serialize defaults without requesting a named preset."""
    root, defaults = catalogs
    assert convert(json_file=root / inputs.CATALOG_FILES[kind]) == serialize(defaults[kind])


def test_save_from_builtin_default_and_protect_its_name(catalogs):
    """Save edits from the baseline; replacing its built-in name must not write files."""
    root, defaults = catalogs
    save_tunable_config("default", {"tunable": {"C8": "0.4"}}, "saved", config_root=root)
    assert inputs.resolve_input_settings("parameters", config="saved")["model"]["C8"] == 0.4
    assert inputs.resolve_input_settings("flags", config="saved") == defaults["flags"]
    paths = [root / filename for filename in inputs.CATALOG_FILES.values()]
    before = [path.read_text() for path in paths]
    with pytest.raises(ValueError, match="reserved"):
        save_tunable_config("saved", {}, "default", force=True, config_root=root)
    assert [path.read_text() for path in paths] == before
