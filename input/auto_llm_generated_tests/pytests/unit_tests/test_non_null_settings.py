"""Provisional input contracts for explicit values and the legacy tuner boundary."""

import copy
import json
from pathlib import Path

import pytest

from input.scripts.create_case_namelist import load_catalog, resolve_input_settings
from input.scripts.namelist import merge_settings, namelists_to_text, parse_value, settings_patch
from input.scripts.stats_json_to_namelist import select_stats_entries


@pytest.mark.parametrize("kind,patch", [
    ("cases", {"model_setting": {"dt_main": None}}),
    ("cases", {"model_setting": {"sclr_tol_nl": [0.01, None]}}),
    ("parameters", {"model": {"C8": None}}),
    ("flags", {"model": {"l_use_cloud_cover": None}}),
    ("stats", [None]),
])
def test_explicit_input_files_reject_null_before_emission(tmp_path, kind, patch):
    """A custom override cannot delete a field or hide nulls inside a list."""
    path = tmp_path / "invalid.jsonc"
    path.write_text(json.dumps(patch))
    with pytest.raises(ValueError, match="Null is not supported"):
        if kind == "stats":
            select_stats_entries(path)
        else:
            resolve_input_settings(kind, path)
    assert not list(tmp_path.glob("*.in"))


@pytest.mark.parametrize("section", ["defaults", "presets"])
def test_catalog_rejects_null_even_in_an_unselected_preset(tmp_path, section):
    """Reject invalid shipped definitions rather than accepting latent deletions."""
    catalog = {"defaults": {"model_setting": {"dt_main": 60}},
               "presets": {"good": {}, "bad": {}}}
    if section == "defaults":
        catalog["defaults"]["model_setting"]["dt_main"] = None
    else:
        catalog["presets"]["bad"] = {"model_setting": {"dt_main": None}}
    path = tmp_path / "case_definitions.jsonc"
    path.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="Null is not supported"):
        load_catalog("cases", path)


def test_sparse_save_rejects_deleting_an_inherited_field_without_mutation():
    """Saving cannot manufacture the removed null deletion representation."""
    defaults = {"group": {"enabled": False, "value": 2}}
    before = copy.deepcopy(defaults)
    with pytest.raises(ValueError, match=r"Cannot omit inherited setting settings.group.value"):
        settings_patch(defaults, {"group": {"enabled": True}})
    assert defaults == before
    settings = {"group": {"enabled": True, "value": 0, "extra": []}}
    assert merge_settings(defaults, settings_patch(defaults, settings)) == settings


@pytest.mark.parametrize("settings", [
    {"group": {"value": None}},
    {"group": {"value": [1, None]}},
    {"group": {"_start_index": {"value": [None]}}},
])
def test_direct_namelist_emission_rejects_null_settings(settings):
    """Direct callers receive the same error instead of empty namelist assignments."""
    with pytest.raises(ValueError, match="Null is not supported"):
        namelists_to_text(settings)


def test_literal_null_is_rejected_but_quoted_text_is_preserved():
    """CLI null cannot become an implicit reset; a literal string remains usable."""
    with pytest.raises(ValueError, match="Null is not supported"):
        parse_value("null")
    assert parse_value('"null"') == "null"


def test_optional_provenance_is_separate_from_input_settings(tmp_path):
    """Nullable provenance neither enables deletion nor reaches emitted model data."""
    path = tmp_path / "case.json"
    path.write_text(json.dumps({"_metadata": {"benchmark_case": None},
                                "model_setting": {"dt_main": 30}}))
    settings = resolve_input_settings("cases", path)
    assert settings["model_setting"]["dt_main"] == 30
    assert "_metadata" not in settings


def test_saved_shared_config_can_have_only_a_parameter_preset(tmp_path):
    """A missing flag preset inherits defaults when saving a shared configuration."""
    from input.scripts.save_tunable_config import save_tunable_config
    from input.scripts.flags_json_to_namelist import flag_settings
    from input.scripts.parameters_json_to_namelist import parameter_settings

    (tmp_path / "clubb_params.jsonc").write_text(json.dumps({
        "defaults": {"model": {"C8": 0.5}, "silhs": {}},
        "presets": {"source": {"model": {"C8": 0.2}}},
    }))
    (tmp_path / "clubb_config_flags.jsonc").write_text(json.dumps({
        "defaults": {"model": {"enabled": False}, "silhs": {}}, "presets": {},
    }))
    save_tunable_config("source", {"tunable": {"C8": "0.4"}}, "saved", config_root=tmp_path)
    assert parameter_settings("saved", tmp_path / "clubb_params.jsonc")["model"]["C8"] == 0.4
    assert flag_settings("saved", tmp_path / "clubb_config_flags.jsonc")["model"]["enabled"] is False
