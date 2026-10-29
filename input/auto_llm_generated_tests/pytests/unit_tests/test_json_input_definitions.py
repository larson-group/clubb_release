"""Provisional input-boundary contracts, owned by the shared Python Jenkins unit stage.

Synthetic catalogs test ordered sparse overrides, copyable files, strict JSON
parsing, and final emission without running any model cases.
"""
import copy
import json
import re
from pathlib import Path

import pytest

from input.scripts import create_case_namelist as builder
from input.scripts.parameters_json_to_namelist import parameter_settings
from input.scripts.flags_json_to_namelist import flag_settings
from input.scripts.stats_json_to_namelist import select_stats_entries


@pytest.fixture
def catalogs(tmp_path, monkeypatch):
    """Provide case, parameter and flag catalogs with overlapping sparse presets."""
    values = {
        "parameters": {"defaults": {"model": {"C8": 0.5, "C11": 0.4}, "silhs": {"vert_decorr_coef": 0.03}},
                       "presets": {"experiment": {"model": {"C8": 0.2}, "silhs": {"vert_decorr_coef": 0.1}},
                                   "specific": {"model": {"C11": 0.8}}}},
        "flags": {"defaults": {"model": {"first": False, "second": False}, "silhs": {"sample": False}},
                  "presets": {"experiment": {"model": {"first": True}}, "specific": {"silhs": {"sample": True}}}},
        "cases": {"defaults": {"model_setting": {"dt_main": 60, "time_initial": 0, "time_final": 600},
                               "stats_setting": {"l_stats": True}},
                  "presets": {"experiment": {"model_setting": {"dt_rad": 30}},
                              "arm": {"model_setting": {"runtype": "arm", "dt_main": 10}}}},
    }
    for kind, value in values.items():
        (tmp_path / builder.CATALOG_FILES[kind]).write_text(json.dumps(value))
    monkeypatch.setattr(builder, "INPUT_ROOT", tmp_path)
    return tmp_path


@pytest.mark.parametrize("kind", ["parameters", "flags", "cases"])
def test_custom_file_and_named_preset_are_identical_sparse_overrides(catalogs, kind):
    """Copy a preset into a file and require identical config-before-section precedence."""
    catalog = json.loads((catalogs / builder.CATALOG_FILES[kind]).read_text())
    name = "arm" if kind == "cases" else "specific"
    custom = catalogs / "custom.json"
    custom.write_text(json.dumps(catalog["presets"][name]))
    before = sorted(catalogs.iterdir())
    named = builder.resolve_input_settings(kind, name, config="experiment")
    saved = builder.resolve_input_settings(kind, str(custom), config="experiment")
    assert named == saved
    assert sorted(catalogs.iterdir()) == before
    if kind == "parameters":
        assert saved == {"model": {"C8": 0.2, "C11": 0.8}, "silhs": {"vert_decorr_coef": 0.1}}
    elif kind == "flags":
        assert saved == {"model": {"first": True, "second": False}, "silhs": {"sample": True}}
    elif kind == "cases":
        assert saved["model_setting"]["dt_main"] == 10
        assert saved["model_setting"]["dt_rad"] == 30


def test_config_missing_from_one_section_preserves_its_defaults(catalogs):
    """A shared name may exist only in cases; explicit unknown section names still fail."""
    settings = builder.resolve_input_settings("parameters", config="arm")
    assert settings["model"]["C8"] == 0.5
    with pytest.raises(ValueError, match="Unknown parameters preset"):
        builder.resolve_input_settings("parameters", "arm")
    with pytest.raises(ValueError, match="Unknown config"):
        builder.resolve_input_settings("parameters", config="typo")


@pytest.mark.parametrize("body", [
    '{"model":{"C8":0.2,"C8":0.8}}',
    '{"model":{"C8":NaN}}',
    '{"model":{"C8":true}}',
    '{"model":{"C8":[0.2,0.8]}}',
    '{"model":{"C8":[0.2,0.8]},"multicol":{"ngrdcol":3}}',
])
def test_invalid_parameter_overrides_fail_before_emission(tmp_path, body):
    """Reject duplicate keys, non-finite values, wrong types and inconsistent column arrays."""
    path = tmp_path / "invalid.json"
    path.write_text(body)
    with pytest.raises(ValueError):
        parameter_settings(path)


def test_stats_custom_list_uses_catalog_metadata_and_preserves_list_order(tmp_path):
    """Only names need copying; deduplication keeps the first explicit occurrence."""
    path = tmp_path / "selected.json"
    path.write_text('["rtm", "thlm", "rtm"]')
    entries = select_stats_entries(str(path))
    assert [entry["name"] for entry in entries] == ["rtm", "thlm"]
    assert entries[0]["units"] == "kg/kg"
    path.write_text('["typo"]')
    with pytest.raises(ValueError, match="Unknown stats variable"):
        select_stats_entries(str(path))


def test_only_final_builder_emits_a_namelist(tmp_path):
    """Resolve all sections in memory, then emit exactly one complete runner input."""
    settings = builder.create_case_settings("bomex", tmp_path, stats="core", override="dt_main=10")
    assert not list(tmp_path.rglob("*.in"))
    assert settings["model_setting"]["dt_main"] == 10
    path = builder.write_case_namelist(settings, tmp_path / "bomex.in")
    assert list(tmp_path.rglob("*.in")) == [path]
    text = path.read_text()
    for group in ("clubb_params_nl", "silhs_params_nl", "configurable_clubb_flags_nl", "model_setting", "clubb_stats_nl"):
        assert "&" + group in text


def test_case_file_can_be_pasted_directly_into_catalog(tmp_path):
    """A case file and a named case use precisely the same override object."""
    from input.scripts.case_json_to_namelist import case_definition_record
    catalog = {"defaults": {"model_setting": {"dt_main": 60}}, "presets": {}}
    patch = {"model_setting": {"runtype": "bomex", "dt_main": 30}, "stats_setting": {"l_stats": False}}
    path = tmp_path / "case_definitions.jsonc"
    path.write_text(json.dumps(catalog))
    custom = tmp_path / "my_case.json"
    custom.write_text(json.dumps(patch))
    saved = case_definition_record(str(custom), path)
    catalog["presets"]["copied"] = patch
    path.write_text(json.dumps(catalog))
    assert saved["namelists"] == case_definition_record("copied", path)["namelists"]


def test_custom_paths_expand_home_and_remain_distinct_from_preset_names(tmp_path, monkeypatch):
    """Resolve explicit home-relative, relative and absolute override files consistently."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    custom = tmp_path / "clubb_config_flags.jsonc"
    custom.write_text('{"silhs":{"l_lh_straight_mc":true}}')
    expected = flag_settings(str(custom))
    assert flag_settings("~/clubb_config_flags.jsonc") == expected
    assert flag_settings("./clubb_config_flags.jsonc") == expected
    with pytest.raises(FileNotFoundError):
        flag_settings(str(tmp_path / "missing.json"))


def test_native_tuner_bounds_cover_the_current_fortran_reader():
    """Every active bounds group names and initializes exactly the native tunables."""
    source = (builder.INPUT_ROOT.parent / "src/CLUBB_core/parameters_tunable.F90").read_text()
    declaration = source.split("namelist /init_minmax/", 1)[1].split("\n\n", 1)[0]
    supported = set(re.findall(r"\b\w+_minmax\b", declaration.lower()))
    from clubb_python_driver.io.namelist import _read_namelist_groups
    controls = builder.INPUT_ROOT.parent / "input_misc/tuner"
    paths = list(controls.glob("error_*.in"))
    assert paths, "Native tuner controls must be present for bounds validation"
    for path in paths:
        settings = _read_namelist_groups(path)
        bounds = settings.get("init_minmax")
        if bounds is not None:
            assert set(bounds) == supported, path.name
