"""Provisional coverage for editable case catalogs and immutable copies."""
import copy
import json

import pytest

from utilities import case_json_to_namelist as cases


@pytest.fixture
def catalog(tmp_path, monkeypatch):
    path = tmp_path / "case_definitions.json"
    defaults = {"model_setting": {"runtype": "bomex", "time_initial": 0.0, "time_final": 21600.0,
                                 "dt_main": 60.0, "dt_rad": 60.0,
                                 "thlm_sponge_damp_settings": {"l_sponge_damping": False, "tau_sponge_damp_min": 60.0}},
                "stats_setting": {"fname_prefix": "source", "l_stats": True}}
    path.write_text(json.dumps({"defaults": defaults, "cases": {"source": {}, "other": {}}}))
    monkeypatch.setattr(cases, "DEFAULT_CASES_JSON", path)
    return path


def test_merge_preserves_omitted_fields_and_replaces_arrays():
    defaults = {"model_setting": {"sclr_tol_nl": [0.01, 1e-8], "sfctype": 1,
                                 "thlm_sponge_damp_settings": {"l_sponge_damping": False, "tau_sponge_damp_min": 60.0}},
                "setfields": {"input_type": "SAM"}}
    patch = {"model_setting": {"sclr_tol_nl": [1e-9], "sfctype": None,
                              "thlm_sponge_damp_settings": {"l_sponge_damping": True}}, "setfields": None}
    before = copy.deepcopy(defaults)
    actual = cases.merge_case_settings(defaults, patch)
    assert actual == {"model_setting": {"sclr_tol_nl": [1e-9],
                                       "thlm_sponge_damp_settings": {"l_sponge_damping": True, "tau_sponge_damp_min": 60.0}}}
    assert defaults == before
    assert cases.merge_case_settings(defaults, cases.case_settings_patch(defaults, actual)) == actual


def test_serializer_keeps_quotes_comment_characters_and_fortran_array_indices():
    text = cases.namelists_to_text({"model_setting": {"label": 'a "quote" ! /', "l_restart": False,
        "sponge": {"enabled": True}}, "sounding": {"sclr": [[1.0, 2.0], [3.0, 4.0]],
        "_start_index": {"sclr": [None, 0]}}})
    assert 'label = "a ""quote"" ! /"' in text
    assert "l_restart = .false." in text and "sponge%enabled = .true." in text
    assert "sclr(:,0) = 1.0, 2.0" in text and "sclr(:,1) = 3.0, 4.0" in text


def test_saved_copy_survives_catalog_edits_and_uses_its_own_output_name(catalog):
    settings = cases.resolved_case_settings("source")
    settings["model_setting"]["dt_main"] = 30.0
    saved = cases.save_case_definition("copied", settings, copied_from="source")
    assert saved["benchmark_case"] == "source"
    assert saved["namelists"]["model_setting"]["runtype"] == "bomex"
    assert saved["namelists"]["stats_setting"]["fname_prefix"] == "copied"
    contents = json.loads(catalog.read_text())
    contents["defaults"]["model_setting"]["dt_main"] = 10.0
    catalog.write_text(json.dumps(contents))
    assert cases.resolved_case_settings("copied")["model_setting"]["dt_main"] == 30.0
    assert "copied" in cases.available_case_names()


def test_overwrite_requires_matching_confirmation_and_preserves_other_case(catalog):
    old = cases.case_definition_record("source")
    settings = copy.deepcopy(old["namelists"])
    settings["model_setting"]["dt_main"] = 30.0
    other = cases.resolved_case_settings("other")
    with pytest.raises(FileExistsError):
        cases.save_case_definition("source", settings, copied_from="source")
    with pytest.raises(ValueError, match="changed since confirmation"):
        cases.save_case_definition("source", settings, copied_from="source", expected_sha256="stale")
    cases.save_case_definition("source", settings, copied_from="source", expected_sha256=old["sha256"])
    assert cases.resolved_case_settings("source") == settings
    assert cases.resolved_case_settings("other") == other


@pytest.mark.parametrize("name", ["../bomex", "/tmp/case", "a/b", "invalid-name"])
def test_case_names_cannot_escape_saved_directory(name, catalog):
    with pytest.raises(ValueError):
        cases.custom_case_path(name)


def test_replacing_physical_case_does_not_silently_reuse_the_original_benchmark(catalog):
    settings = cases.resolved_case_settings("source")
    settings["model_setting"]["runtype"] = "rico"
    result = cases.save_case_definition("different", settings, copied_from="source")
    assert result["benchmark_case"] is None


def test_model_values_are_saved_and_emitted_without_policy_checks(catalog):
    settings = cases.resolved_case_settings("source")
    settings["model_setting"].update(runtype="unknown scenario", dt_main=-5.0, dt_rad=7.0,
        time_final=-1.0, grid_type=99, zt_grid_fname="missing-zt", zm_grid_fname="missing-zm")
    settings["microphysics_setting"] = {"microphys_scheme": "unknown_scheme"}
    settings["stats_setting"]["fname_prefix"] = "unvalidated"
    saved = cases.save_case_definition("unvalidated", settings, copied_from="source")
    assert saved["namelists"] == settings
    text = cases.namelists_to_text(saved["namelists"])
    assert "dt_main = -5.0" in text and "grid_type = 99" in text
    assert 'runtype = "unknown scenario"' in text and 'microphys_scheme = "unknown_scheme"' in text
