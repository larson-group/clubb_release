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


def test_native_tuner_preserves_legacy_controls_at_the_entry_boundary(tmp_path):
    """Legacy array holes stay in namelist text; parameter overrides come from JSON."""
    from run_scripts.run_tuner import create_native_tuner_inputs
    from clubb_python_driver.io.namelist import _read_namelist_groups

    control = tmp_path / "legacy.in"
    text = "&stats\n tune_type = 2\n/\n&cases\n timestep_intvls_nl(3,1) = , , 121\n/\n"
    control.write_text(text)
    params = tmp_path / "params.json"
    params.write_text('{"model": {"C8": 0.2}}')
    destination = tmp_path / "prepared"
    error, cases, tune_type = create_native_tuner_inputs(
        ["bomex"], control, destination, params=params, override="dt_main=30")
    assert error.read_text().startswith(text)
    assert tune_type == 2 and cases == [destination / "bomex_hoc.in"]
    assert _read_namelist_groups(error)["clubb_params_nl"]["c8"] == 0.2
    assert _read_namelist_groups(cases[0])["model_setting"]["dt_main"] == 30


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


def test_native_reruns_inherit_defaults_without_requesting_an_empty_preset(monkeypatch, tmp_path):
    """The SCM boundary receives no default selector that could be absent from a catalog."""
    from run_scripts import run_tuner

    commands = []
    monkeypatch.setattr(run_tuner, "run_and_tee", lambda command, **kwargs: commands.append(command) or 0)
    assert run_tuner.run_case_with_run_scm("bomex", None, None, None, "core", tmp_path / "log") == 0
    assert commands[0][-1] == "bomex"
    assert not {"-params", "-flags", "-silhs_params_file"}.intersection(commands[0])


@pytest.mark.parametrize("tune_type,result_prefix,option", [
    (2, "tunable_parameters", "-params"), (3, "configurable_model_flags", "-flags"),
])
@pytest.mark.parametrize("initial_params", [False, True])
def test_native_execution_keeps_json_results_at_the_modern_boundary(
    monkeypatch, tmp_path, tune_type, result_prefix, option, initial_params,
):
    """A bounded fake tuner exercises preparation, JSON reruns and temporary-file cleanup."""
    from run_scripts import run_tuner
    launch = tmp_path / "launch"
    launch.mkdir()
    results = tmp_path / "results"
    results.mkdir()
    controls = tmp_path / "legacy.in"
    controls.write_text(f"&stats\n tune_type = {tune_type}\n/\n")
    (tmp_path / "rand_seed.dat").write_text("1\n")
    executable = tmp_path / "clubb_tuner"
    executable.touch()
    monkeypatch.setattr(run_tuner, "RUN_SCRIPTS", launch)
    monkeypatch.setattr(run_tuner, "INPUT_DIR", results)
    monkeypatch.setattr(run_tuner, "TUNABLE_DIR", results)
    monkeypatch.setattr(run_tuner, "TUNER_MISC_DIR", tmp_path)
    monkeypatch.setattr(run_tuner, "OUTPUT_DIR", tmp_path / "output")
    monkeypatch.setattr(run_tuner, "TUNER_EXE", executable)
    monkeypatch.setattr(run_tuner, "move_case_outputs", lambda *args: None)
    params = tmp_path / "initial.json"
    params.write_text('{"model": {"C8": 0.2}}')
    arguments = ["run_tuner.py", "-run_case", "bomex",
                 "-controls", str(controls), "-initial_output"]
    if initial_params:
        arguments.extend(["-params", str(params)])
    monkeypatch.setattr("sys.argv", arguments)
    scm_commands = []
    result = results / f"{result_prefix}_test.json"

    def run(command, **kwargs):
        """Model only the external execution boundary, keeping real input preparation."""
        if command == [str(executable)]:
            assert (launch / "error.in").is_file()
            assert (launch / "bomex_hoc.in").is_file()
            result.write_text('{"model": {}}')
        else:
            scm_commands.append(command)
        return 0

    monkeypatch.setattr(run_tuner, "run_and_tee", run)
    assert run_tuner.main() == 0
    assert len(scm_commands) == 2  # Initial output and optimized output.
    assert "-flags" not in scm_commands[0]
    assert ("-params" in scm_commands[0]) is initial_params
    if initial_params:
        assert scm_commands[0][scm_commands[0].index("-params") + 1] == str(params)
    assert scm_commands[1][scm_commands[1].index(option) + 1] == str(result)
    if tune_type == 3 and initial_params:
        assert scm_commands[1][scm_commands[1].index("-params") + 1] == str(params)
    else:
        assert ("-flags" if option == "-params" else "-params") not in scm_commands[1]
    assert not list(launch.iterdir())
