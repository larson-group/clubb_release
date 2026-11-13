"""Check explicit input selection through catalog and namelist boundaries.

Manual route: python -m pytest input/auto_llm_generated_tests/pytests/unit_tests/test_case_input_paths.py
The existing input pytest suite owns these checks.
"""
from pathlib import Path

from input.scripts.case_json_to_namelist import case_definition_record, case_to_namelist
from input.scripts.namelist import resolve_case_input_paths
from clubb_python_driver.io.namelist import read_namelist


def test_unrelated_filename_survives_namelist_round_trip(tmp_path):
    """A case can select a file with no relationship to its runtype or namelist location."""
    run_dir = tmp_path / "run_scripts"
    run_dir.mkdir()
    selected = tmp_path / "observations" / "measured profile.dat"
    selected.parent.mkdir()
    selected.write_text("z[m] thm[K]\n0 300\n")
    settings = case_definition_record("bomex")["namelists"]
    settings["model_setting"]["sounding_file"] = "../observations/measured profile.dat"
    nml = tmp_path / "elsewhere" / "custom.in"
    nml.parent.mkdir()
    nml.write_text(case_to_namelist(settings))
    cfg = resolve_case_input_paths(read_namelist(str(nml)), run_dir)
    assert cfg["runtype"] == "bomex"
    assert Path(cfg["sounding_file"]) == selected
    assert Path(cfg["sounding_file"]).read_text() == selected.read_text()


def test_absolute_and_unused_paths_are_preserved(tmp_path):
    """Absolute selections pass through and omitted roles stay empty without guessing names."""
    chosen = tmp_path / "forcing table.dat"
    cfg = resolve_case_input_paths({"runtype": "invented", "forcing_file": str(chosen)}, tmp_path / "run_scripts")
    assert cfg["forcing_file"] == str(chosen)
    assert cfg["sounding_file"] == ""
    assert cfg["sfc_file"] == ""


def test_scalar_roles_can_share_one_input():
    """The catalog explicitly assigns one shared profile to both scalar readers."""
    cfg = case_definition_record("bomex")["namelists"]["model_setting"]
    assert cfg["sclr_sounding_file"] == cfg["edsclr_sounding_file"]
    assert cfg["sclr_sounding_file"]
