"""Provisional JSONC reader/save checks owned by the shared Python unit stage.

Synthetic inputs check comments, strict values and comment-preserving edits.
Catalog checks exercise case and tunable-config save paths without running
the model.
"""
import copy
from pathlib import Path
import shutil

import pytest

from input.scripts.namelist import load_json, updated_json_text


def test_jsonc_comments_do_not_change_strings_or_strict_validation(tmp_path):
    """Read both comment styles without consuming comment markers inside strings."""
    path = tmp_path / "settings.jsonc"
    path.write_text('/* header */ {"url": "https://host/a/*b*/", // explanation\n'
                    '"quote": "\\\" // literal", "enabled": true /* note */}\n')
    assert load_json(path) == {"url": "https://host/a/*b*/", "quote": '\" // literal', "enabled": True}
    for text, message in [('{"a":1, /* note */ "a":2}', "Duplicate JSON key"),
                          ('{"a": /* note */ NaN}', "Non-finite"),
                          ('{"a":1 /* unfinished', "Unterminated")]:
        path.write_text(text)
        with pytest.raises(ValueError, match=message):
            load_json(path)
    path.write_text('{"a":1, // note\n}')
    with pytest.raises(ValueError):
        load_json(path)


@pytest.mark.parametrize("mask", range(8))
def test_object_membership_edits_keep_comments_and_valid_commas(tmp_path, mask):
    """Remove adjacent/last members and insert a new member in a compact object."""
    path = tmp_path / "settings.jsonc"
    path.write_text('{"a":1,/* b note */"b":2,"c":3 // last note\n}')
    expected = {key: value for index, (key, value) in enumerate({"a": 1, "b": 2, "c": 3}.items())
                if mask & (1 << index)}
    expected["added"] = {"enabled": True}
    updated = updated_json_text(path, expected)
    assert "/* b note */" in updated and "// last note" in updated
    path.write_text(updated)
    assert load_json(path) == expected
    assert updated_json_text(path, expected) == updated


def test_nested_edits_preserve_inline_block_and_container_comments(tmp_path):
    """Change scalar types and arrays while retaining comments in replaced values."""
    path = tmp_path / "settings.jsonc"
    path.write_text('// shared\n{"defaults":{"x":1 /* units */},'
                    '"presets":{"case":{"x":2 // case choice\n}},'
                    '"array":[1, /* element */ 2]}\n')
    expected = load_json(path)
    expected["defaults"]["x"] = True
    expected["presets"]["case"]["x"] = 3
    expected["array"] = [4]
    updated = updated_json_text(path, expected)
    for comment in ("// shared", "/* units */", "// case choice", "/* element */"):
        assert comment in updated
    path.write_text(updated)
    actual = load_json(path)
    assert actual == expected and actual["defaults"]["x"] is True


def test_case_save_keeps_shared_and_case_comments(tmp_path):
    """Replace a confirmed catalog case without removing either scope's notes."""
    from input.scripts.case_json_to_namelist import case_definition_record, save_case_definition
    path = tmp_path / "case_definitions.jsonc"
    path.write_text('{"defaults":{"model_setting":{"runtype":"bomex", "dt_main":60 // timestep\n}},'
                    '"presets":{"case":{"model_setting":{"dt_main":30 /* case choice */}}}}')
    record = case_definition_record("case", path)
    settings = copy.deepcopy(record["namelists"])
    settings["model_setting"]["dt_main"] = 15
    saved = save_case_definition("case", settings, json_file=path, expected_sha256=record["sha256"])
    assert saved["namelists"]["model_setting"]["dt_main"] == 15
    assert "// timestep" in path.read_text() and "/* case choice */" in path.read_text()


def test_case_reset_to_defaults_keeps_nested_comment_and_valid_jsonc(tmp_path):
    """Remove an entire override group without its line comment hiding closing braces."""
    from input.scripts.case_json_to_namelist import case_definition_record, save_case_definition
    path = tmp_path / "case_definitions.jsonc"
    path.write_text('{"defaults":{"model_setting":{"runtype":"bomex","dt_main":60}},'
                    '"presets":{"case":{"model_setting":{"dt_main":30 // experiment\n}}}}')
    record = case_definition_record("case", path)
    settings = copy.deepcopy(record["namelists"])
    settings["model_setting"]["dt_main"] = 60
    saved = save_case_definition("case", settings, json_file=path, expected_sha256=record["sha256"])
    assert saved["namelists"]["model_setting"]["dt_main"] == 60
    assert load_json(path)["presets"]["case"] == {}
    assert "// experiment" in path.read_text()


def test_named_custom_jsonc_case_is_saved_back_to_its_source(tmp_path):
    """Discover a named JSONC override and retain its comments when replacing it."""
    from input.scripts.case_json_to_namelist import available_case_names, case_definition_record, save_case_definition
    catalog = tmp_path / "case_definitions.jsonc"
    catalog.write_text('{"defaults":{"model_setting":{"dt_main":60}},"presets":{}}')
    custom = tmp_path / "experiment.case.jsonc"
    custom.write_text('{"model_setting":{"dt_main":30 // experiment timestep\n}}')
    assert "experiment" in available_case_names(catalog)
    record = case_definition_record("experiment", catalog)
    settings = copy.deepcopy(record["namelists"])
    settings["model_setting"]["dt_main"] = 15
    saved = save_case_definition("experiment", settings, json_file=catalog, expected_sha256=record["sha256"])
    assert saved["namelists"]["model_setting"]["dt_main"] == 15
    assert "// experiment timestep" in custom.read_text()
    assert not custom.with_suffix(".json").exists()


def test_tunable_save_preserves_comments_when_adding_and_replacing_presets(tmp_path):
    """Use the actual commented catalogs through the parameter/flag save owner."""
    from input.scripts.save_tunable_config import save_tunable_config
    from input.scripts.create_case_namelist import CATALOG_FILES
    root = Path(__file__).resolve().parents[2]
    originals = {}
    for kind in ("parameters", "flags"):
        file = CATALOG_FILES[kind]
        shutil.copyfile(root / file, tmp_path / file)
        originals[kind] = (tmp_path / file).read_text()
    for value, force in ((2.5, False), (3.0, True)):
        save_tunable_config("default", {"tunable": {"C4": str(value)}}, "experiment",
                            force=force, config_root=tmp_path)
        for kind in originals:
            text = (tmp_path / CATALOG_FILES[kind]).read_text()
            for line in originals[kind].splitlines():
                if "//" in line:
                    assert line in text
        assert load_json(tmp_path / "clubb_params.jsonc")["presets"]["experiment"]["model"]["C4"] == value
