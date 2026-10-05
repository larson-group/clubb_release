"""Category selection must follow editable JSON, not fixed registry sizes."""
import json

import pytest

from utilities.stats_json_to_namelist import load_stats_categories, select_stats_entries, stats_json_to_namelist


def definition(name):
    return dict(name=name, grid="zt", units="m", long_name=name)


def test_catalog_edits_need_no_consumer_changes(tmp_path):
    path = tmp_path / "stats.json"
    tree = {"all": {"standard": {"radiation": [definition("first")]},
                    "non-standard": {"radiation": [definition("second")]}}}
    path.write_text(json.dumps(tree))
    assert [entry["name"] for entry in select_stats_entries("radiation", path)] == ["first", "second"]
    assert [entry["name"] for entry in select_stats_entries("standard/radiation", path)] == ["first"]
    tree["all"]["standard"]["radiation"].append(definition("added"))
    moved = tree["all"]["non-standard"].pop("radiation")
    tree["all"]["new_branch"] = {"new_category": moved}
    tree["all"]["standard"]["radiation"].pop(0)
    path.write_text(json.dumps(tree))
    assert [entry["name"] for entry in select_stats_entries("all", path)] == ["added", "second"]
    assert [entry["name"] for entry in select_stats_entries("new_branch/new_category", path)] == ["second"]
    assert "all/new_branch/new_category" in load_stats_categories(path)


def test_unions_deduplicate_and_preserve_catalog_order(tmp_path):
    path = tmp_path / "stats.json"
    path.write_text(json.dumps({"all": {"a": [definition("x")], "b": [definition("y")]}}))
    assert stats_json_to_namelist("a,b", path) == stats_json_to_namelist("b+a+a", path)
    assert stats_json_to_namelist("all,a", path) == stats_json_to_namelist("all", path)


def test_individual_variables_union_with_categories_and_keep_metadata(tmp_path):
    path = tmp_path / "stats.json"
    entries = [definition("x"), definition("a"), definition("y")]
    path.write_text(json.dumps({"all": {"a": entries[:2], "b": entries[2:]}}))
    # The prefix distinguishes a variable from a category with the same name.
    assert select_stats_entries("var:a", path) == [entries[1]]
    assert select_stats_entries("var:y+a+var:x", path) == entries
    assert select_stats_entries("var:y,var:x,var:y", path) == [entries[0], entries[2]]
    with pytest.raises(ValueError, match="Unknown stats variable"):
        select_stats_entries("var:missing", path)


@pytest.mark.parametrize("selection", ["", "core,", "none+core", "missing", "../all"])
def test_invalid_selections_are_explicit(selection):
    with pytest.raises(ValueError):
        select_stats_entries(selection)


def test_empty_categories_and_quoted_metadata(tmp_path):
    path = tmp_path / "stats.json"
    entry = definition("x") | {"long_name": 'A "quoted" label!'}
    path.write_text(json.dumps({"all": {"empty": [], "values": [entry]}}))
    assert stats_json_to_namelist("empty", path) == "&clubb_stats_nl\n/\n"
    assert 'A ""quoted"" label!' in stats_json_to_namelist("values", path)
    assert stats_json_to_namelist("none", path) == ""
