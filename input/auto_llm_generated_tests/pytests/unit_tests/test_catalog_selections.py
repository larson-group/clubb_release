"""Catalog inheritance and selection contracts for the shared Python unit stage.

Synthetic definitions exercise nested settings, independent run lists and local
flag names. Runner checks replace only model execution, retaining real parsing
and namelist preparation.
"""
import copy
import json
from pathlib import Path
import sys

import pytest

from input.scripts.input_catalog import load_catalog, resolve_selection, resolve_input_settings
from input.scripts.case_json_to_namelist import case_definition_record, save_case_definition


@pytest.fixture
def catalog_path(tmp_path):
    """Define a two-level family and a run list with repeated membership."""
    path = tmp_path / "cases.jsonc"
    path.write_text(json.dumps({
        "defaults": {"model_setting": {"dt_main": 60, "time_final": 600, "runtype": "bomex"}},
        "presets": {"config": {"model_setting": {"dt_main": 40}},
                    "unrelated": {"model_setting": {"dt_main": 20}}},
        "groups": {"family_group": {
            "settings": {"model_setting": {"dt_main": 30, "numbers": [1, 2]}},
            "presets": {"cold": {"model_setting": {"temperature": 270}}},
            "groups": {"warm_group": {
                "settings": {"model_setting": {"dt_main": 15}},
                "presets": {"hot": {"model_setting": {"temperature": 300}},
                            "warm": {"model_setting": {"numbers": [3]}}}}}}},
        "lists": {"smoke_cases": ["cold", "warm_group", "cold", "unrelated"],
                  "nested_cases": ["smoke_cases", "family_group"]},
    }))
    return path


def test_nested_inheritance_and_selection_lists_are_independent(catalog_path):
    """Groups override config settings; lists select leaves without changing values."""
    assert resolve_selection("cases", "nested_cases", json_file=catalog_path) == ["cold", "hot", "warm", "unrelated"]
    hot = resolve_input_settings("cases", "hot", config="config", json_file=catalog_path)
    warm = resolve_input_settings("cases", "warm", config="config", json_file=catalog_path)
    assert hot["model_setting"]["dt_main"] == 15
    assert hot["model_setting"]["numbers"] == [1, 2]
    assert warm["model_setting"]["numbers"] == [3]
    assert resolve_input_settings("cases", "cold", config="config", json_file=catalog_path)["model_setting"]["dt_main"] == 30


@pytest.mark.parametrize("lists,selection,match", [
    ({"a": ["b"], "b": ["a"]}, "a", "cycle"),
    ({"a": ["missing"]}, "a", "Unknown"),
    ({"cold": ["hot"]}, "cold", "ambiguous"),
    ({"a": []}, "a", "Invalid"),
])
def test_invalid_lists_fail_before_execution(catalog_path, lists, selection, match):
    """Reject cycles, unknown references and ambiguous or empty list definitions."""
    catalog = json.loads(catalog_path.read_text())
    catalog["lists"] = lists
    catalog_path.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match=match):
        resolve_selection("cases", selection, json_file=catalog_path)


def test_grouped_case_save_edits_only_its_leaf(catalog_path):
    """Replacing a grouped case retains inherited fields and its sibling definitions."""
    before = json.loads(catalog_path.read_text())
    record = case_definition_record("hot", catalog_path)
    settings = copy.deepcopy(record["namelists"])
    settings["model_setting"]["dt_main"] = 12
    saved = save_case_definition("hot", settings, expected_sha256=record["sha256"], json_file=catalog_path)
    assert saved["namelists"] == settings
    after = json.loads(catalog_path.read_text())
    group = after["groups"]["family_group"]["groups"]["warm_group"]
    assert group["settings"] == before["groups"]["family_group"]["groups"]["warm_group"]["settings"]
    assert group["presets"]["hot"] == {"model_setting": {"dt_main": 12, "temperature": 300}}
    assert case_definition_record("warm", catalog_path)["namelists"]["model_setting"]["dt_main"] == 15


def test_case_names_are_global_and_group_suffix_is_required(catalog_path):
    """Keep case selectors unambiguous and make inheritance groups recognizable."""
    catalog = json.loads(catalog_path.read_text())
    catalog["presets"]["cold"] = {}
    catalog_path.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="Duplicate"):
        load_catalog("cases", catalog_path)
    del catalog["presets"]["cold"]
    catalog["groups"]["family"] = catalog["groups"].pop("family_group")
    catalog_path.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="_group"):
        load_catalog("cases", catalog_path)


def test_flags_keep_local_member_names_and_independent_bases(tmp_path):
    """Two collections can define flag1 without one inheriting the other's patch."""
    path = tmp_path / "flags.jsonc"
    path.write_text(json.dumps({
        "defaults": {"model": {"first": False, "second": False}, "silhs": {}}, "presets": {},
        "groups": {"core_group": {"presets": {"flag1": {"model": {"first": True}},
                                                "flag2": {"model": {"second": True}}}},
                   "host_group": {"presets": {"flag1": {"model": {"second": True}}}}},
        "lists": {"both": ["core_group", "host_group"]},
    }))
    assert resolve_selection("flags", "both", json_file=path) == ["core_group/flag1", "core_group/flag2", "host_group/flag1"]
    assert resolve_input_settings("flags", "core_group/flag2", json_file=path)["model"] == {"first": False, "second": True}
    assert resolve_input_settings("flags", "host_group/flag1", json_file=path)["model"] == {"first": False, "second": True}


def test_single_case_avoids_preparation_subprocess_and_worker_pool(tmp_path, monkeypatch):
    """A one-case selection prepares a real namelist and launches only the model runner."""
    from concurrent import futures
    from run_scripts import run_scm
    calls = []
    class Runner:
        """Record the prepared input in place of launching the executable."""
        def __init__(self, args, **kwargs):
            self.args = args
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass
        def run_and_wait(self, path, **kwargs):
            calls.append(Path(path).read_text())
            return {"state": "completed"}
    def forbidden(*args, **kwargs):
        """Fail on orchestration that should be unnecessary for one run."""
        pytest.fail("single-case preparation spawned a process or pool")
    monkeypatch.setattr(run_scm, "clubb_case_runner", Runner)
    monkeypatch.setattr(run_scm.subprocess, "run", forbidden)
    monkeypatch.setattr(futures, "ThreadPoolExecutor", forbidden)
    monkeypatch.setattr(sys, "argv", ["run_scm.py", "bomex", "-stats", "none", "-output_dir", str(tmp_path)])
    assert run_scm.main() == 0
    assert len(calls) == 1 and "&model_setting" in calls[0]


def test_multi_case_continues_after_failure_and_preserves_selection(tmp_path, monkeypatch):
    """A failed member does not prevent remaining cases from running sequentially."""
    from run_scripts import run_scm
    calls = []
    def run_case(args, path, **kwargs):
        """Record each expanded case and simulate one model failure."""
        calls.append(args.case_name)
        return int(args.case_name == "clex9_nov02")
    monkeypatch.setattr(run_scm, "run_case", run_case)
    monkeypatch.setattr(sys, "argv", ["run_scm.py", "-cases", "clex_group", "bomex", "-workers", "1"])
    assert run_scm.main() == 1
    assert calls == ["clex9_nov02", "clex9_oct14", "bomex"]


def test_parallel_runs_respect_worker_limit_and_report_exceptions(monkeypatch):
    """Exercise the real scheduler with two overlapping tasks and one failed run."""
    import threading
    from run_scripts import run_scm
    lock = threading.Lock()
    barrier = threading.Barrier(2)
    active = maximum = started = 0
    finished = []
    def run_case(args, path, **kwargs):
        """Overlap the first pair and raise on the final case without launching models."""
        nonlocal active, maximum, started
        with lock:
            active += 1
            started += 1
            position = started
            maximum = max(maximum, active)
        try:
            if position <= 2:
                barrier.wait(timeout=5)
            if args.case_name == "atex":
                raise RuntimeError("synthetic model failure")
            return 0
        finally:
            with lock:
                active -= 1
                finished.append(args.case_name)
    monkeypatch.setattr(run_scm, "run_case", run_case)
    monkeypatch.setattr(sys, "argv", ["run_scm.py", "-cases", "arm,bomex,atex", "-workers", "2"])
    assert run_scm.main() == 1
    assert maximum == 2
    assert sorted(finished) == ["arm", "atex", "bomex"]


def test_flag_group_expands_independent_run_tasks(monkeypatch):
    """Selecting a collection runs its members with distinct outputs and no injected baseline."""
    from run_scripts import run_scm
    calls = []
    def run_case(args, path, **kwargs):
        """Capture member identity and final overrides at the runner boundary."""
        calls.append((args.flags, path.name, args.override))
        return 0
    monkeypatch.setattr(run_scm, "run_case", run_case)
    group = "run_bindiff_w_flags_config_example_group"
    monkeypatch.setattr(sys, "argv", ["run_scm.py", "bomex", "-flags", group,
                                    "-workers", "1", "-override", "dt_main=30"])
    assert run_scm.main() == 0
    assert calls == [(f"{group}/flag{i}", f"{group}__flag{i}", ["dt_main=30"]) for i in (1, 2)]


def test_combined_flag_collections_keep_both_local_flag1_definitions():
    """Flattening a multi-collection sweep cannot overwrite repeated local names."""
    from input.scripts.flag_sets import read_flag_settings
    core = "run_bindiff_w_flags_config_core_flags_group"
    example = "run_bindiff_w_flags_config_example_group"
    combined = read_flag_settings(f"{core},{example}")
    assert len(combined) == 20
    assert combined[f"{core}__flag1"] == read_flag_settings(core)["flag1"]
    assert combined[f"{example}__flag1"] == read_flag_settings(example)["flag1"]


def test_jax_comparison_keeps_existing_collection_output_names():
    """Migrated collection identities retain output names consumed by bindiff."""
    from human_reviewed_tests.run_jax_vs_fortran_cases import FlagData
    group = "run_bindiff_w_flags_config_core_flags_group"
    assert FlagData("flag1", {"saturation_formula": 2}, group).dir_name == "run_bindiff_w_flags_config_core_flags_flag1"
    assert FlagData("default", None, group).dir_name == "default"
