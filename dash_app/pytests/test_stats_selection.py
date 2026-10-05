"""Shared Dash/MCP selection, exact saved lists and category checkbox behavior."""
import json
import shlex
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from dash import Input
from dash_app.run_tab import stats
from dash_app.shared import stats as stats_service
from dash_app.run_tab.runtime import build_case_command
from dash_app.services.models import ScmRunBatchRequest, ScmRunRequest
from utilities.stats_json_to_namelist import load_stats_categories, select_stats_entries, stats_json_to_namelist


@pytest.mark.parametrize("selection", ["all", "standard", "core", "standard/radiation", "core+wp2_budgets", "var:rtm,var:wp2", "none"])
def test_ui_command_and_mcp_keep_category_expression(selection):
    command = shlex.split(build_case_command("bomex", selection))
    assert command[command.index("-stats") + 1] == selection
    assert ScmRunRequest(request_id="stats-test-request", case="bomex", stats_file=selection).stats_file == selection
    assert ScmRunBatchRequest(request_id="stats-test-request", cases=["bomex"], stats_file=selection).stats_file == selection


def test_saved_lists_are_discoverable_portable_and_never_overwritten(tmp_path, monkeypatch):
    monkeypatch.setattr(stats_service, "STATS_DIR", str(tmp_path))
    selection = "core+standard/radiation"
    saved = stats_service.save_stats_config("core_radiation", selection)
    assert saved == "custom/core_radiation.in"
    assert saved in stats_service.list_stats_files()
    assert (tmp_path / saved).read_text() == stats_json_to_namelist(selection)
    assert stats_service.stats_cli_argument(saved) == "input/stats/custom/core_radiation.in"
    assert stats_service.stats_selection_names(saved) == {entry["name"] for entry in select_stats_entries(selection)}
    assert ScmRunRequest(request_id="saved-stats-test", case="bomex", stats_file=saved).stats_file == saved
    with pytest.raises(ValueError, match="already exists"):
        stats_service.save_stats_config("core_radiation", "all")
    assert (tmp_path / saved).read_text() == stats_json_to_namelist(selection)
    for name in ("../escape", "/absolute", "has spaces", ""):
        with pytest.raises(ValueError):
            stats_service.save_stats_config(name, selection)


@pytest.mark.parametrize("selection", ["../secret.in", "/tmp/secret.in", "custom/missing.in", "core,,all"])
def test_mcp_rejects_unknown_lists_and_traversal(selection):
    with pytest.raises(ValueError):
        ScmRunRequest(request_id="invalid-stats-test", case="bomex", stats_file=selection)


@pytest.mark.parametrize("selection", ["all", "standard", "core", "none", "var:rtm,var:wp2", "var:rtm+wp2_budgets"])
def test_variable_selection_round_trip(selection):
    names = stats_service.stats_selection_names(selection)
    expression = stats_service.stats_selection_expression(names)
    assert stats_service.stats_selection_names(expression) == names
    with pytest.raises(ValueError, match="missing from the stats catalog"):
        stats_service.stats_selection_expression(names | {"unknown_variable"})


def registered_callbacks():
    callbacks = {"client_names": [], "server_inputs": []}
    class App:
        def callback(self, *args, **kwargs):
            callbacks["server_inputs"].extend(arg.component_id for arg in args if isinstance(arg, Input))
            def register(fn):
                callbacks[fn.__name__] = fn
                return fn
            return register
        def clientside_callback(self, *args, **kwargs):
            callbacks["client_names"].append(args[0].function_name)
    stats.register_stats_callbacks(App())
    return callbacks


def checkbox_ids(categories):
    return [{"type": "run-stats-category", "path": path}
            for path in categories]


def edit_tree(callbacks, monkeypatch, trigger, selection, previous=None, *, category=None, variables=None):
    categories = load_stats_categories()
    ids = checkbox_ids(categories)
    variable_ids = [{"type": "run-stats-variables", "path": path} for path in stats.leaf_paths(categories)]
    values = [list(value) for value in previous.checks] if previous else [[] for _ in ids]
    variable_values = [list(value) for value in previous.variables] if previous else [[] for _ in variable_ids]
    if category is not None:
        values[ids.index(trigger)] = ["selected"] if category else []
    if variables is not None:
        variable_values[variable_ids.index(trigger)] = variables
    seed = previous.seed if previous else {}
    if trigger == "run-stats-custom-open":
        monkeypatch.setattr(stats, "callback_context", SimpleNamespace(triggered_id=None))
        seed = callbacks["describe_selection"]([], 0, None, selection, [])[3]
    payload = {"trigger": trigger, "args": [1, 0, 0, 0, values, variable_values, ids, variable_ids,
                                           previous.draft if previous else [],
                                           {path: list(dict.fromkeys(entry["name"] for entry in entries))
                                            for path, entries in categories.items()}, seed]}
    result = run_clientside(payload)
    updates = result["changes"]
    merge = lambda changed, old: [before if after == "NO_UPDATE" else after for after, before in zip(changed, old)]
    return SimpleNamespace(
        modal=updates[0], draft=previous.draft if updates[1] == "NO_UPDATE" else updates[1],
        checks=merge(updates[2], values), variables=merge(updates[3], variable_values),
        preview=previous.preview if updates[4] == "NO_UPDATE" else updates[4],
        mixed=result["mixed"], updates=updates, seed=seed,
    )


def run_clientside(payload):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is needed to test the actual browser selection callbacks")
    script = r"""
const fs = require('node:fs');
const vm = require('node:vm');
const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const window = {dash_clientside: {no_update: 'NO_UPDATE', callback_context: {triggered_id: payload.trigger}},
                requestAnimationFrame: fn => fn()};
const inputs = payload.args[6].map(id => ({
  path: id.path, attributes: {}, indeterminate: false,
  closest: () => ({id: JSON.stringify(id)}),
  setAttribute(name, value) {this.attributes[name] = value;},
  removeAttribute(name) {delete this.attributes[name];}
}));
const document = {querySelectorAll: () => inputs};
vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window, document});
const actions = window.dash_clientside.statsSelection;
const changes = actions.editSelection(...payload.args);
actions.syncMixedCheckboxes(changes[1] === 'NO_UPDATE' ? payload.args[8] : changes[1], payload.args[9]);
process.stdout.write(JSON.stringify({changes, mixed: inputs.filter(input => input.indeterminate).map(input => input.path)}));
"""
    path = Path(stats.__file__).parents[1] / "assets/33_stats_tree.js"
    result = subprocess.run([node, "-e", script, str(path)], input=json.dumps(payload), text=True,
                            capture_output=True, check=True, timeout=10)
    return json.loads(result.stdout)


def test_checkbox_parents_propagate_and_cancel_preserves_selection(monkeypatch):
    callbacks = registered_callbacks()
    categories = load_stats_categories()
    ids = checkbox_ids(categories)
    result = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "core")
    core_names = stats_service.stats_selection_names("core")
    assert set(result.draft) == core_names
    parent = {"type": "run-stats-category", "path": "all/standard/budgets"}
    result = edit_tree(callbacks, monkeypatch, parent, "core", result, category=True)
    assert {entry["name"] for entry in categories[parent["path"]]} <= set(result.draft)
    assert "all" in result.mixed
    # The large budget update does not rerender unrelated leaf controls.
    assert result.updates[3][stats.leaf_paths(categories).index("all/standard/core")] == "NO_UPDATE"
    result = edit_tree(callbacks, monkeypatch, parent, "core", result, category=False)
    assert set(result.draft) == core_names
    cancelled = edit_tree(callbacks, monkeypatch, "run-stats-custom-cancel", "core", result)
    assert cancelled.updates[1] == "NO_UPDATE"
    assert cancelled.updates[2] == ["NO_UPDATE"] * len(ids)
    assert cancelled.updates[3] == ["NO_UPDATE"] * len(stats.leaf_paths(categories))
    applied = edit_tree(callbacks, monkeypatch, "run-stats-apply", "core", result)
    assert applied.modal == stats.HIDDEN
    assert applied.updates[2] == ["NO_UPDATE"] * len(ids)
    monkeypatch.setattr(stats, "callback_context", SimpleNamespace(triggered_id="run-stats-apply"))
    selected = callbacks["describe_selection"]([], 1, None, "core", result.draft)
    assert stats_service.stats_selection_names(selected[2]) == core_names
    for preset in ("all", "standard", "core"):
        monkeypatch.setattr(stats, "callback_context", SimpleNamespace(triggered_id={"type": "run-stats-button", "name": preset}))
        selected = callbacks["describe_selection"]([1, 1, 1], 1, None, selected[2], result.draft)
        assert selected[2] == preset
        assert selected[0] == f"{preset.capitalize()} · {len(select_stats_entries(preset))} stats"


def test_individual_variable_removal_applies_and_saves_exact_subset(tmp_path, monkeypatch):
    callbacks = registered_callbacks()
    categories = load_stats_categories()
    ids = checkbox_ids(categories)
    result = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "core")
    variable = {"type": "run-stats-variables", "path": "all/standard/core"}
    values = result.variables[stats.leaf_paths(categories).index(variable["path"])]
    assert "thlm" in values
    result = edit_tree(callbacks, monkeypatch, variable, "core", result,
                       variables=[name for name in values if name != "thlm"])
    expected = stats_service.stats_selection_names("core") - {"thlm"}
    assert set(result.draft) == expected
    for path in ("all", "all/standard", "all/standard/core"):
        assert path in result.mixed
        assert result.checks[ids.index({"type": "run-stats-category", "path": path})] == []
    monkeypatch.setattr(stats, "callback_context", SimpleNamespace(triggered_id="run-stats-apply"))
    applied = callbacks["describe_selection"]([], 1, None, "core", result.draft)
    assert stats_service.stats_selection_names(applied[2]) == expected
    assert applied[0] == f"Custom selection · {len(expected)} stats"
    monkeypatch.setattr(stats_service, "STATS_DIR", tmp_path)
    monkeypatch.setattr(stats, "callback_context", SimpleNamespace(triggered_id="run-stats-save"))
    feedback, options = callbacks["save_selection"](1, 1, "partial_core", result.draft)
    assert feedback.className == "run-stats-feedback-success"
    assert {option["value"] for option in options} == {"custom/partial_core.in"}
    assert stats_service.stats_selection_names("custom/partial_core.in") == expected
    restored = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "custom/partial_core.in")
    assert set(restored.draft) == expected


def test_legacy_partial_category_restores_exact_variables(monkeypatch):
    callbacks = registered_callbacks()
    categories = load_stats_categories()
    result = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "fire_stats.in")
    names = stats_service.stats_selection_names("fire_stats.in")
    assert set(result.draft) == names
    for path, values in zip(stats.leaf_paths(categories), result.variables):
        assert set(values) == {entry["name"] for entry in categories[path]} & names
    assert not any(component["props"].get("className") == "run-stats-feedback-warning" for component in result.preview)


def test_bulk_selection_clear_and_duplicate_updates(monkeypatch):
    callbacks = registered_callbacks()
    result = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "none")
    parent = {"type": "run-stats-category", "path": "all"}
    result = edit_tree(callbacks, monkeypatch, parent, "none", result, category=True)
    assert set(result.draft) == stats_service.stats_selection_names("all")
    assert not result.mixed
    repeated = edit_tree(callbacks, monkeypatch, parent, "none", result, category=True)
    assert repeated.updates[1] == "NO_UPDATE"
    cleared = edit_tree(callbacks, monkeypatch, "run-stats-clear", "none", result)
    assert not cleared.draft and not cleared.mixed
    assert not any(cleared.checks) and not any(cleared.variables)


def test_toggles_and_preview_have_no_server_callback(monkeypatch):
    callbacks = registered_callbacks()
    assert {"editSelection", "syncMixedCheckboxes", "filterTree"} <= set(callbacks["client_names"])
    assert "run-stats-draft" not in callbacks["server_inputs"]
    assert "run-stats-search" not in callbacks["server_inputs"]
    assert not any(isinstance(item, dict) and item["type"] in {"run-stats-category", "run-stats-variables"}
                   for item in callbacks["server_inputs"])


def test_unknown_legacy_names_warn_when_opening(tmp_path, monkeypatch):
    monkeypatch.setattr(stats_service, "STATS_DIR", tmp_path)
    (tmp_path / "extra.in").write_text('&clubb_stats_nl\nentry(1) = "outside_catalog | zt | m | test"\n/\n')
    callbacks = registered_callbacks()
    result = edit_tree(callbacks, monkeypatch, "run-stats-custom-open", "extra.in")
    assert not result.draft
    assert result.preview[-1]["props"]["className"] == "run-stats-feedback-warning"
    assert "outside_catalog" in result.preview[-1]["props"]["children"]


def test_tree_includes_variables_and_long_names_with_categories_collapsed():
    def components(component):
        yield component
        children = getattr(component, "children", [])
        for child in children if isinstance(children, list) else [children]:
            yield from components(child)
    nodes = list(components(stats.build_stats_editor()))
    details = [node for node in nodes if type(node).__name__ == "Details"]
    assert details and all(node.open is False for node in details)
    checkboxes = [node for node in nodes if type(node).__name__ == "Checklist"]
    categories = load_stats_categories()
    for path in stats.leaf_paths(categories):
        checkbox = next(node for node in checkboxes if node.id == {"type": "run-stats-variables", "path": path})
        assert {option["value"] for option in checkbox.options} == {entry["name"] for entry in categories[path]}
        for entry in categories[path]:
            label = next(option["label"] for option in checkbox.options if option["value"] == entry["name"])
            assert label.children[0].children == entry["name"]
            assert label.children[1].children == entry["units"]
            assert label.children[2].children == stats.grid_label(entry["grid"])
            assert label.children[2].title == entry["grid"]
            assert label.children[3].children == entry["long_name"]


def search_tree(steps, opened=(), selected=()):
    """Exercise the shipped filter against the catalog's actual rendered metadata."""
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is needed to test the actual browser search callback")
    def components(component):
        yield component
        children = getattr(component, "children", [])
        for child in children if isinstance(children, list) else [children]:
            yield from components(child)
    controls = list(components(stats.build_stats_editor()))
    payload = {
        "nodes": [{"path": control.id["path"], "open": control.id["path"] in opened}
                  for control in controls if type(control).__name__ == "Details"],
        "labels": [option["label"].to_plotly_json()["props"] for control in controls
                   if getattr(control, "className", "") == "run-stats-check run-stats-variables"
                   for option in control.options],
        "selected": list(selected), "steps": steps,
    }
    # Component children aren't needed by the DOM fixture, only data-* attributes.
    payload["labels"] = [{key: value for key, value in label.items() if key.startswith("data-")}
                         for label in payload["labels"]]
    script = r"""
const fs = require('node:fs');
const vm = require('node:vm');
const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const window = {dash_clientside: {no_update: 'NO_UPDATE', callback_context: {}}};
const nodes = payload.nodes.map(node => ({dataset: {path: node.path}, open: node.open, hidden: false}));
const labels = payload.labels.map(label => ({
  dataset: {name: label['data-name'], path: label['data-path'], search: label['data-search']},
  row: {hidden: false, checked: payload.selected.includes(label['data-name'])},
  closest() {return this.row;}
}));
const tree = {dataset: {}, querySelectorAll: selector => selector === '.run-stats-node' ? nodes : labels};
const document = {getElementById: () => tree};
vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window, document});
const states = payload.steps.map(step => {
  window.dash_clientside.callback_context.triggered_id = step.trigger || 'run-stats-search';
  const result = window.dash_clientside.statsSelection.filterTree(step.query, 0, 0);
  return {result, names: labels.filter(label => !label.row.hidden).map(label => label.dataset.name),
          nodes: nodes.filter(node => !node.hidden).map(node => node.dataset.path),
          opened: nodes.filter(node => node.open).map(node => node.dataset.path),
          selected: labels.filter(label => label.row.checked).map(label => label.dataset.name)};
});
process.stdout.write(JSON.stringify(states));
"""
    path = Path(stats.__file__).parents[1] / "assets/33_stats_tree.js"
    result = subprocess.run([node, "-e", script, str(path)], input=json.dumps(payload), text=True,
                            capture_output=True, check=True, timeout=10)
    return json.loads(result.stdout)


def test_live_search_refines_matches_without_changing_selection():
    selected = stats_service.stats_selection_names("core") - {"thlm"}
    broad, narrow, empty, cleared = search_tree([
        {"query": "wp2"}, {"query": "WP2 budget"}, {"query": "no-such-stat-123"},
        {"query": "no-such-stat-123", "trigger": "run-stats-search-clear"},
    ], opened=["all", "all/standard"], selected=selected)
    assert set(narrow["names"]) < set(broad["names"])
    assert "wp2" in broad["names"] and "wp2" not in narrow["names"]
    assert "wp2_bt" in narrow["names"]
    assert "all/standard/budgets/wp2_budgets" in narrow["nodes"]
    assert set(narrow["nodes"]) <= set(narrow["opened"])
    assert "all/standard/core" not in narrow["nodes"]
    assert empty["names"] == empty["nodes"] == []
    assert empty["result"][1] == "No variables match. Try a different search."
    assert cleared["result"] == ["", ""]
    assert set(cleared["names"]) == stats_service.stats_selection_names("all")
    assert cleared["opened"] == ["all", "all/standard"]
    for state in [broad, narrow, empty, cleared]:
        assert set(state["selected"]) == selected


@pytest.mark.parametrize("query, expected", [
    ("liquid water potential", "thlm"),
    ("momentum m^2/s^2", "wp2"),
    ("zt thlm", "thlm"),
    ("standard core", "thlm"),
    ("wp2_budgets", "wp2_bt"),
])
def test_search_matches_descriptions_units_grids_and_category_paths(query, expected):
    filtered, restored = search_tree([{"query": query}, {"query": query, "trigger": "run-stats-custom-open"}])
    assert expected in filtered["names"]
    assert filtered["result"][0] == "NO_UPDATE"
    count = len(set(filtered["names"]))
    assert filtered["result"][1] == f"{count} matching variable{'s' if count != 1 else ''}"
    assert restored["result"] == ["", ""]
    assert restored["opened"] == []
