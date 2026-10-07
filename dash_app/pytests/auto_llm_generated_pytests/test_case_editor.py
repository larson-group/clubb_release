"""Provisional checks of the actual case form and browser draft callback."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from dash_app.run_tab import cases
from utilities.case_json_to_namelist import load_case_catalog, resolved_case_settings


def test_logical_values_never_render_as_text_or_number_inputs():
    baseline = resolved_case_settings("bomex")
    schema = cases.field_schema(baseline)
    def walk(component):
        yield component
        children = getattr(component, "children", [])
        for child in children if isinstance(children, list) else [children]:
            yield from walk(child)
    for group in cases.render_fields(baseline, schema):
        for control in walk(group):
            if type(control).__name__ == "Input":
                assert not isinstance(control.value, bool), control.id
    assert schema["model_setting.l_t_dependent"]["kind"] == "logical"


def run_draft(trigger, *, value="30", logical=0, array="[1e-9]", name="bomex", repeat=False,
              missing_bool=False, source_change=False, measure_paints=False):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to exercise the shipped case editor callback")
    base = {"key": "bomex", "namelists": {"model_setting": {"dt_main": 60.0, "l_restart": False,
                                                                       "sclr_tol_nl": [0.01, 1e-8]},
                                           "stats_setting": {"fname_prefix": "bomex"}, "gfdl_activation_setting": {}}}
    paths = ["model_setting.dt_main", "model_setting.l_restart", "model_setting.sclr_tol_nl", "stats_setting.fname_prefix"]
    schema = {path: {"kind": kind, "readonly": index == 3}
              for index, (path, kind) in enumerate(zip(paths, ["number", "logical", "array", "text"]))}
    if missing_bool:
        base["namelists"]["model_setting"].pop("l_restart")
    args = [[value, logical, array, "bomex"], base, name,
            [{"type": "run-case-value", "path": path} for path in paths], schema, None]
    script = r"""
const fs = require('node:fs');
const vm = require('node:vm');
const data = JSON.parse(fs.readFileSync(0, 'utf8'));
let paints = 0;
const window = {dash_clientside: {no_update:'NO_UPDATE', callback_context:{triggered_id:data.trigger}}, requestAnimationFrame:fn=>{paints++; fn();}};
const document = {querySelectorAll:()=>[], getElementById:()=>null};
vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window,document});
const edit = window.dash_clientside.caseEditor.editDraft;
let result = edit(...data.args);
if (data.repeat) {
  // Another keystroke changes the draft, but none of the unaffected controls.
  data.args[0][0] = "31";
  data.args[5] = result[0];
  paints = 0;
  result = edit(...data.args);
}
if (data.source_change) {
  data.args[5] = result[0];
  data.args[1] = JSON.parse(JSON.stringify(data.args[1]));
  data.args[1].key = data.args[2] = 'rico';
  data.args[1].namelists.model_setting.dt_main = 90;
  delete data.args[1].namelists.model_setting.l_restart;
  data.args[1].namelists.model_setting.sclr_tol_nl = [0.001];
  data.args[1].namelists.stats_setting.fname_prefix = 'rico';
  window.dash_clientside.callback_context.triggered_id = 'run-case-baseline';
  result = edit(...data.args);
}
process.stdout.write(JSON.stringify({result, paints}));
"""
    path = Path(cases.__file__).parents[1] / "assets/34_case_editor.js"
    result = subprocess.run([node, "-e", script, str(path)], input=json.dumps({"trigger": trigger, "args": args, "repeat": repeat, "source_change": source_change}),
                            text=True, capture_output=True, check=True, timeout=10)
    output = json.loads(result.stdout)
    return output if measure_paints else output["result"]


def test_browser_draft_keeps_false_and_replaces_array_values():
    result = run_draft({"type": "run-case-value", "path": "model_setting.dt_main"})
    model = result[0]["namelists"]["model_setting"]
    assert model == {"dt_main": 30.0, "l_restart": False, "sclr_tol_nl": [1e-9]}
    assert result[0]["errors"] == [] and result[3] is False


def test_blank_cells_are_omitted_without_field_specific_requirements():
    result = run_draft({"type": "run-case-value", "path": "model_setting.sclr_tol_nl"}, array="")
    assert "sclr_tol_nl" not in result[0]["namelists"]["model_setting"]
    empty = run_draft({"type": "run-case-value", "path": "model_setting.dt_main"}, value="")
    assert "dt_main" not in empty[0]["namelists"]["model_setting"]
    assert empty[0]["errors"] == [] and empty[3] is False


@pytest.mark.parametrize("missing_bool", [False, True])
@pytest.mark.parametrize("position, expected", [(0, False), (1, None), (2, True)])
def test_boolean_slider_sets_false_unset_or_true(position, expected, missing_bool):
    result = run_draft({"type": "run-case-value", "path": "model_setting.l_restart"},
                       missing_bool=missing_bool, logical=position)
    model = result[0]["namelists"]["model_setting"]
    if expected is None:
        assert "l_restart" not in model
    else:
        assert model["l_restart"] is expected
    assert result[0]["errors"] == []


def test_new_case_name_updates_output_identity_without_changing_physical_case():
    result = run_draft("run-case-name", name="bomex_30s")
    assert result[0]["namelists"]["stats_setting"]["fname_prefix"] == "bomex_30s"
    assert result[0]["namelists"]["model_setting"]["dt_main"] == 30.0


def test_invalid_array_blocks_save_and_reports_the_entry():
    result = run_draft({"type": "run-case-value", "path": "model_setting.sclr_tol_nl"}, array="[invalid]")
    assert result[3] is True
    assert "model_setting.sclr_tol_nl" in result[2]


def test_typing_keeps_unchanged_controls_and_summaries_mounted():
    output = run_draft({"type": "run-case-value", "path": "model_setting.dt_main"}, repeat=True, measure_paints=True)
    result = output["result"]
    assert result[0]["namelists"]["model_setting"]["dt_main"] == 31
    assert result[1:4] == ["NO_UPDATE"] * 3
    assert result[4] == "NO_UPDATE"
    assert output["paints"] == 0


def test_source_change_replaces_cached_cells_and_clears_previous_errors():
    result = run_draft({"type": "run-case-value", "path": "model_setting.sclr_tol_nl"},
                       array="[invalid]", source_change=True)
    model = result[0]["namelists"]["model_setting"]
    assert model == {"dt_main": 90, "sclr_tol_nl": [0.001]}
    assert result[0]["errors"] == [] and result[0]["changed"] == []
    assert result[4] == ["90", 1, "[0.001]", "rico"]


@pytest.mark.parametrize("new_field", [False, True])
def test_source_changes_reuse_controls_until_the_json_fields_change(monkeypatch, new_field):
    from types import SimpleNamespace
    from dash import Dash, no_update
    from utilities.case_json_to_namelist import case_definition_record
    app = Dash(__name__, suppress_callback_exceptions=True)
    cases.register_case_callbacks(app)
    load = next(entry["callback"].__wrapped__ for entry in app.callback_map.values()
                if entry.get("callback") and entry["callback"].__wrapped__.__name__ == "load_editor")
    previous = {"key": "bomex", **case_definition_record("bomex")}
    record = case_definition_record("rico")
    if new_field:
        record["namelists"]["model_setting"]["editor_added_parameter"] = 4.2
    monkeypatch.setattr(cases, "case_definition_record", lambda _name: record)
    monkeypatch.setattr(cases, "callback_context", SimpleNamespace(triggered_id="run-case-source"))
    result = load("rico", 1, None, [], previous, cases.field_schema(previous["namelists"]))
    assert result[4]["key"] == "rico" and result[4]["namelists"] == record["namelists"]
    if new_field:
        assert result[3] is not no_update
        assert result[5]["model_setting.editor_added_parameter"]["kind"] == "number"
    else:
        assert result[3] is no_update and result[5] is no_update


def test_saved_case_buttons_stay_selectable_below_catalog_cases(monkeypatch):
    from dash_app.run_tab import layout
    monkeypatch.setattr(layout, "load_case_catalog", lambda: {"cases": {"bomex": {}}})
    buttons = layout.build_case_buttons(["bomex_copy", "bomex"])
    assert buttons[0].id == {"type": "run-case-button", "name": "bomex"}
    custom = layout.build_case_buttons(["bomex_copy", "bomex"], custom=True)
    assert [button.id for button in custom] == [{"type": "run-case-button", "name": "bomex_copy"}]


def test_every_catalog_case_round_trips_through_the_actual_cell_parser():
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to exercise the shipped case editor callback")
    catalog = load_case_catalog()
    samples = []
    for name in catalog["cases"]:
        settings = resolved_case_settings(name)
        schema = cases.field_schema(settings, catalog)
        actual = dict(cases.setting_leaves(settings))
        ids = [{"type": "run-case-value", "path": path} for path in schema]
        values = [cases.field_value(actual.get(path), schema[path]) for path in schema]
        samples.append([values, {"key": name, "namelists": settings}, name, ids, schema, None])
    script = r"""
const fs = require('node:fs'), vm = require('node:vm');
const window = {dash_clientside:{no_update:'NO_UPDATE', callback_context:{triggered_id:'run-case-baseline'}}, requestAnimationFrame:fn=>fn()};
const document = {getElementById:()=>null};
vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window,document});
const samples = JSON.parse(fs.readFileSync(0,'utf8'));
process.stdout.write(JSON.stringify(samples.map(args => window.dash_clientside.caseEditor.editDraft(...args)[0])));
"""
    path = Path(cases.__file__).parents[1] / "assets/34_case_editor.js"
    result = subprocess.run([node, "-e", script, str(path)], input=json.dumps(samples),
                            text=True, capture_output=True, check=True, timeout=10)
    for expected, actual in zip(samples, json.loads(result.stdout), strict=True):
        assert actual["namelists"] == expected[1]["namelists"], expected[2]
        assert actual["errors"] == [] and actual["changed"] == [], expected[2]
