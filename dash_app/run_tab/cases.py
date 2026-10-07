"""Custom case form: rendering only; the shared case helper owns saved data."""
import copy
import json

from dash import ALL, ClientsideFunction, Input, Output, State, callback_context, dcc, html, no_update
from dash.exceptions import PreventUpdate

from utilities.case_json_to_namelist import (
    available_case_names, case_definition_record, load_case_catalog,
    namelists_to_text, save_case_definition,
)
from dash_app.shared.components import styled_dropdown

HIDDEN = "shared-notecard-overlay run-case-modal--hidden"


def setting_leaves(settings, prefix=()):
    for key, value in settings.items():
        if key.startswith("_"):
            continue
        path = (*prefix, key)
        if isinstance(value, dict):
            yield from setting_leaves(value, path)
        else:
            yield ".".join(path), value


def field_schema(baseline, catalog=None):
    catalog = catalog if catalog is not None else load_case_catalog()
    definitions = [catalog["defaults"], baseline]
    definitions.extend(catalog["cases"].values())
    values = {}
    for definition in definitions:
        for path, value in setting_leaves(definition):
            if value is not None:
                values.setdefault(path, []).append(value)
    schema = {}
    for path, samples in values.items():
        kind = ("array" if any(isinstance(value, list) for value in samples) else
                "logical" if all(isinstance(value, bool) for value in samples) else
                "number" if all(isinstance(value, (int, float)) for value in samples) else "text")
        schema[path] = {"kind": kind,
                        "readonly": path in {"stats_setting.fname_prefix", "stats_setting.stats_output_filename"}}
    return schema


def field_value(value, spec):
    if spec["kind"] == "logical":
        return 1 if value is None else 2 if value else 0
    if value is None:
        return ""
    if isinstance(value, list) or value == "":
        return json.dumps(value)
    return str(value)


def render_fields(baseline, schema):
    actual = dict(setting_leaves(baseline))
    groups = {}
    paths = sorted(schema, key=lambda path: (path.split('.')[0] != 'model_setting', path != 'model_setting.runtype'))
    for path in paths:
        spec = schema[path]
        group, label = path.split(".", 1)
        props = {"id": {"type": "run-case-value", "path": path},
                 "value": field_value(actual.get(path), spec)}
        if spec["kind"] == "logical":
            control = html.Div([
                html.Div([html.Span("false"), html.Span("unset"), html.Span("true")],
                         className="run-case-bool-labels", **{"aria-hidden": "true"}),
                dcc.Input(**props, type="range", min=0, max=2, step=1, debounce=False),
            ], className="run-case-bool", **{"data-state": ("false", "unset", "true")[props["value"]]})
        else:
            control = dcc.Input(**props, type="text", disabled=spec["readonly"],
                                debounce=False, placeholder="Not specified")
        groups.setdefault(group, []).append(html.Div([
            html.Label(html.Code(label.replace(".", "%")), htmlFor=json.dumps(props["id"], sort_keys=True, separators=(",", ":"))),
            html.Div(control, className="run-case-control"),
        ], id=f"run-case-row-{path}", className="run-case-field", **{"data-path": path}))
    for group in baseline:
        groups.setdefault(group, [])
    return [html.Section([
        html.H3([html.Span(group.replace("_", " ").capitalize()),
                 html.Span("", id=f"run-case-group-changes-{group}", className="run-case-group-changes")]),
        html.Div(rows or "Uses model defaults.", className="run-case-group-fields"),
    ], className="run-case-group") for group, rows in groups.items()]


def build_case_editor():
    return html.Div([
        dcc.Store(id="run-case-baseline"), dcc.Store(id="run-case-schema"),
        dcc.Store(id="run-case-draft"), dcc.Store(id="run-case-saved"), dcc.Store(id="run-case-pending"),
        html.Div([
            html.Div([html.Div([html.Div("Custom case", className="shared-notecard-title"),
                                 html.Div("Copy a case or start from common defaults. Changes stay in this draft until saved.", className="run-stats-editor-note")]),
                      html.Button("Close", id="run-case-cancel", n_clicks=0, className="shared-notecard-close")],
                     className="shared-notecard-header"),
            html.Div([
                html.Div([html.Label("Start from", htmlFor="run-case-source"),
                          styled_dropdown(id="run-case-source", options=[], value="__defaults__", clearable=False)]),
                html.Div([html.Label("Case name", htmlFor="run-case-name"),
                          dcc.Input(id="run-case-name", value="", placeholder="e.g. bomex_30s", type="text", debounce=False)]),
                html.Div("Choosing another starting point replaces this draft. The physical runtype supplies soundings and forcings; the case name identifies saved settings and output.",
                         className="run-case-source-note run-stats-editor-note"),
            ], className="run-case-top"),
            html.Div([
                html.Div('Blank values use model defaults; enter "" for an empty string.', className="run-case-table-note"),
                html.Div(id="run-case-fields", className="run-case-fields"),
                html.Details([
                    html.Summary("Namelist preview"),
                    html.Button("Update preview", id="run-case-preview", n_clicks=0),
                    html.Pre("Preview your resolved model settings here.", id="run-case-preview-text"),
                ], className="run-case-preview"),
            ], className="run-case-body"),
            html.Div([
                html.Div([html.Div(id="run-case-change-count", role="status"),
                          html.Div(id="run-case-draft-errors", className="run-stats-feedback-error", role="status"),
                          html.Div(id="run-case-feedback", role="status")]),
                html.Button("Save case", id="run-case-save", n_clicks=0, disabled=True, className="run-stats-apply"),
            ], className="run-stats-editor-footer"),
            html.Div([
                html.Div([html.H3("Replace existing case?"), html.Div(id="run-case-overwrite-details"),
                          html.Div([html.Button("Keep existing", id="run-case-overwrite-cancel", n_clicks=0),
                                    html.Button("Replace case", id="run-case-overwrite-confirm", n_clicks=0, className="run-stats-apply")],
                                   className="run-case-confirm-actions")], className="run-case-confirm-panel", role="alertdialog"),
            ], id="run-case-confirm", hidden=True, className="run-case-confirm"),
        ], className="shared-notecard-panel run-case-editor-panel", role="dialog", **{"aria-modal": "true", "aria-label": "Custom case"}),
    ], id="run-case-modal", className=HIDDEN)


def register_case_callbacks(app):
    app.clientside_callback(
        ClientsideFunction(namespace="caseEditor", function_name="showEditor"),
        Output("run-case-modal", "className"), Input("run-case-custom-open", "n_clicks"),
        Input("run-case-cancel", "n_clicks"), prevent_initial_call=True,
    )

    @app.callback(Output("run-case-source", "options"), Output("run-case-source", "value"),
                  Output("run-case-name", "value"), Output("run-case-fields", "children"),
                  Output("run-case-baseline", "data"), Output("run-case-schema", "data"),
                  Input("run-case-source", "value"), Input("run-case-custom-open", "n_clicks"),
                  Input("run-case-saved", "data"), State("run-selected-cases", "data"),
                  State("run-case-baseline", "data"), State("run-case-schema", "data"),
                  prevent_initial_call=True)
    def load_editor(source, opened, saved, selected, previous, previous_schema):
        if not opened:
            raise PreventUpdate
        trigger = callback_context.triggered_id
        if trigger == "run-case-saved" and saved:
            source = saved["name"]
        elif trigger == "run-case-custom-open" and len(selected or []) == 1:
            source = selected[0]
        elif trigger == "run-case-source" and previous and previous["key"] == source:
            raise PreventUpdate
        source = source or "__defaults__"
        if trigger == "run-case-custom-open" and previous and previous["key"] == source:
            # Keep the mounted controls and browser-local draft on reopening.
            raise PreventUpdate
        catalog = load_case_catalog()
        names = available_case_names()
        options = [{"label": "Defaults", "value": "__defaults__"}] + [
            {"label": f"{'Case' if name in catalog['cases'] else 'Saved'} · {name}", "value": name} for name in names]
        record = (case_definition_record(source) if source != "__defaults__" else
                  {"namelists": copy.deepcopy(catalog["defaults"]), "benchmark_case": None})
        baseline = {"key": source, **record}
        schema = field_schema(record["namelists"], catalog)
        reuse_fields = previous and previous_schema == schema and previous["namelists"].keys() == record["namelists"].keys()
        fields = no_update if reuse_fields else render_fields(record["namelists"], schema)
        return options, source, "" if source == "__defaults__" else source, fields, baseline, no_update if reuse_fields else schema

    app.clientside_callback(
        ClientsideFunction(namespace="caseEditor", function_name="editDraft"),
        Output("run-case-draft", "data"), Output("run-case-change-count", "children"),
        Output("run-case-draft-errors", "children"), Output("run-case-save", "disabled"),
        Output({"type": "run-case-value", "path": ALL}, "value"),
        Input({"type": "run-case-value", "path": ALL}, "value"),
        Input("run-case-baseline", "data"), Input("run-case-name", "value"),
        State({"type": "run-case-value", "path": ALL}, "id"),
        State("run-case-schema", "data"), State("run-case-draft", "data"),
    )

    @app.callback(Output("run-case-preview-text", "children"), Input("run-case-preview", "n_clicks"),
                  State("run-case-draft", "data"), prevent_initial_call=True)
    def preview_case(_clicks, draft):
        try:
            if not draft or draft.get("errors"):
                raise ValueError("Complete the highlighted entries first")
            return namelists_to_text(draft["namelists"])
        except (ValueError, OSError) as exc:
            return str(exc)

    @app.callback(Output("run-case-feedback", "children"), Output("run-case-confirm", "hidden"),
                  Output("run-case-overwrite-details", "children"), Output("run-case-pending", "data"),
                  Output("run-case-saved", "data"), Output("run-case-buttons", "children"),
                  Output("run-custom-case-buttons", "children"),
                  Output("run-selected-cases", "data", allow_duplicate=True),
                  Input("run-case-save", "n_clicks"), Input("run-case-overwrite-confirm", "n_clicks"),
                  Input("run-case-overwrite-cancel", "n_clicks"), Input("run-case-cancel", "n_clicks"),
                  State("run-case-name", "value"), State("run-case-draft", "data"),
                  State("run-case-baseline", "data"), State("run-case-pending", "data"),
                  State("run-selected-cases", "data"), prevent_initial_call=True)
    def save_case(_save, _confirm, _cancel_overwrite, _cancel, name, draft, baseline, pending, selected):
        from .layout import build_case_buttons
        trigger = callback_context.triggered_id
        if trigger in {"run-case-overwrite-cancel", "run-case-cancel"}:
            return "", True, "", None, no_update, no_update, no_update, no_update
        try:
            if trigger == "run-case-overwrite-confirm":
                if not pending:
                    raise PreventUpdate
                payload = pending
            else:
                if not draft or draft.get("errors"):
                    raise ValueError("Complete the highlighted entries first")
                payload = {"name": str(name or "").strip(), "settings": draft["namelists"],
                           "copied_from": None if baseline["key"] == "__defaults__" else baseline["key"],
                           "source_record": {key: value for key, value in baseline.items() if key != "key"}}
                if payload["name"] in available_case_names():
                    current = case_definition_record(payload["name"])
                    payload["expected_sha256"] = current["sha256"]
                    old = dict(setting_leaves(current["namelists"]))
                    new = dict(setting_leaves(payload["settings"]))
                    rows = [html.Tr([html.Td(path), html.Td(json.dumps(old.get(path, "Unspecified"))),
                                     html.Td(json.dumps(new.get(path, "Unspecified")))])
                            for path in dict.fromkeys([*old, *new]) if old.get(path) != new.get(path)]
                    details = [html.P(f"Replace {payload['name']} in {'the main case catalog' if current['kind'] == 'catalog' else 'its saved custom file'}?"),
                               html.Table([html.Thead(html.Tr([html.Th("Entry"), html.Th("Before"), html.Th("After")])), html.Tbody(rows)])]
                    return "", False, details, payload, no_update, no_update, no_update, no_update
            saved = save_case_definition(**payload)
        except (ValueError, OSError) as exc:
            return html.Span(str(exc), className="run-stats-feedback-error"), True, "", None, no_update, no_update, no_update, no_update
        names = available_case_names()
        selection = list(dict.fromkeys([*(selected or []), saved["name"]]))
        return f"Saved {saved['name']}. It is selected for the next run.", True, "", None, {
            "name": saved["name"], "sha256": saved["sha256"]}, build_case_buttons(names), build_case_buttons(names, custom=True), selection
