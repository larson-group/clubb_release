"""The Run tab's category tree and reusable stats-list editor."""

from dash import ALL, ClientsideFunction, Input, Output, State, callback_context, dcc, html, no_update
from dash.exceptions import PreventUpdate

from utilities.stats_json_to_namelist import load_stats_categories
from dash_app.shared.stats import list_stats_files, save_stats_config, stats_list_options, stats_selection_expression, stats_selection_names

HIDDEN = "shared-notecard-overlay run-stats-modal--hidden"
VISIBLE = "shared-notecard-overlay run-stats-modal"


def leaf_paths(categories):
    """Find terminal categories from the actual JSON hierarchy."""
    parents = {path.rsplit("/", 1)[0] for path in categories if "/" in path}
    return [path for path in categories if path not in parents]


def category_label(path):
    name = path.rsplit("/", 1)[-1].replace("_", " ")
    return name[0].upper() + name[1:]


def grid_label(grid):
    """Use level names in the UI while retaining the catalog's grid IDs."""
    return {"zt": "thermo", "zm": "momentum", "sfc": "surface",
            "lh_zt": "SILHS thermo", "lh_sfc": "SILHS surface",
            "rad_zt": "radiation thermo", "rad_zm": "radiation momentum"}.get(grid, grid)


def build_stats_editor():
    categories = load_stats_categories()

    def variable_option(entry, path):
        return {"label": html.Span([
                html.Code(entry["name"], className="run-stats-variable-name"),
                html.Span(entry["units"], className="run-stats-variable-units"),
                html.Span(grid_label(entry["grid"]), title=entry["grid"], className="run-stats-variable-grid"),
                html.Span(entry["long_name"], className="run-stats-variable-description"),
            ], className="run-stats-variable-label", **{
                "data-name": entry["name"], "data-path": path,
                "data-search": " ".join([
                    path, path.replace("_", " "), entry["name"], entry["units"],
                    entry["grid"], grid_label(entry["grid"]), entry["long_name"],
                ]).lower(),
            }), "value": entry["name"],
                "title": f"{entry['grid']} · {entry['units']}"}

    def node(path):
        names = {entry["name"] for entry in categories[path]}
        check = dcc.Checklist(
            id={"type": "run-stats-category", "path": path},
            options=[{"label": html.Span([
                html.Span(category_label(path)),
                html.Span(f"{len(names)}", className="run-stats-count"),
            ], className="run-stats-node-label"), "value": "selected", "title": path}],
            value=[], className="run-stats-check", labelStyle={"display": "flex"},
        )
        children = [child for child in categories
                    if "/" in child and child.rsplit("/", 1)[0] == path]
        props = {"id": {"type": "run-stats-node", "path": path},
                 "className": "run-stats-node", "data-path": path}
        rows = ([node(child) for child in children] if children else
                dcc.Checklist(
                    id={"type": "run-stats-variables", "path": path},
                    options=[variable_option(entry, path) for entry in
                             {entry["name"]: entry for entry in categories[path]}.values()],
                    value=[], className="run-stats-check run-stats-variables", labelStyle={"display": "flex"},
                ))
        return html.Details([
            html.Summary(check),
            html.Div(rows, className="run-stats-children"),
        ], open=False, **props)

    return html.Div([
        dcc.Store(id="run-stats-draft", data=[]),
        dcc.Store(id="run-stats-seed"),
        dcc.Store(id="run-stats-catalog", data={
            path: list(dict.fromkeys(entry["name"] for entry in entries))
            for path, entries in categories.items()
        }),
        dcc.Store(id="run-stats-checkbox-signal"),
        html.Div([
            html.Div([
                html.Div([
                    html.Div("Custom statistics", id="run-stats-editor-title", className="shared-notecard-title"),
                    html.Div("Expand categories to choose individual variables, or check a category to include everything below it.", className="run-stats-editor-note"),
                ]),
                html.Button("Cancel", id="run-stats-custom-cancel", className="shared-notecard-close", n_clicks=0),
            ], className="shared-notecard-header"),
            html.Div([
                html.Div([
                    html.Div([html.Span("Categories"), html.Button("Clear all", id="run-stats-clear", n_clicks=0)],
                             className="run-stats-tree-heading"),
                    html.Div([
                        html.Label("Search statistics", htmlFor="run-stats-search"),
                        html.Div([
                            dcc.Input(id="run-stats-search", value="", type="search", debounce=False,
                                      placeholder="Names, descriptions, units or categories…"),
                            html.Button("Clear search", id="run-stats-search-clear", n_clicks=0),
                        ], className="run-stats-search-controls"),
                        html.Div(id="run-stats-search-results", role="status", **{"aria-live": "polite"}),
                        html.Div("Category checkboxes select every variable below them, including hidden variables.",
                                 className="run-stats-editor-note"),
                    ], className="run-stats-search"),
                    html.Div(node("all"), id="run-stats-tree", className="run-stats-tree"),
                ], className="run-stats-tree-pane"),
                html.Div([
                    html.Div(id="run-stats-preview", className="run-stats-preview", role="status"),
                    html.Div("Counts are registry definitions; species and scalar templates expand during a run.",
                             className="run-stats-editor-note"),
                    html.Div([
                        html.Label("Save for later", htmlFor="run-stats-save-name"),
                        dcc.Input(id="run-stats-save-name", value="", placeholder="e.g. core_radiation", type="text"),
                        html.Button("Save stats config", id="run-stats-save", n_clicks=0),
                        html.Div("Saved lists appear alongside the legacy lists and can be reused from the CLI.",
                                 className="run-stats-editor-note"),
                        html.Div(id="run-stats-save-feedback", role="status"),
                    ], className="run-stats-save-form"),
                ], className="run-stats-editor-side"),
            ], className="run-stats-editor-body"),
            html.Div([
                html.Span("Changes apply to your next run.", className="run-stats-editor-note"),
                html.Button("Use selection", id="run-stats-apply", n_clicks=0, className="run-stats-apply"),
            ], className="run-stats-editor-footer"),
        ], className="shared-notecard-panel run-stats-editor-panel", role="dialog",
           **{"aria-modal": "true", "aria-labelledby": "run-stats-editor-title"}),
    ], id="run-stats-modal", className=HIDDEN)


def register_stats_callbacks(app):
    @app.callback(
        Output("run-stats-summary", "children"),
        Output("run-stats-list", "value"),
        Output("run-selected-stats-file", "data"),
        Output("run-stats-seed", "data"),
        Output("run-stats-custom-open", "disabled"),
        Input({"type": "run-stats-button", "name": ALL}, "n_clicks"),
        Input("run-stats-apply", "n_clicks"),
        Input("run-stats-list", "value"),
        State("run-selected-stats-file", "data"), State("run-stats-draft", "data"),
    )
    def describe_selection(_presets, _apply, list_value, selection, draft):
        trigger = callback_context.triggered_id
        if isinstance(trigger, dict):
            selection = trigger["name"]
        elif trigger == "run-stats-apply":
            selection = stats_selection_expression(draft or [])
        elif trigger == "run-stats-list" and list_value:
            selection = list_value
        elif trigger == "run-stats-list":
            raise PreventUpdate
        selection = selection or "standard"
        try:
            names = stats_selection_names(selection)
            available = {entry["name"] for entry in load_stats_categories()["all"]}
        except (ValueError, OSError) as exc:
            return str(exc), None, no_update, no_update, True
        saved = selection in list_stats_files()
        label = selection if selection in {"all", "standard", "core"} else "Custom selection"
        text = "Statistics disabled" if selection.lower() == "none" else f"{label.capitalize()} · {len(names)} stats"
        seed = {"names": sorted(names & available), "missing": sorted(names - available)}
        return text, selection if saved else None, selection, seed, False

    @app.callback(Output("run-stats-list", "options"), Input("run-stats-custom-open", "n_clicks"))
    def refresh_lists(_clicks):
        return stats_list_options()

    app.clientside_callback(
        ClientsideFunction(namespace="statsSelection", function_name="editSelection"),
        Output("run-stats-modal", "className"),
        Output("run-stats-draft", "data"),
        Output({"type": "run-stats-category", "path": ALL}, "value"),
        Output({"type": "run-stats-variables", "path": ALL}, "value"),
        Output("run-stats-preview", "children"),
        Input("run-stats-custom-open", "n_clicks"),
        Input("run-stats-custom-cancel", "n_clicks"),
        Input("run-stats-apply", "n_clicks"),
        Input("run-stats-clear", "n_clicks"),
        Input({"type": "run-stats-category", "path": ALL}, "value"),
        Input({"type": "run-stats-variables", "path": ALL}, "value"),
        State({"type": "run-stats-category", "path": ALL}, "id"),
        State({"type": "run-stats-variables", "path": ALL}, "id"),
        State("run-stats-draft", "data"),
        State("run-stats-catalog", "data"),
        State("run-stats-seed", "data"),
        prevent_initial_call=True,
    )

    @app.callback(
        Output("run-stats-save-feedback", "children"),
        Output("run-stats-list", "options", allow_duplicate=True),
        Input("run-stats-save", "n_clicks"),
        Input("run-stats-custom-open", "n_clicks"),
        State("run-stats-save-name", "value"), State("run-stats-draft", "data"),
        prevent_initial_call=True,
    )
    def save_selection(_clicks, _open, name, names):
        if callback_context.triggered_id == "run-stats-custom-open":
            return "", no_update
        try:
            saved = save_stats_config(name, stats_selection_expression(names or []))
        except (ValueError, OSError) as exc:
            return html.Div(str(exc), className="run-stats-feedback-error"), no_update
        return html.Div(f"Saved {saved.rsplit('/', 1)[-1]}. Choose it in Saved & legacy lists.",
                        className="run-stats-feedback-success"), stats_list_options()

    app.clientside_callback(
        ClientsideFunction(namespace="statsSelection", function_name="filterTree"),
        Output("run-stats-search", "value"),
        Output("run-stats-search-results", "children"),
        Input("run-stats-search", "value"),
        Input("run-stats-search-clear", "n_clicks"),
        Input("run-stats-custom-open", "n_clicks"),
    )

    app.clientside_callback(
        ClientsideFunction(namespace="statsSelection", function_name="syncMixedCheckboxes"),
        Output("run-stats-checkbox-signal", "data"),
        Input("run-stats-draft", "data"), State("run-stats-catalog", "data"),
    )
