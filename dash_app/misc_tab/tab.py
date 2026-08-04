"""Entrypoint for the consolidated Misc dashboard tab."""

from dash import Input, Output, dcc, no_update

from dash_app.misc_tab.registry import discover_subtabs, register_subtab_callbacks

from .layout import build_layout, subtab_page_value


def build_tab(app):
    """Register collected tool callbacks and build the final top-level tab."""
    subtabs = discover_subtabs()
    register_subtab_callbacks(app, subtabs)

    @app.callback(
        Output("misc-pages", "value", allow_duplicate=True),
        Input("dashboard-request", "data"),
        prevent_initial_call=True,
    )
    def select_agent_misc_page(request):
        if (request or {}).get("tab") != "misc":
            return no_update
        operation = str((request or {}).get("operation") or "")
        if operation == "open_report":
            slug = str((request or {}).get("report_slug") or "")
            selected = next((subtab for subtab in subtabs if subtab.slug == slug), None)
            return subtab_page_value(selected) if selected is not None else no_update
        return no_update

    return dcc.Tab(label="Misc", value="misc", children=build_layout(subtabs))
