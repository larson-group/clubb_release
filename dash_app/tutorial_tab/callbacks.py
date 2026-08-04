"""Navigation callbacks for the Tutorial tab."""

from dash import Input, Output, ctx, no_update


def register_callbacks(app):
    @app.callback(
        Output("tutorial-pages", "value"),
        Input("tutorial-start-adg1-explorer", "n_clicks"),
        Input("tutorial-start-equations", "n_clicks"),
        Input("dashboard-request", "data"),
        prevent_initial_call=True,
    )
    def _open_suggested_tutorial(_adg1_explorer_clicks, _equations_clicks, agent_request):
        if ctx.triggered_id == "dashboard-request":
            if (agent_request or {}).get("tab") != "tutorial" or (agent_request or {}).get("operation") != "open_lesson":
                return no_update
            lesson = str((agent_request or {}).get("lesson") or "")
            if lesson in {"tutorial-welcome", "tutorial-equations", "tutorial-adg1-explorer"}:
                return lesson
            return no_update
        if ctx.triggered_id == "tutorial-start-equations":
            return "tutorial-equations"
        return "tutorial-adg1-explorer"
