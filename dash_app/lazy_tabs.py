"""Load tab bodies once per browser page, after their first selection.

Register callbacks up front, but keep unvisited controls out of Dash's layout:
children of an unselected dcc.Tab otherwise still fire initial callbacks.
Related tabs can load together when callbacks share their controls. Loaded
bodies stay in the layout so navigation does not reset runtime or form state.
"""

from __future__ import annotations

import json

from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate

from dash_app.persistence import enable_workspace_persistence


class LazyTabs:
    def __init__(self, navigation_id, *, parent=None, groups=(), handoff_id=None):
        self.navigation_id = navigation_id
        self.parent = parent
        self.groups = groups
        self.handoff_id = handoff_id
        self.pages = {}

    def tab(self, *, value, label, build, **props):
        tab_id = props.pop("id", f"{self.navigation_id}-{value}-tab")
        self.pages[value] = (tab_id, build)
        return dcc.Tab(
            id=tab_id, value=value, label=label,
            children=html.Div(
                f"Loading {label}…", role="status", style={"padding": "24px"},
            ),
            **props,
        )

    def register(self, app, *, shared=None):
        """Return loader stores and optional overlays mounted with each group.

        ``shared`` maps a group's first value to an (element id, builder) pair.
        """
        shared = shared or {}
        grouped = {value for group in self.groups for value in group}
        groups = [*self.groups, *((value,) for value in self.pages if value not in grouped)]
        stores = []
        for group in groups:
            prefix = f"{self.navigation_id}-lazy-{group[0]}"
            request_id, loaded_id = f"{prefix}-request", f"{prefix}-loaded"
            stores.extend([dcc.Store(id=request_id), dcc.Store(id=loaded_id, data=False)])
            inputs = [Input(self.navigation_id, "value")]
            arguments = ["value"]
            if self.parent:
                inputs.append(Input(self.parent[0], "value"))
                arguments.append("parent")
            if self.handoff_id:
                inputs.append(Input(self.handoff_id, "data"))
                arguments.append("handoff")
            arguments.extend(["loaded", "requested"])
            active_parent = f"parent === {json.dumps(self.parent[1])}" if self.parent else "true"
            handoff = " || pages.includes((handoff || {}).tab)" if self.handoff_id else ""
            # Gate in the browser: unvisited pages make no HTTP requests, and
            # rapid navigation cannot discard another page's in-flight load.
            app.clientside_callback(
                f"""function({', '.join(arguments)}) {{
                    const pages = {json.dumps(list(group))};
                    const active = (pages.includes(value) && {active_parent}){handoff};
                    return active && !loaded && !requested ? true : window.dash_clientside.no_update;
                }}""",
                Output(request_id, "data"), *inputs,
                State(loaded_id, "data"), State(request_id, "data"),
            )
            targets = [self.pages[value] for value in group]
            if group[0] in shared:
                targets.append(shared[group[0]])
            self._register_loader(app, targets, request_id, loaded_id)
        # Shared overlays must live outside dcc.Tab children so their DOM stays
        # available when any member of the group is visible.
        return html.Div([
            html.Div(stores, style={"display": "none"}),
            *(html.Div(id=component_id) for component_id, _ in shared.values()),
        ])

    def _register_loader(self, app, targets, request_id, loaded_id):
        @app.callback(
            *(Output(component_id, "children") for component_id, _ in targets),
            Output(loaded_id, "data"),
            Input(request_id, "data"),
            prevent_initial_call=True,
        )
        def load(requested):
            if not requested:
                raise PreventUpdate
            return [
                *(enable_workspace_persistence(build()) for _, build in targets),
                True,
            ]
