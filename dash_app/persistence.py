"""Versioned browser-local persistence for dashboard user workspace state."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from dash import dcc, html
from dash.development.base_component import Component


WORKSPACE_SCHEMA_VERSION = 1
REPO_NAME = Path(__file__).resolve().parents[1].name
WORKSPACE_TOKEN = f"{REPO_NAME}:dashboard-workspace:v{WORKSPACE_SCHEMA_VERSION}"
WORKSPACE_SCHEMA_KEY = f"{REPO_NAME}:dashboard-workspace:schema"

# Only user intent belongs in durable Stores. Runtime jobs, logs, derived
# NetCDF metadata, generated results, and process state must be reconstructed.
PERSISTENT_STORE_IDS = frozenset(
    {
        "theme-store",
        # Run configuration.
        "run-selected-config",
        "run-multicol-rows-state",
        "run-multicol-next-id",
        "run-multicol-row-order",
        "run-selected-cases",
        "run-selected-stats-file",
        # Plot workspace. plots-case-data is intentionally included because it
        # owns case selection; directory callbacks refresh its derived fields.
        "plots-output-dirs",
        "plots-case-data",
        "plots-enabled-benchmarks",
        "plots-plot-order",
        "plots-plot-state",
        "plots-next-id",
        "plots-column-filters",
        "plots-selected-column",
        "plots-time-override",
        # A saved Tune selection is the only persisted Tune state.  Its
        # request, controls, results, status, and active process must be
        # reloaded from that durable workspace; persisting individual form
        # controls can overwrite freshly hydrated case rows after a refresh.
        "tune-workspace-selection",
    }
)

EXTRA_LOCAL_STORAGE_KEYS = (
    "plots-ui-right-width",
    "run-ui-right-width",
    "tune-ui-right-width",
)


def workspace_metadata() -> dict[str, Any]:
    return {
        "schema_version": WORKSPACE_SCHEMA_VERSION,
        "token": WORKSPACE_TOKEN,
        "repo_name": REPO_NAME,
        "schema_key": WORKSPACE_SCHEMA_KEY,
        "store_ids": sorted(PERSISTENT_STORE_IDS),
        "extra_keys": list(EXTRA_LOCAL_STORAGE_KEYS),
    }


def _children(component: Component) -> list[Any]:
    children = getattr(component, "children", None)
    if children is None:
        return []
    return list(children) if isinstance(children, (list, tuple)) else [children]


def enable_workspace_persistence(root: Component) -> Component:
    """Enable local persistence on user controls and selected configuration stores."""

    stack: list[Any] = [root]
    while stack:
        component = stack.pop()
        if not isinstance(component, Component):
            continue
        component_id = getattr(component, "id", None)
        properties = set(getattr(component, "_prop_names", ()))
        is_tune_control = (
            isinstance(component_id, str)
            and component_id.startswith("tune-")
            or isinstance(component_id, dict)
            and str(component_id.get("type") or "").startswith("tune-")
        )
        # Tune is workspace-backed, unlike the free-form Run/Plot surfaces.
        # Do not let Dash's automatic component persistence restore stale
        # values into a saved revision after its loader rebuilds the controls.
        if component_id is not None and "persistence" in properties and not is_tune_control:
            component.persistence = WORKSPACE_TOKEN
            if "persistence_type" in properties:
                component.persistence_type = "local"
        if (
            isinstance(component, dcc.Store)
            and isinstance(component_id, str)
            and component_id in PERSISTENT_STORE_IDS
        ):
            component.storage_type = "local"
        stack.extend(_children(component))
    return root


def workspace_toolbar() -> Component:
    return html.Div(
        [
            dcc.Store(id="dashboard-workspace-meta", data=workspace_metadata()),
            html.Span("Refresh-safe workspace", className="dashboard-workspace-label"),
            html.Button(
                "Export workspace",
                id="dashboard-workspace-export",
                n_clicks=0,
                className="dashboard-workspace-button",
            ),
            dcc.Upload(
                id="dashboard-workspace-import",
                children=html.Button(
                    "Import workspace",
                    className="dashboard-workspace-button",
                ),
                accept="application/json,.json",
                multiple=False,
            ),
            html.Button(
                "Reset saved workspace",
                id="dashboard-workspace-reset",
                n_clicks=0,
                className="dashboard-workspace-button dashboard-workspace-reset",
            ),
            html.Span(
                "Cases, settings, flags, and plot layouts survive refreshes.",
                className="dashboard-workspace-detail",
            ),
            html.Span(id="dashboard-workspace-ready-status", className="dashboard-workspace-status"),
            html.Span(id="dashboard-workspace-export-status", className="dashboard-workspace-status"),
            html.Span(id="dashboard-workspace-import-status", className="dashboard-workspace-status"),
            html.Span(id="dashboard-workspace-reset-status", className="dashboard-workspace-status"),
        ],
        className="dashboard-workspace-toolbar",
    )


def mcp_endpoint_panel(details: dict[str, Any] | None) -> Component:
    """Show the stable broker MCP URL and bearer token for manual setup."""
    details = dict(details or {})
    endpoint_url = str(details.get("endpoint_url") or "")
    endpoint_token = str(details.get("endpoint_token") or "")
    if endpoint_url and endpoint_token:
        content = [
            html.Div("Add this dashboard as a Streamable HTTP MCP server in the chat you want to use.", className="dashboard-mcp-endpoint-help"),
            html.Label("Server URL", className="dashboard-mcp-endpoint-label"),
            dcc.Input(id="dashboard-mcp-endpoint-url", value=endpoint_url, readOnly=True, className="dashboard-mcp-endpoint-field", type="text"),
            html.Label("Bearer token", className="dashboard-mcp-endpoint-label"),
            dcc.Input(id="dashboard-mcp-endpoint-token", value=endpoint_token, readOnly=True, className="dashboard-mcp-endpoint-field", type="text"),
            html.Div(f"Instance: {details.get('instance_id', '')}", className="dashboard-mcp-endpoint-instance"),
            html.Div("The endpoint is loopback-only and survives ordinary Dash restarts. Browser actions require one live dashboard; broker-owned jobs do not.", className="dashboard-mcp-endpoint-help"),
        ]
    else:
        content = [
            html.Div("The stable broker MCP endpoint will appear after the local broker migration is ready.", className="dashboard-mcp-endpoint-help")
        ]
    return html.Div(
        [html.Div("MCP connection", className="dashboard-utility-section-title"), *content],
        id="dashboard-mcp-endpoint",
        className="dashboard-mcp-endpoint",
    )


def utility_drawer(endpoint_details: dict[str, Any] | None = None) -> Component:
    """Return the compact bottom-left menu for global dashboard utilities.

    The workspace controls retain their existing component IDs and callbacks;
    only their presentation moves out of the main document flow.
    """

    return html.Div(
        [
            dcc.Store(id="dashboard-utility-open", data=False),
            html.Button(
                "☰",
                id="dashboard-utility-toggle",
                n_clicks=0,
                className="dashboard-utility-toggle",
                title="Dashboard utilities",
                **{"aria-label": "Open dashboard utilities"},
            ),
            html.Div(
                [
                    html.Div("Workspace", className="dashboard-utility-section-title"),
                    workspace_toolbar(),
                    mcp_endpoint_panel(endpoint_details),
                    html.Button(
                        "Theme: Dark",
                        id="theme-toggle-button",
                        n_clicks=0,
                        className="dashboard-utility-item",
                    ),
                ],
                id="dashboard-utility-drawer",
                className="dashboard-utility-drawer dashboard-utility-drawer-closed",
            ),
        ]
    )


def metadata_json() -> str:
    """Stable representation used by focused tests and future migrations."""

    return json.dumps(workspace_metadata(), sort_keys=True)
