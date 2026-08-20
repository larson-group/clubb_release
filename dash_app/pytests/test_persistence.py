"""Focused persistence policy checks for workspace-backed Tune controls."""

from dash import dcc, html

from dash_app.persistence import WORKSPACE_TOKEN, enable_workspace_persistence, mcp_endpoint_panel


def test_tune_form_values_are_not_browser_persisted():
    """A loaded Tune revision must win over stale browser form values."""
    root = html.Div(
        [
            dcc.Input(id="tune-random-max-samples", value=2000),
            dcc.Input(id={"type": "tune-case-name", "index": 0}, value="arm"),
            dcc.Input(id="run-some-control", value="kept"),
            dcc.Store(id="tune-workspace-selection", data={}),
        ]
    )

    enable_workspace_persistence(root)
    tune_scalar, tune_pattern, run_control, selection = root.children

    assert getattr(tune_scalar, "persistence", None) is None
    assert getattr(tune_pattern, "persistence", None) is None
    assert getattr(run_control, "persistence", None) == WORKSPACE_TOKEN
    assert selection.storage_type == "local"


def test_run_config_persists_only_the_user_selection():
    """Derived config metadata and rendered-control state must be rebuilt."""
    root = html.Div(
        [
            dcc.Store(id="run-selected-config", data="before"),
            dcc.Store(id="run-tunable-configs", data=[]),
            dcc.Store(id="run-rendered-config", data="default"),
        ]
    )

    enable_workspace_persistence(root)
    selected, configs, rendered = root.children

    assert selected.storage_type == "local"
    assert getattr(configs, "storage_type", "memory") == "memory"
    assert getattr(rendered, "storage_type", "memory") == "memory"


def test_mcp_endpoint_panel_exposes_manual_url_and_bearer_token():
    panel = mcp_endpoint_panel(
        {
            "instance_id": "a" * 32,
            "endpoint_url": "http://127.0.0.1:41236/mcp",
            "endpoint_token": "endpoint-secret",
        }
    )
    fields = {child.id: child.value for child in panel.children if getattr(child, "id", None)}
    assert fields["dashboard-mcp-endpoint-url"] == "http://127.0.0.1:41236/mcp"
    assert fields["dashboard-mcp-endpoint-token"] == "endpoint-secret"
