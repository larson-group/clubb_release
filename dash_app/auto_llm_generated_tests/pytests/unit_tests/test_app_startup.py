from dash import Dash, dcc, html
from dash._utils import to_json

from dash_app import app
from dash_app.lazy_tabs import LazyTabs
from dash_app.persistence import WORKSPACE_TOKEN


def test_dashboard_generation_route_is_process_specific_and_not_cached():
    dash = Dash(__name__)
    dash.layout = html.Div()
    app._register_dashboard_generation_route(dash)

    response = dash.server.test_client().get("/_clubb-dashboard-generation")

    assert response.status_code == 200
    assert response.get_json() == {"generation": app.DASH_GENERATION}
    assert response.headers["Cache-Control"] == "no-store"
    assert f'window.__CLUBB_DASH_GENERATION__ = "{app.DASH_GENERATION}"' in dash.index_string


def test_main_reopens_existing_dashboard_without_starting_second_dash(monkeypatch):
    opened = {}
    monkeypatch.setattr(app.sys, "argv", ["app.py"])
    monkeypatch.setattr(app, "_resolve_port", lambda _host, _requested: (23404, True, False))
    monkeypatch.setattr(app, "ensure_broker", lambda **_kwargs: None)
    monkeypatch.setattr(app, "_reuse_existing_dashboard", lambda _host: opened.update(reused=True) or True)
    monkeypatch.setattr(app, "_open_browser", lambda url, *, new: opened.update(url=url, new=new))

    def unexpected_dash_start(*_args, **_kwargs):
        raise AssertionError("a second Dash app should not be constructed")

    monkeypatch.setattr(app, "Dash", unexpected_dash_start)

    app.main()

    assert opened == {"reused": True}


def test_restart_runtime_refuses_to_replace_live_dashboard(monkeypatch):
    import pytest

    monkeypatch.setattr(app.sys, "argv", ["app.py", "--restart-runtime"])
    monkeypatch.setattr(app, "_resolve_port", lambda _host, _requested: (23404, True, False))
    monkeypatch.setattr(app, "_existing_dashboard_url", lambda _host: "http://127.0.0.1:23404")

    with pytest.raises(SystemExit, match="dashboard is still live"):
        app.main()


def test_reuse_existing_dashboard_opens_registered_url(monkeypatch):
    opened = {}
    monkeypatch.setattr(app, "_existing_dashboard_url", lambda _host: "http://127.0.0.1:23407")
    monkeypatch.setattr(app, "_open_browser", lambda url, *, new: opened.update(url=url, new=new))

    assert app._reuse_existing_dashboard("127.0.0.1") is True
    assert opened == {"url": "http://127.0.0.1:23407", "new": 0}


def test_existing_dashboard_url_uses_registered_dashboard(monkeypatch):
    monkeypatch.setattr(
        app.client,
        "connect",
        lambda: {"dashboard": {"status": "available", "port": 23407}},
    )

    assert app._existing_dashboard_url("127.0.0.1") == "http://127.0.0.1:23407"


def test_existing_dashboard_url_returns_none_without_live_registration(monkeypatch):
    monkeypatch.setattr(
        app.client,
        "connect",
        lambda: {"dashboard": {"status": "unavailable"}},
    )

    assert app._existing_dashboard_url("127.0.0.1") is None


def test_existing_dashboard_url_handles_broker_connection_failure(monkeypatch):
    monkeypatch.setattr(app.client, "connect", lambda: (_ for _ in ()).throw(RuntimeError("offline")))

    assert app._existing_dashboard_url("127.0.0.1") is None


def test_initial_dashboard_defers_discovery_and_scientific_controls(monkeypatch):
    from dash_app.compile_tab import tab as compile_tab
    from dash_app.profile_tab import tab as profile_tab
    from dash_app.run_tab import tab as run_tab
    from dash_app.tune_tab import tab as tune_tab
    from dash_app.tutorial_tab import layout as tutorial_layout
    from dash_app.misc_tab.registry import SubtabSpec
    from dash_app.misc_tab import tab as misc_tab

    def unexpected():
        raise AssertionError("unvisited page was built at startup")

    for module, name in (
        (compile_tab, "discover_compile_state"),
        (profile_tab, "discover_profile_state"),
        (run_tab, "build_initial_run_state"),
        (tune_tab, "build_initial_tune_state"),
        (tutorial_layout, "build_equations_page"),
        (tutorial_layout, "build_adg1_explorer_layout"),
    ):
        monkeypatch.setattr(module, name, unexpected)
    monkeypatch.setattr(misc_tab, "discover_subtabs", lambda: (
        SubtabSpec(slug="diagnostic", title="Diagnostic", summary="", build_layout=unexpected),
    ))
    dash = Dash(__name__, suppress_callback_exceptions=True)
    dash.layout = app.build_dashboard_tabs(dash)
    payload = to_json(dash.layout)

    assert "tutorial-start-equations" in payload
    assert "compile-env-select" not in payload
    assert "compile-build-selector-anchor" not in payload
    assert "run-selected-cases" not in payload
    assert "plots-case-data" not in payload
    assert "notes-adg1-gaussian-figure" not in payload
    assert len(payload) < 20000
    # All callbacks must be known before the browser fetches dependencies.
    assert "compile-discovery.data" in dash.callback_map
    assert "plots-case-data.data" in str(dash.callback_map)


def test_lazy_group_loads_atomically_with_workspace_persistence():
    dash = Dash(__name__, suppress_callback_exceptions=True)
    lazy = LazyTabs("pages", groups=(("one", "two"),))
    builds = []

    def layout(name):
        builds.append(name)
        return html.Div(dcc.Input(id=f"{name}-value", value="default"))

    tabs = [lazy.tab(label=name, value=name, build=lambda name=name: layout(name))
            for name in ("one", "two")]
    loaders = lazy.register(dash, shared={
        "one": ("shared-overlay", lambda: dcc.Store(id="plots-output-dirs", data=[])),
    })
    dash.layout = html.Div([loaders, dcc.Tabs(tabs, id="pages", value="one")])
    assert builds == []
    key, entry = next((key, entry) for key, entry in dash.callback_map.items()
                      if "pages-one-tab.children" in key)
    client = dash.server.test_client()
    response = client.post("/_dash-update-component", json={
        "output": key,
        "outputs": [output.to_dict() for output in entry["output"]],
        "inputs": [{"id": "pages-lazy-one-request", "property": "data", "value": True}],
        "state": [],
        "changedPropIds": ["pages-lazy-one-request.data"],
    })
    assert response.status_code == 200
    assert builds == ["one", "two"]
    payload = response.get_json()["response"]
    assert payload["pages-lazy-one-loaded"]["data"] is True
    for name in ("one", "two"):
        control = payload[f"pages-{name}-tab"]["children"]["props"]["children"]
        assert control["props"]["persistence"] == WORKSPACE_TOKEN
        assert control["props"]["persistence_type"] == "local"
    assert payload["shared-overlay"]["children"]["props"]["storage_type"] == "local"


def test_browser_gate_loads_restored_page_once_and_requires_active_parent():
    import shutil
    import subprocess

    import pytest

    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is needed to check the browser gate")
    dash = Dash(__name__)
    lazy = LazyTabs("lessons", parent=("dashboard-tabs", "tutorial"))
    lazy.tab(value="equations", label="Equations", build=lambda: html.Div())
    lazy.register(dash)
    top = LazyTabs("tabs", handoff_id="dashboard-request")
    top.tab(value="plots", label="Plots", build=lambda: html.Div())
    top.register(dash)
    script = "const window = {dash_clientside: {no_update: Symbol('no-update')}};\n"
    script += "\n".join(dash._inline_scripts)
    script += """
        const gate = Object.values(window.dash_clientside._dashprivate_clientside_funcs)[0];
        const noUpdate = window.dash_clientside.no_update;
        const assert = require('node:assert/strict');
        assert.equal(gate('equations', 'misc', false, null), noUpdate);
        assert.equal(gate('equations', 'tutorial', false, null), true);
        assert.equal(gate('equations', 'tutorial', false, true), noUpdate);
        assert.equal(gate('equations', 'tutorial', true, true), noUpdate);
        assert.equal(gate('welcome', 'tutorial', false, null), noUpdate);
        const topGate = Object.values(window.dash_clientside._dashprivate_clientside_funcs)[1];
        assert.equal(topGate('tutorial', null, false, null), noUpdate);
        assert.equal(topGate('tutorial', {tab: 'plots', preserve_tab: true}, false, null), true);
        assert.equal(topGate('tutorial', {tab: 'plots'}, false, true), noUpdate);
    """
    subprocess.run([node, "-e", script], check=True, capture_output=True, text=True)
