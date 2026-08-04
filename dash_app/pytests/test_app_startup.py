from dash_app import app


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
