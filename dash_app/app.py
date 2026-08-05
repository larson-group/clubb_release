#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import importlib.util
import logging
import os
from pathlib import Path
import pkgutil
import socket
import sys
import threading
import webbrowser

from dash import ClientsideFunction, Dash, dcc, html, Input, Output, State

DASH_APP_ROOT = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(DASH_APP_ROOT, ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from dash_app.agent_integration import client
from dash_app.agent_integration.handoff import dashboard_handoff, register_dashboard_handoff_callbacks
from dash_app.agent_integration.dashboard_lifecycle import DashboardRegistration
from dash_app.compile_tab.tab import build_tab as build_compile_tab
from dash_app.misc_tab.tab import build_tab as build_misc_tab
from dash_app.persistence import enable_workspace_persistence, utility_drawer
from dash_app.plot_tab.tab import build_tab as build_plots_tab
from dash_app.reports_tab.static import register_static_report_routes
from dash_app.reports_tab.tab import build_tab as build_reports_tab
from dash_app.run_tab.tab import build_tab as build_run_tab
from dash_app.shared.background import create_background_manager
from dash_app.shared.broker import ensure_broker
from dash_app.shared.gateway import BROKER_LOG_PATH, read_connection, update_connection_logs
from dash_app.shared.runtime import private_path
from dash_app.shared.runtime_logging import dashboard_log_path
from dash_app.shared.provenance import runtime_source_fingerprint
from dash_app.tune_tab.tab import build_tab as build_tune_tab
from dash_app.tutorial_tab.tab import build_tab as build_tutorial_tab


if not hasattr(pkgutil, "find_loader"):
    # Dash 3.2.0 still calls pkgutil.find_loader in debug mode; Python 3.14 removed it
    def _find_loader(name):
        spec = importlib.util.find_spec(name)
        return None if spec is None else spec.loader

    pkgutil.find_loader = _find_loader


DEFAULT_PORT = 23404
SELECTED_PORT_ENV = "CLUBB_DASH_SELECTED_PORT"
BROWSER_OPENED_ENV = "CLUBB_DASH_BROWSER_OPENED"
def _port_is_available(host: str, port: int) -> bool:
    addrinfos = socket.getaddrinfo(
        host,
        port,
        family=socket.AF_UNSPEC,
        type=socket.SOCK_STREAM,
        flags=socket.AI_PASSIVE,
    )
    for family, socktype, proto, _canonname, sockaddr in addrinfos:
        try:
            with socket.socket(family, socktype, proto) as sock:
                sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                sock.bind(sockaddr)
        except OSError:
            continue
        return True
    return False


def _find_available_port(host: str, start_port: int) -> int:
    port = start_port
    while not _port_is_available(host, port):
        port += 1
    return port


def _browser_host(host: str) -> str:
    if host in {"0.0.0.0", "::", ""}:
        return "127.0.0.1"
    return host


def _app_title() -> str:
    """Return the browser title for the current CLUBB checkout."""
    return os.path.basename(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def _should_open_browser(debug: bool) -> bool:
    if os.environ.get(BROWSER_OPENED_ENV) == "1":
        return False
    if not debug:
        return True
    return os.environ.get("WERKZEUG_RUN_MAIN") != "true"


def _resolve_port(host: str, requested_port: int | None) -> tuple[int, bool, bool]:
    """Resolve the Dash port once and preserve it across Werkzeug reloads."""
    if requested_port is not None:
        os.environ[SELECTED_PORT_ENV] = str(requested_port)
        return requested_port, False, False

    inherited_port = os.environ.get(SELECTED_PORT_ENV)
    if inherited_port:
        try:
            return int(inherited_port), True, True
        except ValueError:
            pass

    port = _find_available_port(host, DEFAULT_PORT)
    os.environ[SELECTED_PORT_ENV] = str(port)
    return port, True, False


def _open_browser(url: str, *, new: int = 2) -> None:
    def _launch():
        try:
            webbrowser.open(url, new=new)
        except Exception:
            pass

    threading.Timer(1.0, _launch).start()


def _existing_dashboard_url(host: str) -> str | None:
    """Return the live checkout dashboard URL, if one is already registered."""
    try:
        status = client.connect()
    except (OSError, RuntimeError, ValueError):
        return None
    dashboard = status.get("dashboard") or {}
    if dashboard.get("status") != "available":
        return None
    try:
        port = int(dashboard["port"])
    except (KeyError, TypeError, ValueError):
        return None
    return f"http://{_browser_host(host)}:{port}"


def _reuse_existing_dashboard(host: str) -> bool:
    """Open the current dashboard, returning whether startup should stop."""
    dashboard_url = _existing_dashboard_url(host)
    if not dashboard_url:
        return False
    print(f"CLUBB Dash is already running at {dashboard_url}; opening that dashboard.")
    _open_browser(dashboard_url, new=0)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Dash app with multiple tabs for CLUBB NetCDF analysis."
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=None)
    parser.add_argument(
        "-debug",
        action="store_true",
        help="Run Dash in debug mode.",
    )
    parser.add_argument(
        "--threaded",
        action="store_true",
        help=(
            "Enable concurrent request handling. This is unsafe with this dashboard's "
            "NetCDF/HDF5 readers and is intended only for short diagnostic sessions."
        ),
    )
    parser.add_argument(
        "--single-threaded",
        action="store_false",
        dest="threaded",
        help="Serialize Dash request handling (the safe default).",
    )
    parser.add_argument(
        "--restart-runtime",
        action="store_true",
        help="Safely replace an idle stale broker/MCP endpoint; requires no live dashboard or active jobs.",
    )
    args = parser.parse_args()
    port, auto_port, inherited_port = _resolve_port(args.host, args.port)

    if args.restart_runtime:
        if _existing_dashboard_url(args.host):
            raise SystemExit("A dashboard is still live; close it before retrying --restart-runtime.")
    elif _reuse_existing_dashboard(args.host):
        return

    # Starting the launcher again should focus/reopen the current checkout's
    # dashboard instead of creating a second Dash process that the broker must
    # reject.  ``new=0`` lets the browser reuse an existing window or tab; if
    # none exists, the browser decides how to create one.
    connection_path = ensure_broker(
        expected_fingerprint=runtime_source_fingerprint(REPO_ROOT),
        require_fresh=args.restart_runtime,
    )

    # HDF5-backed NetCDF access is not reliably thread-safe in the scientific
    # stack used here.  Keep Flask single-threaded and move explicitly marked
    # expensive callbacks into Diskcache worker *processes* instead.
    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        eager_loading=True,
        title=_app_title(),
        update_title=None,
        background_callback_manager=create_background_manager(REPO_ROOT),
    )
    register_static_report_routes(app)

    tabs = [
        build_tutorial_tab(app),
        build_compile_tab(app),
        build_run_tab(app),
        build_tune_tab(app),
        build_plots_tab(app),
        build_reports_tab(app),
        build_misc_tab(app),
    ]

    broker_connection = read_connection()
    update_connection_logs(
        {
            **dict(broker_connection.get("log_paths") or {}),
            "broker": str(BROKER_LOG_PATH),
            "endpoint": str(private_path(REPO_ROOT, "mcp-endpoint/endpoint.log")),
            "app": str(dashboard_log_path(REPO_ROOT)),
        }
    )
    dashboard_registration = None
    if not args.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        endpoint_details = dict(broker_connection.get("mcp_endpoint") or {})
        if endpoint_details.get("endpoint_url") and endpoint_details.get("endpoint_token"):
            try:
                dashboard_registration = DashboardRegistration.register(port=port)
            except client.DashboardAlreadyRegisteredError:
                # A second process can race the preflight above.  Reuse the
                # winning dashboard if it registered between the two checks.
                if not _reuse_existing_dashboard(args.host):
                    raise
                return
            atexit.register(dashboard_registration.stop)
        else:
            endpoint_details = {
                "status": "waiting for the broker endpoint migration to finish",
            }
    else:
        endpoint_details = {"status": "waiting for the debug reloader"}

    app.layout = html.Div(
        [
            dcc.Store(id="theme-store", data="dark"),
            dashboard_handoff(_app_title()),
            dcc.Tabs(tabs, id="dashboard-tabs", value="tutorial"),
            utility_drawer(endpoint_details),
        ],
        id="app-root",
        className="theme-dark",
        style={"padding": "0px", "minHeight": "100vh"},
    )

    @app.callback(
        Output("theme-store", "data"),
        Input("theme-toggle-button", "n_clicks"),
        State("theme-store", "data"),
        prevent_initial_call=True,
    )
    def _toggle_theme(_n_clicks, current_theme):
        return "light" if (current_theme or "dark") == "dark" else "dark"

    @app.callback(
        Output("dashboard-utility-drawer", "className"),
        Output("dashboard-utility-open", "data"),
        Input("dashboard-utility-toggle", "n_clicks"),
        State("dashboard-utility-open", "data"),
        prevent_initial_call=True,
    )
    def _toggle_utility_drawer(_clicks, is_open):
        next_open = not bool(is_open)
        return (
            "dashboard-utility-drawer dashboard-utility-drawer-open"
            if next_open
            else "dashboard-utility-drawer dashboard-utility-drawer-closed",
            next_open,
        )

    @app.callback(
        Output("app-root", "className"),
        Output("theme-toggle-button", "children"),
        Input("theme-store", "data"),
    )
    def _apply_theme(theme_name):
        theme = (theme_name or "dark").lower()
        if theme == "light":
            return "theme-light", "Theme: Light"
        return "theme-dark", "Theme: Dark"

    register_dashboard_handoff_callbacks(app)

    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace", function_name="initializeWorkspace"
        ),
        Output("dashboard-workspace-ready-status", "children"),
        Input("dashboard-workspace-meta", "data"),
    )
    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace", function_name="exportWorkspace"
        ),
        Output("dashboard-workspace-export-status", "children"),
        Input("dashboard-workspace-export", "n_clicks"),
        State("dashboard-workspace-meta", "data"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace", function_name="importWorkspace"
        ),
        Output("dashboard-workspace-import-status", "children"),
        Input("dashboard-workspace-import", "contents"),
        State("dashboard-workspace-import", "filename"),
        State("dashboard-workspace-meta", "data"),
        prevent_initial_call=True,
    )
    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardWorkspace", function_name="resetWorkspace"
        ),
        Output("dashboard-workspace-reset-status", "children"),
        Input("dashboard-workspace-reset", "n_clicks"),
        State("dashboard-workspace-meta", "data"),
        prevent_initial_call=True,
    )

    app.layout = enable_workspace_persistence(app.layout)

    logging.getLogger("werkzeug").setLevel(logging.WARNING)

    if auto_port and not inherited_port and port != DEFAULT_PORT:
        print(f"Port {DEFAULT_PORT} is in use; starting Dash app on port {port} instead.")
    if _should_open_browser(args.debug):
        os.environ[BROWSER_OPENED_ENV] = "1"
        _open_browser(f"http://{_browser_host(args.host)}:{port}")

    print(f"Local dashboard runtime broker ready: {connection_path}")
    if args.threaded:
        print(
            "WARNING: --threaded bypasses the safe NetCDF/HDF5 process-isolation "
            "model and can cause intermittent crashes."
        )

    try:
        app.run(
            host=args.host,
            port=port,
            debug=args.debug,
            threaded=args.threaded,
        )
    finally:
        if dashboard_registration is not None:
            dashboard_registration.stop()


if __name__ == "__main__":
    main()
