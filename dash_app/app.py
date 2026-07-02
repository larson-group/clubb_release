#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import logging
import os
import pkgutil
import socket
import sys
import threading
import webbrowser

from dash import Dash, dcc, html, Input, Output, State

DASH_APP_ROOT = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(DASH_APP_ROOT, ".."))
for import_root in (REPO_ROOT, DASH_APP_ROOT):
    if import_root not in sys.path:
        sys.path.insert(0, import_root)

from run_tab.tab import build_tab as build_run_tab
from compile_tab.tab import build_tab as build_compile_tab
from tune_tab.tab import build_tab as build_tune_tab
from plot_tab.tab import build_tab as build_plots_tab


if not hasattr(pkgutil, "find_loader"):
    # Dash 3.2.0 still calls pkgutil.find_loader in debug mode; Python 3.14 removed it
    def _find_loader(name):
        spec = importlib.util.find_spec(name)
        return None if spec is None else spec.loader

    pkgutil.find_loader = _find_loader


DEFAULT_PORT = 8050
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


def _open_browser(url: str) -> None:
    def _launch():
        try:
            webbrowser.open(url, new=2)
        except Exception:
            pass

    threading.Timer(1.0, _launch).start()


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
        help="Enable threaded request handling. Disabled by default because netCDF/HDF5 access is more stable without it.",
    )
    args = parser.parse_args()
    port, auto_port, inherited_port = _resolve_port(args.host, args.port)

    app = Dash(__name__, suppress_callback_exceptions=True, eager_loading=True, title=_app_title())

    tabs = [
        build_compile_tab(app),
        build_run_tab(app),
        build_tune_tab(app),
        build_plots_tab(app),
    ]

    app.layout = html.Div(
        [
            dcc.Store(id="theme-store", data="dark"),
            dcc.Tabs(tabs, value="compile"),
            html.Button(
                "Theme: Dark",
                id="theme-toggle-button",
                n_clicks=0,
                className="theme-toggle-button",
            ),
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
        Output("app-root", "className"),
        Output("theme-toggle-button", "children"),
        Input("theme-store", "data"),
    )
    def _apply_theme(theme_name):
        theme = (theme_name or "dark").lower()
        if theme == "light":
            return "theme-light", "Theme: Light"
        return "theme-dark", "Theme: Dark"

    logging.getLogger("werkzeug").setLevel(logging.WARNING)

    if auto_port and not inherited_port and port != DEFAULT_PORT:
        print(f"Port {DEFAULT_PORT} is in use; starting Dash app on port {port} instead.")
    if _should_open_browser(args.debug):
        os.environ[BROWSER_OPENED_ENV] = "1"
        _open_browser(f"http://{_browser_host(args.host)}:{port}")

    app.run(
        host=args.host,
        port=port,
        debug=args.debug,
        threaded=args.threaded,
    )


if __name__ == "__main__":
    main()
