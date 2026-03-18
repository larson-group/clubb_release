#!/usr/bin/env python3
import argparse
import importlib.util
import logging
import os
import pkgutil
import socket
import threading
import webbrowser

from dash import Dash, dcc, html, Input, Output, State

from run_tab.tab import build_tab as build_run_tab
from plot_tab.tab import build_tab as build_plots_tab


if not hasattr(pkgutil, "find_loader"):
    # Dash 3.2.0 still calls pkgutil.find_loader in debug mode; Python 3.14 removed it.
    def _find_loader(name):
        spec = importlib.util.find_spec(name)
        return None if spec is None else spec.loader

    pkgutil.find_loader = _find_loader


DEFAULT_PORT = 8050


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


def _should_open_browser(debug: bool) -> bool:
    if not debug:
        return True
    return os.environ.get("WERKZEUG_RUN_MAIN") == "true"


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
    auto_port = args.port is None
    port = args.port if args.port is not None else _find_available_port(args.host, DEFAULT_PORT)

    app = Dash(__name__, suppress_callback_exceptions=True)

    tabs = [
        build_run_tab(app),
        build_plots_tab(app),
    ]

    app.layout = html.Div(
        [
            dcc.Store(id="theme-store", data="dark"),
            dcc.Tabs(tabs, value="run"),
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

    if auto_port and port != DEFAULT_PORT:
        print(f"Port {DEFAULT_PORT} is in use; starting Dash app on port {port} instead.")
    if _should_open_browser(args.debug):
        _open_browser(f"http://{_browser_host(args.host)}:{port}")

    app.run(
        host=args.host,
        port=port,
        debug=args.debug,
        threaded=args.threaded,
    )


if __name__ == "__main__":
    main()
