#!/usr/bin/env python3
import argparse

from dash import Dash, dcc, html, Input, Output, State

from tab_run import build_tab as build_run_tab
from tab_compare import build_tab as build_compare_tab
from tab_multicol import build_tab as build_multicol_tab


def main():
    parser = argparse.ArgumentParser(
        description="Dash app with multiple tabs for CLUBB NetCDF analysis."
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument(
        "-debug",
        action="store_true",
        help="Run Dash in debug mode (disables threaded mode).",
    )
    args = parser.parse_args()

    app = Dash(__name__)

    tabs = [
        build_run_tab(app),
        build_compare_tab(app),
        build_multicol_tab(app),
    ]

    app.layout = html.Div(
        [
            dcc.Store(id="theme-store", data="dark"),
            dcc.Tabs(tabs),
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

    app.run(
        host=args.host,
        port=args.port,
        debug=args.debug,
        threaded=not args.debug,
    )


if __name__ == "__main__":
    main()
