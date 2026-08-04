"""Vertical page layout for the Tutorial tab."""

from dash import dcc, html

from .adg1_gaussian_explorer import build_layout as build_adg1_explorer_layout
from .clubb_equations import build_page as build_equations_page
from .welcome import build_page as build_welcome_page


def build_layout():
    """Build the tutorial page rail and its initial lessons."""
    return html.Div(
        dcc.Tabs(
            id="tutorial-pages",
            value="tutorial-welcome",
            vertical=True,
            persistence=True,
            persistence_type="local",
            parent_className="tutorial-tabs-shell",
            className="tutorial-page-rail",
            content_className="tutorial-page-content",
            children=[
                dcc.Tab(
                    label="Welcome",
                    value="tutorial-welcome",
                    className="tutorial-page-tab",
                    selected_className="tutorial-page-tab-selected",
                    children=build_welcome_page(),
                ),
                dcc.Tab(
                    label="CLUBB Equations",
                    value="tutorial-equations",
                    className="tutorial-page-tab",
                    selected_className="tutorial-page-tab-selected",
                    children=build_equations_page(),
                ),
                dcc.Tab(
                    label="ADG1 two-Gaussian explorer",
                    value="tutorial-adg1-explorer",
                    className="tutorial-page-tab",
                    selected_className="tutorial-page-tab-selected",
                    children=build_adg1_explorer_layout(),
                ),
            ],
        ),
        className="tutorial-tab-root",
    )
