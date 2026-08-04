"""Entrypoint for the dashboard Tutorial tab."""

from dash import dcc

from .adg1_gaussian_explorer import register_callbacks as register_adg1_callbacks
from .clubb_equations_demo.app import register_callbacks as register_equation_callbacks

from .callbacks import register_callbacks
from .layout import build_layout


def build_tab(app):
    register_callbacks(app)
    register_equation_callbacks(app)
    register_adg1_callbacks(app)
    return dcc.Tab(label="Tutorial", value="tutorial", children=build_layout())
