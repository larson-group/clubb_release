"""Entrypoint for the dashboard Tutorial tab."""

from dash import dcc

from .adg1_gaussian_explorer import register_callbacks as register_adg1_callbacks
from .clubb_equations_demo.app import register_callbacks as register_equation_callbacks

from .callbacks import register_callbacks
from .layout import build_layout
from dash_app.lazy_tabs import LazyTabs


def build_tab(app, *, defer=False):
    register_callbacks(app)
    register_equation_callbacks(app)
    register_adg1_callbacks(app)
    lazy = LazyTabs("tutorial-pages", parent=("dashboard-tabs", "tutorial")) if defer else None
    layout = build_layout(lazy=lazy)
    children = [lazy.register(app), layout] if lazy is not None else layout
    return dcc.Tab(label="Tutorial", value="tutorial", children=children)
