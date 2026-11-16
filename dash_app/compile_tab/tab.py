"""Thin entrypoint for the compile tab."""

from __future__ import annotations

from dash import dcc

from .callbacks import register_compile_callbacks
from .layout import build_layout


def build_tab(app, *, lazy=None):
    """Build the Compile tab and register callbacks."""
    register_compile_callbacks(app)

    def layout():
        return build_layout()

    if lazy is not None:
        return lazy.tab(id="dashboard-tab-compile", label="Compile", value="compile", build=layout)
    return dcc.Tab(id="dashboard-tab-compile", label="Compile", value="compile", children=layout())
