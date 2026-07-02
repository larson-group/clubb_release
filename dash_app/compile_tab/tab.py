"""Thin entrypoint for the compile tab."""

from __future__ import annotations

from dash import dcc

from .callbacks import register_compile_callbacks
from .discovery import discover_compile_state
from .layout import build_layout


def build_tab(app):
    """Build the Compile tab and register callbacks."""
    initial_state = discover_compile_state()
    register_compile_callbacks(app)
    return dcc.Tab(label="Compile", value="compile", children=build_layout(initial_state))
