"""Thin entrypoint for the Profile tab."""

from __future__ import annotations

from dash import dcc

from .callbacks import register_profile_callbacks
from .discovery import discover_profile_state
from .layout import build_layout


def build_tab(app, *, lazy=None):
    register_profile_callbacks(app)

    def layout():
        return build_layout(discover_profile_state())

    if lazy is not None:
        return lazy.tab(id="dashboard-tab-profile", label="Profile", value="profile", build=layout)
    return dcc.Tab(id="dashboard-tab-profile", label="Profile", value="profile", children=layout())
