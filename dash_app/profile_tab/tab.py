"""Thin entrypoint for the Profile tab."""

from __future__ import annotations

from dash import dcc

from .callbacks import register_profile_callbacks
from .discovery import discover_profile_state
from .layout import build_layout


def build_tab(app):
    initial_state = discover_profile_state()
    register_profile_callbacks(app)
    return dcc.Tab(id="dashboard-tab-profile", label="Profile", value="profile", children=build_layout(initial_state))
