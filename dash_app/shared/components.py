"""Small reusable Dash components with application-wide visual fixes."""

from __future__ import annotations

from typing import Any

from dash import dcc


def styled_dropdown(*, className: str = "", **kwargs: Any):
    """Build a Dash dropdown with the shared focus/overflow repair class.

    Keep this intentionally small: tabs retain control over options, ids, and
    callbacks while recurring React-Select rendering fixes live in one place.
    """

    classes = " ".join(part for part in ("clubb-dropdown", className) if part)
    return dcc.Dropdown(className=classes, **kwargs)
