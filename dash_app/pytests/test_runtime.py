"""Focused tests for private dashboard runtime-directory discovery."""

from pathlib import Path

from dash_app.shared import runtime


def test_runtime_candidates_find_standard_uid_path_without_xdg(monkeypatch):
    """A separately launched agent must find a Dash broker's normal runtime root."""

    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)

    assert runtime._runtime_base_candidates("1000") == [
        Path("/run/user/1000"),
        Path("/tmp"),
    ]


def test_runtime_candidates_keep_explicit_xdg_before_standard_path(monkeypatch):
    monkeypatch.setenv("XDG_RUNTIME_DIR", "/custom/runtime")

    assert runtime._runtime_base_candidates("1000") == [
        Path("/custom/runtime"),
        Path("/run/user/1000"),
        Path("/tmp"),
    ]
