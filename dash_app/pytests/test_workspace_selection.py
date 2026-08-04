"""Regression checks for Tune workspace selection ownership."""

from dash_app.tune_tab.callbacks_workspaces import (
    _is_positive_click_count,
    _sealed_selection_for_job,
    _workspace_button_class,
)


def test_old_active_job_cannot_reselect_workspace_browser():
    """Opening Workspaces must not be overwritten by a previous live job."""
    active_job = {"workspace_id": "old", "revision_id": "original"}

    assert _sealed_selection_for_job(active_job, {}) is None
    assert _sealed_selection_for_job(active_job, {"mode": "new", "draft_token": "new-draft"}) is None


def test_new_draft_token_allows_only_its_own_start_to_be_sealed():
    active_job = {"workspace_id": "new", "revision_id": "original", "draft_token": "token-1"}

    assert _sealed_selection_for_job(active_job, {"mode": "new", "draft_token": "token-1"}) == {
        "workspace_id": "new",
        "revision_id": "original",
        "mode": "readonly",
    }
    assert _sealed_selection_for_job(active_job, {"mode": "new", "draft_token": "token-2"}) is None


def test_workspace_activity_stripes_only_for_a_different_revision():
    activity = [
        {"workspace_id": "selected", "revision_id": "original", "state": "running"},
        {"workspace_id": "other", "revision_id": "rev2", "state": "stopped"},
    ]
    selection = {"workspace_id": "selected", "revision_id": "original"}

    assert _workspace_button_class(selection, activity) == "tune-workspace-change"
    activity[1]["state"] = "running"
    assert _workspace_button_class(selection, activity).endswith("tune-workspace-change-active")


def test_workspace_actions_ignore_rehydrated_zero_click_counts():
    """Dynamic workspace rows must not navigate when Dash rebuilds them."""
    assert not _is_positive_click_count(None)
    assert not _is_positive_click_count(0)
    assert not _is_positive_click_count("0")
    assert _is_positive_click_count(1)
