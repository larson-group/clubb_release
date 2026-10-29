"""Focused lifecycle checks for durable Tune workspace lineages."""

from __future__ import annotations

import os

from tuner import workspaces


def _request():
    return {
        "cases": ["bomex"],
        "case_configs": [
            {
                "case_name": "bomex",
                "time_average_range": [0, 60],
                "average_time_seconds": 60,
                "altitude_comparison_range": [0.0, 100.0],
            }
        ],
        "selected_fields": ["cloud_frac"],
        "parameter_ranges": [{"name": "C_uu_shr", "min": 0.2, "max": 0.22}],
        "strategy": {"name": "random", "options": {"max_samples": 1}},
    }


def test_workspace_lineage_uses_original_revisions_and_restart_names(tmp_path, monkeypatch):
    """A workspace preserves immutable execution directory names and metadata."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    workspace_id, _original = workspaces.create_workspace(_request())
    revision = workspaces.create_revision(workspace_id, "original")
    restart = workspaces.create_revision(workspace_id, "original", restart=True)

    workspaces.rename_workspace(workspace_id, "PDF smoke")
    summary = workspaces.list_workspaces()

    assert summary[0]["display_name"] == "PDF smoke"
    assert {item["revision_id"] for item in summary[0]["revisions"]} == {
        "original",
        revision.job_dir.name,
        restart.job_dir.name,
    }
    assert restart.job_dir.name == "original_restart1"


def test_workspace_rejects_delete_while_a_revision_is_active(tmp_path, monkeypatch):
    """Deletion cannot remove an active job's durable control and results files."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    workspace_id, original = workspaces.create_workspace(_request())
    workspaces.atomic_write_json(original.status_path, {"state": "running"})

    try:
        workspaces.delete_workspace(workspace_id)
    except RuntimeError as exc:
        assert "active revision" in str(exc)
    else:
        raise AssertionError("active workspace deletion unexpectedly succeeded")


def test_reset_execution_erases_results_and_reopens_the_same_draft(tmp_path, monkeypatch):
    """Reset preserves the request but makes one inactive revision editable again."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    workspace_id, original = workspaces.create_workspace(_request())
    workspaces.atomic_write_json(original.status_path, {"state": "finished", "samples_evaluated": 1})
    original.results_path.write_text('{"state": "finished", "best_results": [{"loss": 1.0}]}', encoding="utf-8")
    original.log_path.write_text("old worker output\n", encoding="utf-8")

    reset = workspaces.reset_execution(workspace_id, "original")
    loaded = workspaces.load_execution(workspace_id, "original")

    assert reset.job_dir == original.job_dir
    assert loaded["request"] == _request()
    assert loaded["execution"]["state"] == "draft"
    assert loaded["results"]["best_results"] == []
    assert not reset.log_path.exists()


def test_workspace_listing_uses_actual_latest_data_modification(tmp_path, monkeypatch):
    """The browser order follows worker/result activity, not stale metadata."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    older_id, older = workspaces.create_workspace(_request(), display_name="older")
    newer_id, newer = workspaces.create_workspace(_request(), display_name="newer")
    os.utime(older.results_path, (1_900_000_000.0, 1_900_000_000.0))
    os.utime(newer.results_path, (2_000_000_000.0, 2_000_000_000.0))

    listed = workspaces.list_workspaces()

    assert [item["workspace_id"] for item in listed] == [newer_id, older_id]
    assert listed[0]["modified_at"].startswith("2033-05-18T03:33:20")


def test_workspace_display_names_default_to_unique_timecodes_and_reject_duplicates(tmp_path, monkeypatch):
    """The browser-facing workspace labels remain unambiguous."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    first_id, _first = workspaces.create_workspace(_request())
    second_id, _second = workspaces.create_workspace(_request())
    labels = {item["workspace_id"]: item["display_name"] for item in workspaces.list_workspaces()}

    assert labels[first_id].startswith("unnamed-")
    assert labels[second_id].startswith("unnamed-")
    assert labels[first_id] != labels[second_id]

    workspaces.rename_workspace(first_id, "PDF comparison")
    try:
        workspaces.rename_workspace(second_id, "pdf COMPARISON")
    except ValueError as exc:
        assert "already exists" in str(exc)
    else:
        raise AssertionError("duplicate workspace display name unexpectedly accepted")


def test_workspace_activity_reads_states_without_rich_browser_metadata(tmp_path, monkeypatch):
    """The activity indicator needs only direct revision-state files."""
    monkeypatch.setattr(workspaces, "WORKSPACE_ROOT", tmp_path / "tuner")
    workspace_id, original = workspaces.create_workspace(_request())
    workspaces.atomic_write_json(original.status_path, {"state": "running"})

    def unexpected_recursive_scan(_path):
        raise AssertionError("activity polling must not scan revision artifacts")

    monkeypatch.setattr(workspaces, "_latest_modified_epoch", unexpected_recursive_scan)

    assert workspaces.list_workspace_activity() == [
        {"workspace_id": workspace_id, "revision_id": "original", "state": "running"}
    ]
