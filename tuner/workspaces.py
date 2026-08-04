"""Durable Tune workspace and revision lineage helpers.

Each workspace is a human-named tuning lineage under ``output/tuner``.  An
execution directory is immutable once started: ``original`` is the first
request, ``revN`` is an editable-clone revision, and ``<revision>_restartN``
is a fresh retry of that exact request.
"""

from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path
import re
import shutil
from typing import Any

from tuner.job_runtime import TunerJob
from tuner.paths import REPO_ROOT
from tuner.status import atomic_write_json, read_json_or_default


WORKSPACE_ROOT = REPO_ROOT / "output" / "tuner"
WORKSPACE_METADATA = "workspace.json"
WORKSPACE_SCHEMA_VERSION = 1
_SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _safe_id(value: str, *, label: str) -> str:
    value = str(value or "").strip()
    if not _SAFE_ID_RE.fullmatch(value):
        raise ValueError(f"{label} must contain only letters, digits, '_' or '-'")
    return value


def _workspace_path(workspace_id: str) -> Path:
    return WORKSPACE_ROOT / _safe_id(workspace_id, label="workspace id")


def _execution_path(workspace_id: str, revision_id: str) -> Path:
    return _workspace_path(workspace_id) / _safe_id(revision_id, label="revision id")


def _workspace_metadata(path: Path) -> dict[str, Any]:
    return read_json_or_default(path / WORKSPACE_METADATA, {})


def _write_workspace_metadata(path: Path, payload: dict[str, Any]) -> None:
    payload = dict(payload)
    payload["updated_at"] = _now()
    atomic_write_json(path / WORKSPACE_METADATA, payload)


def _directory_size(path: Path) -> int:
    total = 0
    try:
        for child in path.rglob("*"):
            if child.is_file():
                try:
                    total += child.stat().st_size
                except OSError:
                    continue
    except OSError:
        pass
    return total


def _latest_modified_epoch(path: Path) -> float:
    """Return the newest local file timestamp below a workspace/execution."""
    try:
        newest = path.stat().st_mtime
    except OSError:
        return 0.0
    try:
        for child in path.rglob("*"):
            try:
                newest = max(newest, child.stat().st_mtime)
            except OSError:
                continue
    except OSError:
        pass
    return newest


def _timestamp_from_epoch(value: float) -> str:
    return datetime.fromtimestamp(value, timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _status_for(path: Path) -> dict[str, Any]:
    return read_json_or_default(path / "status.json", {"state": "draft"})


def _execution_summary(path: Path) -> dict[str, Any]:
    status = _status_for(path)
    return {
        "revision_id": path.name,
        "path": str(path),
        "state": str(status.get("state") or "draft"),
        "samples_evaluated": int(status.get("samples_evaluated") or 0),
        "total_samples": status.get("total_samples"),
        "created_at": _timestamp_from_epoch(path.stat().st_ctime),
        "modified_at": _timestamp_from_epoch(_latest_modified_epoch(path)),
        "size_bytes": _directory_size(path),
        "request_path": str(path / "request.json"),
        "status_path": str(path / "status.json"),
        "results_path": str(path / "results.json"),
    }


def _revision_paths(workspace_path: Path, *, order_by_activity: bool = True) -> list[Path]:
    """Return execution directories for one workspace.

    The rich saved-workspace browser orders revisions by their newest artifact,
    which requires a recursive scan.  Lightweight status polling must not pay
    that cost, so it can request the stable directory-name order instead.
    """
    if not workspace_path.is_dir():
        return []
    paths = [path for path in workspace_path.iterdir() if path.is_dir() and (path / "request.json").is_file()]
    if order_by_activity:
        return sorted(paths, key=_latest_modified_epoch)
    return sorted(paths, key=lambda path: path.name)


def _next_revision_id(workspace_path: Path) -> str:
    highest = 0
    for path in _revision_paths(workspace_path):
        match = re.fullmatch(r"rev(\d+)", path.name)
        if match:
            highest = max(highest, int(match.group(1)))
    return f"rev{highest + 1}"


def _next_restart_id(workspace_path: Path, base_revision_id: str) -> str:
    highest = 0
    pattern = re.compile(re.escape(base_revision_id) + r"_restart(\d+)$")
    for path in _revision_paths(workspace_path):
        match = pattern.fullmatch(path.name)
        if match:
            highest = max(highest, int(match.group(1)))
    return f"{base_revision_id}_restart{highest + 1}"


def _display_name_for_request(request: dict[str, Any]) -> str:
    """Return a neutral, collision-resistant default workspace label."""
    # A generated label is deliberately independent of the selected cases.
    # It is visible in the Tune browser before a user chooses a descriptive
    # name, and remains unique even for several otherwise identical launches.
    return suggested_workspace_display_name()


def suggested_workspace_display_name() -> str:
    """Return the next available ``unnamed-<timecode>`` display label."""
    WORKSPACE_ROOT.mkdir(parents=True, exist_ok=True)
    base = f"unnamed-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    existing = {
        str(_workspace_metadata(path).get("display_name") or "").strip().casefold()
        for path in WORKSPACE_ROOT.iterdir()
        if path.is_dir() and (path / WORKSPACE_METADATA).is_file()
    }
    if base.casefold() not in existing:
        return base
    suffix = 2
    while f"{base}-{suffix}".casefold() in existing:
        suffix += 1
    return f"{base}-{suffix}"


def _validate_unique_display_name(display_name: str, *, excluding_workspace_id: str | None = None) -> str:
    """Validate a human-visible workspace name against all saved workspaces."""
    label = str(display_name or "").strip()
    if not label:
        raise ValueError("workspace name cannot be blank")
    if len(label) > 120:
        raise ValueError("workspace name must be at most 120 characters")
    wanted = label.casefold()
    if WORKSPACE_ROOT.is_dir():
        for path in WORKSPACE_ROOT.iterdir():
            if not path.is_dir() or not (path / WORKSPACE_METADATA).is_file():
                continue
            if path.name == excluding_workspace_id:
                continue
            existing = str(_workspace_metadata(path).get("display_name") or "").strip()
            if existing.casefold() == wanted:
                raise ValueError(f"workspace name already exists: {label}")
    return label


def _new_workspace_id(request: dict[str, Any]) -> str:
    case_text = "-".join(str(case).strip() for case in request.get("cases", []) if str(case).strip()) or "tune"
    case_text = re.sub(r"[^A-Za-z0-9_-]+", "-", case_text).strip("-") or "tune"
    return f"{case_text}-{datetime.now().strftime('%Y%m%d-%H%M%S-%f')}"


def create_workspace(request: dict[str, Any], *, display_name: str | None = None) -> tuple[str, TunerJob]:
    """Create a workspace and its unstarted ``original`` execution."""
    WORKSPACE_ROOT.mkdir(parents=True, exist_ok=True)
    workspace_id = _new_workspace_id(request)
    path = _workspace_path(workspace_id)
    path.mkdir(parents=True, exist_ok=False)
    created_at = _now()
    label = str(display_name or "").strip() or _display_name_for_request(request)
    label = _validate_unique_display_name(label)
    _write_workspace_metadata(
        path,
        {
            "schema_version": WORKSPACE_SCHEMA_VERSION,
            "workspace_id": workspace_id,
            "display_name": label,
            "created_at": created_at,
            "updated_at": created_at,
        },
    )
    job = TunerJob.create(request, job_dir=path / "original", initial_state="draft")
    return workspace_id, job


def create_revision(workspace_id: str, source_revision_id: str, *, restart: bool = False) -> TunerJob:
    """Clone one saved request into a new revision or an exact restart attempt."""
    workspace_path = _workspace_path(workspace_id)
    source_path = _execution_path(workspace_id, source_revision_id)
    if not (source_path / "request.json").is_file():
        raise ValueError(f"unknown saved revision: {source_revision_id}")
    source_status = _status_for(source_path)
    if str(source_status.get("state") or "") in {"running", "initializing", "stopping"}:
        raise RuntimeError("stop the selected revision before creating a revision or restart")
    request = read_json_or_default(source_path / "request.json", {})
    if not request:
        raise RuntimeError("saved revision has no readable request")
    revision_id = _next_restart_id(workspace_path, source_revision_id) if restart else _next_revision_id(workspace_path)
    job = TunerJob.create(request, job_dir=workspace_path / revision_id, initial_state="draft")
    metadata = _workspace_metadata(workspace_path)
    _write_workspace_metadata(workspace_path, metadata)
    return job


def replace_draft_request(workspace_id: str, revision_id: str, request: dict[str, Any]) -> TunerJob:
    """Replace an unstarted draft request before its first execution."""
    path = _execution_path(workspace_id, revision_id)
    status = _status_for(path)
    if str(status.get("state") or "draft") not in {"draft", "idle"}:
        raise RuntimeError("only an unstarted revision can be edited")
    job = TunerJob.from_dir(path)
    atomic_write_json(job.request_path, dict(request))
    metadata = _workspace_metadata(_workspace_path(workspace_id))
    _write_workspace_metadata(_workspace_path(workspace_id), metadata)
    return job


def reset_execution(workspace_id: str, revision_id: str) -> TunerJob:
    """Erase one inactive execution's results and restore its editable draft.

    The workspace identity and the exact prior request are retained, while all
    generated samples, logs, status, and result data are deliberately removed.
    """
    workspace_path = _workspace_path(workspace_id)
    path = _execution_path(workspace_id, revision_id)
    status = _status_for(path)
    if str(status.get("state") or "") in {"initializing", "running", "stopping"}:
        raise RuntimeError("stop the selected revision before resetting its data")
    request = read_json_or_default(path / "request.json", {})
    if not request:
        raise RuntimeError("saved revision has no readable request")
    shutil.rmtree(path)
    job = TunerJob.create(request, job_dir=path, initial_state="draft")
    _write_workspace_metadata(workspace_path, _workspace_metadata(workspace_path))
    return job


def list_workspaces() -> list[dict[str, Any]]:
    """Return lightweight metadata for each new-style saved Tune workspace."""
    if not WORKSPACE_ROOT.is_dir():
        return []
    entries = []
    for path in WORKSPACE_ROOT.iterdir():
        if not path.is_dir() or not (path / WORKSPACE_METADATA).is_file():
            continue
        metadata = _workspace_metadata(path)
        revisions = [_execution_summary(revision) for revision in _revision_paths(path)]
        modified_epoch = _latest_modified_epoch(path)
        entries.append(
            {
                "workspace_id": path.name,
                "display_name": str(metadata.get("display_name") or path.name),
                "created_at": metadata.get("created_at"),
                "modified_at": _timestamp_from_epoch(modified_epoch),
                "_modified_epoch": modified_epoch,
                "size_bytes": _directory_size(path),
                "revision_count": len(revisions),
                "latest_revision_id": revisions[-1]["revision_id"] if revisions else None,
                "latest_state": revisions[-1]["state"] if revisions else "draft",
                "revisions": revisions,
            }
        )
    entries.sort(key=lambda entry: float(entry.get("_modified_epoch") or 0.0), reverse=True)
    for entry in entries:
        entry.pop("_modified_epoch", None)
    return entries


def list_workspace_activity() -> list[dict[str, str]]:
    """Return only revision states needed for lightweight dashboard activity UI.

    This intentionally avoids recursive size and modification-time scans. A
    Tune workspace can contain thousands of worker artifacts, and polling the
    rich workspace browser merely to animate an activity indicator makes the
    dashboard feel sluggish and can rebuild controls while the user is
    interacting with them.
    """
    if not WORKSPACE_ROOT.is_dir():
        return []
    activity: list[dict[str, str]] = []
    for workspace_path in WORKSPACE_ROOT.iterdir():
        if not workspace_path.is_dir() or not (workspace_path / WORKSPACE_METADATA).is_file():
            continue
        for revision_path in _revision_paths(workspace_path, order_by_activity=False):
            activity.append(
                {
                    "workspace_id": workspace_path.name,
                    "revision_id": revision_path.name,
                    "state": str(_status_for(revision_path).get("state") or "draft"),
                }
            )
    return activity


def load_execution(workspace_id: str, revision_id: str) -> dict[str, Any]:
    """Read one immutable execution's request, live status, and retained results."""
    path = _execution_path(workspace_id, revision_id)
    if not (path / "request.json").is_file():
        raise ValueError(f"unknown saved revision: {revision_id}")
    return {
        "workspace_id": workspace_id,
        "revision_id": revision_id,
        "workspace": _workspace_metadata(_workspace_path(workspace_id)),
        "execution": _execution_summary(path),
        "request": read_json_or_default(path / "request.json", {}),
        "status": _status_for(path),
        "results": read_json_or_default(path / "results.json", {}),
        "job": TunerJob.from_dir(path).to_dict(),
    }


def workspace_display_name(workspace_id: str) -> str:
    """Return the current human label for a durable Tune workspace."""
    path = _workspace_path(workspace_id)
    metadata = _workspace_metadata(path)
    return str(metadata.get("display_name") or workspace_id)


def rename_workspace(workspace_id: str, display_name: str) -> dict[str, Any]:
    """Rename a workspace display label without moving immutable execution paths."""
    path = _workspace_path(workspace_id)
    if not path.is_dir():
        raise ValueError(f"unknown workspace: {workspace_id}")
    label = _validate_unique_display_name(display_name, excluding_workspace_id=workspace_id)
    metadata = _workspace_metadata(path)
    metadata["display_name"] = label
    _write_workspace_metadata(path, metadata)
    return next(entry for entry in list_workspaces() if entry["workspace_id"] == workspace_id)


def delete_workspace(workspace_id: str) -> None:
    """Delete an inactive workspace lineage after the caller has confirmed it."""
    path = _workspace_path(workspace_id)
    if not path.is_dir():
        raise ValueError(f"unknown workspace: {workspace_id}")
    active = [entry["revision_id"] for entry in (_execution_summary(item) for item in _revision_paths(path)) if entry["state"] in {"running", "initializing", "stopping"}]
    if active:
        raise RuntimeError("cannot delete a workspace with active revision(s): " + ", ".join(active))
    shutil.rmtree(path)
