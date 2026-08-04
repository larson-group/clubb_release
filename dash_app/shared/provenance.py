"""Compact, deterministic provenance records for dashboard-created artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from pathlib import Path
from typing import Any


RUNTIME_FINGERPRINT_SCOPE = "dash_app,tuner,utilities:non-test-python"
_RUNTIME_ROOTS = ("dash_app", "tuner", "utilities")
_RUNTIME_EXCLUDED_PARTS = {"__pycache__", "pytests", "tests"}


def sha256_file(path: Path) -> str | None:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def canonical_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def runtime_source_fingerprint(repo_root: Path | str) -> str:
    """Hash the Python runtime source used by the broker and Dash actions.

    This intentionally hashes file contents rather than only Git metadata, so
    tracked edits, uncommitted edits, and untracked runtime modules all cause
    a new fingerprint. Test files, caches, and generated output are excluded
    to avoid restarting the durable broker for test-only changes.
    """
    root = Path(repo_root).resolve()
    digest = hashlib.sha256(RUNTIME_FINGERPRINT_SCOPE.encode("utf-8"))
    files: list[Path] = []
    for relative_root in _RUNTIME_ROOTS:
        source_root = root / relative_root
        if not source_root.is_dir():
            continue
        for path in source_root.rglob("*.py"):
            if path.is_file() and not any(part in _RUNTIME_EXCLUDED_PARTS for part in path.relative_to(root).parts):
                files.append(path)
    for path in sorted(files, key=lambda item: item.relative_to(root).as_posix()):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        try:
            content = path.read_bytes()
        except OSError:
            content = b"<unreadable>"
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


def _git(repo_root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args], cwd=repo_root, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=3, check=False
        )
    except OSError:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def source_provenance(repo_root: Path | str) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    status = _git(root, "status", "--porcelain")
    return {
        "repository": str(root),
        "commit": _git(root, "rev-parse", "HEAD"),
        "dirty": bool(status),
        "dirty_hash": hashlib.sha256((status or "").encode("utf-8")).hexdigest(),
        "captured_at_unix_seconds": time.time(),
    }
