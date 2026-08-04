"""Durable idempotent job and artifact records independent of Dash layout."""

from __future__ import annotations

import json
import os
import re
import secrets
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dash_app.shared.provenance import canonical_hash, source_provenance
from dash_app.shared.runtime import atomic_write_json, private_path


REPO_ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_ROOT = REPO_ROOT / "output" / "agent_artifacts"
# MCP bundles are staging material: they make a running job inspectable and
# let an agent collect compact results, but they are not a durable evidence
# store.  Twenty completed bundles is ample for a live investigation while
# keeping this ignored directory safe to clear between investigations.
DEFAULT_ARTIFACT_RETENTION_COUNT = 20
SERVICE_PROVENANCE = {
    "name": "CLUBB Dash local MCP service",
    "version": 3,
    "transport_scope": "local-only",
}


class SubmissionConflict(ValueError):
    """A request id was reused with a different immutable request payload."""


class JobStore:
    """Small JSON-backed registry for broker-owned job identity and retries."""

    def __init__(self, path: Path | None = None):
        self.path = path or private_path(REPO_ROOT, "jobs.json")
        self.lock_path = self.path.with_suffix(".lock")

    @contextmanager
    def _locked(self):
        """Serialize retries from independent MCP stdio processes.

        The dashboard and each MCP client are separate processes.  Atomic
        replacement alone prevents partial JSON, but not two simultaneous
        readers both deciding that a request id is new.
        """
        import fcntl

        self.lock_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(self.lock_path, os.O_RDWR | os.O_CREAT, 0o600)
        try:
            os.fchmod(descriptor, 0o600)
            fcntl.flock(descriptor, fcntl.LOCK_EX)
            yield
        finally:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
            os.close(descriptor)

    def _read(self) -> dict[str, Any]:
        try:
            payload = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            payload = {}
        return payload if isinstance(payload, dict) else {}

    def _write(self, payload: dict[str, Any]) -> None:
        atomic_write_json(self.path, payload)

    def submit(self, kind: str, request_id: str, request: dict[str, Any]) -> tuple[dict[str, Any], bool]:
        with self._locked():
            state = self._read()
            key = f"{kind}:{request_id}"
            request_hash = canonical_hash(request)
            existing = state.get(key)
            if isinstance(existing, dict):
                if existing.get("request_hash") != request_hash:
                    raise SubmissionConflict("request_id is already associated with a different request")
                return dict(existing), False
            record = {
                "job_id": f"{kind}_{secrets.token_urlsafe(12)}",
                "run_id": f"run_{secrets.token_urlsafe(12)}" if kind == "scm" else None,
                "batch_id": f"batch_{secrets.token_urlsafe(12)}" if kind == "scm_batch" else None,
                "kind": kind,
                "request_id": request_id,
                "request": request,
                "request_hash": request_hash,
                "state": "submitting",
                "created_at_unix_seconds": time.time(),
                "provenance": {
                    "source": source_provenance(REPO_ROOT),
                    "service": dict(SERVICE_PROVENANCE),
                },
            }
            state[key] = record
            self._write(state)
            return dict(record), True

    def update(self, job_id: str, **updates: Any) -> dict[str, Any] | None:
        with self._locked():
            state = self._read()
            for key, record in state.items():
                if isinstance(record, dict) and record.get("job_id") == job_id:
                    record.update(updates)
                    record["updated_at_unix_seconds"] = time.time()
                    state[key] = record
                    self._write(state)
                    return dict(record)
        return None

    def get(self, job_id: str) -> dict[str, Any] | None:
        for record in self._read().values():
            if isinstance(record, dict) and record.get("job_id") == job_id:
                return dict(record)
        return None

    def get_submission(self, kind: str, request_id: str) -> dict[str, Any] | None:
        """Read one idempotent submission by its typed kind/request key."""
        record = self._read().get(f"{kind}:{request_id}")
        return dict(record) if isinstance(record, dict) else None

    def list_kind(self, kind: str) -> list[dict[str, Any]]:
        """Return compact copies of durable records of one job kind."""
        return [dict(record) for record in self._read().values() if isinstance(record, dict) and record.get("kind") == kind]

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        for record in self._read().values():
            if isinstance(record, dict) and record.get("run_id") == run_id:
                return dict(record)
        return None


class ArtifactStore:
    """Own disposable MCP staging bundles without making NetCDF copies.

    ``output/agent_artifacts`` is intentionally not a report or experiment
    archive.  A bundle may contain an isolated SCM output while its broker job
    is alive, then compact manifests, logs, and PNGs after completion.  Any
    result worth retaining must be copied into a durable output directory or
    a ``doc/reports/<report-id>/`` bundle.  The root can therefore be removed
    between jobs; it is recreated lazily on the next request.
    """

    _ACTIVE_MARKER = ".active"

    def __init__(self, root: Path | None = None):
        self.root = root or ARTIFACT_ROOT

    @staticmethod
    def _safe_id(artifact_id: str) -> str:
        value = str(artifact_id or "")
        if not re.fullmatch(r"[A-Za-z0-9_-]+", value):
            raise ValueError("artifact id must be a simple bundle name")
        return value

    def bundle_dir(self, artifact_id: str) -> Path:
        """Return one owner-private service bundle below the fixed artifact root."""
        safe_id = self._safe_id(artifact_id)
        self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
        try:
            self.root.chmod(0o700)
        except OSError:
            pass
        directory = self.root / safe_id
        directory.mkdir(mode=0o700, parents=False, exist_ok=True)
        try:
            directory.chmod(0o700)
        except OSError:
            pass
        return directory

    @staticmethod
    def retention_count() -> int:
        raw = os.environ.get("CLUBB_AGENT_ARTIFACT_RETENTION_COUNT", str(DEFAULT_ARTIFACT_RETENTION_COUNT))
        try:
            return max(1, min(int(raw), 1000))
        except (TypeError, ValueError):
            return DEFAULT_ARTIFACT_RETENTION_COUNT

    def activate(self, artifact_id: str) -> Path:
        """Protect a live bundle from automatic completed-bundle cleanup."""
        marker = self.bundle_dir(artifact_id) / self._ACTIVE_MARKER
        atomic_write_json(marker, {"artifact_id": self._safe_id(artifact_id), "started_at_unix_seconds": time.time()})
        return marker

    def release(self, artifact_id: str) -> None:
        """Mark a bundle as completed and eligible for bounded retention."""
        try:
            (self.root / self._safe_id(artifact_id) / self._ACTIVE_MARKER).unlink()
        except OSError:
            pass

    def prune(self, *, keep: set[str] | None = None) -> list[str]:
        """Remove oldest inactive artifact bundles beyond the configured cap.

        Only direct service-owned children with manifests are considered.  A
        live bundle has an explicit marker, so a later submission cannot prune
        output still being read or written by an active agent job.  Unmarked
        legacy directories are deliberately ignored rather than guessed at.
        """
        keep = set(keep or set())
        try:
            bundles = [path for path in self.root.iterdir() if path.is_dir() and (path / "manifest.json").is_file()]
        except OSError:
            return []
        limit = self.retention_count()
        candidates = sorted(
            (
                path
                for path in bundles
                if path.name not in keep and not (path / self._ACTIVE_MARKER).exists()
            ),
            # Directory mtimes change when the active marker is released.
            # Prefer the immutable manifest creation time so releasing a
            # bundle cannot make an older artifact appear newer than a later
            # completed bundle.  Legacy manifests fall back to directory mtime.
            key=self._retention_timestamp,
            reverse=True,
        )
        remaining_limit = max(0, limit - len(keep))
        removed: list[str] = []
        for path in candidates[remaining_limit:]:
            # Artifact bundles are service-owned and ignored by Git.  Delete
            # only validated direct children; never follow a supplied path.
            import shutil

            shutil.rmtree(path)
            removed.append(path.name)
        return removed

    @staticmethod
    def _retention_timestamp(path: Path) -> float:
        """Return stable bundle age metadata for retention ordering."""
        try:
            payload = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
            created_at = float(payload.get("created_at_unix_seconds")) if isinstance(payload, dict) else float("nan")
            if created_at == created_at and created_at not in {float("inf"), float("-inf")}:
                return created_at
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            pass
        try:
            return path.stat().st_mtime
        except OSError:
            return 0.0

    def create_manifest(self, artifact_id: str, payload: dict[str, Any]) -> Path:
        self.prune(keep={artifact_id})
        directory = self.bundle_dir(artifact_id)
        manifest = {
            "manifest_version": 1,
            "artifact_id": artifact_id,
            "created_at_unix_seconds": time.time(),
            "lifecycle": {
                "storage": "ephemeral_staging",
                "retention": "bounded; completed bundles may be removed at any time",
                "preservation": "copy final evidence to doc/reports or a named output directory",
            },
            "service": dict(SERVICE_PROVENANCE),
            **payload,
        }
        path = directory / "manifest.json"
        if path.exists():
            raise FileExistsError(f"artifact manifest already exists: {artifact_id}")
        atomic_write_json(path, manifest)
        return path

    def get_manifest(self, artifact_id: str) -> dict[str, Any] | None:
        try:
            path = self.root / self._safe_id(artifact_id) / "manifest.json"
        except ValueError:
            return None
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return None
        return payload if isinstance(payload, dict) else None

    def write_summary(self, artifact_id: str, name: str, payload: dict[str, Any]) -> Path:
        """Add one immutable compact post-execution summary to a bundle."""
        if not name.endswith(".json") or Path(name).name != name:
            raise ValueError("artifact summary name must be a simple .json filename")
        path = self.bundle_dir(artifact_id) / name
        if not (self.root / artifact_id / "manifest.json").is_file():
            raise ValueError("artifact manifest must exist before a summary")
        if path.exists():
            raise FileExistsError(f"artifact summary already exists: {name}")
        atomic_write_json(path, payload)
        return path
