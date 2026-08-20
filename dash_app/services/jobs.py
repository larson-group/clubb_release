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
from dash_app.shared.runtime import atomic_write_json, exclusive_file_lock, private_path


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


class JobTransaction:
    """Mutate several job records with one registry read and atomic write."""

    def __init__(self, state: dict[str, Any]):
        self.state = state
        self.keys_by_job_id = {
            str(record.get("job_id")): key
            for key, record in state.items()
            if isinstance(record, dict) and record.get("job_id")
        }

    def _bump_scm_revision(self) -> None:
        metadata = dict(self.state.get("__metadata__") or {})
        metadata["scm_revision"] = int(metadata.get("scm_revision") or 0) + 1
        self.state["__metadata__"] = metadata

    def submit(self, kind: str, request_id: str, request: dict[str, Any]) -> tuple[dict[str, Any], bool]:
        key = f"{kind}:{request_id}"
        request_hash = canonical_hash(request)
        existing = self.state.get(key)
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
        self.state[key] = record
        self.keys_by_job_id[record["job_id"]] = key
        if kind in {"scm", "scm_batch"}:
            self._bump_scm_revision()
        return dict(record), True

    def update(self, job_id: str, **updates: Any) -> dict[str, Any] | None:
        key = self.keys_by_job_id.get(str(job_id))
        record = self.state.get(key) if key is not None else None
        if not isinstance(record, dict):
            return None
        record.update(updates)
        record["updated_at_unix_seconds"] = time.time()
        self.state[key] = record
        if record.get("kind") in {"scm", "scm_batch"}:
            self._bump_scm_revision()
        return dict(record)

    def get(self, job_id: str) -> dict[str, Any] | None:
        key = self.keys_by_job_id.get(str(job_id))
        record = self.state.get(key) if key is not None else None
        return dict(record) if isinstance(record, dict) else None

    def delete(self, job_ids: set[str]) -> None:
        deleted = {str(job_id) for job_id in job_ids}
        deleted_scm = any(
            isinstance(record, dict)
            and record.get("kind") in {"scm", "scm_batch"}
            and str(record.get("job_id") or "") in deleted
            for record in self.state.values()
        )
        self.state = {
            key: record
            for key, record in self.state.items()
            if not isinstance(record, dict)
            or str(record.get("job_id") or "") not in deleted
        }
        for job_id in deleted:
            self.keys_by_job_id.pop(job_id, None)
        if deleted_scm:
            self._bump_scm_revision()


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

    def transaction(self, operation):
        """Run one registry operation with one read and one atomic write."""
        with self._locked():
            transaction = JobTransaction(self._read())
            result = operation(transaction)
            self._write(transaction.state)
            return result

    def snapshot(self) -> dict[str, dict[str, Any]]:
        """Return the complete small registry for indexed batch reconciliation."""
        return {
            str(key): dict(record)
            for key, record in self._read().items()
            if isinstance(record, dict) and record.get("job_id")
        }

    def scm_revision(self) -> int:
        """Return the durable revision for user-facing SCM lifecycle changes."""
        metadata = dict(self._read().get("__metadata__") or {})
        return int(metadata.get("scm_revision") or 0)

    def scm_snapshot(self) -> tuple[int, list[dict[str, Any]]]:
        """Read the canonical SCM records and revision from one file snapshot."""
        payload = self._read()
        metadata = dict(payload.get("__metadata__") or {})
        records = [
            dict(record)
            for record in payload.values()
            if isinstance(record, dict) and record.get("kind") == "scm"
        ]
        return int(metadata.get("scm_revision") or 0), records

    def scm_summary(self) -> dict[str, Any]:
        """Return compact user-visible SCM counts for global status surfaces."""
        revision, records = self.scm_snapshot()
        counts = {"queued": 0, "running": 0, "stopping": 0}
        for record in records:
            if str(record.get("visibility") or "user") != "user":
                continue
            state = str(record.get("state") or "")
            if state in {"submitting", "starting"}:
                state = "queued"
            if state in counts:
                counts[state] += 1
        active_count = sum(counts.values())
        return {
            "state": "running" if active_count else "idle",
            "revision": revision,
            "active_count": active_count,
            **counts,
        }

    def submit(self, kind: str, request_id: str, request: dict[str, Any]) -> tuple[dict[str, Any], bool]:
        return self.transaction(lambda transaction: transaction.submit(kind, request_id, request))

    def update(self, job_id: str, **updates: Any) -> dict[str, Any] | None:
        return self.transaction(lambda transaction: transaction.update(job_id, **updates))

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

    def _prune_unlocked(self, *, keep: set[str]) -> list[str]:
        """Remove excess inactive bundles while the retention lock is held."""
        try:
            bundles = [
                path
                for path in self.root.iterdir()
                if path.is_dir() and (path / "manifest.json").is_file()
            ]
        except OSError:
            return []
        limit = self.retention_count()
        candidates = sorted(
            (
                path
                for path in bundles
                if path.name not in keep
                and not (path / self._ACTIVE_MARKER).exists()
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

            try:
                shutil.rmtree(path)
            except FileNotFoundError:
                # Manual cleanup or a pre-lock service process may have won
                # the deletion race.  The requested end state already holds.
                pass
            removed.append(path.name)
        return removed

    def prune(self, *, keep: set[str] | None = None) -> list[str]:
        """Remove oldest inactive artifact bundles beyond the configured cap.

        Only direct service-owned children with manifests are considered.  A
        live bundle has an explicit marker, so a later submission cannot prune
        output still being read or written by an active agent job.  Unmarked
        legacy directories are deliberately ignored rather than guessed at.
        """
        with exclusive_file_lock(self.root / ".prune.lock"):
            return self._prune_unlocked(keep=set(keep or set()))

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

    def create_manifest(
        self,
        artifact_id: str,
        payload: dict[str, Any],
        *,
        active: bool = False,
    ) -> Path:
        """Create a manifest and, when requested, protect it atomically."""
        with exclusive_file_lock(self.root / ".prune.lock"):
            self._prune_unlocked(keep={artifact_id})
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
                raise FileExistsError(
                    f"artifact manifest already exists: {artifact_id}"
                )
            atomic_write_json(path, manifest)
            if active:
                marker = directory / self._ACTIVE_MARKER
                atomic_write_json(
                    marker,
                    {
                        "artifact_id": self._safe_id(artifact_id),
                        "started_at_unix_seconds": time.time(),
                    },
                )
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
