"""Private runtime-state helpers shared by the local dashboard broker."""

from __future__ import annotations

import hashlib
import json
import os
import re
import stat
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any


def process_is_alive(pid: Any) -> bool:
    """Return whether a PID exists and is not a zombie."""
    import psutil

    try:
        process = psutil.Process(int(pid))
        return process.status() != psutil.STATUS_ZOMBIE
    except (psutil.Error, TypeError, ValueError, OSError):
        return False


def read_file_chunk(
    path: Path | str, cursor: int, max_bytes: int
) -> tuple[bytes, int, bool]:
    """Read one bounded byte range, resetting a cursor beyond a replaced file."""
    offset = max(0, int(cursor))
    try:
        file_path = Path(path)
        size = file_path.stat().st_size
        if offset > size:
            offset = 0
        with file_path.open("rb") as handle:
            handle.seek(offset)
            chunk = handle.read(max(1, int(max_bytes)))
    except (OSError, TypeError, ValueError):
        return b"", offset, True
    next_cursor = offset + len(chunk)
    return chunk, next_cursor, next_cursor >= size


def read_file_tail(path: Path | str, max_lines: int) -> tuple[bytes, int, bool]:
    """Read the newest complete log lines and position the cursor at EOF."""
    try:
        file_path = Path(path)
        size = file_path.stat().st_size
        position = size
        chunks = []
        newline_count = 0
        with file_path.open("rb") as handle:
            while position > 0 and newline_count <= max_lines:
                length = min(65536, position)
                position -= length
                handle.seek(position)
                chunk = handle.read(length)
                chunks.append(chunk)
                newline_count += chunk.count(b"\n")
    except (OSError, TypeError, ValueError):
        return b"", 0, True
    lines = b"".join(reversed(chunks)).splitlines(keepends=True)
    return b"".join(lines[-max_lines:]), size, True


_RUN_PROGRESS_RE = re.compile(rb"iteration:\s*(\d+)\s*/\s*(\d+)", re.IGNORECASE)


def read_latest_run_progress(path: Path | str, max_bytes: int = 8192) -> tuple[int, int] | None:
    """Read only a log's tail and return its latest CLUBB iteration pair."""
    try:
        file_path = Path(path)
        size = file_path.stat().st_size
        with file_path.open("rb") as handle:
            handle.seek(max(0, size - max(1, int(max_bytes))))
            matches = list(_RUN_PROGRESS_RE.finditer(handle.read()))
    except (OSError, TypeError, ValueError):
        return None
    if not matches:
        return None
    return int(matches[-1].group(1)), int(matches[-1].group(2))


def _runtime_base_candidates(uid: str) -> list[Path]:
    """Return private-runtime roots in discovery order for one local user."""

    candidates = [Path(os.environ["XDG_RUNTIME_DIR"])] if os.environ.get("XDG_RUNTIME_DIR") else []
    # A dashboard launched from a graphical/login shell commonly receives
    # XDG_RUNTIME_DIR=/run/user/<uid>, while a separately launched Codex/MCP
    # client may not inherit that environment variable.  Both processes still
    # need to discover the *same* private broker record.  Prefer the standard
    # per-user runtime location before the intentionally less-preferred /tmp
    # fallback when the environment is absent.
    standard_runtime = Path("/run/user") / uid
    if standard_runtime not in candidates:
        candidates.append(standard_runtime)
    candidates.append(Path("/tmp"))
    return candidates


def private_runtime_candidates(repo_root: Path | str, name: str) -> list[Path]:
    """Return every plausible private runtime path for one checkout.

    A Dash process launched from a graphical shell and an MCP adapter launched
    by an editor do not always inherit the same ``XDG_RUNTIME_DIR``.  They must
    nevertheless discover the same broker record.  Writers still use
    :func:`private_runtime_dir`; readers use this ordered list to find an
    already-published record before treating the broker as absent.
    """

    root = Path(repo_root).resolve()
    token = hashlib.sha256(str(root).encode("utf-8")).hexdigest()[:16]
    uid = str(os.getuid()) if hasattr(os, "getuid") else "user"
    directory_name = f"clubb_dash_{uid}_{token}"
    seen: set[Path] = set()
    paths: list[Path] = []
    for base in _runtime_base_candidates(uid):
        candidate = (base / directory_name / name)
        if candidate not in seen:
            seen.add(candidate)
            paths.append(candidate)
    return paths


def readable_private_paths(repo_root: Path | str, name: str) -> list[Path]:
    """Return existing, owner-private candidate files safe to read.

    This deliberately does not create a directory.  It is for broker discovery
    only, where creating an empty XDG directory must not hide a valid record
    published under another per-user runtime root.
    """

    uid = os.getuid() if hasattr(os, "getuid") else None
    result: list[Path] = []
    for path in private_runtime_candidates(repo_root, name):
        try:
            info = path.stat()
        except OSError:
            continue
        if not stat.S_ISREG(info.st_mode):
            continue
        if uid is not None and info.st_uid != uid:
            continue
        if stat.S_IMODE(info.st_mode) & 0o077:
            continue
        result.append(path)
    return result


def private_runtime_dir(repo_root: Path | str) -> Path:
    """Return a user-private runtime directory for one checkout.

    The broker is deliberately local-only, but `/tmp` is often shared on HPC
    systems.  State therefore lives in a mode-0700 per-user directory instead
    of directly in `/tmp`.
    """

    root = Path(repo_root).resolve()
    token = hashlib.sha256(str(root).encode("utf-8")).hexdigest()[:16]
    uid = str(os.getuid()) if hasattr(os, "getuid") else "user"
    name = f"clubb_dash_{uid}_{token}"
    for base in _runtime_base_candidates(uid):
        directory = base / name
        try:
            directory.mkdir(mode=0o700, parents=True, exist_ok=True)
            directory.chmod(0o700)
            return directory
        except OSError:
            # Containers sometimes expose an XDG path that is read-only to
            # worker processes.  The UID-scoped /tmp fallback is still 0700.
            continue
    raise RuntimeError("could not create a private dashboard runtime directory")


def private_path(repo_root: Path | str, name: str) -> Path:
    return private_runtime_dir(repo_root) / name


@contextmanager
def exclusive_file_lock(path: Path):
    """Hold a mode-0600 advisory lock shared by independent local clients."""
    import fcntl

    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write JSON without leaving world-readable temporary files."""

    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        os.fchmod(fd, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        path.chmod(0o600)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def restrict_existing(path: Path) -> None:
    """Best-effort migration for a legacy state file created under `/tmp`."""

    try:
        if path.exists():
            path.chmod(0o600)
    except OSError:
        pass
