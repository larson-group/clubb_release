"""Process-isolated Dash callback support for expensive local analysis."""

from __future__ import annotations

import os
from pathlib import Path

from dash import DiskcacheManager

from .runtime import private_runtime_dir


def create_background_manager(repo_root: Path | str) -> DiskcacheManager:
    """Build the local process manager used by explicitly expensive callbacks.

    Diskcache starts a separate process per active callback.  That keeps the
    main Dash server responsive while a case is opened and, critically, avoids
    concurrent HDF5 access from Dash request threads.  Its state is placed in
    the same private runtime area as the local agent broker rather than in the
    checkout or a shared ``/tmp`` directory.
    """

    try:
        import diskcache
    except ImportError as exc:  # pragma: no cover - launcher installs runtime deps
        raise RuntimeError(
            "Dash background callbacks require diskcache. Run "
            "./launch_dashboard.sh so dash_app/requirements.txt is installed."
        ) from exc

    cache_path = private_runtime_dir(repo_root) / "background-callbacks"
    cache_path.mkdir(mode=0o700, parents=True, exist_ok=True)
    cache_path.chmod(0o700)
    # Diskcache creates its SQLite files on first construction.  Do that under
    # a private umask so a fresh local worker cache cannot leak callback inputs
    # or results even if a parent runtime directory is later misconfigured.
    previous_umask = os.umask(0o077)
    try:
        cache = diskcache.Cache(str(cache_path))
    finally:
        os.umask(previous_umask)
    return DiskcacheManager(cache, expire=60 * 60)
