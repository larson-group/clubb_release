"""Keep derived Plot metadata on the server, shared with case-loading workers.

Browser stores carry a small descriptor, not time coordinates or variable
catalogs. Snapshots are immutable and content-addressed, so a late callback
cannot accidentally resolve another case's metadata. Evicted snapshots can be
rebuilt from the descriptor; user selections never depend on cache lifetime.
"""

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path

from dash_app.shared.runtime import private_runtime_dir


_DERIVED_FIELDS = frozenset({
    "time_seconds", "time_elapsed_seconds", "time_bounds_seconds",
    "time_bounds_elapsed_seconds", "profile_vars", "budget_groups",
    "timeheight_vars", "subcolumn_vars", "timeseries_vars",
})


@lru_cache(maxsize=1)
def _cache():
    import diskcache

    directory = private_runtime_dir(Path(__file__).resolve().parents[2]) / "plot-metadata-v1"
    directory.mkdir(mode=0o700, parents=True, exist_ok=True)
    previous = os.umask(0o077)
    try:
        return diskcache.Cache(str(directory), size_limit=64 * 1024 * 1024,
                              eviction_policy="least-recently-used")
    finally:
        os.umask(previous)


@lru_cache(maxsize=32)
def _load_metadata(key):
    encoded = _cache().get(key)
    return json.loads(encoded) if encoded is not None else None


def compact_case_data(case_data):
    """Publish a full metadata snapshot and return its small browser descriptor."""
    if not case_data:
        return case_data
    # Accept an existing descriptor without replacing its immutable snapshot.
    if "metadata_key" in case_data:
        return dict(case_data)
    derived = {key: case_data[key] for key in _DERIVED_FIELDS if key in case_data}
    encoded = json.dumps(derived, sort_keys=True, separators=(",", ":"))
    identity = {field: case_data.get(field) for field in ("name", "files", "stats_fingerprints")}
    digest = hashlib.sha256(encoded.encode())
    digest.update(json.dumps(identity, sort_keys=True).encode())
    key = digest.hexdigest()
    _cache().set(key, encoded)
    # A previous cache miss must not hide a snapshot just republished by us.
    _load_metadata.cache_clear()
    return {**{key: value for key, value in case_data.items() if key not in _DERIVED_FIELDS},
            "metadata_key": key}


def resolve_case_data(case_data):
    """Resolve browser descriptors; ordinary service/test metadata passes through."""
    if not case_data or "metadata_key" not in case_data:
        return case_data
    derived = _load_metadata(case_data["metadata_key"])
    if derived is None:
        from dash_app.services.profiles import build_case_metadata

        rebuilt = build_case_metadata(case_data["name"], case_data.get("files") or [],
                                     case_data.get("output_dirs") or [])
        derived = {key: rebuilt[key] for key in _DERIVED_FIELDS if key in rebuilt}
        _cache().set(case_data["metadata_key"], json.dumps(derived))
        _load_metadata.cache_clear()
    return {**derived, **{key: value for key, value in case_data.items() if key != "metadata_key"}}
