"""Small, reusable helpers for reading Dash NetCDF data consistently.

The plotting and analysis tabs use file metadata as cache keys and need the
same dimension and time-coordinate conventions.  Keeping those primitives in
one module prevents subtle differences between tabs without coupling their
higher-level data models.
"""

from __future__ import annotations

import os
import re
from datetime import datetime
from typing import Any, Iterable

import numpy as np


_TIME_UNITS_RE = re.compile(r"^\s*([A-Za-z]+)\s+since\s+(.+?)\s*$")


def file_signature(path: str | os.PathLike[str]) -> tuple[Any, int | None, int | None]:
    """Return a cheap, hashable cache signature for *path*.

    Missing and temporarily inaccessible files retain a stable sentinel
    signature.  This matches the plotting caches' historical behavior.
    """

    try:
        stat = os.stat(path)
    except OSError:
        return (path, None, None)
    return (path, int(stat.st_mtime_ns), int(stat.st_size))


def file_metadata_signature(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Return serializable file metadata for analysis results.

    Unlike :func:`file_signature`, this helper raises when the file cannot be
    inspected.  Callers use it only after validating that an input exists.
    """

    stat = os.stat(path)
    return {
        "path": path,
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def find_dimension(
    dimension_names: Iterable[str], candidates: Iterable[str]
) -> str | None:
    """Return the first dimension whose name matches *candidates*.

    Matching is case-insensitive while the original dataset spelling is
    preserved in the returned name.
    """

    lowered_candidates = {str(candidate).lower() for candidate in candidates}
    for dimension_name in dimension_names:
        if str(dimension_name).lower() in lowered_candidates:
            return dimension_name
    return None


def coordinate_values(
    dataset: Any,
    dimension_name: str | None,
    length: int,
    *,
    dtype: Any = np.float64,
    masked_value: float = np.nan,
) -> np.ndarray:
    """Read a one-dimensional coordinate or return an index coordinate.

    NetCDF dimensions do not have to have a matching coordinate variable.  A
    numeric ``0..length-1`` fallback gives every caller the same behavior for
    those datasets.
    """

    if dimension_name and dimension_name in dataset.variables:
        values = dataset.variables[dimension_name][:]
        if np.ma.isMaskedArray(values):
            values = np.ma.filled(values, masked_value)
        return np.asarray(values, dtype=dtype).reshape(-1)
    return np.arange(int(length), dtype=dtype)


def time_units_factor(units: str | None) -> float:
    """Return the multiplier that converts common NetCDF time units to seconds."""

    if not units:
        return 1.0
    unit = units.split()[0].lower()
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return 1.0
    if unit in {"min", "mins", "minute", "minutes"}:
        return 60.0
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return 3600.0
    return 1.0


def _parse_time_origin(origin_text: str | None) -> datetime | None:
    text = str(origin_text or "").strip().replace("T", " ")
    if text.endswith("Z"):
        text = text[:-1].strip()
    for date_format in (
        "%Y-%m-%d %H:%M:%S.%f",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
    ):
        try:
            return datetime.strptime(text, date_format)
        except ValueError:
            continue
    return None


def time_values_to_seconds(values: Any, units: str | None) -> np.ndarray:
    """Convert time coordinates to seconds using existing Dash conventions.

    For ``"units since DATE"`` coordinates, the returned values include the
    origin's offset from midnight.  This lets benchmark files with different
    same-day origins align while preserving the established plot behavior.
    """

    array = np.asarray(values, dtype=float)
    match = _TIME_UNITS_RE.match(str(units or ""))
    if match:
        origin = _parse_time_origin(match.group(2))
        if origin is not None:
            midnight = origin.replace(hour=0, minute=0, second=0, microsecond=0)
            return (
                array * float(time_units_factor(match.group(1)))
                + (origin - midnight).total_seconds()
            )
    return array * float(time_units_factor(units))
