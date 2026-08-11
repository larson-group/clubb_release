"""NetCDF loading and summaries for the fill-holes viewer."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


MOMENTS = {
    "wp2": {"label": "w′²", "grid": "zm"},
    "rtp2": {"label": "rₜ′²", "grid": "zm"},
    "thlp2": {"label": "θₗ′²", "grid": "zm"},
    "up2": {"label": "u′²", "grid": "zm"},
    "vp2": {"label": "v′²", "grid": "zm"},
    "rtm": {"label": "rₜ", "grid": "zt"},
    "thlm": {"label": "θₗ", "grid": "zt"},
}

THRESHOLDS = {
    "wp2": 4.0e-4,
    "rtp2": 1.0e-16,
    "thlp2": 1.0e-4,
    "up2": 4.0e-4,
    "vp2": 4.0e-4,
    "rtm": 1.0e-8,
    "thlm": 1.0e-2,
}

THRESHOLD_NAMES = {
    "wp2": "w_tol_sqd",
    "rtp2": "rt_tol**2",
    "thlp2": "thl_tol**2",
    "up2": "w_tol_sqd",
    "vp2": "w_tol_sqd",
    "rtm": "rt_tol",
    "thlm": "thl_tol",
}


def diagnostic_names(moment: str) -> tuple[str, str]:
    return (
        f"{moment}_before_hf",
        f"{moment}_after_hf",
    )


def _as_float(values) -> np.ndarray:
    return np.asarray(np.ma.filled(values, np.nan), dtype=float)


def _profile_array(variable) -> np.ndarray:
    """Return a profile variable in (time, level, column) order."""
    values = _as_float(variable[:])
    dimensions = list(variable.dimensions)
    try:
        order = [dimensions.index("time")]
        order.append(next(i for i, name in enumerate(dimensions) if name in {"zt", "zm"}))
        order.append(dimensions.index("col"))
    except (StopIteration, ValueError) as error:
        raise ValueError(
            f"{variable.name} does not have time, vertical-grid, and col dimensions"
        ) from error
    return np.transpose(values, order)


def inspect_dataset(path: str | Path) -> dict:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ValueError(f"Stats file is unavailable: {source}")

    with Dataset(source) as dataset:
        available = [
            moment
            for moment in MOMENTS
            if all(name in dataset.variables for name in diagnostic_names(moment))
        ]
        if not available:
            raise ValueError(
                "This file has no dedicated CLUBB hole-filling snapshots. "
                "Rebuild CLUBB and rerun the case with the updated all_stats.in."
            )
        times = _as_float(dataset.variables["time"][:]).reshape(-1)
        column_count = len(dataset.dimensions.get("col", ()))
        units = {
            moment: str(getattr(dataset.variables[diagnostic_names(moment)[0]], "units", ""))
            for moment in available
        }

    return {
        "path": str(source),
        "moments": available,
        "times": times.tolist(),
        "record_count": len(times),
        "column_count": column_count,
        "units": units,
    }


@lru_cache(maxsize=32)
def _load_moment_cached(
    path: str, modified_ns: int, size_bytes: int, moment: str
) -> dict:
    del modified_ns, size_bytes
    with Dataset(path) as dataset:
        before_name, after_name = diagnostic_names(moment)
        before = _profile_array(dataset.variables[before_name])
        after = _profile_array(dataset.variables[after_name])
        grid_name = MOMENTS[moment]["grid"]
        height = _as_float(dataset.variables[grid_name][:]).reshape(-1)
        times = _as_float(dataset.variables["time"][:]).reshape(-1)
        column_count = len(dataset.dimensions["col"])
        threshold = np.full(
            (len(times), column_count), THRESHOLDS[moment], dtype=float
        )
        units = str(getattr(dataset.variables[before_name], "units", ""))
    return {
        "before": before,
        "after": after,
        "threshold": threshold,
        "height": height,
        "times": times,
        "units": units,
    }


def load_moment(path: str | Path, moment: str) -> dict:
    if moment not in MOMENTS:
        raise ValueError(f"Unknown CLUBB moment: {moment}")
    source = Path(path).expanduser().resolve()
    stat = source.stat()
    return _load_moment_cached(str(source), stat.st_mtime_ns, stat.st_size, moment)


def below_threshold_counts(data: dict, column: int) -> tuple[np.ndarray, np.ndarray]:
    threshold = data["threshold"][:, column, np.newaxis]
    before = data["before"][:, :, column]
    after = data["after"][:, :, column]
    valid_before = np.isfinite(before) & np.isfinite(threshold)
    valid_after = np.isfinite(after) & np.isfinite(threshold)
    return (
        np.sum(valid_before & (before < threshold), axis=1),
        np.sum(valid_after & (after < threshold), axis=1),
    )
