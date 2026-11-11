"""Static metadata shared by JAX stats accumulation and host NetCDF output."""

from __future__ import annotations

from dataclasses import dataclass


GRID_ORDER = ("zt", "zm", "sfc", "lh_zt", "lh_sfc", "rad_zt", "rad_zm")


@dataclass(frozen=True)
class StatsLayout:
    """Immutable registry layout used to initialize fixed-shape JAX banks."""

    names: tuple[str, ...]
    grids: tuple[str, ...]
    grid_nlev: tuple[int, ...]
    samples_per_write: int
