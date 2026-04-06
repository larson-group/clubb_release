"""NumPy port of CLUBB grid_class utilities.

This file keeps the same public API shape as `jax_examples/grid_class.py` but
uses pure NumPy and standard Python loops (no JAX dependency).
"""

from typing import NamedTuple, Optional, Tuple, Union

import numpy as np

try:
    from clubb_precision import core_rknd  # type: ignore
    from constants_clubb import one, zero, one_half  # type: ignore
except ImportError:
    core_rknd = np.float64
    one = 1.0
    zero = 0.0
    one_half = 0.5

Array = np.ndarray

# Interpolation weight indices
T_ABOVE = 0
T_BELOW = 1
M_ABOVE = 0
M_BELOW = 1

# Grid type constants
GRID_TYPE_EVEN = 1
GRID_TYPE_STRETCHED_ZT = 2
GRID_TYPE_STRETCHED_ZM = 3

LSCALE_SFCLYR_DEPTH = 500.0
NWARNING = 500


class Grid(NamedTuple):
    """Grid structure for CLUBB vertical discretization."""

    nzm: int
    nzt: int
    ngrdcol: int

    zm: Array
    zt: Array

    dzm: Array
    dzt: Array

    invrs_dzm: Array
    invrs_dzt: Array

    weights_zt2zm: Array
    weights_zm2zt: Array

    k_lb_zm: int
    k_ub_zm: int
    k_lb_zt: int
    k_ub_zt: int

    grid_dir_indx: int
    grid_dir: float


def setup_grid(
    ngrdcol: int,
    deltaz: Union[float, Array],
    zm_init: Union[float, Array],
    zm_top: Union[float, Array],
    l_ascending_grid: bool = True,
    grid_type: int = GRID_TYPE_EVEN,
    momentum_heights: Optional[Array] = None,
    thermodynamic_heights: Optional[Array] = None,
) -> Grid:
    """Initialize a CLUBB vertical grid structure.

    Notes:
      - `nzm`/`nzt` are computed following `grid_class.F90`.
      - For stretched grids in standalone mode, level selection is based on the
        first column's `zm_init`/`zm_top`, matching Fortran behavior.
        `nzm = floor((zm_top(1) - zm_init(1) + deltaz(1)) / deltaz(1))`
        `nzt = nzm - 1`
    """
    if isinstance(deltaz, (int, float)):
        deltaz_arr = np.full(ngrdcol, deltaz, dtype=core_rknd)
    else:
        deltaz_arr = np.asarray(deltaz, dtype=core_rknd)

    if isinstance(zm_init, (int, float)):
        zm_init_arr = np.full(ngrdcol, zm_init, dtype=core_rknd)
    else:
        zm_init_arr = np.asarray(zm_init, dtype=core_rknd)

    if isinstance(zm_top, (int, float)):
        zm_top_arr = np.full(ngrdcol, zm_top, dtype=core_rknd)
    else:
        zm_top_arr = np.asarray(zm_top, dtype=core_rknd)

    if grid_type not in (GRID_TYPE_EVEN, GRID_TYPE_STRETCHED_ZT, GRID_TYPE_STRETCHED_ZM):
        raise ValueError(f"Unsupported grid_type={grid_type}.")

    if grid_type == GRID_TYPE_EVEN:
        # Mirror Fortran standalone calculation (using first column):
        # gr%nzm = floor( (zm_top(1) - zm_init(1) + deltaz(1)) / deltaz(1) )
        nzm = int(np.floor((zm_top_arr[0] - zm_init_arr[0] + deltaz_arr[0]) / deltaz_arr[0]))
        nzt = nzm - 1
        if nzm < 2 or nzt < 1:
            raise ValueError(
                f"Invalid derived grid dimensions: nzm={nzm}, nzt={nzt}. "
                "Check zm_top, zm_init, and deltaz."
            )
        zm, zt = _setup_even_grid(nzm, nzt, ngrdcol, deltaz_arr, zm_init_arr, l_ascending_grid)
    elif grid_type == GRID_TYPE_STRETCHED_ZT:
        zt_in = _prepare_height_array("thermodynamic_heights", thermodynamic_heights, ngrdcol)
        _validate_monotonic_increasing("thermodynamic_heights", zt_in)

        if np.any(zt_in[:, 0] <= zm_init_arr):
            raise ValueError(
                "Stretched zt grid cannot fulfill zm_init requirement: "
                "lowest thermodynamic level must be above zm_init."
            )

        # Fortran standalone determines the active range from column 1.
        end_idx = int(np.searchsorted(zt_in[0], zm_top_arr[0], side="right") - 1)
        if end_idx < 0:
            raise ValueError(
                "Stretched zt grid cannot fulfill zm_top requirement: "
                "no thermodynamic level is <= zm_top."
            )
        zt_used = zt_in[:, : end_idx + 1]
        nzt = zt_used.shape[1]
        nzm = nzt + 1
        if nzt < 1:
            raise ValueError("Stretched zt grid produced no thermodynamic levels.")
        zm, zt = _setup_stretched_zt_grid(
            nzm, nzt, ngrdcol, zt_used, zm_init_arr, l_ascending_grid
        )
    else:  # grid_type == GRID_TYPE_STRETCHED_ZM
        zm_in = _prepare_height_array("momentum_heights", momentum_heights, ngrdcol)
        _validate_monotonic_increasing("momentum_heights", zm_in)

        begin_idx = int(np.searchsorted(zm_in[0], zm_init_arr[0], side="left"))
        if begin_idx >= zm_in.shape[1]:
            raise ValueError(
                "Stretched zm grid cannot fulfill zm_init requirement: "
                "no momentum level is >= zm_init."
            )
        end_idx = int(np.searchsorted(zm_in[0], zm_top_arr[0], side="right") - 1)
        if end_idx < begin_idx:
            raise ValueError(
                "Stretched zm grid cannot fulfill zm_top requirement: "
                "no valid momentum levels in [zm_init, zm_top]."
            )

        zm_used = zm_in[:, begin_idx : end_idx + 1]
        nzm = zm_used.shape[1]
        nzt = nzm - 1
        if nzm < 2 or nzt < 1:
            raise ValueError(
                f"Stretched zm grid produced invalid dimensions: nzm={nzm}, nzt={nzt}."
            )
        zm, zt = _setup_stretched_zm_grid(nzm, nzt, ngrdcol, zm_used)

    dzm, dzt, invrs_dzm, invrs_dzt = _calc_grid_spacings(
        nzm, nzt, ngrdcol, zm, zt, l_ascending_grid
    )

    weights_zt2zm = _calc_zt2zm_weights(nzm, nzt, ngrdcol, zm, zt)
    weights_zm2zt = _calc_zm2zt_weights(nzm, nzt, ngrdcol, zm, zt, dzt)

    if l_ascending_grid:
        # Keep Python-native 0-based bounds inside Grid.
        k_lb_zm, k_ub_zm = 0, nzm - 1
        k_lb_zt, k_ub_zt = 0, nzt - 1
        grid_dir_indx = 1
        grid_dir = 1.0
    else:
        # Descending grid, still in Python 0-based indexing.
        k_lb_zm, k_ub_zm = nzm - 1, 0
        k_lb_zt, k_ub_zt = nzt - 1, 0
        grid_dir_indx = -1
        grid_dir = -1.0

    return Grid(
        nzm=nzm,
        nzt=nzt,
        ngrdcol=ngrdcol,
        zm=zm,
        zt=zt,
        dzm=dzm,
        dzt=dzt,
        invrs_dzm=invrs_dzm,
        invrs_dzt=invrs_dzt,
        weights_zt2zm=weights_zt2zm,
        weights_zm2zt=weights_zm2zt,
        k_lb_zm=k_lb_zm,
        k_ub_zm=k_ub_zm,
        k_lb_zt=k_lb_zt,
        k_ub_zt=k_ub_zt,
        grid_dir_indx=grid_dir_indx,
        grid_dir=grid_dir,
    )


def _setup_even_grid(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    deltaz: Array,
    zm_init: Array,
    l_ascending_grid: bool,
) -> Tuple[Array, Array]:
    if l_ascending_grid:
        k_indices = np.arange(nzm, dtype=core_rknd)
    else:
        k_indices = np.arange(nzm - 1, -1, -1, dtype=core_rknd)

    zm = zm_init[:, None] + deltaz[:, None] * k_indices[None, :]
    zt = 0.5 * (zm[:, :-1] + zm[:, 1:])
    return zm, zt


def _setup_stretched_zt_grid(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    thermodynamic_heights: Array,
    zm_init: Array,
    l_ascending_grid: bool,
) -> Tuple[Array, Array]:
    zt = np.asarray(thermodynamic_heights, dtype=core_rknd)
    zm = _calc_zm_from_zt(nzm, nzt, ngrdcol, zt, zm_init, l_ascending_grid)
    return zm, zt


def _setup_stretched_zm_grid(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    momentum_heights: Array,
) -> Tuple[Array, Array]:
    zm = np.asarray(momentum_heights, dtype=core_rknd)
    zt = _calc_zt_from_zm(nzm, nzt, ngrdcol, zm)
    return zm, zt


def _calc_zm_from_zt(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    zt: Array,
    zm_init: Array,
    l_ascending_grid: bool,
) -> Array:
    zm = np.zeros((ngrdcol, nzm), dtype=core_rknd)
    zm[:, 1:-1] = 0.5 * (zt[:, :-1] + zt[:, 1:])
    if l_ascending_grid:
        zm[:, 0] = np.asarray(zm_init, dtype=core_rknd)
        if nzt > 1:
            zm[:, -1] = zt[:, -1] + 0.5 * (zt[:, -1] - zt[:, -2])
        else:
            zm[:, -1] = zt[:, -1] + (zt[:, -1] - zm[:, 0])
    else:
        zm[:, -1] = np.asarray(zm_init, dtype=core_rknd)
        if nzt > 1:
            zm[:, 0] = zt[:, 0] + 0.5 * (zt[:, 0] - zt[:, 1])
        else:
            zm[:, 0] = zt[:, 0] + (zt[:, 0] - zm[:, -1])
    return zm


def _calc_zt_from_zm(nzm: int, nzt: int, ngrdcol: int, zm: Array) -> Array:
    return 0.5 * (zm[:, :-1] + zm[:, 1:])


def _prepare_height_array(name: str, heights: Optional[Array], ngrdcol: int) -> Array:
    if heights is None:
        raise ValueError(f"{name} must be provided for stretched grid setup.")
    arr = np.asarray(heights, dtype=core_rknd)
    if arr.ndim == 1:
        arr = np.tile(arr[None, :], (ngrdcol, 1))
    if arr.ndim != 2:
        raise ValueError(f"{name} must be 1D or 2D, got shape {arr.shape}.")
    if arr.shape[0] != ngrdcol:
        raise ValueError(
            f"{name} first dimension must match ngrdcol={ngrdcol}, got {arr.shape[0]}."
        )
    return arr


def _validate_monotonic_increasing(name: str, arr: Array) -> None:
    if arr.shape[1] < 1:
        raise ValueError(f"{name} must contain at least one level.")
    if np.any(np.diff(arr, axis=1) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing with level index.")


def _calc_grid_spacings(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    zm: Array,
    zt: Array,
    l_ascending_grid: bool,
) -> Tuple[Array, Array, Array, Array]:
    dzm = np.zeros((ngrdcol, nzm), dtype=core_rknd)
    dzt = np.zeros((ngrdcol, nzt), dtype=core_rknd)

    # Match grid_class.F90 setup_grid_heights:
    # interior dzm uses signed zt differences (negative on descending grids).
    dzm[:, 1:-1] = zt[:, 1:] - zt[:, :-1]

    if l_ascending_grid:
        # Lower boundary uses twice the zt-zm half-distance.
        dzm[:, 0] = 2.0 * (zt[:, 0] - zm[:, 0])
        # Upper boundary copies adjacent interior spacing.
        dzm[:, -1] = dzm[:, -2]
    else:
        # Descending-grid lower boundary at index nzm (Python -1).
        dzm[:, -1] = 2.0 * (zm[:, -1] - zt[:, -1])
        # Descending-grid upper boundary at index 1 (Python 0).
        dzm[:, 0] = dzm[:, 1]

    # Match grid_class.F90: signed spacing between momentum levels.
    dzt = zm[:, 1:] - zm[:, :-1]

    invrs_dzm = np.where(np.abs(dzm) > 1e-30, 1.0 / dzm, 0.0)
    invrs_dzt = np.where(np.abs(dzt) > 1e-30, 1.0 / dzt, 0.0)

    return dzm, dzt, invrs_dzm, invrs_dzt


def _calc_zt2zm_weights(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    zm: Array,
    zt: Array,
) -> Array:
    weights = np.zeros((ngrdcol, nzm, 2), dtype=core_rknd)

    # Interior momentum levels (Fortran k=2..nzm-1):
    # interpolate between zt(k-1) and zt(k) onto zm(k).
    for k in range(1, nzm - 1):
        denom = zt[:, k] - zt[:, k - 1]
        weights[:, k, T_ABOVE] = (zm[:, k] - zt[:, k - 1]) / (denom + 1e-30)
        weights[:, k, T_BELOW] = (zt[:, k] - zm[:, k]) / (denom + 1e-30)

    if nzt >= 2:
        # Boundary momentum level k=1 (Python 0): linear extension using zt(1), zt(2).
        denom0 = zt[:, 1] - zt[:, 0]
        weights[:, 0, T_ABOVE] = (zm[:, 0] - zt[:, 0]) / (denom0 + 1e-30)
        weights[:, 0, T_BELOW] = (zt[:, 1] - zm[:, 0]) / (denom0 + 1e-30)

        # Boundary momentum level k=nzm (Python -1): extension using zt(nzt-1), zt(nzt).
        denomn = zt[:, -1] - zt[:, -2]
        weights[:, -1, T_ABOVE] = (zm[:, -1] - zt[:, -2]) / (denomn + 1e-30)
        weights[:, -1, T_BELOW] = (zt[:, -1] - zm[:, -1]) / (denomn + 1e-30)
    else:
        # Degenerate 2-level grid fallback; not used in standard CLUBB setups.
        weights[:, 0, T_ABOVE] = 1.0
        weights[:, 0, T_BELOW] = 0.0
        weights[:, -1, T_ABOVE] = 1.0
        weights[:, -1, T_BELOW] = 0.0

    return weights


def _calc_zm2zt_weights(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    zm: Array,
    zt: Array,
    dzt: Array,
) -> Array:
    weights = np.zeros((ngrdcol, nzt, 2), dtype=core_rknd)

    for k in range(nzt):
        dist_lower = np.abs(zt[:, k] - zm[:, k])
        total_dist = dzt[:, k] + 1e-30
        w_upper = 1.0 - dist_lower / total_dist
        w_lower = dist_lower / total_dist
        weights[:, k, M_ABOVE] = w_upper
        weights[:, k, M_BELOW] = w_lower

    return weights


def zt2zm_api(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azt: Array,
    zm_min: Optional[float] = None,
) -> Array:
    azm = np.zeros((ngrdcol, nzm), dtype=core_rknd)

    for k in range(1, nzm - 1):
        w_upper = gr.weights_zt2zm[:, k, T_ABOVE]
        w_lower = gr.weights_zt2zm[:, k, T_BELOW]
        k_upper = int(np.clip(k - 1, 0, nzt - 1))
        k_lower = int(np.clip(k, 0, nzt - 1))
        azm[:, k] = w_upper * azt[:, k_upper] + w_lower * azt[:, k_lower]

    if gr.grid_dir_indx == 1:
        factor = gr.dzm[:, 0] / (gr.dzt[:, 0] + 1e-30)
        azm[:, 0] = azt[:, 0] - (azt[:, 1] - azt[:, 0]) * factor
    else:
        factor = gr.dzm[:, 0] / (gr.dzt[:, 0] + 1e-30)
        azm[:, 0] = azt[:, 0] + (azt[:, 0] - azt[:, 1]) * factor

    if gr.grid_dir_indx == 1:
        factor = gr.dzm[:, -1] / (gr.dzt[:, -1] + 1e-30)
        azm[:, -1] = azt[:, -1] + (azt[:, -1] - azt[:, -2]) * factor
    else:
        factor = gr.dzm[:, -1] / (gr.dzt[:, -1] + 1e-30)
        azm[:, -1] = azt[:, -1] - (azt[:, -2] - azt[:, -1]) * factor

    if zm_min is not None:
        azm = np.maximum(azm, zm_min)

    return azm


def zm2zt_api(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azm: Array,
    zt_min: Optional[float] = None,
) -> Array:
    azt = np.zeros((ngrdcol, nzt), dtype=core_rknd)

    for k in range(nzt):
        w_upper = gr.weights_zm2zt[:, k, M_ABOVE]
        w_lower = gr.weights_zm2zt[:, k, M_BELOW]
        k_upper = k
        k_lower = int(np.clip(k + 1, 0, nzm - 1))
        azt[:, k] = w_upper * azm[:, k_upper] + w_lower * azm[:, k_lower]

    if zt_min is not None:
        azt = np.maximum(azt, zt_min)

    return azt


def zt2zm2zt(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azt: Array,
    zt_min: Optional[float] = None,
) -> Array:
    azm = zt2zm_api(nzm, nzt, ngrdcol, gr, azt)
    return zm2zt_api(nzm, nzt, ngrdcol, gr, azm, zt_min)


def zm2zt2zm(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azm: Array,
    zm_min: Optional[float] = None,
) -> Array:
    azt = zm2zt_api(nzm, nzt, ngrdcol, gr, azm)
    return zt2zm_api(nzm, nzt, ngrdcol, gr, azt, zm_min)


def ddzm(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azm: Array,
) -> Array:
    return (azm[:, 1:] - azm[:, :-1]) * gr.invrs_dzt


def ddzt(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    azt: Array,
) -> Array:
    dazt_dz = np.zeros((ngrdcol, nzm), dtype=core_rknd)
    for k in range(1, nzm - 1):
        dazt_dz[:, k] = (azt[:, k] - azt[:, k - 1]) * gr.invrs_dzm[:, k]
    return dazt_dz


def flip_vertical(x: Array, axis: int = -1) -> Array:
    return np.flip(x, axis=axis)


__all__ = [
    'Grid',
    'setup_grid',
    'zt2zm_api',
    'zm2zt_api',
    'zt2zm2zt',
    'zm2zt2zm',
    'ddzm',
    'ddzt',
    'flip_vertical',
]
