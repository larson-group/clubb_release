"""JAX implementations of selected ``grid_class.F90`` helpers.

Description:

Definition of a grid class and associated functions

The grid specification is as follows for an ASCENDING grid:

    +                ================== zm(nzm) =================== top
    |
    |
1/dzt(nzt)   +       ------------------ zt(nzt) -------------------
    |        |
    |        |
    +  1/dzm(nzm-1)  ================== zm(nzm-1) =================
             |
             |
             +       ------------------ zt(nzt-1) -----------------

                                          .
                                          .
                                          .
                                          .

                     ================== zm(k+1) ===================


                     ------------------ zt(k+1) -------------------


    +                ================== zm(k+1) ===================
    |
    |
1/dzt(k)     +       ------------------ zt(k) ---------------------
    |        |
    |        |
    +    1/dzm(k)    ================== zm(k) =====================
             |
             |
             +       ------------------ zt(k-1) -------------------


                     ================== zm(k-1) ===================


                     ------------------ zt(k-2) -------------------

                                          .
                                          .
                                          .
                                          .

             +       ------------------ zt(2) ---------------------
             |
             |
    +    1/dzm(2)    ================== zm(2) =====================
    |        |
    |        |
1/dzt(1)     +       ------------------ zt(1) ---------------------
    |
    |
    +                ================== zm(1) =====================  zm_init
                     //////////////////////////////////////////////  surface


The variable zm(k) stands for the momentum level altitude at momentum
level k; the variable zt(k) stands for the thermodynamic level altitude at
thermodynamic level k; the variable invrs_dzt(k) is the inverse distance
between momentum levels (over a central thermodynamic level k); and the
variable invrs_dzm(k) is the inverse distance between thermodynamic levels
(over a central momentum level k).  Please note that in the above diagram,
"invrs_dzt" is denoted "dzt", and "invrs_dzm" is denoted "dzm", such that
1/dzt is the distance between successive momentum levels k and k+1 (over a
central thermodynamic level k), and 1/dzm is the distance between successive
thermodynamic levels k-1 and k (over a central momentum level k).

The grid setup is compatible with a stretched (unevely-spaced) grid.  Thus,
the distance between successive grid levels may not always be constant.

NOTE:  Any future code written for use in the CLUBB parameterization should
       use interpolation formulas consistent with a stretched grid.  The
       simplest way to do so is to call the appropriate interpolation
       function from this module.  Interpolations should *not* be handled in
       the form of:  ( var_zm(k+1) + var_zm(k) ) / 2; *nor* in the form of:
       0.5*( var_zt(k) + var_zt(k-1) ).

References:

https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_grid

Section 3c, p. 3548 /Numerical discretization/ of:
 ``A PDF-Based Model for Boundary Layer Clouds. Part I:
   Method and Model Description'' Golaz, et al. (2002)
   JAS, Vol. 59, pp. 3540--3551.

Porting deviations:
  * Fortran generic interfaces include scalar/1D/2D overloads and optional
    cubic interpolation through ``l_cubic_interp``.  JAX implements the active
    2D linear operator path used by the core.
  * Fortran band selector constants are 1-based; JAX uses Python 0-based
    ``T_ABOVE``, ``T_BELOW``, ``M_ABOVE``, and ``M_BELOW``.
  * Fortran mutates output arrays.  JAX returns arrays.
"""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

core_rknd = np.float64

T_ABOVE = 0
T_BELOW = 1
M_ABOVE = 0
M_BELOW = 1

GRID_TYPE_EVEN = 1
GRID_TYPE_STRETCHED_ZT = 2
GRID_TYPE_STRETCHED_ZM = 3


_GRID_ARRAY_FIELDS = (
    "zm",
    "zt",
    "dzm",
    "dzt",
    "invrs_dzm",
    "invrs_dzt",
    "weights_zt2zm",
    "weights_zm2zt",
)

_GRID_STATIC_FIELDS = (
    "nzm",
    "nzt",
    "ngrdcol",
    "k_lb_zm",
    "k_ub_zm",
    "k_lb_zt",
    "k_ub_zt",
    "grid_dir_indx",
    "grid_dir",
)


@jax.tree_util.register_pytree_node_class
class Grid(NamedTuple):
    """Grid structure for CLUBB vertical discretization."""

    nzm: int
    nzt: int
    ngrdcol: int

    zm: jnp.ndarray
    zt: jnp.ndarray

    dzm: jnp.ndarray
    dzt: jnp.ndarray

    invrs_dzm: jnp.ndarray
    invrs_dzt: jnp.ndarray

    weights_zt2zm: jnp.ndarray
    weights_zm2zt: jnp.ndarray

    k_lb_zm: int
    k_ub_zm: int
    k_lb_zt: int
    k_ub_zt: int

    grid_dir_indx: int
    grid_dir: float

    def replace(self, **kwargs):
        return self._replace(**kwargs)

    def tree_flatten(self):
        children = tuple(getattr(self, name) for name in _GRID_ARRAY_FIELDS)
        aux_data = tuple(getattr(self, name) for name in _GRID_STATIC_FIELDS)
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        data = dict(zip(_GRID_STATIC_FIELDS, aux_data))
        data.update(dict(zip(_GRID_ARRAY_FIELDS, children)))
        return cls(**data)


def setup_grid(
    ngrdcol: int,
    deltaz,
    zm_init,
    zm_top,
    l_ascending_grid: bool = True,
    grid_type: int = GRID_TYPE_EVEN,
    momentum_heights=None,
    thermodynamic_heights=None,
) -> Grid:
    """Construct a JAX `Grid` with CLUBB-style heights, spacings, and weights."""
    deltaz_arr = _as_column_input(deltaz, ngrdcol, "deltaz")
    zm_init_arr = _as_column_input(zm_init, ngrdcol, "zm_init")
    zm_top_arr = _as_column_input(zm_top, ngrdcol, "zm_top")

    if grid_type not in (
        GRID_TYPE_EVEN,
        GRID_TYPE_STRETCHED_ZT,
        GRID_TYPE_STRETCHED_ZM,
    ):
        raise ValueError(f"Unsupported grid_type={grid_type}.")
    if not l_ascending_grid:
        raise ValueError(
            "l_ascending_grid=False is not supported by the JAX grid operators."
        )

    if grid_type == GRID_TYPE_EVEN:
        nzm = int(
            np.floor(
                (zm_top_arr[0] - zm_init_arr[0] + deltaz_arr[0])
                / deltaz_arr[0]
            )
        )
        nzt = nzm - 1
        if nzm < 2 or nzt < 1:
            raise ValueError(
                f"Invalid derived grid dimensions: nzm={nzm}, nzt={nzt}."
            )
        zm, zt = _setup_even_grid(nzm, ngrdcol, deltaz_arr, zm_init_arr)
    elif grid_type == GRID_TYPE_STRETCHED_ZT:
        zt_in = _prepare_height_array(
            "thermodynamic_heights", thermodynamic_heights, ngrdcol
        )
        _validate_monotonic_increasing("thermodynamic_heights", zt_in)
        if np.any(zt_in[:, 0] <= zm_init_arr):
            raise ValueError(
                "Stretched zt grid lowest thermodynamic level must be above zm_init."
            )
        end_idx = int(np.searchsorted(zt_in[0], zm_top_arr[0], side="right") - 1)
        if end_idx < 0:
            raise ValueError("Stretched zt grid cannot fulfill zm_top requirement.")
        zt = zt_in[:, : end_idx + 1]
        nzt = zt.shape[1]
        nzm = nzt + 1
        zm = _calc_zm_from_zt(nzm, nzt, ngrdcol, zt, zm_init_arr)
    else:
        zm_in = _prepare_height_array("momentum_heights", momentum_heights, ngrdcol)
        _validate_monotonic_increasing("momentum_heights", zm_in)
        begin_idx = int(np.searchsorted(zm_in[0], zm_init_arr[0], side="left"))
        end_idx = int(np.searchsorted(zm_in[0], zm_top_arr[0], side="right") - 1)
        if begin_idx >= zm_in.shape[1]:
            raise ValueError("Stretched zm grid cannot fulfill zm_init requirement.")
        if end_idx < begin_idx:
            raise ValueError("Stretched zm grid cannot fulfill zm_top requirement.")
        zm = zm_in[:, begin_idx : end_idx + 1]
        nzm = zm.shape[1]
        nzt = nzm - 1
        if nzm < 2 or nzt < 1:
            raise ValueError(
                f"Stretched zm grid produced invalid dimensions: nzm={nzm}, nzt={nzt}."
            )
        zt = 0.5 * (zm[:, :-1] + zm[:, 1:])

    dzm, dzt, invrs_dzm, invrs_dzt = _calc_grid_spacings(
        nzm, ngrdcol, zm, zt
    )
    weights_zt2zm = _calc_zt2zm_weights(nzm, nzt, ngrdcol, zm, zt)
    weights_zm2zt = _calc_zm2zt_weights(nzt, ngrdcol, zm, zt, dzt)

    return Grid(
        nzm=nzm,
        nzt=nzt,
        ngrdcol=ngrdcol,
        zm=jnp.asarray(zm),
        zt=jnp.asarray(zt),
        dzm=jnp.asarray(dzm),
        dzt=jnp.asarray(dzt),
        invrs_dzm=jnp.asarray(invrs_dzm),
        invrs_dzt=jnp.asarray(invrs_dzt),
        weights_zt2zm=jnp.asarray(weights_zt2zm),
        weights_zm2zt=jnp.asarray(weights_zm2zt),
        k_lb_zm=0,
        k_ub_zm=nzm - 1,
        k_lb_zt=0,
        k_ub_zt=nzt - 1,
        grid_dir_indx=1,
        grid_dir=1.0,
    )


def _as_column_input(value, ngrdcol: int, name: str):
    if isinstance(value, (int, float)):
        return np.full(ngrdcol, value, dtype=core_rknd)
    arr = np.asarray(value, dtype=core_rknd)
    if arr.shape != (ngrdcol,):
        raise ValueError(f"{name} must have shape ({ngrdcol},), got {arr.shape}.")
    return arr


def _prepare_height_array(name: str, heights, ngrdcol: int):
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


def _validate_monotonic_increasing(name: str, arr) -> None:
    if arr.shape[1] < 1:
        raise ValueError(f"{name} must contain at least one level.")
    if np.any(np.diff(arr, axis=1) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing with level index.")


def _setup_even_grid(nzm: int, ngrdcol: int, deltaz, zm_init):
    k_indices = np.arange(nzm, dtype=core_rknd)
    zm = zm_init[:, None] + deltaz[:, None] * k_indices[None, :]
    zt = 0.5 * (zm[:, :-1] + zm[:, 1:])
    return zm, zt


def _calc_zm_from_zt(nzm: int, nzt: int, ngrdcol: int, zt, zm_init):
    zm = np.zeros((ngrdcol, nzm), dtype=core_rknd)
    zm[:, 1:-1] = 0.5 * (zt[:, :-1] + zt[:, 1:])
    zm[:, 0] = np.asarray(zm_init, dtype=core_rknd)
    if nzt > 1:
        zm[:, -1] = zt[:, -1] + 0.5 * (zt[:, -1] - zt[:, -2])
    else:
        zm[:, -1] = zt[:, -1] + (zt[:, -1] - zm[:, 0])
    return zm


def _calc_grid_spacings(nzm: int, ngrdcol: int, zm, zt):
    dzm = np.zeros((ngrdcol, nzm), dtype=core_rknd)
    dzt = zm[:, 1:] - zm[:, :-1]

    dzm[:, 1:-1] = zt[:, 1:] - zt[:, :-1]
    dzm[:, 0] = 2.0 * (zt[:, 0] - zm[:, 0])
    dzm[:, -1] = dzm[:, -2]

    invrs_dzm = np.where(np.abs(dzm) > 1e-30, 1.0 / dzm, 0.0)
    invrs_dzt = np.where(np.abs(dzt) > 1e-30, 1.0 / dzt, 0.0)
    return dzm, dzt, invrs_dzm, invrs_dzt


def _calc_zt2zm_weights(nzm: int, nzt: int, ngrdcol: int, zm, zt):
    weights = np.zeros((ngrdcol, nzm, 2), dtype=core_rknd)

    for k in range(1, nzm - 1):
        denom = zt[:, k] - zt[:, k - 1]
        weights[:, k, T_ABOVE] = (zm[:, k] - zt[:, k - 1]) / (denom + 1e-30)
        weights[:, k, T_BELOW] = (zt[:, k] - zm[:, k]) / (denom + 1e-30)

    if nzt >= 2:
        denom0 = zt[:, 1] - zt[:, 0]
        weights[:, 0, T_ABOVE] = (zm[:, 0] - zt[:, 0]) / (denom0 + 1e-30)
        weights[:, 0, T_BELOW] = (zt[:, 1] - zm[:, 0]) / (denom0 + 1e-30)

        denomn = zt[:, -1] - zt[:, -2]
        weights[:, -1, T_ABOVE] = (zm[:, -1] - zt[:, -2]) / (denomn + 1e-30)
        weights[:, -1, T_BELOW] = (zt[:, -1] - zm[:, -1]) / (denomn + 1e-30)
    else:
        weights[:, 0, T_ABOVE] = 1.0
        weights[:, -1, T_ABOVE] = 1.0

    return weights


def _calc_zm2zt_weights(nzt: int, ngrdcol: int, zm, zt, dzt):
    weights = np.zeros((ngrdcol, nzt, 2), dtype=core_rknd)
    for k in range(nzt):
        total_dist = dzt[:, k] + 1e-30
        weights[:, k, M_ABOVE] = (zt[:, k] - zm[:, k]) / total_dist
        weights[:, k, M_BELOW] = (zm[:, k + 1] - zt[:, k]) / total_dist
    return weights


def zt2zm(nzm: int, nzt: int, ngrdcol: int, gr, azt, zm_min=None):
    """Function to interpolate a variable located on the thermodynamic grid
    levels (azt) to the momentum grid levels (azm).  This function inputs the
    entire azt array and outputs the results as an azm array.  The
    formulation used is compatible with a stretched (unevenly-spaced) grid.
    """
    azt = jnp.asarray(azt, dtype=jnp.float64)

    # Interpolate the value of a thermodynamic-level variable to the central
    # momentum level, k, between two successive thermodynamic levels using
    # linear interpolation.
    interior = (
        gr.weights_zt2zm[:, 1:nzm - 1, T_ABOVE] * azt[:, 1:nzt]
        + gr.weights_zt2zm[:, 1:nzm - 1, T_BELOW] * azt[:, :nzt - 1]
    )

    # Set the value of the thermodynamic-level variable, azt, at momentum
    # level 1.  The name of the variable when interpolated/extended to momentum
    # levels is azm.  This is the lower boundary for an ascending grid and the
    # upper boundary for a descending grid.
    lower_ascending = azt[:, :1]
    upper_ascending = (
        gr.weights_zt2zm[:, nzm - 1:nzm, T_ABOVE] * azt[:, nzt - 1:nzt]
        + gr.weights_zt2zm[:, nzm - 1:nzm, T_BELOW] * azt[:, nzt - 2:nzt - 1]
    )
    # Use a linear extension based on the values of azt at levels 1 and 2 to
    # find the value of azm at level 1.
    lower_descending = (
        gr.weights_zt2zm[:, :1, T_ABOVE] * azt[:, 1:2]
        + gr.weights_zt2zm[:, :1, T_BELOW] * azt[:, :1]
    )
    upper_descending = azt[:, nzt - 1:nzt]

    is_ascending = gr.grid_dir_indx == 1
    lower = jnp.where(is_ascending, lower_ascending, lower_descending)
    upper = jnp.where(is_ascending, upper_ascending, upper_descending)

    azm = jnp.concatenate([lower, interior, upper], axis=1)
    if zm_min is not None:
        azm = jnp.maximum(azm, zm_min)
    return azm


def zm2zt(nzm: int, nzt: int, ngrdcol: int, gr, azm, zt_min=None):
    """Function to interpolate a variable located on the momentum grid levels
    (azm) to the thermodynamic grid levels (azt).  This function inputs the
    entire azm array and outputs the results as an azt array.  The formulation
    used is compatible with a stretched (unevenly-spaced) grid.
    """
    azm = jnp.asarray(azm, dtype=jnp.float64)
    # Interpolate the value of a momentum-level variable to the central
    # thermodynamic level, k, between two successive momentum levels using
    # linear interpolation.
    azt = (
        gr.weights_zm2zt[:, :, M_ABOVE] * azm[:, 1:nzm]
        + gr.weights_zm2zt[:, :, M_BELOW] * azm[:, :nzt]
    )
    if zt_min is not None:
        azt = jnp.maximum(azt, zt_min)
    return azt


def zt2zm2zt(nzm: int, nzt: int, ngrdcol: int, gr, azt, zt_min=None):
    """Function to interpolate a variable located on the thermodynamic grid
    levels (azt) to the momentum grid levels (azm), then interpolate back
    to thermodynamic grid levels (azt).

    Note:
      This is intended for smoothing variables.
    """
    # Interpolate azt to momentum levels
    # Interpolate back to thermodynamic levels
    return zm2zt(nzm, nzt, ngrdcol, gr, zt2zm(nzm, nzt, ngrdcol, gr, azt), zt_min)


def zm2zt2zm(nzm: int, nzt: int, ngrdcol: int, gr, azm, zm_min=None):
    """Function to interpolate a variable located on the momentum grid
    levels(azm) to thermodynamic grid levels (azt), then interpolate
    back to momentum grid levels (azm).

    Note:
      This is intended for smoothing variables.
    """
    # Interpolate azt to termodynamic levels
    # Interpolate back to momentum levels
    return zt2zm(nzm, nzt, ngrdcol, gr, zm2zt(nzm, nzt, ngrdcol, gr, azm), zm_min)


def ddzm(nzm: int, nzt: int, ngrdcol: int, gr, azm):
    """2D version of gradzm."""
    azm = jnp.asarray(azm, dtype=jnp.float64)
    # Vertical derivative of azm (thermo. levs.) [units vary / m]
    return (azm[:, 1:nzm] - azm[:, :nzt]) * gr.invrs_dzt


def ddzt(nzm: int, nzt: int, ngrdcol: int, gr, azt):
    """2D version of gradzt."""
    azt = jnp.asarray(azt, dtype=jnp.float64)
    # Vertical derivative of azt (mom.levs.) [units vary / m]
    interior = (azt[:, 1:nzt] - azt[:, :nzt - 1]) * gr.invrs_dzm[:, 1:nzm - 1]
    return jnp.concatenate([interior[:, :1], interior, interior[:, -1:]], axis=1)


__all__ = [
    "Grid",
    "setup_grid",
    "GRID_TYPE_EVEN",
    "GRID_TYPE_STRETCHED_ZT",
    "GRID_TYPE_STRETCHED_ZM",
    "zt2zm",
    "zm2zt",
    "zt2zm2zt",
    "zm2zt2zm",
    "ddzm",
    "ddzt",
]
