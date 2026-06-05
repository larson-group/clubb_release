"""User-facing wrappers for routines from CLUBB_core/sponge_layer_damping.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def _fa_2d(arr):
    out = f_arr(arr, dtype=np.float64)
    if out.ndim != 2:
        raise ValueError(f"Expected a rank-2 array, got rank-{out.ndim}.")
    return out


def _stack_column_outputs(column_outputs):
    return f_arr(np.stack(column_outputs, axis=0))


def sponge_damp_xm(
    gr: Grid, nzm: int, nzt: int,
    dt: float, zt, zm, xm_ref, xm, damping_profile=None,
):
    """Damp a mean-field profile toward a reference profile in sponge layer."""
    ngrdcol = np.asarray(zt).shape[0]
    ngrdcol = int(ngrdcol)
    tau_sponge_damp = getattr(damping_profile, "tau_sponge_damp", None)
    sponge_layer_depth = getattr(damping_profile, "sponge_layer_depth", None)
    if tau_sponge_damp is None or sponge_layer_depth is None:
        raise ValueError("sponge_damp_xm requires damping_profile.")
    set_fortran_grid(gr)
    zt_arr = _fa_2d(zt)
    zm_arr = _fa_2d(zm)
    xm_ref_arr = _fa_2d(xm_ref)
    xm_arr = _fa_2d(xm)
    tau_arr = _fa_2d(tau_sponge_damp)
    depth_arr = np.asarray(sponge_layer_depth, dtype=np.float64).reshape(ngrdcol)
    return _stack_column_outputs([
        clubb_f2py.f2py_sponge_damp_xm(
            float(dt), zt_arr[i, :], zm_arr[i, :], xm_ref_arr[i, :], xm_arr[i, :],
            tau_arr[i, :], float(depth_arr[i]), nzm=int(nzm), nzt=int(nzt)
        )
        for i in range(ngrdcol)
    ])


def sponge_damp_xp2(
    gr: Grid, nzm: int,
    dt: float, zm, xp2, x_tol_sqd, damping_profile=None,
):
    """Damp variance profile in sponge layer with a lower floor."""
    ngrdcol = np.asarray(zm).shape[0]
    ngrdcol = int(ngrdcol)
    tau_sponge_damp = getattr(damping_profile, "tau_sponge_damp", None)
    sponge_layer_depth = getattr(damping_profile, "sponge_layer_depth", None)
    if tau_sponge_damp is None or sponge_layer_depth is None:
        raise ValueError("sponge_damp_xp2 requires damping_profile.")
    set_fortran_grid(gr)
    zm_arr = _fa_2d(zm)
    xp2_arr = _fa_2d(xp2)
    tau_arr = _fa_2d(tau_sponge_damp)
    x_tol_arr = np.asarray(x_tol_sqd, dtype=np.float64).reshape(ngrdcol)
    depth_arr = np.asarray(sponge_layer_depth, dtype=np.float64).reshape(ngrdcol)
    return _stack_column_outputs([
        clubb_f2py.f2py_sponge_damp_xp2(
            float(dt), zm_arr[i, :], xp2_arr[i, :], float(x_tol_arr[i]),
            tau_arr[i, :], float(depth_arr[i]), nzm=int(nzm)
        )
        for i in range(ngrdcol)
    ])


def sponge_damp_xp3(
    gr: Grid, nzm: int, nzt: int,
    dt: float, z, zm, xp3, damping_profile=None,
):
    """Damp third-moment profile in sponge layer."""
    ngrdcol = np.asarray(z).shape[0]
    ngrdcol = int(ngrdcol)
    tau_sponge_damp = getattr(damping_profile, "tau_sponge_damp", None)
    sponge_layer_depth = getattr(damping_profile, "sponge_layer_depth", None)
    if tau_sponge_damp is None or sponge_layer_depth is None:
        raise ValueError("sponge_damp_xp3 requires damping_profile.")
    set_fortran_grid(gr)
    z_arr = _fa_2d(z)
    zm_arr = _fa_2d(zm)
    xp3_arr = _fa_2d(xp3)
    tau_arr = _fa_2d(tau_sponge_damp)
    depth_arr = np.asarray(sponge_layer_depth, dtype=np.float64).reshape(ngrdcol)
    return _stack_column_outputs([
        clubb_f2py.f2py_sponge_damp_xp3(
            float(dt), z_arr[i, :], zm_arr[i, :], xp3_arr[i, :],
            tau_arr[i, :], float(depth_arr[i]), nzm=int(nzm), nzt=int(nzt)
        )
        for i in range(ngrdcol)
    ])
