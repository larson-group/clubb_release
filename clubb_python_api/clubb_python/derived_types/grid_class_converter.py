"""Push/pull converters for Grid <-> Fortran module storage."""

import numpy as np

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid


def set_fortran_grid(gr: Grid):
    """Push a Python Grid (0-based bounds) into Fortran module storage (1-based)."""
    clubb_f2py.set_grid(
        np.asfortranarray(gr.zm), np.asfortranarray(gr.zt),
        np.asfortranarray(gr.dzm), np.asfortranarray(gr.dzt),
        np.asfortranarray(gr.invrs_dzm), np.asfortranarray(gr.invrs_dzt),
        np.asfortranarray(gr.weights_zt2zm), np.asfortranarray(gr.weights_zm2zt),
        # Grid bounds are stored as Python 0-based indices in Grid, but
        # the Fortran grid type expects 1-based indices.
        int(gr.k_lb_zm) + 1,
        int(gr.k_ub_zm) + 1,
        int(gr.k_lb_zt) + 1,
        int(gr.k_ub_zt) + 1,
        gr.grid_dir_indx, gr.grid_dir,
    )


def get_fortran_grid() -> Grid:
    """Pull the current Fortran grid data from module storage into a Grid."""
    ngrdcol, nzm, nzt = clubb_f2py.f2py_get_grid_shape()
    ngrdcol = int(ngrdcol)
    nzm = int(nzm)
    nzt = int(nzt)
    if ngrdcol <= 0 or nzm <= 0 or nzt <= 0:
        raise RuntimeError(
            "Fortran grid data is not initialized in module storage; call "
            "setup_grid or set_fortran_grid first."
        )

    (
        zm,
        zt,
        dzm,
        dzt,
        invrs_dzm,
        invrs_dzt,
        weights_zt2zm,
        weights_zm2zt,
        k_lb_zm,
        k_ub_zm,
        k_lb_zt,
        k_ub_zt,
        grid_dir_indx,
        grid_dir,
    ) = clubb_f2py.get_grid(ngrdcol, nzm, nzt)
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
        # Fortran returns 1-based bounds; convert back to Python 0-based.
        k_lb_zm=int(k_lb_zm) - 1,
        k_ub_zm=int(k_ub_zm) - 1,
        k_lb_zt=int(k_lb_zt) - 1,
        k_ub_zt=int(k_ub_zt) - 1,
        grid_dir_indx=int(grid_dir_indx),
        grid_dir=float(grid_dir),
    )
