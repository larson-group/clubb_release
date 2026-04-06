"""User-facing wrappers for routines from CLUBB_core/diffusion.F90."""

import numpy as np
from numpy import asfortranarray as f_arr
import warnings

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def _diffusion_zm_lhs_legacy_call(k_zm, k_zt, nu, invrs_rho_ds_zm, rho_ds_zt, nzm, nzt, ngrdcol):
    """Compatibility call path for pre-reorder compiled wrappers."""
    return clubb_f2py.f2py_diffusion_zm_lhs(
        f_arr(k_zm), f_arr(k_zt), f_arr(nu),
        f_arr(invrs_rho_ds_zm), f_arr(rho_ds_zt),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def diffusion_zt_lhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, k_zm, k_zt, nu, invrs_rho_ds_zt, rho_ds_zm
):
    """Eddy-diffusion lhs contribution for zt-grid variables."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_diffusion_zt_lhs(
        f_arr(k_zm), f_arr(k_zt), f_arr(nu),
        f_arr(invrs_rho_ds_zt), f_arr(rho_ds_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def diffusion_zm_lhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, k_zt, k_zm, nu, invrs_rho_ds_zm, rho_ds_zt
):
    """Eddy-diffusion lhs contribution for zm-grid variables."""
    set_fortran_grid(gr)
    try:
        return clubb_f2py.f2py_diffusion_zm_lhs(
            f_arr(k_zt), f_arr(k_zm), f_arr(nu),
            f_arr(invrs_rho_ds_zm), f_arr(rho_ds_zt),
            nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))
    except Exception as exc:
        # Compatibility path for older compiled wrappers that still use (k_zm, k_zt).
        msg = str(exc)
        if "shape(k_zm, 1)" not in msg:
            raise
        warnings.warn(
            "Using legacy diffusion_zm_lhs argument order (k_zm, k_zt); rebuild "
            "clubb_f2py to use canonical source order (k_zt, k_zm).",
            DeprecationWarning,
            stacklevel=2,
        )
        return _diffusion_zm_lhs_legacy_call(
            k_zm, k_zt, nu, invrs_rho_ds_zm, rho_ds_zt, nzm, nzt, ngrdcol)
