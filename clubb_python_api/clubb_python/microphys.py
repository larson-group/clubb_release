"""User-facing wrappers for microphysics routines (src/Microphys)."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def cloud_drop_sed(nzm: int, nzt: int, ngrdcol: int, gr: Grid, rcm, Ncm, rho_zm, rho, exner,
                   sigma_g: float, rcm_mc, thlm_mc):
    """Add cloud droplet sedimentation to the rcm and thlm microphysics tendencies."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_cloud_drop_sed(
        f_arr(rcm), f_arr(Ncm), f_arr(rho_zm), f_arr(rho), f_arr(exner), float(sigma_g),
        f_arr(rcm_mc), f_arr(thlm_mc), nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))
