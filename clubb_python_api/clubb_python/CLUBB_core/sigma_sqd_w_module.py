"""User-facing wrappers for routines from CLUBB_core/sigma_sqd_w_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def compute_sigma_sqd_w(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int,
    wp3, wp2, thlp2, rtp2, up2, vp2,
    wpthlp, wprtp, upwp, vpwp, clubb_params,
    l_predict_upwp_vpwp: bool,
):
    """Compute the PDF width parameter sigma_sqd_w."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_compute_sigma_sqd_w(
        f_arr(wp3), f_arr(wp2), f_arr(thlp2), f_arr(rtp2),
        f_arr(up2), f_arr(vp2), f_arr(wpthlp), f_arr(wprtp),
        f_arr(upwp), f_arr(vpwp), f_arr(clubb_params), bool(l_predict_upwp_vpwp),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))
