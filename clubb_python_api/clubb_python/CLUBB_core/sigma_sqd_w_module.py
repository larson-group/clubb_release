"""User-facing wrappers for routines from CLUBB_core/sigma_sqd_w_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def compute_sigma_sqd_w(
    nzm: int, ngrdcol: int,
    gamma_Skw_fnc, wp2, thlp2, rtp2, up2, vp2,
    wpthlp, wprtp, upwp, vpwp,
    l_predict_upwp_vpwp: bool,
):
    """Compute the PDF width parameter sigma_sqd_w."""
    return clubb_f2py.f2py_compute_sigma_sqd_w(
        f_arr(gamma_Skw_fnc), f_arr(wp2), f_arr(thlp2), f_arr(rtp2),
        f_arr(up2), f_arr(vp2), f_arr(wpthlp), f_arr(wprtp),
        f_arr(upwp), f_arr(vpwp), bool(l_predict_upwp_vpwp),
        nzm=int(nzm), ngrdcol=int(ngrdcol))
