"""User-facing wrappers for routines from CLUBB_core/sfc_varnce_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def calc_sfc_varnce(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, dt: float,
    sfc_elevation, upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc, um, vm,
    lscale_up, wpsclrp_sfc, lhs_splat_wp2, tau_zm, l_vary_convect_depth: bool, t0: float,
    up2_sfc_coef, a_const, wp2, up2, vp2, thlp2, rtp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp,
    sclr_idx: SclrIdx,
    err_info: ErrInfo,
):
    """Compute/update surface variances using grid/sclr_idx/stats/err_info data in module storage."""
    set_fortran_grid(gr)
    set_fortran_sclr_idx(sclr_idx)
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_calc_sfc_varnce(
        int(sclr_dim), float(dt), f_arr(sfc_elevation), f_arr(upwp_sfc), f_arr(vpwp_sfc), f_arr(wpthlp), f_arr(wprtp_sfc),
        f_arr(um), f_arr(vm), f_arr(lscale_up), f_arr(wpsclrp_sfc), f_arr(lhs_splat_wp2), f_arr(tau_zm),
        l_vary_convect_depth, float(t0), f_arr(up2_sfc_coef), f_arr(a_const), f_arr(wp2), f_arr(up2), f_arr(vp2),
        f_arr(thlp2), f_arr(rtp2), f_arr(rtpthlp), f_arr(sclrp2), f_arr(sclrprtp), f_arr(sclrpthlp),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()
