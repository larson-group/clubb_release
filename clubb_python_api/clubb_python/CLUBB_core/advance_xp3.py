"""User-facing wrappers for routines from CLUBB_core/advance_xp3_module.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def advance_xp3(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, dt: float,
    rtm, thlm, rtp2, thlp2, wprtp, wpthlp, wprtp2, wpthlp2,
    rho_ds_zm, invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,
    sclrm, sclrp2, wpsclrp, wpsclrp2,
    wp2, wp3, upwp, vpwp, up2, vp2, thvm, sigma_sqd_w, clubb_params,
    l_lmm_stepping: bool,
    rtp3, thlp3, sclrp3, up3, vp3,
):
    """Advance third-order moments one model timestep."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_advance_xp3(
        int(sclr_dim), f_arr(sclr_tol), float(dt),
        f_arr(rtm), f_arr(thlm), f_arr(rtp2), f_arr(thlp2), f_arr(wprtp), f_arr(wpthlp),
        f_arr(wprtp2), f_arr(wpthlp2), f_arr(rho_ds_zm), f_arr(invrs_rho_ds_zt),
        f_arr(invrs_tau_zt), f_arr(tau_max_zt),
        f_arr(sclrm), f_arr(sclrp2), f_arr(wpsclrp), f_arr(wpsclrp2),
        f_arr(wp2), f_arr(wp3), f_arr(upwp), f_arr(vpwp),
        f_arr(up2), f_arr(vp2), f_arr(thvm), f_arr(sigma_sqd_w),
        f_arr(clubb_params),
        l_lmm_stepping,
        f_arr(rtp3), f_arr(thlp3), f_arr(sclrp3), f_arr(up3), f_arr(vp3),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )


def compute_xp3(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol,
    iipdf_type: int, clubb_params,
    wp2, wp3, thvm,
    wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2,
    sigma_sqd_w, wpsclrp, sclrp2,
    rtp3, thlp3, up3, vp3, sclrp3,
):
    """Diagnose third-order moments without the predictive xp3 solve."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_compute_xp3(
        int(sclr_dim), f_arr(sclr_tol), int(iipdf_type), f_arr(clubb_params),
        f_arr(wp2), f_arr(wp3), f_arr(thvm),
        f_arr(wprtp), f_arr(wpthlp), f_arr(rtp2), f_arr(thlp2),
        f_arr(upwp), f_arr(vpwp), f_arr(up2), f_arr(vp2),
        f_arr(sigma_sqd_w), f_arr(wpsclrp), f_arr(sclrp2),
        f_arr(rtp3), f_arr(thlp3), f_arr(up3), f_arr(vp3), f_arr(sclrp3),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
