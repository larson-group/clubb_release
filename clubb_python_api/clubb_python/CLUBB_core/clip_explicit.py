"""User-facing wrappers for routines from CLUBB_core/clip_explicit.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def clip_rcm(nzt: int, ngrdcol: int, rtm, message: str, rcm):
    """Clip cloud water mixing ratio so rcm <= rtm and rcm >= 0."""
    return clubb_f2py.f2py_clip_rcm(f_arr(rtm), message, f_arr(rcm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def clip_covar(
    nzm: int, ngrdcol: int,
    solve_type: int,
    l_first_clip_ts: bool, l_last_clip_ts: bool,
    dt: float, xp2, yp2,
    l_predict_upwp_vpwp: bool, xpyp,
):
    """Clip covariance x'y' to satisfy Cauchy-Schwarz."""
    return clubb_f2py.f2py_clip_covar(
        solve_type, l_first_clip_ts, l_last_clip_ts,
        dt, f_arr(xp2), f_arr(yp2), l_predict_upwp_vpwp, f_arr(xpyp),
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def clip_variance(
    gr: Grid, nzm: int, ngrdcol: int, solve_type: int, dt: float,
    threshold_lo, xp2,
):
    """Clip variance x'^2 to lower threshold."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_clip_variance(
        solve_type, dt, f_arr(threshold_lo), f_arr(xp2), nzm=int(nzm), ngrdcol=int(ngrdcol))


def clip_covars_denom(
    nzm: int, ngrdcol: int, sclr_dim: int, dt: float,
    rtp2, thlp2, up2, vp2, wp2, sclrp2,
    wprtp_cl_num: int, wpthlp_cl_num: int, wpsclrp_cl_num: int,
    upwp_cl_num: int, vpwp_cl_num: int,
    l_predict_upwp_vpwp: bool, l_tke_aniso: bool,
    l_linearize_pbl_winds: bool,
    wprtp, wpthlp, upwp, vpwp, wpsclrp, upwp_pert, vpwp_pert,
):
    """Clip all covariances using denominator (variance) bounds."""
    return clubb_f2py.f2py_clip_covars_denom(
        int(sclr_dim), float(dt),
        f_arr(rtp2), f_arr(thlp2), f_arr(up2), f_arr(vp2), f_arr(wp2), f_arr(sclrp2),
        wprtp_cl_num, wpthlp_cl_num, wpsclrp_cl_num,
        upwp_cl_num, vpwp_cl_num,
        l_predict_upwp_vpwp, l_tke_aniso, l_linearize_pbl_winds,
        f_arr(wprtp), f_arr(wpthlp), f_arr(upwp), f_arr(vpwp),
        f_arr(wpsclrp), f_arr(upwp_pert), f_arr(vpwp_pert),
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def clip_skewness(
    gr: Grid, nzt: int, ngrdcol: int, dt: float,
    sfc_elevation, skw_max_mag, wp2_zt,
    l_use_wp3_lim_with_smth_heaviside: bool, wp3,
):
    """Clip wp3 using the skewness-based limiter from clip_explicit.F90."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_clip_skewness(
        float(dt), f_arr(sfc_elevation), f_arr(skw_max_mag), f_arr(wp2_zt),
        l_use_wp3_lim_with_smth_heaviside, f_arr(wp3),
        nzt=int(nzt), ngrdcol=int(ngrdcol))
