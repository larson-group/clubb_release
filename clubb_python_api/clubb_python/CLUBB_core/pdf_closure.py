"""User-facing wrappers for routines from CLUBB_core pdf_closure_module."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import (
    get_fortran_implicit_coefs,
    get_fortran_pdf_params,
    get_fortran_pdf_params_zm,
    set_fortran_implicit_coefs,
)
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params_converter import (
    set_fortran_pdf_params,
    set_fortran_pdf_params_zm,
)
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def _transport_3d(arr, ncol: int, nz: int, logical_dim: int, dtype=np.float64):
    """Pad 3D transported arrays to at least size 1 in the trailing dim."""
    transport_dim = max(int(logical_dim), 1)
    out = np.zeros((int(ncol), int(nz), transport_dim), dtype=dtype, order="F")
    if int(logical_dim) > 0:
        out[:, :, : int(logical_dim)] = f_arr(arr)
    return out


def _transport_1d(arr, logical_dim: int, dtype=np.float64):
    """Pad 1D transported arrays to at least size 1."""
    transport_dim = max(int(logical_dim), 1)
    out = np.zeros((transport_dim,), dtype=dtype, order="F")
    if int(logical_dim) > 0:
        out[: int(logical_dim)] = f_arr(arr)
    return out

def pdf_closure_driver(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, dt: float, hydromet_dim: int, sclr_dim: int,
    sclr_tol, wprtp, thlm, wpthlp, rtp2, rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, wm_zm, wm_zt,
    um, up2, upwp, up3, vm, vp2, vpwp, vp3, p_in_pa, exner, thv_ds_zm, thv_ds_zt, rtm_ref,
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt,
    sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, sclrp3, p_sfc,
    l_samp_stats_in_pdf_call: bool, mixt_frac_max_mag: float, ts_nudge: float,
    rtm_min: float, rtm_nudge_max_altitude: float, clubb_params,
    iiPDF_type: int, saturation_formula: int, l_predict_upwp_vpwp: bool, l_rtm_nudge: bool,
    l_trapezoidal_rule_zt: bool, l_trapezoidal_rule_zm: bool, l_call_pdf_closure_twice: bool,
    l_use_cloud_cover: bool, l_rcm_supersat_adj: bool, l_mix_rat_hm,
    pdf_params: pdf_parameter,
    pdf_params_zm: pdf_parameter,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    err_info: ErrInfo,
    rtm,
):
    """Run pdf_closure_driver with strict direct argument mapping to Fortran."""
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    set_fortran_pdf_params_zm(pdf_params_zm)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_pdf_closure_driver(
        float(dt), int(hydromet_dim), int(sclr_dim), _transport_1d(sclr_tol, sclr_dim),
        f_arr(wprtp), f_arr(thlm), f_arr(wpthlp),
        f_arr(rtp2), f_arr(rtp3), f_arr(thlp2), f_arr(thlp3),
        f_arr(rtpthlp), f_arr(wp2), f_arr(wp3),
        f_arr(wm_zm), f_arr(wm_zt),
        f_arr(um), f_arr(up2), f_arr(upwp), f_arr(up3),
        f_arr(vm), f_arr(vp2), f_arr(vpwp), f_arr(vp3),
        f_arr(p_in_pa), f_arr(exner), f_arr(thv_ds_zm),
        f_arr(thv_ds_zt), f_arr(rtm_ref),
        _transport_3d(wphydrometp, ngrdcol, nzm, hydromet_dim),
        _transport_3d(wp2hmp, ngrdcol, nzt, hydromet_dim),
        _transport_3d(rtphmp_zt, ngrdcol, nzt, hydromet_dim),
        _transport_3d(thlphmp_zt, ngrdcol, nzt, hydromet_dim),
        _transport_3d(sclrm, ngrdcol, nzt, sclr_dim),
        _transport_3d(wpsclrp, ngrdcol, nzm, sclr_dim),
        _transport_3d(sclrp2, ngrdcol, nzm, sclr_dim),
        _transport_3d(sclrprtp, ngrdcol, nzm, sclr_dim),
        _transport_3d(sclrpthlp, ngrdcol, nzm, sclr_dim),
        _transport_3d(sclrp3, ngrdcol, nzt, sclr_dim),
        f_arr(p_sfc), bool(l_samp_stats_in_pdf_call), float(mixt_frac_max_mag), float(ts_nudge),
        float(rtm_min), float(rtm_nudge_max_altitude), f_arr(clubb_params), int(iiPDF_type),
        int(saturation_formula), bool(l_predict_upwp_vpwp), bool(l_rtm_nudge),
        bool(l_trapezoidal_rule_zt), bool(l_trapezoidal_rule_zm), bool(l_call_pdf_closure_twice),
        bool(l_use_cloud_cover), bool(l_rcm_supersat_adj),
        _transport_1d(l_mix_rat_hm, hydromet_dim, dtype=np.bool_), f_arr(rtm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return (
        result[0],
        get_fortran_implicit_coefs(),
        get_fortran_pdf_params(),
        get_fortran_pdf_params_zm(),
        get_fortran_err_info(),
        *result[1:],
    )


def calc_wp4_pdf(nz: int, ngrdcol: int, wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac):
    """Compute <w'^4> from double-Gaussian PDF component moments."""
    return clubb_f2py.f2py_calc_wp4_pdf(
        f_arr(wm), f_arr(w_1), f_arr(w_2),
        f_arr(varnce_w_1), f_arr(varnce_w_2), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_wp2xp_pdf(nz: int, ngrdcol: int, wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2,
                   varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac):
    """Compute <w'^2 x'> from bivariate PDF component moments."""
    return clubb_f2py.f2py_calc_wp2xp_pdf(
        f_arr(wm), f_arr(xm), f_arr(w_1), f_arr(w_2),
        f_arr(x_1), f_arr(x_2), f_arr(varnce_w_1),
        f_arr(varnce_w_2), f_arr(varnce_x_1), f_arr(varnce_x_2),
        f_arr(corr_w_x_1), f_arr(corr_w_x_2), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_wpxp2_pdf(nz: int, ngrdcol: int, wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, varnce_w_2,
                   varnce_x_1, varnce_x_2, corr_w_x_1, corr_w_x_2, mixt_frac):
    """Compute <w' x'^2> from bivariate PDF component moments."""
    return clubb_f2py.f2py_calc_wpxp2_pdf(
        f_arr(wm), f_arr(xm), f_arr(w_1), f_arr(w_2),
        f_arr(x_1), f_arr(x_2), f_arr(varnce_w_1),
        f_arr(varnce_w_2), f_arr(varnce_x_1), f_arr(varnce_x_2),
        f_arr(corr_w_x_1), f_arr(corr_w_x_2), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_wpxpyp_pdf(nz: int, ngrdcol: int, wm, xm, ym, w_1, w_2, x_1, x_2, y_1, y_2, varnce_w_1, varnce_w_2,
                    varnce_x_1, varnce_x_2, varnce_y_1, varnce_y_2, corr_w_x_1, corr_w_x_2,
                    corr_w_y_1, corr_w_y_2, corr_x_y_1, corr_x_y_2, mixt_frac):
    """Compute <w' x' y'> from trivariate PDF component moments."""
    return clubb_f2py.f2py_calc_wpxpyp_pdf(
        f_arr(wm), f_arr(xm), f_arr(ym), f_arr(w_1), f_arr(w_2),
        f_arr(x_1), f_arr(x_2), f_arr(y_1), f_arr(y_2),
        f_arr(varnce_w_1), f_arr(varnce_w_2),
        f_arr(varnce_x_1), f_arr(varnce_x_2),
        f_arr(varnce_y_1), f_arr(varnce_y_2),
        f_arr(corr_w_x_1), f_arr(corr_w_x_2),
        f_arr(corr_w_y_1), f_arr(corr_w_y_2),
        f_arr(corr_x_y_1), f_arr(corr_x_y_2),
        f_arr(mixt_frac), nz=int(nz), ngrdcol=int(ngrdcol))


def calc_w_up_in_cloud(
    nz: int, ngrdcol: int, mixt_frac, cloud_frac_1, cloud_frac_2, w_1, w_2, varnce_w_1, varnce_w_2
):
    """Compute cloudy updraft/downdraft velocity and fractions."""
    return clubb_f2py.f2py_calc_w_up_in_cloud(
        f_arr(mixt_frac), f_arr(cloud_frac_1), f_arr(cloud_frac_2),
        f_arr(w_1), f_arr(w_2), f_arr(varnce_w_1), f_arr(varnce_w_2),
        nz=int(nz), ngrdcol=int(ngrdcol))
