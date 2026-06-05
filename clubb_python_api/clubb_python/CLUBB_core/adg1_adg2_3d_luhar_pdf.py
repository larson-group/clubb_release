"""User-facing wrappers for routines from CLUBB_core/adg1_adg2_3d_luhar_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import set_fortran_err_info
from clubb_python.CLUBB_core.err_info_type_module import init_err_info


def _fa_2d(arr):
    out = f_arr(arr, dtype=np.float64)
    if out.ndim != 2:
        raise ValueError(f"Expected a rank-2 array, got rank-{out.ndim}.")
    return out


def _stack_tuple_outputs(column_outputs):
    return tuple(f_arr(np.stack(items, axis=0)) for items in zip(*column_outputs, strict=True))


def _ensure_err_info(ngrdcol: int, err_info):
    if err_info is None:
        return init_err_info(int(ngrdcol))
    if isinstance(err_info, ErrInfo):
        set_fortran_err_info(err_info)
    return err_info


def calc_luhar_params(nz: int, skx, wpxp, xp2, x_tol_sqd: float):
    """Compute Luhar closure mixture fraction and shape parameters (M, m)."""
    arrays = [_fa_2d(arg) for arg in (skx, wpxp, xp2)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_luhar_params(*(arr[i, :] for arr in arrays), float(x_tol_sqd), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def close_luhar_pdf(nz: int, xm, xp2, mixt_frac, small_m, wpxp, x_tol_sqd: float):
    """Compute closed Luhar component moments and variances from (a, m)."""
    arrays = [_fa_2d(arg) for arg in (xm, xp2, mixt_frac, small_m, wpxp)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_close_luhar_pdf(*(arr[i, :] for arr in arrays), float(x_tol_sqd), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def adg1_w_closure(
    nz: int, ngrdcol: int, wm, wp2, skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag: float,
    err_info=None,
):
    """Compute ADG1 w-mixture fraction and component means/variances."""
    _ensure_err_info(ngrdcol, err_info)
    return clubb_f2py.f2py_adg1_w_closure(
        f_arr(wm), f_arr(wp2), f_arr(skw), f_arr(sigma_sqd_w), f_arr(sqrt_wp2), float(mixt_frac_max_mag),
        nz=int(nz), ngrdcol=int(ngrdcol))


def adg2_pdf_driver(
    nz: int, ngrdcol: int, sclr_dim: int,
    sclr_tol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, sqrt_wp2,
    beta, sclrm, sclrp2, wpsclrp, l_scalar_calc: bool, err_info=None,
):
    """Compute ADG2 PDF component means/variances for w/rt/thl/(sclr)."""
    _ensure_err_info(ngrdcol, err_info)
    return clubb_f2py.f2py_adg2_pdf_driver(
        int(sclr_dim), f_arr(sclr_tol), f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2),
        f_arr(thlp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp), f_arr(sqrt_wp2), f_arr(beta), f_arr(sclrm),
        f_arr(sclrp2), f_arr(wpsclrp), l_scalar_calc, nz=int(nz), ngrdcol=int(ngrdcol),
    )


def adg1_pdf_driver(
    nz: int, ngrdcol: int, sclr_dim: int,
    sclr_tol, wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp,
    upwp, vpwp, sqrt_wp2, sigma_sqd_w, beta, mixt_frac_max_mag: float,
    sclrm, sclrp2, wpsclrp, l_scalar_calc: bool, err_info=None,
):
    """Compute ADG1 PDF component means/variances for w/rt/thl/u/v/(sclr)."""
    _ensure_err_info(ngrdcol, err_info)
    return clubb_f2py.f2py_adg1_pdf_driver(
        int(sclr_dim), f_arr(sclr_tol), f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(um), f_arr(vm), f_arr(wp2),
        f_arr(rtp2), f_arr(thlp2), f_arr(up2), f_arr(vp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp),
        f_arr(upwp), f_arr(vpwp), f_arr(sqrt_wp2), f_arr(sigma_sqd_w), f_arr(beta), float(mixt_frac_max_mag),
        f_arr(sclrm), f_arr(sclrp2), f_arr(wpsclrp), l_scalar_calc, nz=int(nz), ngrdcol=int(ngrdcol),
    )


def luhar_3d_pdf_driver(nz: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp):
    """Compute 3D Luhar PDF component means/variances for w/rt/thl."""
    arrays = [_fa_2d(arg) for arg in (wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_luhar_3d_pdf_driver(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)
