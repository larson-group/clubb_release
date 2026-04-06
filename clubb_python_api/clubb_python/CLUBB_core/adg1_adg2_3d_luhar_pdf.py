"""User-facing wrappers for routines from CLUBB_core/adg1_adg2_3d_luhar_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def calc_luhar_params(nz: int, ngrdcol: int, skx, wpxp, xp2, x_tol_sqd: float):
    """Compute Luhar closure mixture fraction and shape parameters (M, m)."""
    return clubb_f2py.f2py_calc_luhar_params(
        f_arr(skx), f_arr(wpxp), f_arr(xp2), float(x_tol_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol))


def close_luhar_pdf(nz: int, ngrdcol: int, xm, xp2, mixt_frac, small_m, wpxp, x_tol_sqd: float):
    """Compute closed Luhar component moments and variances from (a, m)."""
    return clubb_f2py.f2py_close_luhar_pdf(
        f_arr(xm), f_arr(xp2), f_arr(mixt_frac), f_arr(small_m), f_arr(wpxp), float(x_tol_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol))


def adg1_w_closure(
    nz: int, ngrdcol: int, wm, wp2, skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag: float
):
    """Compute ADG1 w-mixture fraction and component means/variances."""
    return clubb_f2py.f2py_adg1_w_closure(
        f_arr(wm), f_arr(wp2), f_arr(skw), f_arr(sigma_sqd_w), f_arr(sqrt_wp2), float(mixt_frac_max_mag),
        nz=int(nz), ngrdcol=int(ngrdcol))


def adg2_pdf_driver(
    nz: int, ngrdcol: int, sclr_dim: int,
    sclr_tol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, sqrt_wp2,
    beta, sclrm, sclrp2, wpsclrp, l_scalar_calc: bool,
):
    """Compute ADG2 PDF component means/variances for w/rt/thl/(sclr)."""
    return clubb_f2py.f2py_adg2_pdf_driver(
        int(sclr_dim), f_arr(sclr_tol), f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2),
        f_arr(thlp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp), f_arr(sqrt_wp2), f_arr(beta), f_arr(sclrm),
        f_arr(sclrp2), f_arr(wpsclrp), l_scalar_calc, nz=int(nz), ngrdcol=int(ngrdcol),
    )


def adg1_pdf_driver(
    nz: int, ngrdcol: int, sclr_dim: int,
    sclr_tol, wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp,
    upwp, vpwp, sqrt_wp2, sigma_sqd_w, beta, mixt_frac_max_mag: float,
    sclrm, sclrp2, wpsclrp, l_scalar_calc: bool,
):
    """Compute ADG1 PDF component means/variances for w/rt/thl/u/v/(sclr)."""
    return clubb_f2py.f2py_adg1_pdf_driver(
        int(sclr_dim), f_arr(sclr_tol), f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(um), f_arr(vm), f_arr(wp2),
        f_arr(rtp2), f_arr(thlp2), f_arr(up2), f_arr(vp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp),
        f_arr(upwp), f_arr(vpwp), f_arr(sqrt_wp2), f_arr(sigma_sqd_w), f_arr(beta), float(mixt_frac_max_mag),
        f_arr(sclrm), f_arr(sclrp2), f_arr(wpsclrp), l_scalar_calc, nz=int(nz), ngrdcol=int(ngrdcol),
    )


def luhar_3d_pdf_driver(
    nz: int, ngrdcol: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp
):
    """Compute 3D Luhar PDF component means/variances for w/rt/thl."""
    return clubb_f2py.f2py_luhar_3d_pdf_driver(
        f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2),
        f_arr(skw), f_arr(skrt), f_arr(skthl), f_arr(wprtp), f_arr(wpthlp),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )
