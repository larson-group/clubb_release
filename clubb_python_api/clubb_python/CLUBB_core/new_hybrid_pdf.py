"""User-facing wrappers for routines from CLUBB_core/new_hybrid_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def calculate_w_params(nz: int, ngrdcol: int, wm, wp2, skw, f_w, zeta_w):
    """Compute the setting-variable component means/stddevs and mixture fraction."""
    return clubb_f2py.f2py_calculate_w_params(
        f_arr(wm), f_arr(wp2), f_arr(skw), f_arr(f_w), f_arr(zeta_w),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calculate_responder_params(
    nz: int, ngrdcol: int, xm, xp2, skx, wpxp, wp2, f_w, mixt_frac,
):
    """Compute responder component means/variances and the variance coefficients."""
    return clubb_f2py.f2py_calculate_responder_params(
        f_arr(xm), f_arr(xp2), f_arr(skx), f_arr(wpxp), f_arr(wp2), f_arr(f_w), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calculate_coef_wp4_implicit(
    nz: int, ngrdcol: int, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
):
    """Compute coefficient such that <w'^4> = coef_wp4_implicit * <w'^2>^2."""
    return clubb_f2py.f2py_calculate_coef_wp4_implicit(
        f_arr(mixt_frac), f_arr(f_w), f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calc_coef_wp2xp_implicit(
    nz: int, ngrdcol: int, wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
):
    """Compute coefficient such that <w'^2 x'> = coef_wp2xp_implicit * <w'x'>."""
    return clubb_f2py.f2py_calc_coef_wp2xp_implicit(
        f_arr(wp2), f_arr(mixt_frac), f_arr(f_w), f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calc_coefs_wpxp2_semiimpl(
    nz: int, ngrdcol: int, wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
):
    """Compute (<coef_wpxp2_implicit>, <term_wpxp2_explicit>) for semi-implicit transport."""
    return clubb_f2py.f2py_calc_coefs_wpxp2_semiimpl(
        f_arr(wp2), f_arr(wpxp), f_arr(mixt_frac), f_arr(f_w), f_arr(coef_sigma_x_1_sqd), f_arr(coef_sigma_x_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )
