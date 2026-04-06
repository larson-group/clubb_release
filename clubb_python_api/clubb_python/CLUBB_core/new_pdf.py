"""User-facing wrappers for routines from CLUBB_core/new_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def calc_setter_var_params(nz: int, ngrdcol: int, xm, xp2, skx, sgn_wpxp, f_x, zeta_x):
    """Compute new-PDF setter component means/stddevs, mixture fraction, and coefficients."""
    return clubb_f2py.f2py_calc_setter_var_params(
        f_arr(xm), f_arr(xp2), f_arr(skx), f_arr(sgn_wpxp), f_arr(f_x), f_arr(zeta_x),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_responder_params(nz: int, ngrdcol: int, xm, xp2, skx, sgn_wpxp, f_x, mixt_frac):
    """Compute responder-PDF component means/variances and the variance coefficients."""
    return clubb_f2py.f2py_calc_responder_params(
        f_arr(xm), f_arr(xp2), f_arr(skx), f_arr(sgn_wpxp), f_arr(f_x), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_limits_f_x_responder(
    nz: int, ngrdcol: int, mixt_frac, skx, sgn_wpxp,
    max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp,
):
    """Compute allowable lower/upper bounds for responder F_x."""
    return clubb_f2py.f2py_calc_limits_f_x_responder(
        f_arr(mixt_frac), f_arr(skx), f_arr(sgn_wpxp),
        f_arr(max_skx2_pos_skx_sgn_wpxp), f_arr(max_skx2_neg_skx_sgn_wpxp),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calc_coef_wp4_implicit(nz: int, ngrdcol: int, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Compute coefficient such that <w'^4> = coef_wp4_implicit * <w'^2>^2."""
    return clubb_f2py.f2py_calc_coef_wp4_implicit(
        f_arr(mixt_frac), f_arr(f_w), f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol))


def calc_coef_wpxp2_implicit(
    nz: int, ngrdcol: int,
    wp2, xp2, wpxp, sgn_wpxp, mixt_frac, f_w, f_x,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
):
    """Compute coefficient such that <w'x'^2> = coef_wpxp2_implicit * <x'^2>."""
    return clubb_f2py.f2py_calc_coef_wpxp2_implicit(
        f_arr(wp2), f_arr(xp2), f_arr(wpxp), f_arr(sgn_wpxp), f_arr(mixt_frac), f_arr(f_w), f_arr(f_x),
        f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd), f_arr(coef_sigma_x_1_sqd), f_arr(coef_sigma_x_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calc_coefs_wp2xp_semiimpl(
    nz: int, ngrdcol: int,
    wp2, xp2, sgn_wpxp, mixt_frac, f_w, f_x,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
):
    """Compute (<coef_wp2xp_implicit>, <term_wp2xp_explicit>) for semi-implicit transport."""
    return clubb_f2py.f2py_calc_coefs_wp2xp_semiimpl(
        f_arr(wp2), f_arr(xp2), f_arr(sgn_wpxp), f_arr(mixt_frac), f_arr(f_w), f_arr(f_x),
        f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd), f_arr(coef_sigma_x_1_sqd), f_arr(coef_sigma_x_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )


def calc_coefs_wpxpyp_semiimpl(
    nz: int, ngrdcol: int,
    wp2, xp2, yp2, wpxp, wpyp, sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, f_y,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    coef_sigma_y_1_sqd, coef_sigma_y_2_sqd,
):
    """Compute (<coef_wpxpyp_implicit>, <term_wpxpyp_explicit>) for semi-implicit transport."""
    return clubb_f2py.f2py_calc_coefs_wpxpyp_semiimpl(
        f_arr(wp2), f_arr(xp2), f_arr(yp2), f_arr(wpxp), f_arr(wpyp), f_arr(sgn_wpxp), f_arr(sgn_wpyp),
        f_arr(mixt_frac), f_arr(f_w), f_arr(f_x), f_arr(f_y), f_arr(coef_sigma_w_1_sqd), f_arr(coef_sigma_w_2_sqd),
        f_arr(coef_sigma_x_1_sqd), f_arr(coef_sigma_x_2_sqd), f_arr(coef_sigma_y_1_sqd), f_arr(coef_sigma_y_2_sqd),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )
