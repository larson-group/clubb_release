"""User-facing wrappers for routines from CLUBB_core/new_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def _fa_2d(arr):
    out = f_arr(arr, dtype=np.float64)
    if out.ndim != 2:
        raise ValueError(f"Expected a rank-2 array, got rank-{out.ndim}.")
    return out


def _stack_tuple_outputs(column_outputs):
    return tuple(f_arr(np.stack(items, axis=0)) for items in zip(*column_outputs, strict=True))


def _stack_single_output(column_outputs):
    return f_arr(np.stack(column_outputs, axis=0))


def calc_setter_var_params(nz: int, xm, xp2, skx, sgn_wpxp, f_x, zeta_x, **compat_kwargs):
    """Compute new-PDF setter component means/stddevs, mixture fraction, and coefficients."""
    arrays = [_fa_2d(arg) for arg in (xm, xp2, skx, sgn_wpxp, f_x, zeta_x)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_setter_var_params(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def calc_responder_params(nz: int, xm, xp2, skx, sgn_wpxp, f_x, mixt_frac, **compat_kwargs):
    """Compute responder-PDF component means/variances and the variance coefficients."""
    arrays = [_fa_2d(arg) for arg in (xm, xp2, skx, sgn_wpxp, f_x, mixt_frac)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_responder_params(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def calc_limits_f_x_responder(
    nz: int, mixt_frac, skx, sgn_wpxp,
    max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp,
    **compat_kwargs,
):
    """Compute allowable lower/upper bounds for responder F_x."""
    arrays = [_fa_2d(arg) for arg in (
        mixt_frac, skx, sgn_wpxp, max_skx2_pos_skx_sgn_wpxp, max_skx2_neg_skx_sgn_wpxp
    )]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_limits_f_x_responder(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def calc_coef_wp4_implicit(nz: int, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, **compat_kwargs):
    """Compute coefficient such that <w'^4> = coef_wp4_implicit * <w'^2>^2."""
    arrays = [_fa_2d(arg) for arg in (mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_coef_wp4_implicit(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_single_output(column_outputs)


def calc_coef_wpxp2_implicit(
    nz: int,
    wp2, xp2, wpxp, sgn_wpxp, mixt_frac, f_w, f_x,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    **compat_kwargs,
):
    """Compute coefficient such that <w'x'^2> = coef_wpxp2_implicit * <x'^2>."""
    arrays = [_fa_2d(arg) for arg in (
        wp2, xp2, wpxp, sgn_wpxp, mixt_frac, f_w, f_x,
        coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
    )]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_coef_wpxp2_implicit(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_single_output(column_outputs)


def calc_coefs_wp2xp_semiimpl(
    nz: int,
    wp2, xp2, sgn_wpxp, mixt_frac, f_w, f_x,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    **compat_kwargs,
):
    """Compute (<coef_wp2xp_implicit>, <term_wp2xp_explicit>) for semi-implicit transport."""
    arrays = [_fa_2d(arg) for arg in (
        wp2, xp2, sgn_wpxp, mixt_frac, f_w, f_x,
        coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd
    )]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_coefs_wp2xp_semiimpl(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def calc_coefs_wpxpyp_semiimpl(
    nz: int,
    wp2, xp2, yp2, wpxp, wpyp, sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, f_y,
    coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
    coef_sigma_y_1_sqd, coef_sigma_y_2_sqd,
    **compat_kwargs,
):
    """Compute (<coef_wpxpyp_implicit>, <term_wpxpyp_explicit>) for semi-implicit transport."""
    arrays = [_fa_2d(arg) for arg in (
        wp2, xp2, yp2, wpxp, wpyp, sgn_wpxp, sgn_wpyp, mixt_frac, f_w, f_x, f_y,
        coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
        coef_sigma_y_1_sqd, coef_sigma_y_2_sqd
    )]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_coefs_wpxpyp_semiimpl(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)
