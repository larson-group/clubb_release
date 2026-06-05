"""User-facing wrappers for routines from CLUBB_core/new_hybrid_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def _fa_2d(arr):
    out = f_arr(arr, dtype=np.float64)
    if out.ndim != 2:
        raise ValueError(f"Expected a rank-2 array, got rank-{out.ndim}.")
    return out


def _loop_scalar_tuple(*arrays, func):
    arrays = [_fa_2d(arr) for arr in arrays]
    ngrdcol, nz = arrays[0].shape
    sample = func(*(float(arrays[idx][0, 0]) for idx in range(len(arrays))))
    outputs = [np.empty((ngrdcol, nz), dtype=np.float64) for _ in range(len(sample))]
    for i in range(ngrdcol):
        for k in range(nz):
            result = func(*(float(arr[i, k]) for arr in arrays))
            for out, value in zip(outputs, result, strict=True):
                out[i, k] = value
    return tuple(f_arr(out) for out in outputs)


def _loop_scalar_array(*arrays, func):
    arrays = [_fa_2d(arr) for arr in arrays]
    ngrdcol, nz = arrays[0].shape
    out = np.empty((ngrdcol, nz), dtype=np.float64)
    for i in range(ngrdcol):
        for k in range(nz):
            out[i, k] = func(*(float(arr[i, k]) for arr in arrays))
    return f_arr(out)


def calculate_w_params(nz: int, ngrdcol: int, wm, wp2, skw, f_w, zeta_w):
    """Compute the setting-variable component means/stddevs and mixture fraction."""
    return _loop_scalar_tuple(
        wm, wp2, skw, f_w, zeta_w,
        func=clubb_f2py.f2py_calculate_w_params,
    )


def calculate_responder_params(
    nz: int, ngrdcol: int, xm, xp2, skx, wpxp, wp2, f_w, mixt_frac,
):
    """Compute responder component means/variances and the variance coefficients."""
    return _loop_scalar_tuple(
        xm, xp2, skx, wpxp, wp2, f_w, mixt_frac,
        func=clubb_f2py.f2py_calculate_responder_params,
    )


def calculate_coef_wp4_implicit(mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Compute coefficient such that <w'^4> = coef_wp4_implicit * <w'^2>^2."""
    return _loop_scalar_array(
        mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
        func=clubb_f2py.f2py_calculate_coef_wp4_implicit,
    )


def calc_coef_wp2xp_implicit(wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Compute coefficient such that <w'^2 x'> = coef_wp2xp_implicit * <w'x'>."""
    return _loop_scalar_array(
        wp2, mixt_frac, f_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
        func=clubb_f2py.f2py_calc_coef_wp2xp_implicit,
    )


def calc_coefs_wpxp2_semiimpl(
    nz: int, ngrdcol: int, wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
):
    """Compute (<coef_wpxp2_implicit>, <term_wpxp2_explicit>) for semi-implicit transport."""
    return _loop_scalar_tuple(
        wp2, wpxp, mixt_frac, f_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
        func=clubb_f2py.f2py_calc_coefs_wpxp2_semiimpl,
    )
