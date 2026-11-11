"""User-facing wrappers for routines from CLUBB_core/new_tsdadg_pdf.F90."""

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


def calc_l_x_skx_fnc(skx: float, sgn_wpxp: float, small_l_x_1: float, small_l_x_2: float):
    """Compute tsdadg skewness-function parameters big_L_x_1 and big_L_x_2."""
    return clubb_f2py.f2py_calc_l_x_skx_fnc(
        float(skx), float(sgn_wpxp), float(small_l_x_1), float(small_l_x_2))


def calc_setter_parameters_tsdadg(
    xm: float, xp2: float, skx: float, sgn_wpxp: float, big_l_x_1: float, big_l_x_2: float,
):
    """Compute tsdadg setter means/variances and mixture fraction from inputs."""
    return clubb_f2py.f2py_calc_setter_parameters_tsdadg(
        float(xm), float(xp2), float(skx), float(sgn_wpxp), float(big_l_x_1), float(big_l_x_2))


def tsdadg_pdf_driver(nz: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp):
    """Compute full tsdadg PDF moments for w/rt/thl."""
    arrays = [_fa_2d(arg) for arg in (wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_tsdadg_pdf_driver(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)
