"""User-facing wrappers for routines from CLUBB_core/ly93_pdf.F90."""

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


def calc_params_ly93(nz: int, xm, xp2, skx, mixt_frac, **compat_kwargs):
    """Compute LY93 component means/variances from mean, variance, and skewness."""
    arrays = [_fa_2d(arg) for arg in (xm, xp2, skx, mixt_frac)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_calc_params_ly93(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)


def ly93_driver(nz: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, **compat_kwargs):
    """Compute full LY93 means/variances and mixture fraction for w/rt/thl."""
    arrays = [_fa_2d(arg) for arg in (wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl)]
    ngrdcol = arrays[0].shape[0]
    column_outputs = [
        clubb_f2py.f2py_ly93_driver(*(arr[i, :] for arr in arrays), nz=int(nz))
        for i in range(ngrdcol)
    ]
    return _stack_tuple_outputs(column_outputs)
