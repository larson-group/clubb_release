"""User-facing wrappers for routines from CLUBB_core/ly93_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def _fa_2d(arr):
    out = f_arr(arr, dtype=np.float64)
    if out.ndim != 2:
        raise ValueError(f"Expected a rank-2 array, got rank-{out.ndim}.")
    return out


def calc_params_ly93(nz: int, ngrdcol: int, xm, xp2, skx, mixt_frac):
    """Compute LY93 component means/variances from mean, variance, and skewness."""
    arrays = [_fa_2d(arg) for arg in (xm, xp2, skx, mixt_frac)]
    return clubb_f2py.f2py_calc_params_ly93(*arrays, nz=int(nz), ngrdcol=int(ngrdcol))


def ly93_driver(nz: int, ngrdcol: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl):
    """Compute full LY93 means/variances and mixture fraction for w/rt/thl."""
    arrays = [_fa_2d(arg) for arg in (wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl)]
    return clubb_f2py.f2py_ly93_driver(*arrays, nz=int(nz), ngrdcol=int(ngrdcol))
