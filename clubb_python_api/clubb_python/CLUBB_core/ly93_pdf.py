"""User-facing wrappers for routines from CLUBB_core/ly93_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def calc_params_ly93(nz: int, ngrdcol: int, xm, xp2, skx, mixt_frac):
    """Compute LY93 component means/variances from mean, variance, and skewness."""
    return clubb_f2py.f2py_calc_params_ly93(
        f_arr(xm), f_arr(xp2), f_arr(skx), f_arr(mixt_frac),
        nz=int(nz), ngrdcol=int(ngrdcol))


def ly93_driver(nz: int, ngrdcol: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl):
    """Compute full LY93 means/variances and mixture fraction for w/rt/thl."""
    return clubb_f2py.f2py_ly93_driver(
        f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2),
        f_arr(skw), f_arr(skrt), f_arr(skthl), nz=int(nz), ngrdcol=int(ngrdcol))
