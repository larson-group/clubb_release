"""User-facing wrappers for routines from CLUBB_core/new_tsdadg_pdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


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


def tsdadg_pdf_driver(nz: int, ngrdcol: int, wm, rtm, thlm, wp2, rtp2, thlp2, skw, skrt, skthl, wprtp, wpthlp):
    """Compute full tsdadg PDF moments for w/rt/thl."""
    return clubb_f2py.f2py_tsdadg_pdf_driver(
        f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2),
        f_arr(skw), f_arr(skrt), f_arr(skthl), f_arr(wprtp), f_arr(wpthlp),
        nz=int(nz), ngrdcol=int(ngrdcol),
    )
