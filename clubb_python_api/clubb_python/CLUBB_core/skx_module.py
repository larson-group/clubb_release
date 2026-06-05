"""User-facing wrappers for routines from CLUBB_core/Skx_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def skx_func(nz: int, ngrdcol: int, xp2, xp3, x_tol: float, clubb_params):
    """Compute skewness from variance and third moment."""
    return clubb_f2py.f2py_skx_func(
        f_arr(xp2), f_arr(xp3), float(x_tol), f_arr(clubb_params),
        nz=int(nz), ngrdcol=int(ngrdcol))


def lg_2005_ansatz(nz: int, ngrdcol: int, skw, wpxp, wp2, xp2, beta, sigma_sqd_w, x_tol: float):
    """Compute skewness using the Larson and Golaz (2005) ansatz."""
    return clubb_f2py.f2py_lg_2005_ansatz(
        f_arr(skw), f_arr(wpxp), f_arr(wp2), f_arr(xp2), f_arr(beta),
        f_arr(sigma_sqd_w), float(x_tol), nz=int(nz), ngrdcol=int(ngrdcol))


def xp3_lg_2005_ansatz(
    nzt: int, ngrdcol: int, skw_zt, wpxp_zt, wp2_zt, xp2_zt, sigma_sqd_w_zt, clubb_params,
    x_tol: float,
):
    """Compute third moment from LG05 ansatz skewness closure."""
    return clubb_f2py.f2py_xp3_lg_2005_ansatz(
        f_arr(skw_zt), f_arr(wpxp_zt), f_arr(wp2_zt), f_arr(xp2_zt),
        f_arr(sigma_sqd_w_zt), f_arr(clubb_params), float(x_tol),
        nzt=int(nzt), ngrdcol=int(ngrdcol))


def compute_gamma_skw(nz: int | None = None, ngrdcol: int | None = None, skw=None, clubb_params=None, **compat_kwargs):
    """Compute gamma as a function of w skewness on zm and zt grids."""
    if skw is None:
        skw = compat_kwargs.pop("skw_zm")
    if clubb_params is None:
        clubb_params = compat_kwargs.pop("clubb_params")
    if nz is None:
        nz = int(np.asarray(skw).shape[1])
    if ngrdcol is None:
        ngrdcol = int(np.asarray(skw).shape[0])
    l_gamma_skw = compat_kwargs.pop("l_gamma_skw", True)
    nzt = int(compat_kwargs.pop("nzt", nz))
    skw_zt = compat_kwargs.pop("skw_zt", skw)
    return clubb_f2py.f2py_compute_gamma_skw(
        l_gamma_skw, f_arr(skw), f_arr(clubb_params), f_arr(skw_zt),
        nzm=int(nz), nzt=int(nzt), ngrdcol=int(ngrdcol))
