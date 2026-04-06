"""User-facing wrappers for routines from CLUBB_core/pdf_utilities.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def smooth_corr_quotient(ngrdcol: int, nz: int, numerator, denominator, denom_thresh: float):
    """Compute smoothed/bounded quotient used for PDF correlations."""
    return clubb_f2py.f2py_smooth_corr_quotient(
        f_arr(numerator), f_arr(denominator), float(denom_thresh),
        ngrdcol=int(ngrdcol), nz=int(nz))


def calc_comp_corrs_binormal(
    nz: int, ngrdcol: int,
    xpyp, xm, ym, mu_x_1, mu_x_2, mu_y_1, mu_y_2,
    sigma_x_1_sqd, sigma_x_2_sqd, sigma_y_1_sqd, sigma_y_2_sqd, mixt_frac,
):
    """Compute binormal component correlations corr_x_y_1 and corr_x_y_2."""
    return clubb_f2py.f2py_calc_comp_corrs_binormal(
        f_arr(xpyp), f_arr(xm), f_arr(ym), f_arr(mu_x_1), f_arr(mu_x_2), f_arr(mu_y_1), f_arr(mu_y_2),
        f_arr(sigma_x_1_sqd), f_arr(sigma_x_2_sqd), f_arr(sigma_y_1_sqd), f_arr(sigma_y_2_sqd),
        f_arr(mixt_frac), nz=int(nz), ngrdcol=int(ngrdcol),
    )


def mean_l2n(mu_x: float, sigma2_on_mu2: float) -> float:
    """Compute the lognormal-space mean from linear-space mean and variance ratio."""
    return float(clubb_f2py.f2py_mean_l2n(float(mu_x), float(sigma2_on_mu2)))


def stdev_l2n(sigma2_on_mu2: float) -> float:
    """Compute the lognormal-space standard deviation from variance ratio."""
    return float(clubb_f2py.f2py_stdev_l2n(float(sigma2_on_mu2)))


def corr_nn2nl(nz: int, ngrdcol: int, corr_x_y_n, sigma_y_n, y_sigma2_on_mu2):
    """Convert corr(x, ln y) to corr(x, y) for 2D fields."""
    return clubb_f2py.f2py_corr_nn2nl(
        f_arr(corr_x_y_n), f_arr(sigma_y_n), f_arr(y_sigma2_on_mu2), nz=int(nz), ngrdcol=int(ngrdcol)
    )


def corr_nn2ll(nz: int, ngrdcol: int, corr_x_y_n, sigma_x_n, sigma_y_n, x_sigma2_on_mu2, y_sigma2_on_mu2):
    """Convert corr(ln x, ln y) to corr(x, y) for 2D fields."""
    return clubb_f2py.f2py_corr_nn2ll(
        f_arr(corr_x_y_n), f_arr(sigma_x_n), f_arr(sigma_y_n), f_arr(x_sigma2_on_mu2), f_arr(y_sigma2_on_mu2),
        nz=int(nz), ngrdcol=int(ngrdcol)
    )


def compute_mean_binormal(mu_x_1: float, mu_x_2: float, mixt_frac: float) -> float:
    """Compute the overall mean of a two-component binormal mixture."""
    return float(clubb_f2py.f2py_compute_mean_binormal(float(mu_x_1), float(mu_x_2), float(mixt_frac)))


def compute_variance_binormal(
    xm: float, mu_x_1: float, mu_x_2: float, stdev_x_1: float, stdev_x_2: float, mixt_frac: float
) -> float:
    """Compute the overall variance of a two-component binormal mixture."""
    return float(
        clubb_f2py.f2py_compute_variance_binormal(
            float(xm), float(mu_x_1), float(mu_x_2), float(stdev_x_1), float(stdev_x_2), float(mixt_frac)
        )
    )


def calc_corr_chi_x(
    crt_i: float, cthl_i: float, sigma_rt_i: float, sigma_thl_i: float, sigma_chi_i: float,
    corr_rt_x_i: float, corr_thl_x_i: float,
) -> float:
    """Compute the ith-component correlation between chi and x."""
    return float(
        clubb_f2py.f2py_calc_corr_chi_x(
            float(crt_i), float(cthl_i), float(sigma_rt_i), float(sigma_thl_i), float(sigma_chi_i),
            float(corr_rt_x_i), float(corr_thl_x_i),
        )
    )


def calc_corr_eta_x(
    crt_i: float, cthl_i: float, sigma_rt_i: float, sigma_thl_i: float, sigma_eta_i: float,
    corr_rt_x_i: float, corr_thl_x_i: float,
) -> float:
    """Compute the ith-component correlation between eta and x."""
    return float(
        clubb_f2py.f2py_calc_corr_eta_x(
            float(crt_i), float(cthl_i), float(sigma_rt_i), float(sigma_thl_i), float(sigma_eta_i),
            float(corr_rt_x_i), float(corr_thl_x_i),
        )
    )
