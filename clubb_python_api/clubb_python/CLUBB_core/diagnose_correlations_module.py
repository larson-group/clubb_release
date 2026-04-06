"""User-facing wrappers for routines from CLUBB_core/diagnose_correlations_module.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def diagnose_correlations(pdf_dim: int, iipdf_w: int, corr_array_pre, l_calc_w_corr: bool):
    """Diagnose the SILHS-ready correlation matrix from a prescribed one."""
    return clubb_f2py.f2py_diagnose_correlations(
        int(iipdf_w), f_arr(corr_array_pre), l_calc_w_corr, pdf_dim=int(pdf_dim)
    )


def calc_cholesky_corr_mtx_approx(n_variables: int, iipdf_w: int, corr_matrix):
    """Compute the transposed Cholesky correlation matrix and its approximation."""
    return clubb_f2py.f2py_calc_cholesky_corr_mtx_approx(
        int(iipdf_w), f_arr(corr_matrix), n_variables=int(n_variables)
    )
