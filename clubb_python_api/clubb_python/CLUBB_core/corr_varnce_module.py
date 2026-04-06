"""User-facing wrappers for routines from CLUBB_core/corr_varnce_module.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def assert_corr_symmetric(pdf_dim: int, corr_array_n, err_info: ErrInfo):
    """Assert that a correlation matrix is symmetric with unit diagonal."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_assert_corr_symmetric(f_arr(corr_array_n), pdf_dim=int(pdf_dim))
    return get_fortran_err_info()
