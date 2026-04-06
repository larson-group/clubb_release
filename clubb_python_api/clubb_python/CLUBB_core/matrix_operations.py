"""User-facing wrappers for routines from CLUBB_core/matrix_operations.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def cholesky_factor(ndim: int, a_input):
    """Compute a lower-triangular Cholesky factorization with CLUBB's scaling path."""
    return clubb_f2py.f2py_cholesky_factor(f_arr(a_input), ndim=int(ndim))


def mirror_lower_triangular_matrix(nvars: int, matrix):
    """Mirror a lower-triangular square matrix into the upper triangle."""
    return clubb_f2py.f2py_mirror_lower_triangular_matrix(f_arr(matrix), nvars=int(nvars))
