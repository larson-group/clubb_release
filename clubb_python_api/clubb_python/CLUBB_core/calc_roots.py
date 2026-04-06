"""User-facing wrappers for routines from CLUBB_core/calc_roots.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def cubic_solve(nz: int, a_coef, b_coef, c_coef, d_coef):
    """Solve cubic equations and return complex roots with shape (nz, 3)."""
    roots_real, roots_imag = clubb_f2py.f2py_cubic_solve(
        f_arr(a_coef), f_arr(b_coef), f_arr(c_coef), f_arr(d_coef), nz=int(nz)
    )
    return roots_real + 1j * roots_imag


def quadratic_solve(nz: int, a_coef, b_coef, c_coef):
    """Solve quadratic equations and return complex roots with shape (nz, 2)."""
    roots_real, roots_imag = clubb_f2py.f2py_quadratic_solve(
        f_arr(a_coef), f_arr(b_coef), f_arr(c_coef), nz=int(nz)
    )
    return roots_real + 1j * roots_imag
