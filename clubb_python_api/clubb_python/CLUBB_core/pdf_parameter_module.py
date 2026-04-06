"""User-facing wrappers for routines from CLUBB_core/pdf_parameter_module.F90."""

import clubb_f2py

from clubb_python.derived_types.pdf_params import (
    implicit_coefs_terms,
    pdf_parameter,
)
from clubb_python.derived_types.pdf_params_converter import (
    get_fortran_implicit_coefs,
    get_fortran_pdf_params,
    get_fortran_pdf_params_zm,
    set_fortran_pdf_params,
    set_fortran_pdf_params_zm,
)

def init_pdf_params(nz: int, ngrdcol: int) -> pdf_parameter:
    """Initialize pdf_params in Fortran module storage and return pulled data."""
    clubb_f2py.f2py_init_pdf_params(nz, ngrdcol)
    return get_fortran_pdf_params()


def init_pdf_params_zm(nz: int, ngrdcol: int) -> pdf_parameter:
    """Initialize pdf_params_zm in Fortran module storage and return pulled data."""
    clubb_f2py.f2py_init_pdf_params_zm(nz, ngrdcol)
    return get_fortran_pdf_params_zm()


def init_pdf_implicit(
    nz: int,
    ngrdcol: int,
    sclr_dim: int,
) -> implicit_coefs_terms:
    """Initialize implicit coefs/terms in Fortran module storage and return pulled data."""
    clubb_f2py.f2py_init_pdf_implicit(nz, ngrdcol, sclr_dim)
    return get_fortran_implicit_coefs()


def init_pdf_implicit_coefs_terms(
    nz: int,
    ngrdcol: int,
    sclr_dim: int,
) -> implicit_coefs_terms:
    """Initialize implicit coefs/terms using the Fortran routine name."""
    return init_pdf_implicit(nz=nz, ngrdcol=ngrdcol, sclr_dim=sclr_dim)


def zero_pdf_params(
    which: int,
    pdf_params: pdf_parameter,
    pdf_params_zm: pdf_parameter,
):
    """Zero out stored PDF params."""
    if which == 1:
        set_fortran_pdf_params(pdf_params)
    if which == 2:
        set_fortran_pdf_params_zm(pdf_params_zm)
    clubb_f2py.f2py_zero_pdf_params(which)
    if which == 1:
        return get_fortran_pdf_params()
    if which == 2:
        return get_fortran_pdf_params_zm()
    return get_fortran_pdf_params(), get_fortran_pdf_params_zm()


def zero_pdf_implicit_coefs_terms():
    """Zero out stored implicit coefficient/term arrays."""
    clubb_f2py.f2py_zero_pdf_implicit_coefs_terms()
