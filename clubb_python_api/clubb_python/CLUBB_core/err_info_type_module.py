"""User-facing wrappers for routines from CLUBB_core/err_info_type_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import (
    init_err_info as init_err_info_udt,
    set_fortran_err_info,
)


def init_err_info(err_info: ErrInfo | int):
    """Initialize the error-info derived-type data.

    Preferred usage is passing an `ErrInfo` object, which is pushed before
    return. Integer `ngrdcol` is retained for backward compatibility.
    """
    if isinstance(err_info, ErrInfo):
        init_err_info_udt(err_info)
        return err_info

    clubb_f2py.f2py_init_err_info(int(err_info))
    return ErrInfo(ngrdcol=int(err_info))


def get_err_code(ngrdcol: int) -> np.ndarray:
    """Get the error code array from Fortran module storage."""
    return clubb_f2py.get_err_code(ngrdcol)


def cleanup_err_info(err_info: ErrInfo):
    """Deallocate stored err_info arrays."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_cleanup_err_info()


def set_err_info_values(ngrdcol: int, chunk_idx: int, mpi_rank: int,
                        lat, lon, err_info: ErrInfo):
    """Set error context info (chunk index, MPI rank, lat/lon)."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_set_err_info_values(chunk_idx, mpi_rank, f_arr(lat), f_arr(lon))
