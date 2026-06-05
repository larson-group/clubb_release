"""User-facing wrappers for routines from CLUBB_core/err_info_type_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import (
    init_err_info as init_err_info_udt,
    set_fortran_err_info,
)


def init_err_info(ngrdcol: ErrInfo | int, err_info: ErrInfo | None = None):
    """Initialize the error-info derived-type data.

    Preferred usage is passing an `ErrInfo` object, which is pushed before
    return. Integer `ngrdcol` is retained for backward compatibility.
    """
    if err_info is not None:
        init_err_info_udt(err_info)
        return err_info
    if isinstance(ngrdcol, ErrInfo):
        init_err_info_udt(ngrdcol)
        return ngrdcol

    clubb_f2py.f2py_init_err_info(int(ngrdcol))
    return ErrInfo(ngrdcol=int(ngrdcol))


def get_err_code(ngrdcol: int) -> np.ndarray:
    """Get the error code array from Fortran module storage."""
    return clubb_f2py.get_err_code(ngrdcol)


def cleanup_err_info(err_info: ErrInfo):
    """Deallocate stored err_info arrays."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_cleanup_err_info()


def set_err_info_values(ngrdcol: int, err_info: ErrInfo, chunk_idx_in: int, mpi_rank_in: int,
                        lat_in, lon_in, **compat_kwargs):
    """Set error context info (chunk index, MPI rank, lat/lon)."""
    chunk_idx_in = compat_kwargs.pop("chunk_idx", chunk_idx_in)
    mpi_rank_in = compat_kwargs.pop("mpi_rank", mpi_rank_in)
    lat_in = compat_kwargs.pop("lat", lat_in)
    lon_in = compat_kwargs.pop("lon", lon_in)
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_set_err_info_values(chunk_idx_in, mpi_rank_in, f_arr(lat_in), f_arr(lon_in))
