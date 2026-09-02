"""Push/pull converters for ErrInfo <-> Fortran module storage."""

from typing import Optional

import numpy as np

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo

_cached_err_info: Optional[ErrInfo] = None


def _as_latlon(arr: Optional[np.ndarray], ngrdcol: int, name: str) -> np.ndarray:
    """Normalize lat/lon arrays to Fortran-contiguous 1D vectors."""
    if arr is None:
        out = np.zeros(ngrdcol, dtype=np.float64)
    else:
        out = np.asfortranarray(arr, dtype=np.float64).reshape(-1)
        if out.size != ngrdcol:
            raise ValueError(
                f"{name} length mismatch: expected {ngrdcol}, got {out.size}"
            )
    return out


def init_err_info(err: ErrInfo):
    """Initialize Fortran err_info storage and push metadata values."""
    clubb_f2py.f2py_init_err_info(int(err.ngrdcol))
    set_fortran_err_info(err)


def set_fortran_err_info(err: ErrInfo):
    """Push ErrInfo metadata and err_code into Fortran module storage."""
    global _cached_err_info
    ncol = int(err.ngrdcol)
    lat = _as_latlon(err.lat, ncol, "lat")
    lon = _as_latlon(err.lon, ncol, "lon")
    clubb_f2py.f2py_set_err_info_values(int(err.chunk_idx), int(err.mpi_rank), lat, lon)
    if err.err_code is not None:
        err_code = np.asfortranarray(err.err_code, dtype=np.int32).reshape(-1)
        if err_code.size != ncol:
            raise ValueError(
                f"err_code length mismatch: expected {ncol}, got {err_code.size}"
            )
        clubb_f2py.f2py_set_err_code(err_code, ngrdcol=ncol)
    else:
        err_code = None
    _cached_err_info = err._replace(lat=lat, lon=lon, err_code=err_code)


def reset_cached_err_code() -> None:
    """Mirror a Fortran err_code reset in the cached Python ErrInfo."""
    global _cached_err_info
    if _cached_err_info is None:
        return
    _cached_err_info = _cached_err_info.reset_code()


def get_fortran_err_info() -> ErrInfo:
    """Pull ErrInfo values from Fortran module storage."""
    global _cached_err_info

    has_native_pull = (
        hasattr(clubb_f2py, "get_err_info_dims")
        and hasattr(clubb_f2py, "f2py_get_err_info_values")
    )

    if has_native_pull:
        ncol = int(clubb_f2py.get_err_info_dims())
        if ncol <= 0:
            pulled = ErrInfo(
                ngrdcol=0,
                chunk_idx=1,
                mpi_rank=0,
                lat=np.zeros(0, dtype=np.float64),
                lon=np.zeros(0, dtype=np.float64),
                err_code=np.zeros(0, dtype=np.int32),
            )
            _cached_err_info = pulled
            return pulled

        chunk_idx, mpi_rank, lat, lon = clubb_f2py.f2py_get_err_info_values(ncol)
        err_code = np.asfortranarray(clubb_f2py.get_err_code(ncol), dtype=np.int32)
        pulled = ErrInfo(
            ngrdcol=ncol,
            chunk_idx=int(chunk_idx),
            mpi_rank=int(mpi_rank),
            lat=np.asfortranarray(lat, dtype=np.float64),
            lon=np.asfortranarray(lon, dtype=np.float64),
            err_code=err_code,
        )
        _cached_err_info = pulled
        return pulled

    if _cached_err_info is None:
        raise RuntimeError(
            "get_fortran_err_info() requires prior init/push to infer ngrdcol when "
            "native err_info pull bindings are unavailable."
        )

    ncol = int(_cached_err_info.ngrdcol)
    err_code = np.asfortranarray(clubb_f2py.get_err_code(ncol), dtype=np.int32)
    pulled = ErrInfo(
        ngrdcol=ncol,
        chunk_idx=int(_cached_err_info.chunk_idx),
        mpi_rank=int(_cached_err_info.mpi_rank),
        lat=_as_latlon(_cached_err_info.lat, ncol, "lat"),
        lon=_as_latlon(_cached_err_info.lon, ncol, "lon"),
        err_code=err_code,
    )
    _cached_err_info = pulled
    return pulled
