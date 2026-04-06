"""User-facing wrappers for routines from CLUBB_core/precipitation_fraction.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def _transport_3d(arr, ncol: int, nz: int, logical_dim: int, dtype=np.float64):
    """Pad 3D transported arrays to at least size 1 in the trailing dim."""
    transport_dim = max(int(logical_dim), 1)
    out = np.zeros((int(ncol), int(nz), transport_dim), dtype=dtype, order="F")
    if int(logical_dim) > 0:
        out[:, :, : int(logical_dim)] = f_arr(arr)
    return out


def _transport_1d(arr, logical_dim: int, dtype=np.float64):
    """Pad 1D transported arrays to at least size 1."""
    transport_dim = max(int(logical_dim), 1)
    out = np.zeros((transport_dim,), dtype=dtype, order="F")
    if int(logical_dim) > 0:
        out[: int(logical_dim)] = f_arr(arr)
    return out


def precip_fraction(
    gr: Grid,
    nzt: int,
    ngrdcol: int,
    hydromet_dim: int,
    hydromet,
    cloud_frac,
    cloud_frac_1,
    l_mix_rat_hm,
    l_frozen_hm,
    hydromet_tol,
    cloud_frac_2,
    ice_supersat_frac,
    ice_supersat_frac_1,
    ice_supersat_frac_2,
    mixt_frac,
    clubb_params,
    err_info: ErrInfo,
):
    """Compute overall and component precipitation fractions."""
    set_fortran_grid(gr)
    set_fortran_err_info(err_info)
    precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol = clubb_f2py.f2py_precip_fraction(
        int(hydromet_dim),
        _transport_3d(hydromet, ngrdcol, nzt, hydromet_dim),
        f_arr(cloud_frac),
        f_arr(cloud_frac_1),
        _transport_1d(l_mix_rat_hm, hydromet_dim, dtype=np.bool_),
        _transport_1d(l_frozen_hm, hydromet_dim, dtype=np.bool_),
        _transport_1d(hydromet_tol, hydromet_dim),
        f_arr(cloud_frac_2),
        f_arr(ice_supersat_frac),
        f_arr(ice_supersat_frac_1),
        f_arr(ice_supersat_frac_2),
        f_arr(mixt_frac),
        f_arr(clubb_params).reshape(-1),
        nzt=int(nzt),
        ngrdcol=int(ngrdcol),
        hydromet_dim_transport=max(int(hydromet_dim), 1),
    )
    return get_fortran_err_info(), precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol
