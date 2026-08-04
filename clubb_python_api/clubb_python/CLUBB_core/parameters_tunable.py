"""User-facing wrappers for routines from CLUBB_core/parameters_tunable.F90."""

from __future__ import annotations

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info
from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.nu_vert_res_dep_converter import (
    get_fortran_nu_vert_res_dep,
    set_fortran_nu_vert_res_dep,
)
from clubb_python.string_conversion import fortran_char_matrix_to_python_strings


def _get_nparams() -> int:
    """Return the compiled tunable-parameter count from Fortran."""
    return int(clubb_f2py.f2py_get_nparams())


def init_clubb_params(ngrdcol: int, iunit: int, filename: str) -> np.ndarray:
    """Initialize CLUBB tunable parameters."""
    return clubb_f2py.f2py_init_clubb_params(ngrdcol, iunit, filename, _get_nparams())


def get_param_names() -> list[str]:
    """Return tunable parameter names in CLUBB's native ordering."""
    return fortran_char_matrix_to_python_strings(clubb_f2py.f2py_get_param_names(_get_nparams()))


def get_parameter_hard_bounds(nparams_in: int) -> list[dict]:
    """Return Fortran-owned scalar validity bounds in parameter-list order.

    The Fortran table is a direct ``[minimum, maximum]`` catalog.  ``None``
    denotes an open endpoint.  Values already lie inside strict physical
    endpoints, so clients need no separate strictness metadata.
    """
    hard_bounds = clubb_f2py.f2py_get_parameter_hard_bounds(int(nparams_in))
    no_bound = np.finfo(hard_bounds.dtype).max
    return [
        {
            "name": name,
            "min": float(hard_bounds[0, idx]) if hard_bounds[0, idx] > -no_bound else None,
            "max": float(hard_bounds[1, idx]) if hard_bounds[1, idx] < no_bound else None,
        }
        for idx, name in enumerate(get_param_names())
    ]


def calc_derrived_params(
    gr: Grid, ngrdcol: int, grid_type: int,
    deltaz: np.ndarray,
    clubb_params: np.ndarray,
    l_prescribed_avg_deltaz: bool,
    nu_vert_res_dep: NuVertResDep | None,
):
    """Calculate derived parameters and return updated nu coefficient data."""
    set_fortran_grid(gr)
    if nu_vert_res_dep is not None:
        set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    lmin, mixt_frac_max_mag = clubb_f2py.f2py_calc_derrived_params(
        grid_type,
        f_arr(deltaz),
        f_arr(clubb_params),
        l_prescribed_avg_deltaz,
        ngrdcol=int(ngrdcol),
    )
    return get_fortran_nu_vert_res_dep(int(ngrdcol)), float(lmin), float(mixt_frac_max_mag)


def check_parameters(ngrdcol: int, clubb_params, err_info: ErrInfo):
    """Validate tunable parameters."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_check_parameters(f_arr(clubb_params), ngrdcol=int(ngrdcol))
    return get_fortran_err_info()
