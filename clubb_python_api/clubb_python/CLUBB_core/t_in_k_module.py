"""User-facing wrappers for routines from CLUBB_core/T_in_K_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py


def thlm2t_in_k(nz: int, thlm: np.ndarray, exner: np.ndarray,
                rcm: np.ndarray) -> np.ndarray:
    """Convert liquid water potential temperature to temperature in Kelvin."""
    del nz
    return clubb_f2py.f2py_thlm2t_in_k_1d(
        f_arr(thlm),
        f_arr(exner),
        f_arr(rcm),
    )
