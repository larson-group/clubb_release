"""Initialize, reuse and close native CLUBB profile-loss requests.

init_clubb_loss(runfile, return_default_params=True) returns field names and the
full default parameter matrix; get_clubb_params_all reads that matrix later.
clubb_get_loss_for_params scores full matrices against the prepared benchmark
and fixed weights. It returns (rmse,correlation,std_ratio,bias_norm,case_loss):
four (window,field,candidate) arrays and one candidate vector. finalize_clubb_loss
releases the process-global session. clubb_get_loss performs this lifecycle once
and prepends field names to those five outputs. Only one session is active per
process; rebuild the native/F2PY binding after interface changes.
"""

from __future__ import annotations

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.CLUBB_core.parameters_tunable import _get_nparams
from clubb_python.string_conversion import fortran_char_matrix_to_python_strings

_loss_initialized = False
_loss_num_variables = 0
_loss_total_param_sets = 0
_loss_num_time_windows = 1
_loss_var_names: list[str] = []


def init_clubb_loss(runfile: str, return_default_params: bool = False):
    """Initialize the reusable CLUBB loss session from an aggregate runfile."""
    global _loss_initialized, _loss_num_variables, _loss_total_param_sets, _loss_num_time_windows, _loss_var_names

    num_variables, total_param_sets, num_time_windows, raw_var_names = clubb_f2py.f2py_init_clubb_loss(runfile)

    _loss_num_variables = int(num_variables)
    _loss_total_param_sets = int(total_param_sets)
    _loss_num_time_windows = int(num_time_windows)
    _loss_var_names = fortran_char_matrix_to_python_strings(raw_var_names)[:_loss_num_variables]
    _loss_initialized = True

    clubb_params_all = None
    if return_default_params:
        clubb_params_all = get_clubb_params_all()

    return list(_loss_var_names), clubb_params_all


def get_clubb_params_all() -> np.ndarray:
    """Return the initialized full default CLUBB parameter matrix."""
    if not _loss_initialized:
        raise RuntimeError("get_clubb_params_all requires init_clubb_loss to be called first.")

    return clubb_f2py.f2py_get_clubb_params_all(int(_loss_total_param_sets), _get_nparams())


def clubb_get_loss_for_params(clubb_params_all: np.ndarray):
    """Score one full parameter matrix against the active prepared loss request."""
    if not _loss_initialized:
        raise RuntimeError("clubb_get_loss_for_params requires init_clubb_loss to be called first.")

    params = np.asarray(clubb_params_all, dtype=np.float64)
    expected_shape = (_loss_total_param_sets, _get_nparams())
    if params.shape != expected_shape:
        raise ValueError(
            f"clubb_params_all must have shape {expected_shape}, got {params.shape}."
        )

    result = clubb_f2py.f2py_clubb_get_loss_for_params(
        int(_loss_num_variables),
        int(_loss_num_time_windows),
        f_arr(params),
    )
    if len(result) != 5:
        raise RuntimeError('Loss API requires four profile metrics and one case-loss vector; rebuild the Python API')
    expected_metrics = (_loss_num_time_windows, _loss_num_variables, _loss_total_param_sets)
    if (any(np.shape(value) != expected_metrics for value in result[:4])
            or np.shape(result[4]) != (_loss_total_param_sets,)):
        raise RuntimeError('Loss API profile/case shapes do not match the prepared request; rebuild the Python API')
    return tuple(result)



def finalize_clubb_loss():
    """Release the active reusable CLUBB loss session."""
    global _loss_initialized, _loss_num_variables, _loss_total_param_sets, _loss_num_time_windows, _loss_var_names

    clubb_f2py.f2py_finalize_clubb_loss()
    _loss_initialized = False
    _loss_num_variables = 0
    _loss_total_param_sets = 0
    _loss_num_time_windows = 1
    _loss_var_names = []


def clubb_get_loss(runfile: str):
    """Run the one-shot CLUBB loss path using the runfile's default parameters."""
    clubb_var_names, clubb_params_all = init_clubb_loss(
        runfile,
        return_default_params=True,
    )
    try:
        rmse, correlation, std_ratio, bias_norm, case_loss = clubb_get_loss_for_params(
            clubb_params_all
        )
    finally:
        finalize_clubb_loss()

    return clubb_var_names, rmse, correlation, std_ratio, bias_norm, case_loss
