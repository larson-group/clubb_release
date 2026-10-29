"""Check the JAX Skx_func port (Skx_module.F90:Skx_func) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.skx_module as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func

NG, NZ = 2, 6
NPARAMS = 102
I_SKW_DENOM = 73   # 1-based iSkw_denom_coef
_DENOM_COEF = 0.15
_X_TOL = 1.0e-2


def _params():
    p = np.zeros((NG, NPARAMS))
    p[:, I_SKW_DENOM - 1] = _DENOM_COEF
    return p


def test_f2py_oracle():
    """Skx_func matches Fortran for the same clubb_params (only iSkw_denom_coef is read)."""
    rng = np.random.default_rng(4)
    params = _params()
    worst = 0.0
    for _ in range(10):
        xp2 = rng.uniform(0.0, 2.0, (NG, NZ))
        xp3 = rng.uniform(-1.0, 1.0, (NG, NZ))
        args = (NZ, NG, xp2, xp3, _X_TOL, params)
        worst = max(worst, max_abs_diff("Skx_func", Skx_func(*args), fortran.skx_func(*args)))
    assert worst < 1e-11, f"Skx_func Fortran mismatch {worst:.2e}"


def test_LG_2005_ansatz_f2py():
    """LG_2005_ansatz and its <x'^3> wrapper xp3_LG_2005_ansatz (beta read from clubb_params) match Fortran.

    The ansatz's clip carries a ~1e-11 relative difference, so a relative tolerance is used."""
    from clubb_jax.src.CLUBB_core.Skx_module import LG_2005_ansatz, xp3_LG_2005_ansatz
    from clubb_jax.src.CLUBB_core.parameter_indices import ibeta
    params = _params(); params[:, ibeta] = 2.0   # ibeta is a zero-based column index
    beta = np.ascontiguousarray(params[:, ibeta])
    rng = np.random.default_rng(2)
    w_lg = w_xp3 = 0.0
    for _ in range(30):
        Skw = rng.uniform(-3, 3, (NG, NZ)); wpxp = rng.uniform(-1, 1, (NG, NZ))
        wp2 = rng.uniform(1e-2, 2, (NG, NZ)); xp2 = rng.uniform(1e-4, 1, (NG, NZ)); ssw = rng.uniform(0, 0.6, (NG, NZ))
        args = (NZ, NG, Skw, wpxp, wp2, xp2, beta, ssw, _X_TOL)
        w_lg = max(w_lg, max_rel_diff("LG_2005_ansatz", LG_2005_ansatz(*args), fortran.lg_2005_ansatz(*args), 1.0))
        args = (NZ, NG, Skw, wpxp, wp2, xp2, ssw, params, _X_TOL)
        w_xp3 = max(w_xp3, max_rel_diff("xp3_LG_2005_ansatz", xp3_LG_2005_ansatz(*args),
                                        fortran.xp3_lg_2005_ansatz(*args), 1.0))
    assert w_lg < 1e-9, f"LG_2005_ansatz Fortran rel mismatch {w_lg:.2e}"
    assert w_xp3 < 1e-9, f"xp3_LG_2005_ansatz Fortran rel mismatch {w_xp3:.2e}"
