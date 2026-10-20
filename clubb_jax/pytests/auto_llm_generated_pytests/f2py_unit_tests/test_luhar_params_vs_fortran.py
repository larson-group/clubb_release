"""Check the JAX calc_Luhar_params port (adg1_adg2_3d_luhar_pdf.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.adg1_adg2_3d_luhar_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import calc_Luhar_params

NG, NZ = 2, 6
_X_TOL_SQD = 1.0e-4


def _inputs(seed):
    rng = np.random.default_rng(seed)
    Skx = rng.uniform(-3, 3, (NG, NZ))
    wpxp = rng.uniform(-1, 1, (NG, NZ))
    xp2 = rng.uniform(0.01, 2.0, (NG, NZ))
    return Skx, wpxp, xp2


def test_f2py_oracle():
    worst = 0.0
    n_const = 0
    for seed in (11, 22, 33):
        Skx, wpxp, xp2 = _inputs(seed)
        if seed == 33:                          # force some constant-x points
            xp2[:, :2] = 1e-8; n_const += 1
        f = fortran.calc_luhar_params(NZ, Skx, wpxp, xp2, _X_TOL_SQD)
        g = calc_Luhar_params(Skx, wpxp, xp2, _X_TOL_SQD)
        worst = max(worst, max_abs_diff("calc_Luhar_params", tuple(g), tuple(f)))
    assert worst < 1e-11, f"calc_luhar_params Fortran mismatch {worst:.2e}"
    assert n_const > 0
