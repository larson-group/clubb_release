"""Check the JAX calc_setter_var_params port (new_pdf.F90, Griffin & Larson 2018)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_pdf import calc_setter_var_params

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-2.5, 2.5, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    F_x = rng.uniform(0.05, 0.9, (NG, NZ))
    zeta_x = rng.uniform(0.0, 2.0, (NG, NZ))
    return xm, xp2, Skx, sgn, F_x, zeta_x


def test_f2py_oracle():
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, Skx, sgn, F_x, zeta = _inputs(seed)
        f = fortran.calc_setter_var_params(NZ, xm, xp2, Skx, sgn, F_x, zeta)
        g = calc_setter_var_params(xm, xp2, Skx, sgn, F_x, zeta)
        worst = max(worst, max_abs_diff("setter_var_params", tuple(g), tuple(f)))
    assert worst < 1e-11, f"calc_setter_var_params Fortran mismatch {worst:.2e}"
