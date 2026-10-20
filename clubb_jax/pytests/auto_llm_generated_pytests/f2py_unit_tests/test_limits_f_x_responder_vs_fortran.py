"""Check the JAX calc_limits_F_x_responder + sort_roots ports (new_pdf.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_pdf import calc_limits_F_x_responder

NG, NZ = 2, 6
MAX_POS, MAX_NEG = 4.0, 4.0


def _inputs(seed):
    rng = np.random.default_rng(seed)
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    Skx = rng.uniform(-3, 3, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    return mixt_frac, Skx, sgn


def test_f2py_oracle():
    worst = 0.0
    for seed in (11, 22, 33):
        mf, Skx, sgn = _inputs(seed)
        mp = np.full((NG, NZ), MAX_POS); mn = np.full((NG, NZ), MAX_NEG)
        f_out = fortran.calc_limits_f_x_responder(NZ, mf, Skx, sgn, mp, mn)
        g_out = calc_limits_F_x_responder(mf, Skx, sgn, mp, mn)
        worst = max(worst, max_abs_diff("calc_limits_F_x_responder", tuple(g_out), tuple(f_out)))
    assert worst < 1e-11, f"calc_limits_f_x_responder Fortran mismatch {worst:.2e}"
