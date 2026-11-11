"""Check the JAX calc_responder_params port (new_pdf.F90, Griffin & Larson 2018)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_pdf as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_pdf import calc_responder_params

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-1.5, 1.5, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    F_x = rng.uniform(0.1, 0.9, (NG, NZ))
    mixt_frac = rng.uniform(0.3, 0.7, (NG, NZ))
    return xm, xp2, Skx, sgn, F_x, mixt_frac


def test_f2py_oracle():
    worst = 0.0
    n_zero = 0
    for seed in (11, 22, 33):
        xm, xp2, Skx, sgn, F_x, mf = _inputs(seed)
        if seed == 33:
            F_x[:, :2] = 0.0; n_zero += 1   # exercise the F_x=0 single-Gaussian branch
        f = fortran.calc_responder_params(NZ, xm, xp2, Skx, sgn, F_x, mf)
        g = calc_responder_params(xm, xp2, Skx, sgn, F_x, mf)
        worst = max(worst, max_abs_diff("responder_params", tuple(g), tuple(f)))
    assert worst < 1e-11, f"calc_responder_params Fortran mismatch {worst:.2e}"
    assert n_zero > 0
