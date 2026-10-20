"""Check the JAX calc_params_LY93 port (LY93_pdf.F90, Lewellen & Yoh 1993)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.ly93_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.LY93_pdf import calc_params_LY93, LY93_driver

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-2.5, 2.5, (NG, NZ))
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    return xm, xp2, Skx, mixt_frac


def test_f2py_oracle():
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, Skx, mf = _inputs(seed)
        f = fortran.calc_params_ly93(NZ, xm, xp2, Skx, mf)
        g = calc_params_LY93(xm, xp2, Skx, mf)
        worst = max(worst, max_abs_diff("ly93_pdf", tuple(g), tuple(f)))
    assert worst < 1e-11, f"calc_params_ly93 Fortran mismatch {worst:.2e}"


def test_driver_f2py():
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(5):
        wm = rng.uniform(-1, 1, (NG, NZ)); rtm = rng.uniform(0, 0.02, (NG, NZ)); thlm = rng.uniform(285, 305, (NG, NZ))
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); rtp2 = rng.uniform(1e-8, 1e-6, (NG, NZ)); thlp2 = rng.uniform(0.01, 1.0, (NG, NZ))
        # Skewness spanning both the 0.75 branch (|Sk|<=0.84) and the iterative branch.
        Skw = rng.uniform(-3, 3, (NG, NZ)); Skrt = rng.uniform(-1, 1, (NG, NZ)); Skthl = rng.uniform(-2, 2, (NG, NZ))
        f = fortran.ly93_driver(NZ, wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl)
        g = LY93_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl)
        worst = max(worst, max_abs_diff("ly93_pdf", tuple(g), tuple(f)))
    assert worst < 1e-9, f"ly93_driver Fortran mismatch {worst:.2e}"
