"""Check the JAX new_tsdadg_pdf.F90 ports (calc_L_x_Skx_fnc, calc_setter_parameters,
tsdadg_pdf_driver) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_tsdadg_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_tsdadg_pdf import (
    calc_L_x_Skx_fnc,
    calc_setter_parameters,
    tsdadg_pdf_driver,
)


def test_f2py_oracle():
    # The f2py wrapper is scalar (the Fortran routine is per-level); call it element-wise.
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(60):
        Skx = float(rng.uniform(-3, 3)); sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        l1 = float(rng.uniform(2 / 3, 1.0)); l2 = float(rng.uniform(0.0, 1.0))
        f1, f2 = fortran.calc_l_x_skx_fnc(Skx, sgn, l1, l2)
        g1, g2 = calc_L_x_Skx_fnc(Skx, sgn, l1, l2)
        worst = max(worst, max_abs_diff("calc_l_x_skx_fnc", (g1, g2), (f1, f2)))
    assert worst < 1e-13, f"calc_l_x_skx_fnc Fortran mismatch {worst:.2e}"


def test_setter_f2py():
    # The f2py wrapper is scalar (per-level); call it element-wise.
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(60):
        xm = float(rng.uniform(-2, 2)); xp2 = float(rng.uniform(0.05, 2.0))
        Skx = float(rng.uniform(-2.5, 2.5)); sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        L1 = float(rng.uniform(0.05, 0.5)); L2 = float(rng.uniform(0.05, 0.5))
        f = fortran.calc_setter_parameters_tsdadg(xm, xp2, Skx, sgn, L1, L2)
        g = calc_setter_parameters(xm, xp2, Skx, sgn, L1, L2)
        worst = max(worst, max_abs_diff("new_tsdadg_pdf", tuple(np.asarray(x) for x in g), tuple(np.asarray(x) for x in f)))
    assert worst < 1e-11, f"calc_setter_parameters_tsdadg Fortran mismatch {worst:.2e}"


def test_driver_f2py():
    rng = np.random.default_rng(13)
    worst = 0.0
    ng, nz = 2, 8
    for _ in range(4):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(0, 0.02, (ng, nz)); thlm = rng.uniform(285, 305, (ng, nz))
        wp2 = rng.uniform(0.05, 2, (ng, nz)); rtp2 = rng.uniform(1e-8, 1e-6, (ng, nz)); thlp2 = rng.uniform(0.01, 1, (ng, nz))
        # Skewness with varying magnitudes so each of w/rt/thl is the setter at some points.
        Skw = rng.uniform(-3, 3, (ng, nz)); Skrt = rng.uniform(-3, 3, (ng, nz)); Skthl = rng.uniform(-3, 3, (ng, nz))
        wprtp = rng.uniform(-1, 1, (ng, nz)); wpthlp = rng.uniform(-1, 1, (ng, nz))
        f = fortran.tsdadg_pdf_driver(nz, wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp)
        g = tsdadg_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp)
        worst = max(worst, max_abs_diff("new_tsdadg_pdf", tuple(np.asarray(x) for x in g), tuple(np.asarray(x) for x in f)))
    assert worst < 1e-11, f"tsdadg_pdf_driver Fortran mismatch {worst:.2e}"
