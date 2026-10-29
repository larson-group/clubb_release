"""Check the JAX calc_coefs_wp2xp_semiimpl port (new_pdf.F90) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_hybrid_pdf as fortran_hybrid
import clubb_python.CLUBB_core.new_pdf as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import calc_coefs_wpxp2_semiimpl
from clubb_jax.src.CLUBB_core.new_pdf import calc_coefs_wp2xp_semiimpl, calc_coefs_wpxpyp_semiimpl

NG, NZ = 2, 6


def test_f2py_oracle():
    rng = np.random.default_rng(4)
    worst = 0.0
    n_red = 0
    for trial in range(20):
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
        sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
        mf = rng.uniform(0.1, 0.9, (NG, NZ)); F_w = rng.uniform(0.05, 1.0, (NG, NZ)); F_x = rng.uniform(0.05, 1.0, (NG, NZ))
        cw1 = rng.uniform(0.1, 2.0, (NG, NZ)); cw2 = rng.uniform(0.1, 2.0, (NG, NZ))
        cx1 = rng.uniform(0.1, 2.0, (NG, NZ)); cx2 = rng.uniform(0.1, 2.0, (NG, NZ))
        if trial % 4 == 0:   # reduced branch (zero sigma products)
            cw1 = np.zeros((NG, NZ)); cw2 = np.zeros((NG, NZ)); n_red += 1
        fc, ft = fortran.calc_coefs_wp2xp_semiimpl(NZ, wp2, xp2, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2)
        gc, gt = calc_coefs_wp2xp_semiimpl(wp2, xp2, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2)
        worst = max(worst, max_abs_diff("coefs_semiimpl", (gc, gt), (fc, ft)))
    assert worst < 1e-11, f"calc_coefs_wp2xp_semiimpl Fortran mismatch {worst:.2e}"
    assert n_red > 0


def test_wpxp2_f2py():
    rng = np.random.default_rng(5)
    worst = 0.0; n_zero = 0
    for trial in range(20):
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); wpxp = rng.uniform(-1, 1, (NG, NZ))
        mf = rng.uniform(0.1, 0.9, (NG, NZ)); F_w = rng.uniform(0.05, 1.0, (NG, NZ))
        cx1 = rng.uniform(0.1, 2.0, (NG, NZ)); cx2 = rng.uniform(0.1, 2.0, (NG, NZ))
        if trial % 4 == 0:
            F_w = np.zeros((NG, NZ)); n_zero += 1   # F_w=0 -> 0,0 branch
        fc, ft = fortran_hybrid.calc_coefs_wpxp2_semiimpl(NZ, NG, wp2, wpxp, mf, F_w, cx1, cx2)
        gc, gt = calc_coefs_wpxp2_semiimpl(wp2, wpxp, mf, F_w, cx1, cx2)
        worst = max(worst, max_abs_diff("coefs_semiimpl", (gc, gt), (fc, ft)))
    assert worst < 1e-11, f"calc_coefs_wpxp2_semiimpl Fortran mismatch {worst:.2e}"
    assert n_zero > 0


def test_wpxpyp_f2py():
    rng = np.random.default_rng(6)
    worst = 0.0; n_red = 0
    for trial in range(20):
        sh = (NG, NZ)
        wp2 = rng.uniform(0.05, 2.0, sh); xp2 = rng.uniform(0.05, 2.0, sh); yp2 = rng.uniform(0.05, 2.0, sh)
        wpxp = rng.uniform(-1, 1, sh); wpyp = rng.uniform(-1, 1, sh)
        sx = np.sign(wpxp); sx[sx == 0] = 1.0; sy = np.sign(wpyp); sy[sy == 0] = 1.0
        mf = rng.uniform(0.1, 0.9, sh); F_w = rng.uniform(0.05, 1.0, sh); F_x = rng.uniform(0.05, 1.0, sh); F_y = rng.uniform(0.05, 1.0, sh)
        cw1 = rng.uniform(0.1, 2, sh); cw2 = rng.uniform(0.1, 2, sh)
        cx1 = rng.uniform(0.1, 2, sh); cx2 = rng.uniform(0.1, 2, sh)
        cy1 = rng.uniform(0.1, 2, sh); cy2 = rng.uniform(0.1, 2, sh)
        if trial % 4 == 0:
            cx1 = np.zeros(sh); cx2 = np.zeros(sh); n_red += 1   # reduced branch
        args = (wp2, xp2, yp2, wpxp, wpyp, sx, sy, mf, F_w, F_x, F_y, cw1, cw2, cx1, cx2, cy1, cy2)
        fc, ft = fortran.calc_coefs_wpxpyp_semiimpl(NZ, *args)
        gc, gt = calc_coefs_wpxpyp_semiimpl(*args)
        worst = max(worst, max_abs_diff("coefs_semiimpl", (gc, gt), (fc, ft)))
    assert worst < 1e-11, f"calc_coefs_wpxpyp_semiimpl Fortran mismatch {worst:.2e}"
    assert n_red > 0
