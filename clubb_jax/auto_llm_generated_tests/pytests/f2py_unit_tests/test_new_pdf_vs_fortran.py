"""Check the JAX new_pdf.py ports (new-hybrid PDF helpers, Griffin & Larson 2018)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.new_hybrid_pdf as fortran_hybrid
import clubb_python.CLUBB_core.new_pdf as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import (
    calc_coef_wp2xp_implicit,
    calculate_coef_wp4_implicit,
    calculate_mixture_fraction,
)
from clubb_jax.src.CLUBB_core.new_pdf import (
    calc_coef_wp4_implicit,
    calc_mixture_fraction,
    calc_coef_wpxp2_implicit,
)

NG, NZ = 2, 6


def test_coef_wp4_f2py():
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(20):
        mf = rng.uniform(0.1, 0.9, (NG, NZ))
        F = rng.uniform(0.0, 1.0, (NG, NZ))
        c1 = rng.uniform(0.0, 2.0, (NG, NZ))
        c2 = rng.uniform(0.0, 2.0, (NG, NZ))
        ref = fortran.calc_coef_wp4_implicit(NZ, mf, F, c1, c2)
        got = np.asarray(calc_coef_wp4_implicit(mf, F, c1, c2))
        worst = max(worst, max_abs_diff("new_pdf", got, ref))
    assert worst < 1e-11, f"calc_coef_wp4_implicit Fortran mismatch {worst:.2e}"


def test_coef_wpxp2_f2py():
    rng = np.random.default_rng(8)
    worst = 0.0
    n_branchB = 0
    for trial in range(30):
        mf = rng.uniform(0.1, 0.9, (NG, NZ))
        F_w = rng.uniform(0.05, 1.0, (NG, NZ)); F_x = rng.uniform(0.05, 1.0, (NG, NZ))
        wp2 = rng.uniform(0.1, 2.0, (NG, NZ)); xp2 = rng.uniform(0.1, 2.0, (NG, NZ))
        wpxp = rng.uniform(-1.0, 1.0, (NG, NZ)); sgn = np.sign(wpxp); sgn[sgn == 0] = 1.0
        cw1 = rng.uniform(0.1, 2.0, (NG, NZ)); cw2 = rng.uniform(0.1, 2.0, (NG, NZ))
        cx1 = rng.uniform(0.1, 2.0, (NG, NZ)); cx2 = rng.uniform(0.1, 2.0, (NG, NZ))
        if trial % 3 == 0:   # trigger the reduced branch (zero sigma products)
            cw1 = np.zeros((NG, NZ)); cw2 = np.zeros((NG, NZ)); n_branchB += 1
        ref = (fortran.calc_coef_wpxp2_implicit(
            NZ, wp2, xp2, wpxp, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2))
        got = np.asarray(calc_coef_wpxp2_implicit(wp2, xp2, wpxp, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2))
        worst = max(worst, max_abs_diff("new_pdf", got, ref))
    assert worst < 1e-11, f"calc_coef_wpxp2_implicit Fortran mismatch {worst:.2e}"
    assert n_branchB > 0, "reduced branch not exercised"


def test_coef_wp2xp_f2py():
    rng = np.random.default_rng(9)
    worst = 0.0
    n_zero = 0
    for trial in range(20):
        mf = rng.uniform(0.1, 0.9, (NG, NZ))
        F = rng.uniform(0.05, 1.0, (NG, NZ))
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ))
        c1 = rng.uniform(0.0, 2.0, (NG, NZ)); c2 = rng.uniform(0.0, 2.0, (NG, NZ))
        if trial % 4 == 0:   # exercise the F_w = 0 -> 0 branch
            F = np.zeros((NG, NZ)); n_zero += 1
        ref = fortran_hybrid.calc_coef_wp2xp_implicit(wp2, mf, F, c1, c2)
        got = np.asarray(calc_coef_wp2xp_implicit(wp2, mf, F, c1, c2))
        worst = max(worst, max_abs_diff("new_pdf", got, ref))
    assert worst < 1e-11, f"calc_coef_wp2xp_implicit Fortran mismatch {worst:.2e}"
    assert n_zero > 0, "F_w=0 branch not exercised"


def test_new_hybrid_aliases():
    """The new_hybrid_pdf.F90 calculate_* forms: calculate_coef_wp4_implicit has its own f2py oracle (identical
    formula); calculate_mixture_fraction is the sgn=+1 specialization (no oracle -> equivalence + analytic)."""
    rng = np.random.default_rng(13)
    # calculate_mixture_fraction == calc_mixture_fraction(., ., ., +1)
    worst_mf = 0.0
    for _ in range(100):
        F, Skw, zeta = rng.uniform(0.05, 1.0), rng.uniform(-3, 3), rng.uniform(0, 3)
        a = float(calculate_mixture_fraction(Skw, F, zeta))
        b = float(calc_mixture_fraction(Skw, F, zeta, 1.0))
        worst_mf = max(worst_mf, max_abs_diff("calculate_mixture_fraction", a, b))
    assert worst_mf < 1e-15, f"calculate_mixture_fraction != sgn=+1 specialization: {worst_mf:.2e}"
    mf = rng.uniform(0.1, 0.9, (NG, NZ)); F = rng.uniform(0.0, 1.0, (NG, NZ))
    c1 = rng.uniform(0.0, 2.0, (NG, NZ)); c2 = rng.uniform(0.0, 2.0, (NG, NZ))
    ref = fortran_hybrid.calculate_coef_wp4_implicit(mf, F, c1, c2)
    got = np.asarray(calculate_coef_wp4_implicit(mf, F, c1, c2))
    d = max_abs_diff("new_pdf", got, ref)
    assert d < 1e-11, f"calculate_coef_wp4_implicit Fortran mismatch {d:.2e}"
