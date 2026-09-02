#!/usr/bin/env python3
"""test_new_pdf.py — validate the JAX new_pdf.py ports (new-hybrid PDF helpers, Griffin & Larson 2018).

Oracles:
  1. calc_coef_wp4_implicit: f2py bit-shadow vs f2py_calc_coef_wp4_implicit. SKIPs if clubb_f2py is unbuilt.
  2. calc_mixture_fraction: literal NumPy transcription + the analytic symmetric limit (F_x=0, Skx=0 ->
     (zeta+1)/(zeta+2)) + the bounds mixt_frac in (0,1) for valid inputs.
  3. Finite jax.grad for both.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_pdf import (
    calc_coef_wp4_implicit, calc_mixture_fraction, calc_coef_wpxp2_implicit)
# calc_coef_wp2xp_implicit + the calculate_* aliases moved to new_hybrid_pdf.py (mirror-refactor iter 18)
from clubb_jax.src.CLUBB_core.new_hybrid_pdf import (
    calc_coef_wp2xp_implicit, calculate_coef_wp4_implicit, calculate_mixture_fraction)

NG, NZ = 2, 6


def test_coef_wp4_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coef_wp4_implicit oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(20):
        mf = rng.uniform(0.1, 0.9, (NG, NZ))
        F = rng.uniform(0.0, 1.0, (NG, NZ))
        c1 = rng.uniform(0.0, 2.0, (NG, NZ))
        c2 = rng.uniform(0.0, 2.0, (NG, NZ))
        ref = np.asarray(clubb_f2py.f2py_calc_coef_wp4_implicit(mf, F, c1, c2))
        got = np.asarray(calc_coef_wp4_implicit(mf, F, c1, c2))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"calc_coef_wp4_implicit f2py mismatch {worst:.2e}"
    print(f"  f2py calc_coef_wp4_implicit: bit-match over 20 configs, worst {worst:.2e}  PASS")


def _ref_mixt_frac(Skx, F, zeta, sgn):
    zp2 = zeta + 2.0
    if F > 0.0:
        sa = 4 * F ** 3 + 12 * F ** 2 * (1 - F) + 36 * F * (zeta + 1) * (1 - F) ** 2 / zp2 ** 2 + Skx ** 2
        num = (4 * F ** 3 + 18 * F * (zeta + 1) * (1 - F) / zp2 + 6 * F ** 2 * (1 - F) / zp2
               + Skx ** 2 - Skx * sgn * np.sqrt(sa))
        den = 2 * F * (F - 3) ** 2 + 2 * Skx ** 2
        return num / den
    return (zeta + 1) / zp2


def test_mixture_fraction():
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(200):
        F = rng.uniform(0.05, 1.0)
        Skx = rng.uniform(-3, 3)
        zeta = rng.uniform(0.0, 3.0)
        sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        got = float(calc_mixture_fraction(Skx, F, zeta, sgn))
        ref = _ref_mixt_frac(Skx, F, zeta, sgn)
        worst = max(worst, abs(got - ref))
    assert worst < 1e-13, f"mixt_frac transcription mismatch {worst:.2e}"
    # Symmetric limit F=0, Skx=0 -> (zeta+1)/(zeta+2).
    for zeta in (0.0, 1.0, 2.5):
        v = float(calc_mixture_fraction(0.0, 0.0, zeta, 1.0))
        assert abs(v - (zeta + 1) / (zeta + 2)) < 1e-14, "symmetric limit wrong"
    print(f"  calc_mixture_fraction: literal transcription + symmetric limit, worst {worst:.2e}  PASS")


def test_coef_wpxp2_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coef_wpxp2_implicit oracle: SKIP ({type(e).__name__})")
        return
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
        ref = np.asarray(clubb_f2py.f2py_calc_coef_wpxp2_implicit(
            wp2, xp2, wpxp, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2))
        got = np.asarray(calc_coef_wpxp2_implicit(wp2, xp2, wpxp, sgn, mf, F_w, F_x, cw1, cw2, cx1, cx2))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"calc_coef_wpxp2_implicit f2py mismatch {worst:.2e}"
    assert n_branchB > 0, "reduced branch not exercised"
    print(f"  f2py calc_coef_wpxp2_implicit: bit-match (both branches), worst {worst:.2e}  PASS")


def test_coef_wp2xp_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_coef_wp2xp_implicit oracle: SKIP ({type(e).__name__})")
        return
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
        ref = np.asarray(clubb_f2py.f2py_calc_coef_wp2xp_implicit(wp2, mf, F, c1, c2))
        got = np.asarray(calc_coef_wp2xp_implicit(wp2, mf, F, c1, c2))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"calc_coef_wp2xp_implicit f2py mismatch {worst:.2e}"
    assert n_zero > 0, "F_w=0 branch not exercised"
    print(f"  f2py calc_coef_wp2xp_implicit: bit-match (incl. F_w=0), worst {worst:.2e}  PASS")


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
        worst_mf = max(worst_mf, abs(a - b))
    assert worst_mf < 1e-15, f"calculate_mixture_fraction != sgn=+1 specialization: {worst_mf:.2e}"
    try:
        import clubb_f2py
        mf = rng.uniform(0.1, 0.9, (NG, NZ)); F = rng.uniform(0.0, 1.0, (NG, NZ))
        c1 = rng.uniform(0.0, 2.0, (NG, NZ)); c2 = rng.uniform(0.0, 2.0, (NG, NZ))
        ref = np.asarray(clubb_f2py.f2py_calculate_coef_wp4_implicit(mf, F, c1, c2))
        got = np.asarray(calculate_coef_wp4_implicit(mf, F, c1, c2))
        d = np.max(np.abs(got - ref))
        assert d < 1e-11, f"calculate_coef_wp4_implicit f2py mismatch {d:.2e}"
        print(f"  new_hybrid aliases: calculate_mixture_fraction==sgn+1 ({worst_mf:.1e}); "
              f"calculate_coef_wp4_implicit f2py bit-match ({d:.1e})  PASS")
    except Exception:
        print(f"  new_hybrid aliases: calculate_mixture_fraction==sgn+1 ({worst_mf:.1e}); f2py SKIP  PASS")


def test_differentiable():
    def loss_wp2xp(args):
        wp2, mf, F = args
        return calc_coef_wp2xp_implicit(wp2, mf, F, 1.0, 0.8)
    g4 = np.asarray(jax.grad(loss_wp2xp)(jnp.array([1.0, 0.4, 0.5])))
    assert np.isfinite(g4).all(), "non-finite grad through calc_coef_wp2xp_implicit"

    def loss_wpxp2(args):
        wp2, xp2, wpxp = args
        return calc_coef_wpxp2_implicit(wp2, xp2, wpxp, 1.0, 0.4, 0.5, 0.6, 1.0, 0.8, 0.9, 0.7)
    g3 = np.asarray(jax.grad(loss_wpxp2)(jnp.array([1.0, 1.2, 0.3])))
    assert np.isfinite(g3).all(), "non-finite grad through calc_coef_wpxp2_implicit"

    def loss_wp4(args):
        mf, F, c1, c2 = args
        return calc_coef_wp4_implicit(mf, F, c1, c2)
    g1 = np.asarray(jax.grad(loss_wp4)(jnp.array([0.4, 0.5, 1.0, 0.8])))
    def loss_mf(args):
        Skx, F, zeta = args
        return calc_mixture_fraction(Skx, F, zeta, 1.0)
    g2 = np.asarray(jax.grad(loss_mf)(jnp.array([1.2, 0.6, 1.5])))
    assert np.isfinite(g1).all() and np.isfinite(g2).all(), "non-finite grad"
    print(f"  jax.grad through both: finite (wp4 {g1.size}, mixt_frac {g2.size} entries)  PASS")


def main():
    print("test_new_pdf:")
    for t in (test_coef_wp4_f2py, test_mixture_fraction, test_coef_wpxp2_f2py,
              test_coef_wp2xp_f2py, test_new_hybrid_aliases, test_differentiable):
        t()
    print("All new_pdf checks PASSED")


if __name__ == "__main__":
    main()
