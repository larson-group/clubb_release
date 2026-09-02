#!/usr/bin/env python3
"""test_luhar_params.py — validate the JAX calc_Luhar_params port (adg1_adg2_3d_luhar_pdf.F90).

Luhar PDF closure (Larson, Golaz & Cotton 2002): mixt_frac, big_m (M), small_m (m) from skewness and the w-x
covariance. Oracles:
  1. f2py bit-shadow vs f2py_calc_luhar_params (mixt_frac, big_m, small_m). SKIPs if clubb_f2py is unbuilt.
  2. Invariants: mixt_frac in [0,1]; small_m floored at 0.05 (varying x); constant-x limit (xp2<=tol ->
     mixt_frac=0.5, m=M=0); unskewed (Skx=0) -> mixt_frac=0.5.
  3. A finite jax.grad.
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
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_luhar_params oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    n_const = 0
    for seed in (11, 22, 33):
        Skx, wpxp, xp2 = _inputs(seed)
        if seed == 33:                          # force some constant-x points
            xp2[:, :2] = 1e-8; n_const += 1
        f = clubb_f2py.f2py_calc_luhar_params(Skx, wpxp, xp2, _X_TOL_SQD)
        g = calc_Luhar_params(Skx, wpxp, xp2, _X_TOL_SQD)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"calc_luhar_params f2py mismatch {worst:.2e}"
    assert n_const > 0
    print(f"  f2py calc_luhar_params: bit-match (mixt_frac/M/m incl. constant-x), worst {worst:.2e}  PASS")


def test_invariants():
    Skx, wpxp, xp2 = _inputs(5)
    mf, big_m, small_m = (np.asarray(x) for x in calc_Luhar_params(Skx, wpxp, xp2, _X_TOL_SQD))
    assert np.all(mf >= -1e-12) and np.all(mf <= 1.0 + 1e-12), "mixt_frac out of [0,1]"
    assert np.all(small_m >= 0.05 - 1e-12), "small_m below the 0.05 floor for varying x"
    # Constant-x limit.
    mf0, M0, m0 = (np.asarray(x) for x in
                   calc_Luhar_params(np.array([[1.0]]), np.array([[0.5]]), np.array([[1e-8]]), _X_TOL_SQD))
    assert abs(mf0[0, 0] - 0.5) < 1e-14 and M0[0, 0] == 0.0 and m0[0, 0] == 0.0, "constant-x limit"
    # Unskewed -> mixt_frac = 0.5.
    mfs = np.asarray(calc_Luhar_params(np.array([[0.0]]), np.array([[0.5]]), np.array([[1.0]]), _X_TOL_SQD)[0])
    assert abs(mfs[0, 0] - 0.5) < 1e-14, "Skx=0 should give mixt_frac=0.5"
    print("  invariants: mixt_frac in [0,1], small_m>=0.05, constant-x & unskewed limits  PASS")


def test_differentiable():
    Skx, wpxp, xp2 = _inputs(7)
    def loss(s):
        mf, M, m = calc_Luhar_params(s, wpxp, xp2, _X_TOL_SQD)
        return jnp.sum(mf ** 2 + M ** 2 + m ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(Skx)))
    assert np.isfinite(g).all(), "non-finite grad through calc_Luhar_params"
    print(f"  jax.grad through calc_Luhar_params: finite ({g.size} entries)  PASS")


def main():
    print("test_luhar_params:")
    for t in (test_f2py_oracle, test_invariants, test_differentiable):
        t()
    print("All calc_Luhar_params checks PASSED")


if __name__ == "__main__":
    main()
