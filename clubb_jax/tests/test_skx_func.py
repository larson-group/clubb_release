#!/usr/bin/env python3
"""test_skx_func.py — validate the JAX Skx_func port (Skx_module.F90:Skx_func).

Skx = xp3 · (xp2 + iSkw_denom_coef·x_tol²)^(-3/2) — the sensitivity-reduced skewness. Ported and used throughout
the PDF closure, but had no dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_skx_func, passing the SAME clubb_params array to both (only iSkw_denom_coef=73 is
     read), so the comparison isolates the JAX transcription of the (xp2 + denom)^(-3/2) algebra. SKIPs if
     clubb_f2py is unbuilt.
  2. Closed-form identity + the symmetric (xp3=0 -> Skx=0) limit.
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

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func

NG, NZ = 2, 6
NPARAMS = 102
I_SKW_DENOM = 73   # 1-based iSkw_denom_coef
_DENOM_COEF = 0.15
_X_TOL = 1.0e-2


def _params():
    p = np.zeros((NG, NPARAMS))
    p[:, I_SKW_DENOM - 1] = _DENOM_COEF
    return p


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py Skx_func oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    params = _params()
    worst = 0.0
    for _ in range(10):
        xp2 = rng.uniform(0.0, 2.0, (NG, NZ))
        xp3 = rng.uniform(-1.0, 1.0, (NG, NZ))
        ref = np.asarray(clubb_f2py.f2py_skx_func(xp2, xp3, _X_TOL, params))
        got = np.asarray(Skx_func(xp2, xp3, _X_TOL, params))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"Skx_func f2py mismatch {worst:.2e}"
    print(f"  f2py Skx_func: bit-match over 10 configs, worst {worst:.2e}  PASS")


def test_closed_form():
    params = _params()
    xp2 = np.array([[0.5, 1.0]]); xp3 = np.array([[0.2, -0.3]])
    got = np.asarray(Skx_func(xp2, xp3, _X_TOL, params))
    denom = _DENOM_COEF * _X_TOL ** 2
    ref = xp3 * (xp2 + denom) ** (-1.5)
    assert np.max(np.abs(got - ref)) < 1e-12, "closed-form mismatch"
    # xp3 = 0 -> Skx = 0.
    z = np.asarray(Skx_func(xp2, np.zeros_like(xp3), _X_TOL, params))
    assert np.all(z == 0.0), "xp3=0 should give Skx=0"
    print("  closed-form identity + xp3=0 limit  PASS")


def test_differentiable():
    params = _params()
    xp2 = jnp.asarray(np.random.default_rng(1).uniform(0.1, 2.0, (NG, NZ)))
    def loss(xp3):
        return jnp.sum(Skx_func(xp2, xp3, _X_TOL, params) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(2).uniform(-1, 1, (NG, NZ)))))
    assert np.isfinite(g).all(), "non-finite grad through Skx_func"
    print(f"  jax.grad through Skx_func: finite ({g.size} entries)  PASS")


def test_LG_2005_ansatz_f2py():
    """LG_2005_ansatz (Skx_module.F90 — the Larson-Golaz 2005 skewness ansatz: Skx from Skw/wpxp/wp2/xp2/sigma_sqd_w/
    beta) and its <x'^3> wrapper xp3_LG_2005_ansatz vs the f2py oracle. Both previously untested. The ansatz's clip
    carries a ~1e-11 relative FP difference; checked with a relative tolerance. (iter 422)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py LG_2005_ansatz oracle: SKIP ({type(e).__name__})")
        return
    from clubb_jax.src.CLUBB_core.Skx_module import LG_2005_ansatz, xp3_LG_2005_ansatz
    from clubb_jax.src.CLUBB_core.parameter_indices import ibeta
    params = _params(); params[:, ibeta - 1] = 2.0       # f2py xp3 extracts beta from clubb_params; JAX takes it explicitly
    beta = np.ascontiguousarray(params[:, ibeta - 1])
    rng = np.random.default_rng(2)
    w_lg = w_xp3 = 0.0
    for _ in range(30):
        Skw = rng.uniform(-3, 3, (NG, NZ)); wpxp = rng.uniform(-1, 1, (NG, NZ))
        wp2 = rng.uniform(1e-2, 2, (NG, NZ)); xp2 = rng.uniform(1e-4, 1, (NG, NZ)); ssw = rng.uniform(0, 0.6, (NG, NZ))
        jl = np.asarray(LG_2005_ansatz(jnp.asarray(Skw), jnp.asarray(wpxp), jnp.asarray(wp2), jnp.asarray(xp2),
                                       jnp.asarray(ssw), jnp.asarray(beta)[:, None], _X_TOL))
        fl = np.asarray(clubb_f2py.f2py_lg_2005_ansatz(Skw, wpxp, wp2, xp2, beta, ssw, _X_TOL))
        w_lg = max(w_lg, float(np.max(np.abs(jl - fl) / np.maximum(np.abs(fl), 1.0))))
        jx = np.asarray(xp3_LG_2005_ansatz(jnp.asarray(Skw), jnp.asarray(wpxp), jnp.asarray(wp2), jnp.asarray(xp2),
                                           jnp.asarray(ssw), jnp.asarray(beta)[:, None], jnp.asarray(params), _X_TOL))
        fx = np.asarray(clubb_f2py.f2py_xp3_lg_2005_ansatz(Skw, wpxp, wp2, xp2, ssw, params, _X_TOL))
        w_xp3 = max(w_xp3, float(np.max(np.abs(jx - fx) / np.maximum(np.abs(fx), 1.0))))
    assert w_lg < 1e-9, f"LG_2005_ansatz f2py rel mismatch {w_lg:.2e}"
    assert w_xp3 < 1e-9, f"xp3_LG_2005_ansatz f2py rel mismatch {w_xp3:.2e}"
    print(f"  f2py LG_2005_ansatz + xp3_LG_2005_ansatz (30 cases): rel-match, worst {max(w_lg, w_xp3):.2e}  PASS")


def main():
    print("test_skx_func:")
    for t in (test_f2py_oracle, test_closed_form, test_LG_2005_ansatz_f2py, test_differentiable):
        t()
    print("All Skx_func checks PASSED")


if __name__ == "__main__":
    main()
