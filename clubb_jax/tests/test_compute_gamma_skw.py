#!/usr/bin/env python3
"""test_compute_gamma_skw.py — validate the JAX compute_gamma_Skw port (Skx_module.F90:compute_gamma_Skw).

gamma(Skw) = gamma_coefb + (gamma_coef - gamma_coefb)·exp(-1/2 (Skw/gamma_coefc)^2) when l_gamma_Skw and the
coefficients differ; else the constant gamma_coef. Implemented inline in advance_clubb_core_module.py but not as
a standalone tested function. Oracles:
  1. f2py bit-shadow vs f2py_compute_gamma_skw (both zm and zt outputs), for l_gamma_Skw on/off and a
     degenerate-coefficient case. SKIPs if clubb_f2py is unbuilt.
  2. Closed-form identity + the Skw=0 -> gamma_coef peak.
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

from clubb_jax.src.CLUBB_core.Skx_module import compute_gamma_Skw

NG, NZM = 2, 9
NZT = NZM - 1
NPARAMS = 102
I_GC, I_GB, I_GCF = 57, 58, 59   # 1-based igamma_coef / coefb / coefc


def _params(gc, gb, gcf):
    p = np.zeros((NG, NPARAMS))
    p[:, I_GC - 1] = gc
    p[:, I_GB - 1] = gb
    p[:, I_GCF - 1] = gcf
    return p


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py compute_gamma_skw oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    skw_zm = rng.uniform(-3, 3, (NG, NZM))
    skw_zt = rng.uniform(-3, 3, (NG, NZT))
    worst = 0.0
    configs = [(True, _params(0.32, 0.08, 1.2)),     # varying
               (False, _params(0.32, 0.08, 1.2)),    # constant
               (True, _params(0.25, 0.25, 1.0))]     # degenerate coefs -> constant gamma_coef
    for l_gamma, params in configs:
        f_zm, f_zt = clubb_f2py.f2py_compute_gamma_skw(l_gamma, skw_zm, params, skw_zt)
        g_zm = np.asarray(compute_gamma_Skw(skw_zm, params, l_gamma))
        g_zt = np.asarray(compute_gamma_Skw(skw_zt, params, l_gamma))
        worst = max(worst, np.max(np.abs(g_zm - np.asarray(f_zm))), np.max(np.abs(g_zt - np.asarray(f_zt))))
    assert worst < 1e-12, f"compute_gamma_skw f2py mismatch {worst:.2e}"
    print(f"  f2py compute_gamma_skw: bit-match (zm+zt, varying/constant/degenerate), worst {worst:.2e}  PASS")


def test_closed_form():
    gc, gb, gcf = 0.32, 0.08, 1.2
    params = _params(gc, gb, gcf)[:1]                 # single column to match skw's ngrdcol=1
    skw = np.array([[0.0, 1.0, -2.0]])
    got = np.asarray(compute_gamma_Skw(skw, params, True))
    ref = gb + (gc - gb) * np.exp(-0.5 * (skw / gcf) ** 2)
    assert np.max(np.abs(got - ref)) < 1e-13, "closed-form mismatch"
    # At Skw=0 the Gaussian is 1 -> gamma = gamma_coef.
    assert abs(got[0, 0] - gc) < 1e-13, "Skw=0 should give gamma_coef"
    print("  closed-form identity + Skw=0 peak (= gamma_coef)  PASS")


def test_differentiable():
    params = _params(0.32, 0.08, 1.2)
    def loss(skw):
        return jnp.sum(compute_gamma_Skw(skw, params, True) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(1).uniform(-3, 3, (NG, NZM)))))
    assert np.isfinite(g).all(), "non-finite grad through compute_gamma_Skw"
    print(f"  jax.grad through compute_gamma_Skw: finite ({g.size} entries)  PASS")


def main():
    print("test_compute_gamma_skw:")
    for t in (test_f2py_oracle, test_closed_form, test_differentiable):
        t()
    print("All compute_gamma_skw checks PASSED")


if __name__ == "__main__":
    main()
