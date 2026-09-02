#!/usr/bin/env python3
"""test_calculate_w_responder_params.py — validate the new-hybrid PDF leaf routines (new_hybrid_pdf.F90).

calculate_w_params (w-setter component params + mixt_frac) and calculate_responder_params (a responding
variable's component params from <w'x'>). Oracles:
  1. f2py bit-shadow vs f2py_calculate_w_params / f2py_calculate_responder_params. SKIPs if unbuilt.
  2. The |<w'x'>|=0 single-Gaussian invariant for the responder, and a finite jax.grad.
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

from clubb_jax.src.CLUBB_core.new_hybrid_pdf import calculate_w_params, calculate_responder_params

NG, NZ = 2, 8


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calculate_w/responder_params oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(13)
    worst_w = worst_r = 0.0
    for _ in range(4):
        sh = (NG, NZ)
        wm = rng.uniform(-1, 1, sh); wp2 = rng.uniform(0.05, 2, sh)
        Skw = rng.uniform(-2, 2, sh)
        F_w = rng.uniform(0.05, 0.95, sh); zeta_w = rng.uniform(0.0, 3.0, sh)

        fw = clubb_f2py.f2py_calculate_w_params(wm, wp2, Skw, F_w, zeta_w)
        gw = calculate_w_params(wm, wp2, Skw, F_w, zeta_w)
        for fi, gi in zip(fw, gw):
            worst_w = max(worst_w, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))

        # Responder: use the setter's mixt_frac, a nonzero covariance, and a responder variable.
        mixt_frac = np.asarray(fw[4])
        xm = rng.uniform(285, 305, sh); xp2 = rng.uniform(0.05, 1, sh)
        Skx = rng.uniform(-2, 2, sh); wpxp = rng.uniform(-0.5, 0.5, sh)
        fr = clubb_f2py.f2py_calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac)
        gr = calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac)
        for fi, gi in zip(fr, gr):
            worst_r = max(worst_r, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst_w < 1e-12, f"calculate_w_params f2py mismatch {worst_w:.2e}"
    assert worst_r < 1e-12, f"calculate_responder_params f2py mismatch {worst_r:.2e}"
    print(f"  f2py calculate_w_params: bit-match, worst {worst_w:.2e}  PASS")
    print(f"  f2py calculate_responder_params: bit-match, worst {worst_r:.2e}  PASS")


def test_responder_zero_cov_gaussian():
    sh = (NG, NZ)
    rng = np.random.default_rng(2)
    xm = rng.uniform(285, 305, sh); xp2 = rng.uniform(0.05, 1, sh)
    Skx = rng.uniform(-2, 2, sh); wp2 = rng.uniform(0.05, 2, sh)
    F_w = rng.uniform(0.05, 0.95, sh); mf = rng.uniform(0.2, 0.8, sh)
    wpxp = np.zeros(sh)                                  # |<w'x'>| = 0 → single Gaussian
    mu1, mu2, s1, s2, c1, c2 = (np.asarray(o) for o in
                                calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mf))
    assert np.allclose(mu1, xm) and np.allclose(mu2, xm), "zero-cov means not collapsed to xm"
    assert np.allclose(s1, xp2) and np.allclose(s2, xp2), "zero-cov variances not xp2"
    assert np.allclose(c1, 1.0) and np.allclose(c2, 1.0), "zero-cov coefs not 1"
    print("  responder |<w'x'>|=0 → single Gaussian  PASS")


def test_differentiable():
    sh = (NG, NZ)
    rng = np.random.default_rng(5)
    wm = rng.uniform(-1, 1, sh); wp2 = rng.uniform(0.05, 2, sh)
    Skw = rng.uniform(-2, 2, sh); F_w = rng.uniform(0.05, 0.95, sh); zeta = rng.uniform(0, 3, sh)
    g = np.asarray(jax.grad(lambda s: sum(jnp.sum(o ** 2) for o in
                   calculate_w_params(wm, wp2, s, F_w, zeta)))(jnp.asarray(Skw)))
    assert np.isfinite(g).all(), "non-finite grad through calculate_w_params"
    print(f"  jax.grad through calculate_w_params: finite ({g.size} entries)  PASS")


def main():
    print("test_calculate_w_responder_params:")
    for t in (test_f2py_oracle, test_responder_zero_cov_gaussian, test_differentiable):
        t()
    print("All calculate_w/responder_params checks PASSED")


if __name__ == "__main__":
    main()
