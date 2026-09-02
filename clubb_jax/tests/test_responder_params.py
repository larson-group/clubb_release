#!/usr/bin/env python3
"""test_responder_params.py — validate the JAX calc_responder_params port (new_pdf.F90, Griffin & Larson 2018).

PDF component means/variances for a variable responding to the PDF set by another variable (uses the setter's
mixt_frac). Oracles:
  1. f2py bit-shadow vs f2py_calc_responder_params (6 outputs), incl. the F_x=0 single-Gaussian branch.
     SKIPs if clubb_f2py is unbuilt.
  2. Moment reconstruction: the binormal reproduces the overall mean xm and variance xp2 (using the signed
     component variances — responder variances can be negative, a known feature).
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
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_responder_params oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    n_zero = 0
    for seed in (11, 22, 33):
        xm, xp2, Skx, sgn, F_x, mf = _inputs(seed)
        if seed == 33:
            F_x[:, :2] = 0.0; n_zero += 1   # exercise the F_x=0 single-Gaussian branch
        f = clubb_f2py.f2py_calc_responder_params(xm, xp2, Skx, sgn, F_x, mf)
        g = calc_responder_params(xm, xp2, Skx, sgn, F_x, mf)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"calc_responder_params f2py mismatch {worst:.2e}"
    assert n_zero > 0
    print(f"  f2py calc_responder_params: bit-match (6 outputs incl. F_x=0), worst {worst:.2e}  PASS")


def test_moment_reconstruction():
    xm, xp2, Skx, sgn, F_x, mf = _inputs(5)
    mu1, mu2, s1sq, s2sq, c1, c2 = (np.asarray(x) for x in calc_responder_params(xm, xp2, Skx, sgn, F_x, mf))
    xm_rec = mf * mu1 + (1 - mf) * mu2
    assert np.max(np.abs(xm_rec - xm)) < 1e-10, "overall mean not reproduced"
    xp2_rec = mf * ((mu1 - xm) ** 2 + s1sq) + (1 - mf) * ((mu2 - xm) ** 2 + s2sq)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    print("  moment reconstruction: binormal reproduces overall mean & variance (signed comp. variances)  PASS")


def test_differentiable():
    xm, xp2, Skx, sgn, F_x, mf = _inputs(7)
    def loss(s):
        outs = calc_responder_params(xm, xp2, s, sgn, F_x, mf)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(Skx)))
    assert np.isfinite(g).all(), "non-finite grad through calc_responder_params"
    print(f"  jax.grad through calc_responder_params: finite ({g.size} entries)  PASS")


def main():
    print("test_responder_params:")
    for t in (test_f2py_oracle, test_moment_reconstruction, test_differentiable):
        t()
    print("All calc_responder_params checks PASSED")


if __name__ == "__main__":
    main()
