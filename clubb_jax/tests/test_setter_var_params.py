#!/usr/bin/env python3
"""test_setter_var_params.py — validate the JAX calc_setter_var_params port (new_pdf.F90, Griffin & Larson 2018).

PDF component means/stdevs + mixture fraction for the setting variable. Oracles:
  1. f2py bit-shadow vs f2py_calc_setter_var_params (7 outputs). SKIPs if clubb_f2py is unbuilt.
  2. Moment reconstruction: the resulting binormal reproduces the overall mean xm and variance xp2.
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

from clubb_jax.src.CLUBB_core.new_pdf import calc_setter_var_params

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-2.5, 2.5, (NG, NZ))
    sgn = np.sign(rng.uniform(-1, 1, (NG, NZ))); sgn[sgn == 0] = 1.0
    F_x = rng.uniform(0.05, 0.9, (NG, NZ))
    zeta_x = rng.uniform(0.0, 2.0, (NG, NZ))
    return xm, xp2, Skx, sgn, F_x, zeta_x


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_setter_var_params oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, Skx, sgn, F_x, zeta = _inputs(seed)
        f = clubb_f2py.f2py_calc_setter_var_params(xm, xp2, Skx, sgn, F_x, zeta)
        g = calc_setter_var_params(xm, xp2, Skx, sgn, F_x, zeta)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"calc_setter_var_params f2py mismatch {worst:.2e}"
    print(f"  f2py calc_setter_var_params: bit-match (7 outputs) over 3 seeds, worst {worst:.2e}  PASS")


def test_moment_reconstruction():
    xm, xp2, Skx, sgn, F_x, zeta = _inputs(5)
    mu1, mu2, s1, s2, mf, c1, c2 = (np.asarray(x) for x in calc_setter_var_params(xm, xp2, Skx, sgn, F_x, zeta))
    # Overall mean.
    xm_rec = mf * mu1 + (1 - mf) * mu2
    assert np.max(np.abs(xm_rec - xm)) < 1e-10, "overall mean not reproduced"
    # Overall variance = a((mu1-xm)^2 + s1^2) + (1-a)((mu2-xm)^2 + s2^2) == xp2.
    xp2_rec = mf * ((mu1 - xm) ** 2 + s1 ** 2) + (1 - mf) * ((mu2 - xm) ** 2 + s2 ** 2)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    # coef * xp2 == sigma^2.
    assert np.max(np.abs(c1 * xp2 - s1 ** 2)) < 1e-12 and np.max(np.abs(c2 * xp2 - s2 ** 2)) < 1e-12
    print("  moment reconstruction: binormal reproduces overall mean & variance; coef·xp2 = sigma^2  PASS")


def test_differentiable():
    xm, xp2, Skx, sgn, F_x, zeta = _inputs(7)
    def loss(F):
        outs = calc_setter_var_params(xm, xp2, Skx, sgn, F, zeta)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(F_x)))
    assert np.isfinite(g).all(), "non-finite grad through calc_setter_var_params"
    print(f"  jax.grad through calc_setter_var_params: finite ({g.size} entries)  PASS")


def main():
    print("test_setter_var_params:")
    for t in (test_f2py_oracle, test_moment_reconstruction, test_differentiable):
        t()
    print("All calc_setter_var_params checks PASSED")


if __name__ == "__main__":
    main()
