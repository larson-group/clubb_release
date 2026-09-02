#!/usr/bin/env python3
"""test_binormal_moments.py — validate the JAX compute_mean_binormal / compute_variance_binormal ports.

compute_mean_binormal (pdf_utilities.F90) and compute_variance_binormal (pdf_utilities.F90:720) give the overall
mean and variance of a two-component (binormal) PDF variable. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_compute_mean_binormal / f2py_compute_variance_binormal. SKIPs if unbuilt.
  2. Monte-Carlo: draw samples from the binormal (mix two normals by mixt_frac) and check the formulas match
     the sample mean/variance to sampling tolerance.
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

from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, compute_variance_binormal


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py binormal-moments oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(3)
    worst_m = worst_v = 0.0
    for _ in range(50):
        mu1, mu2 = rng.uniform(-3, 3, 2)
        s1, s2 = rng.uniform(0.05, 2.0, 2)
        a = rng.uniform(0.0, 1.0)
        xm_ref = float(clubb_f2py.f2py_compute_mean_binormal(mu1, mu2, a))
        xm_j = float(compute_mean_binormal(mu1, mu2, a))
        v_ref = float(clubb_f2py.f2py_compute_variance_binormal(xm_ref, mu1, mu2, s1, s2, a))
        v_j = float(compute_variance_binormal(xm_j, mu1, mu2, s1, s2, a))
        worst_m = max(worst_m, abs(xm_j - xm_ref))
        worst_v = max(worst_v, abs(v_j - v_ref))
    assert worst_m < 1e-13 and worst_v < 1e-12, f"f2py mismatch mean {worst_m:.2e} var {worst_v:.2e}"
    print(f"  f2py mean/variance binormal: bit-match over 50 configs (mean {worst_m:.1e}, var {worst_v:.1e})  PASS")


def test_monte_carlo():
    rng = np.random.default_rng(11)
    mu1, mu2, s1, s2, a = 1.0, -2.0, 0.5, 1.3, 0.35
    n = 4_000_000
    pick1 = rng.random(n) < a
    samples = np.where(pick1, rng.normal(mu1, s1, n), rng.normal(mu2, s2, n))
    xm = float(compute_mean_binormal(mu1, mu2, a))
    xp2 = float(compute_variance_binormal(xm, mu1, mu2, s1, s2, a))
    assert abs(xm - samples.mean()) < 5e-3, f"mean {xm} vs MC {samples.mean()}"
    assert abs(xp2 - samples.var()) < 1e-2, f"variance {xp2} vs MC {samples.var()}"
    print(f"  Monte-Carlo: formula mean/variance match samples (xm={xm:.3f}, xp2={xp2:.3f})  PASS")


def test_differentiable():
    def loss(args):
        mu1, mu2, s1, s2, a = args
        xm = compute_mean_binormal(mu1, mu2, a)
        return compute_variance_binormal(xm, mu1, mu2, s1, s2, a)
    g = np.asarray(jax.grad(loss)(jnp.array([1.0, -2.0, 0.5, 1.3, 0.35])))
    assert np.isfinite(g).all(), "non-finite grad through binormal moments"
    print(f"  jax.grad through compute_variance_binormal: finite ({g.size} entries)  PASS")


def main():
    print("test_binormal_moments:")
    for t in (test_f2py_oracle, test_monte_carlo, test_differentiable):
        t()
    print("All binormal-moments checks PASSED")


if __name__ == "__main__":
    main()
