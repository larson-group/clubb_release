#!/usr/bin/env python3
"""test_max_cubic_root.py — validate the JAX max_cubic_root port (adg1_adg2_3d_luhar_pdf.F90).

Largest real root of a cubic (with quadratic/linear fallbacks). No direct f2py oracle; validated by:
  1. Root property: the returned value satisfies a·r³ + b·r² + c·r + d ~ 0.
  2. Max-real-root: matches max(real(numpy.roots)) for true-cubic cases.
  3. The fallback branches (a~0 quadratic; a,b~0 linear).
  4. A finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import max_cubic_root


def test_root_property_and_max():
    rng = np.random.default_rng(3)
    worst_resid = 0.0
    worst_max = 0.0
    for _ in range(200):
        a, b, c, d = rng.uniform(-3, 3, 4)
        if abs(a) < 0.5:
            a = 0.5 * np.sign(a) if a != 0 else 0.5     # keep it a genuine cubic
        r = float(np.asarray(max_cubic_root(np.array([a]), np.array([b]), np.array([c]), np.array([d])))[0])
        worst_resid = max(worst_resid, abs(a * r ** 3 + b * r ** 2 + c * r + d))
        real_roots = np.real(np.roots([a, b, c, d]))[np.abs(np.imag(np.roots([a, b, c, d]))) < 1e-9]
        if real_roots.size:
            worst_max = max(worst_max, abs(r - real_roots.max()))
    assert worst_resid < 1e-9, f"root residual {worst_resid:.2e}"
    assert worst_max < 1e-7, f"max-real-root mismatch {worst_max:.2e}"
    print(f"  root property (resid {worst_resid:.1e}) + max real root vs np.roots ({worst_max:.1e})  PASS")


def test_fallbacks():
    # Quadratic fallback: a=0, b!=0 -> x^2 - 3x + 2 = 0 -> roots 1,2 -> max=2.
    r = float(np.asarray(max_cubic_root(np.array([0.0]), np.array([1.0]), np.array([-3.0]), np.array([2.0])))[0])
    assert abs(r - 2.0) < 1e-12, f"quadratic fallback: {r}"
    # Linear fallback: a=0, b=0 -> 2x - 4 = 0 -> x=2.
    rl = float(np.asarray(max_cubic_root(np.array([0.0]), np.array([0.0]), np.array([2.0]), np.array([-4.0])))[0])
    assert abs(rl - 2.0) < 1e-12, f"linear fallback: {rl}"
    print("  quadratic + linear fallbacks  PASS")


def test_differentiable():
    a = jnp.array([1.0, 2.0]); b = jnp.array([-1.0, 0.5]); c = jnp.array([-2.0, -3.0])
    g = np.asarray(jax.grad(lambda d: jnp.sum(max_cubic_root(a, b, c, d) ** 2))(jnp.array([0.5, -1.0])))
    assert np.isfinite(g).all(), "non-finite grad through max_cubic_root"
    print(f"  jax.grad through max_cubic_root: finite ({g.size} entries)  PASS")


def main():
    print("test_max_cubic_root:")
    for t in (test_root_property_and_max, test_fallbacks, test_differentiable):
        t()
    print("All max_cubic_root checks PASSED")


if __name__ == "__main__":
    main()
