#!/usr/bin/env python3
"""validate the JAX calc_roots port (cubic/quadratic/cube_root)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.calc_roots import cubic_solve, quadratic_solve, cube_root


def _residual_cubic(a, b, c, d, roots):
    r = roots
    return np.abs(a[..., None] * r ** 3 + b[..., None] * r ** 2 + c[..., None] * r + d[..., None])


def _setmatch(got, ref, tol):
    """Order-independent comparison of two equal-length root multisets (complex)."""
    got = sorted(np.asarray(got).tolist(), key=lambda z: (round(z.real, 9), round(z.imag, 9)))
    ref = sorted(np.asarray(ref).tolist(), key=lambda z: (round(z.real, 9), round(z.imag, 9)))
    return max(abs(g - r) / (abs(r) + 1.0) for g, r in zip(got, ref)) < tol


def test_cube_root():
    x = np.array([-27.0, -8.0, -1.0, -1e-9, 0.0, 1e-9, 1.0, 8.0, 64.0])
    got = np.asarray(cube_root(jnp.asarray(x)))
    ref = np.cbrt(x)
    assert np.allclose(got, ref, atol=1e-14, rtol=1e-13), f"cube_root max err {np.max(np.abs(got-ref)):.2e}"
    print(f"  cube_root vs np.cbrt: max err {np.max(np.abs(got - ref)):.2e}  PASS")


def test_quadratic():
    # columns: two real roots; double root; complex-conjugate pair
    a = np.array([1.0, 1.0,  2.0])
    b = np.array([-3.0, -4.0, 2.0])
    c = np.array([2.0, 4.0,  5.0])    # x^2-3x+2=(x-1)(x-2); x^2-4x+4=(x-2)^2; 2x^2+2x+5 (D<0)
    roots = np.asarray(quadratic_solve(jnp.asarray(a), jnp.asarray(b), jnp.asarray(c)))
    res = np.abs(a[:, None] * roots ** 2 + b[:, None] * roots + c[:, None])
    assert res.max() < 1e-12, f"quadratic residual {res.max():.2e}"
    for i in range(3):
        assert _setmatch(roots[i], np.roots([a[i], b[i], c[i]]), 1e-10), f"quadratic col {i} mismatch"
    print(f"  quadratic_solve: max residual {res.max():.2e}, set-match np.roots  PASS")


def test_cubic():
    # D<0 (3 distinct real): (x-1)(x-2)(x-3) = x^3-6x^2+11x-6
    # D=0 (double root):     (x-1)^2(x-2)   = x^3-4x^2+5x-2
    # D>0 (1 real, 2 cplx):  (x-1)(x^2+x+1) = x^3-0x^2+0x-1  -> x^3-1
    a = np.array([1.0, 1.0, 1.0])
    b = np.array([-6.0, -4.0, 0.0])
    c = np.array([11.0, 5.0, 0.0])
    d = np.array([-6.0, -2.0, -1.0])
    roots = np.asarray(cubic_solve(*[jnp.asarray(v) for v in (a, b, c, d)]))
    res = _residual_cubic(a, b, c, d, roots)
    assert res.max() < 1e-9, f"cubic residual {res.max():.2e}"
    for i in range(3):
        assert _setmatch(roots[i], np.roots([a[i], b[i], c[i], d[i]]), 1e-7), f"cubic col {i} mismatch"
    print(f"  cubic_solve: max residual {res.max():.2e}, set-match np.roots  PASS")
