#!/usr/bin/env python3
"""test_calc_roots.py — validate the JAX calc_roots port (cubic/quadratic/cube_root).

Oracles (no Fortran build required):
  1. The roots ARE the roots — the polynomial residual `|a r^3 + b r^2 + c r + d|` must be ~machine-eps
     at every returned root (covers D<0 three-real, D=0 double, D>0 one-real-two-complex).
  2. Set-match vs `numpy.roots` (order-independent).
  3. cube_root vs `numpy.cbrt`.
Plus a `jax.grad` smoke check (smooth away from the determinant=0 / cube-root cusp), and an optional
bit-to-bit f2py `calc_roots` shadow when `clubb_f2py` is importable (SKIPs cleanly otherwise).
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

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


def test_differentiable():
    # d(real part of root1)/d(c) for a well-separated real-root quadratic — finite + FD-correct.
    a, b = 1.0, -3.0
    def f(c):
        return jnp.real(quadratic_solve(jnp.asarray(a), jnp.asarray(b), c)[..., 0])
    g = float(jax.grad(lambda c: f(c).sum())(jnp.asarray(2.0)))
    eps = 1e-6
    fd = (float(f(jnp.asarray(2.0 + eps))) - float(f(jnp.asarray(2.0 - eps)))) / (2 * eps)
    rel = abs(g - fd) / (abs(fd) + 1e-30)
    assert np.isfinite(g) and rel < 1e-4, f"grad {g} vs fd {fd} rel {rel:.2e}"
    print(f"  jax.grad(quadratic root) = {g:+.6e} (fd {fd:+.6e}, rel {rel:.1e})  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
        fn = clubb_f2py.f2py_cubic_solve
    except Exception as e:
        print(f"  f2py cubic_solve oracle: SKIP ({type(e).__name__})")
        return
    # Cases incl. a NEGATIVE-R cubic ([-1.5,0,1.5,-0.75]) where the Cardano args R±sqrt(D) are negative real —
    # this is the case that exposed the principal-vs-real cube-root branch bug (iter 71).
    cases = [(1.0, -6.0, 11.0, -6.0), (1.0, -4.0, 5.0, -2.0), (1.0, 0.0, 0.0, -1.0),
             (-2.0, 0.0, 3.0, 1.5), (-1.5, 0.0, 1.5, -0.75), (1.0, 0.0, -7.0, 6.0)]
    a = np.array([x[0] for x in cases]); b = np.array([x[1] for x in cases])
    c = np.array([x[2] for x in cases]); d = np.array([x[3] for x in cases])
    f = np.asarray(fn(a, b, c, d))          # (2, n, 3): [real, imag]
    ref = f[0] + 1j * f[1]
    got = np.asarray(cubic_solve(*[jnp.asarray(v) for v in (a, b, c, d)]))
    # gfortran's complex **(1/3) branch cut orders the conjugate pair / real root differently than this port
    # for D>0; both return the same correct root SET, so compare the sorted real parts (what the consumers use).
    err = np.max(np.abs(np.sort(np.real(got), axis=-1) - np.sort(np.real(ref), axis=-1)))
    # ~1e-9 near the [1,-4,5,-2] double root (1,1,2) is Cardano's intrinsic ill-conditioning, shared by both
    # implementations; the residual check below is the strict correctness guard.
    assert err < 1e-8, f"f2py cubic sorted-real mismatch {err:.2e}"
    # And the roots must actually satisfy each cubic.
    resid = np.max(np.abs(a[:, None] * got ** 3 + b[:, None] * got ** 2 + c[:, None] * got + d[:, None]))
    assert resid < 1e-11, f"cubic residual {resid:.2e}"
    print(f"  f2py cubic_solve: sorted-real bit-match {err:.2e} (incl. negative-R), residual {resid:.2e}  PASS")

    # quadratic_solve vs f2py_quadratic_solve — was only checked vs np.roots (a re-implementation); the f2py
    # oracle is the independent Fortran truth, incl. the complex-root branch. The JAX entry assumes a/=0, so keep
    # |a| away from 0 (a=0 is the Fortran's separate linear fallback). (iter 412)
    rng = np.random.default_rng(11)
    qworst = 0.0
    for _ in range(50):
        n = 8
        aa = rng.uniform(-3.0, 3.0, n); aa = np.where(np.abs(aa) < 0.3, 1.0, aa)
        bb = rng.uniform(-5.0, 5.0, n); cc = rng.uniform(-5.0, 5.0, n)   # mix of real & complex roots
        jr = np.asarray(quadratic_solve(jnp.asarray(aa), jnp.asarray(bb), jnp.asarray(cc)))   # complex (n,2)
        qre, qim = clubb_f2py.f2py_quadratic_solve(aa, bb, cc)
        qworst = max(qworst, float(np.max(np.abs(jr.real - np.asarray(qre)))),
                     float(np.max(np.abs(jr.imag - np.asarray(qim)))))
    assert qworst < 1e-12, f"f2py quadratic_solve mismatch {qworst:.2e}"
    print(f"  f2py quadratic_solve: bit-match over 50 cases (real + complex roots), worst {qworst:.2e}  PASS")


def main():
    print("test_calc_roots:")
    for t in (test_cube_root, test_quadratic, test_cubic, test_differentiable, test_f2py_oracle):
        t()
    print("All calc_roots checks PASSED")


if __name__ == "__main__":
    main()
