"""Check the JAX calc_roots port (cubic/quadratic/cube_root) against Fortran."""

import jax.numpy as jnp
import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.calc_roots import cubic_solve, quadratic_solve


def test_f2py_oracle():
    import clubb_f2py
    fn = clubb_f2py.f2py_cubic_solve
    # Cases incl. a NEGATIVE-R cubic ([-1.5,0,1.5,-0.75]) where the Cardano args R±sqrt(D) are negative real —
    # this is the case that exposed the principal-vs-real cube-root branch bug.
    cases = [(1.0, -6.0, 11.0, -6.0), (1.0, -4.0, 5.0, -2.0), (1.0, 0.0, 0.0, -1.0),
             (-2.0, 0.0, 3.0, 1.5), (-1.5, 0.0, 1.5, -0.75), (1.0, 0.0, -7.0, 6.0)]
    a = np.array([x[0] for x in cases]); b = np.array([x[1] for x in cases])
    c = np.array([x[2] for x in cases]); d = np.array([x[3] for x in cases])
    f = np.asarray(fn(a, b, c, d))          # (2, n, 3): [real, imag]
    ref = f[0] + 1j * f[1]
    got = np.asarray(cubic_solve(*[jnp.asarray(v) for v in (a, b, c, d)]))
    # gfortran's complex **(1/3) branch cut orders the conjugate pair / real root differently than this port
    # for D>0; both return the same correct root SET, so compare the sorted real parts (what the consumers use).
    err = max_abs_diff("calc_roots", np.sort(np.real(got), axis=-1), np.sort(np.real(ref), axis=-1))
    # ~1e-9 near the [1,-4,5,-2] double root (1,1,2) is Cardano's intrinsic ill-conditioning, shared by both
    # implementations; the residual check below is the strict correctness guard.
    assert err < 1e-8, f"f2py cubic sorted-real mismatch {err:.2e}"
    # And the roots must actually satisfy each cubic.
    resid = np.max(np.abs(a[:, None] * got ** 3 + b[:, None] * got ** 2 + c[:, None] * got + d[:, None]))
    assert resid < 1e-11, f"cubic residual {resid:.2e}"

    # quadratic_solve vs f2py_quadratic_solve — was only checked vs np.roots (a re-implementation); the f2py
    # oracle is the independent Fortran truth, incl. the complex-root branch. The JAX entry assumes a/=0, so keep
    # |a| away from 0 (a=0 is the Fortran's separate linear fallback).
    rng = np.random.default_rng(11)
    qworst = 0.0
    for _ in range(50):
        n = 8
        aa = rng.uniform(-3.0, 3.0, n); aa = np.where(np.abs(aa) < 0.3, 1.0, aa)
        bb = rng.uniform(-5.0, 5.0, n); cc = rng.uniform(-5.0, 5.0, n)   # mix of real & complex roots
        jr = np.asarray(quadratic_solve(jnp.asarray(aa), jnp.asarray(bb), jnp.asarray(cc)))   # complex (n,2)
        qre, qim = clubb_f2py.f2py_quadratic_solve(aa, bb, cc)
        qworst = max(qworst, float(max_abs_diff("calc_roots", jr.real, np.asarray(qre))),
                     float(max_abs_diff("calc_roots", jr.imag, np.asarray(qim))))
    assert qworst < 1e-12, f"f2py quadratic_solve mismatch {qworst:.2e}"
