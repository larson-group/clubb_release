#!/usr/bin/env python3
"""test_mirror_lower_triangular.py — validate the JAX mirror_lower_triangular_matrix port (matrix_operations.F90).

Mirrors a lower-triangular matrix's lower triangle onto the upper triangle to make it symmetric. Oracles:
  1. f2py bit-shadow vs f2py_mirror_lower_triangular_matrix. SKIPs if clubb_f2py is unbuilt.
  2. Symmetry + lower-triangle/diagonal preservation invariants, and a finite jax.grad.
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

from clubb_jax.src.CLUBB_core.matrix_operations import mirror_lower_triangular_matrix


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py mirror_lower_triangular_matrix oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(9)
    worst = 0.0
    for nvars in (1, 2, 4, 7):
        m = rng.uniform(-3, 3, (nvars, nvars))
        # f2py mutates in-place AND returns; pass a fresh Fortran-ordered copy so the input array is untouched.
        f = np.asarray(clubb_f2py.f2py_mirror_lower_triangular_matrix(np.asfortranarray(m.copy())))
        g = np.asarray(mirror_lower_triangular_matrix(m))
        worst = max(worst, np.max(np.abs(g - f)))
    assert worst < 1e-15, f"mirror_lower_triangular_matrix f2py mismatch {worst:.2e}"
    print(f"  f2py mirror_lower_triangular_matrix: bit-match (nvars 1,2,4,7), worst {worst:.2e}  PASS")


def test_invariants():
    rng = np.random.default_rng(1)
    m = rng.uniform(-3, 3, (5, 5))
    g = np.asarray(mirror_lower_triangular_matrix(m))
    assert np.allclose(g, g.T), "result not symmetric"
    il = np.tril_indices(5)
    assert np.allclose(g[il], m[il]), "lower triangle / diagonal not preserved"
    # Upper triangle must equal the original lower triangle's transpose, independent of m's upper triangle.
    m2 = m.copy(); m2[np.triu_indices(5, 1)] = 999.0
    g2 = np.asarray(mirror_lower_triangular_matrix(m2))
    assert np.allclose(g, g2), "result depends on the input's upper triangle (it must not)"
    print("  symmetry + lower-triangle preservation + upper-triangle independence  PASS")


def test_differentiable():
    m = jnp.asarray(np.random.default_rng(2).uniform(-1, 1, (4, 4)))
    grad = np.asarray(jax.grad(lambda x: jnp.sum(mirror_lower_triangular_matrix(x) ** 2))(m))
    assert np.isfinite(grad).all(), "non-finite grad"
    print(f"  jax.grad through mirror_lower_triangular_matrix: finite ({grad.size} entries)  PASS")


def main():
    print("test_mirror_lower_triangular:")
    for t in (test_f2py_oracle, test_invariants, test_differentiable):
        t()
    print("All mirror_lower_triangular_matrix checks PASSED")


if __name__ == "__main__":
    main()
