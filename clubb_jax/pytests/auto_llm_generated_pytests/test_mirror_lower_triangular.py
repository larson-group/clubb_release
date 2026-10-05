#!/usr/bin/env python3
"""validate the JAX mirror_lower_triangular_matrix port (matrix_operations.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.matrix_operations import mirror_lower_triangular_matrix


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
    # Off-diagonal lower entries appear twice in the loss; upper entries do not contribute.
    expected = 4.0 * np.tril(np.asarray(m), -1) + 2.0 * np.diag(np.diag(np.asarray(m)))
    np.testing.assert_allclose(grad, expected, rtol=1e-14, atol=1e-15)
    print("  mirror gradient matches the analytic lower-triangle multiplicities  PASS")
