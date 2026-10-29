import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.matrix_operations import Cholesky_factor, _portable_cholesky


@pytest.mark.parametrize("dimension", [1, 4, 8, 10])
def test_portable_cholesky_matches_numpy_under_jit_and_batching(dimension):
    rng = np.random.default_rng(dimension)
    matrices = rng.normal(size=(2, 3, dimension, dimension))
    matrices = matrices @ np.swapaxes(matrices, -1, -2) + np.eye(dimension) * dimension
    factorize = jax.jit(jax.vmap(jax.vmap(_portable_cholesky)))
    factors = factorize(jnp.asarray(matrices, dtype=jnp.float64))
    np.testing.assert_allclose(factors, np.linalg.cholesky(matrices), atol=1e-13, rtol=1e-13)
    np.testing.assert_allclose(factors @ factors.mT, matrices, atol=1e-13, rtol=1e-13)
    assert factors.dtype == jnp.float64


def test_portable_cholesky_preserves_failure_and_diagonal_retry(monkeypatch):
    monkeypatch.setenv("CLUBB_JAX_PORTABLE_CHOLESKY", "1")
    matrix = jnp.array([[1.0, 1.05], [1.05, 1.0]], dtype=jnp.float64)
    assert np.isnan(np.asarray(_portable_cholesky(matrix))).all()
    scaling, factor, scaled = Cholesky_factor(matrix)
    np.testing.assert_allclose(jnp.tril(factor) @ jnp.tril(factor).T, matrix + jnp.eye(2) * 0.1,
                               atol=1e-13, rtol=1e-13)
    np.testing.assert_array_equal(scaling, jnp.ones(2))
    assert not bool(scaled)


def test_portable_cholesky_remains_differentiable():
    matrix = jnp.array([[2.0, 0.3], [0.3, 1.0]], dtype=jnp.float64)
    gradient = jax.jit(jax.grad(lambda values: jnp.sum(_portable_cholesky(values))))(matrix)
    expected = jax.grad(lambda values: jnp.sum(jnp.linalg.cholesky(values)))(matrix)
    np.testing.assert_allclose(gradient, expected, atol=1e-13, rtol=1e-13)


@pytest.mark.parametrize("diagonal", [0.0, -1.0])
def test_portable_cholesky_rejects_nonpositive_singleton(diagonal):
    factor = jax.jit(_portable_cholesky)(jnp.array([[diagonal]], dtype=jnp.float64))
    assert np.isnan(np.asarray(factor)).all()
