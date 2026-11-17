"""Check the native Cholesky port, including JIT, batching, retries and gradients."""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.matrix_operations import Cholesky_factor, cholesky_lower


@pytest.mark.parametrize("dimension", [1, 4, 8, 10])
def test_native_cholesky_matches_numpy_under_jit_and_batching(dimension):
    rng = np.random.default_rng(dimension)
    matrices = rng.normal(size=(2, 3, dimension, dimension))
    matrices = matrices @ np.swapaxes(matrices, -1, -2) + np.eye(dimension) * dimension
    factorize = jax.jit(jax.vmap(jax.vmap(lambda matrix: cholesky_lower(dimension, matrix)[0])))
    stored = factorize(jnp.asarray(matrices, dtype=jnp.float64))
    np.testing.assert_array_equal(np.triu(stored, 1), np.triu(matrices, 1))
    factors = np.tril(stored)
    np.testing.assert_allclose(factors, np.linalg.cholesky(matrices), atol=1e-13, rtol=1e-13)
    np.testing.assert_allclose(factors @ factors.mT, matrices, atol=1e-13, rtol=1e-13)
    assert factors.dtype == jnp.float64


def test_native_cholesky_preserves_failure_and_diagonal_retry():
    matrix = jnp.array([[1.0, 1.05], [1.05, 1.0]], dtype=jnp.float64)
    assert int(cholesky_lower(2, matrix)[1]) == 2
    scaling, factor, scaled = jax.jit(Cholesky_factor)(matrix)
    np.testing.assert_allclose(jnp.tril(factor) @ jnp.tril(factor).T, matrix + jnp.eye(2) * 0.1,
                               atol=1e-13, rtol=1e-13)
    np.testing.assert_array_equal(scaling, jnp.ones(2))
    assert not bool(scaled)


def test_native_cholesky_remains_differentiable():
    """The compiled factor derivative agrees with a centered finite difference."""
    matrix = jnp.array([[2.0, 0.3], [0.3, 1.0]], dtype=jnp.float64)
    objective = lambda values: jnp.sum(jnp.tril(Cholesky_factor(values)[1]))
    gradient = jax.jit(jax.grad(objective))(matrix)
    direction = jnp.array([[0.2, 0.0], [-0.3, 0.4]])
    step = 1e-5
    difference = (objective(matrix + step * direction) - objective(matrix - step * direction)) / (2 * step)
    np.testing.assert_allclose(jnp.sum(gradient * direction), difference, atol=1e-9, rtol=1e-9)


@pytest.mark.parametrize("diagonal", [0.0, -1.0])
def test_native_cholesky_rejects_nonpositive_singleton(diagonal):
    factor, info = jax.jit(lambda a: cholesky_lower(1, a))(jnp.array([[diagonal]], dtype=jnp.float64))
    assert int(info) == 1
    assert float(factor[0, 0]) == diagonal


@pytest.mark.parametrize("bad_pivot", [1, 2, 3])
def test_failed_pivot_leaves_later_columns_untouched(bad_pivot):
    """Return the first failed pivot and preserve the unfinished matrix region."""
    matrix = np.diag([4.0, 9.0, 16.0])
    matrix[bad_pivot-1, bad_pivot-1] = -1.0
    factor, info = jax.jit(lambda a: cholesky_lower(3, a))(jnp.asarray(matrix))
    assert int(info) == bad_pivot
    np.testing.assert_array_equal(factor[bad_pivot-1:, bad_pivot-1:], matrix[bad_pivot-1:, bad_pivot-1:])


def test_upper_triangle_is_never_read_or_scaled():
    """A NaN upper triangle cannot contaminate the lower-triangle calculation."""
    matrix = jnp.array([[1e-4, jnp.nan], [0.002, 4.0]])
    scaling, factor, scaled = jax.jit(Cholesky_factor)(matrix)
    assert bool(scaled)
    assert np.isnan(factor[0, 1])
    lower = np.tril(np.asarray(factor))
    expected = np.array([[1.0, 0.1], [0.1, 1.0]])
    np.testing.assert_allclose(lower @ lower.T, expected, atol=1e-14)
    np.testing.assert_allclose(scaling, [100.0, 0.5])


def test_retry_exhaustion_matches_ten_attempt_source_contract(capfd):
    """After ten failures Fortran returns the base plus the tenth diagonal shift."""
    matrix = jnp.array([[1.0, 4.0], [4.0, 1.0]])
    _, result, _ = jax.jit(Cholesky_factor)(matrix)
    result.block_until_ready()
    np.testing.assert_array_equal(result, np.array([[2.0, 4.0], [4.0, 2.0]]))
    assert "iteration = 10 itermax = 10" in capfd.readouterr().err



def test_mixed_batch_retries_leave_healthy_factors_unchanged():
    """A failed batch member retries without disturbing a completed member."""
    matrices = jnp.array([[[2.0, 71.0], [0.3, 1.0]], [[1.0, 19.0], [1.05, 1.0]]])
    _, factors, scaled = jax.jit(jax.vmap(Cholesky_factor))(matrices)
    assert factors.shape == (2, 2, 2)
    assert not np.any(scaled)
    for i, shift in enumerate([0.0, 0.1]):
        lower = np.tril(np.asarray(factors[i]))
        original = np.tril(matrices[i]) + np.tril(matrices[i], -1).T
        np.testing.assert_allclose(lower @ lower.T, original + shift * np.eye(2), atol=1e-14)
        assert factors[i, 0, 1] == matrices[i, 0, 1]
