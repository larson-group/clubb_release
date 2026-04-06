"""Test wrappers for routines from CLUBB_core/matrix_operations.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_mirror_lower_triangular_matrix_reflects_lower_triangle():
    """mirror_lower_triangular_matrix should copy lower-triangle entries upward."""
    matrix = np.array(
        [[1.0, 0.0, 0.0],
         [2.0, 3.0, 0.0],
         [4.0, 5.0, 6.0]],
        dtype=np.float64,
    )

    mirrored = clubb_api.mirror_lower_triangular_matrix(3, matrix.copy())
    expected = np.array(
        [[1.0, 2.0, 4.0],
         [2.0, 3.0, 5.0],
         [4.0, 5.0, 6.0]],
        dtype=np.float64,
    )
    np.testing.assert_allclose(mirrored, expected)


def test_cholesky_factor_reconstructs_spd_matrix():
    """cholesky_factor should return a valid factorization for an SPD matrix."""
    a_input = np.array(
        [[4.0, 1.0, 1.0],
         [1.0, 3.0, 0.5],
         [1.0, 0.5, 2.5]],
        dtype=np.float64,
    )

    a_scaling, a_cholesky, l_scaled = clubb_api.cholesky_factor(3, a_input)

    assert a_scaling.shape == (3,)
    assert a_cholesky.shape == (3, 3)
    assert isinstance(bool(l_scaled), bool)
    recon = np.tril(a_cholesky) @ np.tril(a_cholesky).T
    if l_scaled:
        recon = recon / np.outer(a_scaling, a_scaling)
    np.testing.assert_allclose(recon, a_input, rtol=1.0e-10, atol=1.0e-10)
