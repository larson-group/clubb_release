"""Test wrappers for routines from CLUBB_core/diagnose_correlations_module.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_diagnose_correlations_preserves_identity_matrix():
    """An identity prescribed matrix should remain identity when diagnosed."""
    corr_array_pre = np.eye(3, dtype=np.float64)
    diagnosed = clubb_api.diagnose_correlations(3, 1, corr_array_pre, False)
    np.testing.assert_allclose(diagnosed, np.eye(3, dtype=np.float64))


def test_calc_cholesky_corr_mtx_approx_handles_identity_matrix():
    """Identity correlation matrix should yield identity Cholesky and approximation."""
    corr_matrix = np.eye(3, dtype=np.float64)
    corr_cholesky_mtx, corr_mtx_approx = clubb_api.calc_cholesky_corr_mtx_approx(3, 1, corr_matrix)
    np.testing.assert_allclose(corr_cholesky_mtx, np.eye(3, dtype=np.float64))
    np.testing.assert_allclose(corr_mtx_approx, np.eye(3, dtype=np.float64))
