"""Test wrappers for routines from CLUBB_core/calc_roots.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_quadratic_solve_known_roots():
    """quadratic_solve should recover the roots of x^2 - 3x + 2."""
    roots = clubb_api.quadratic_solve(
        nz=1,
        a_coef=np.array([1.0], dtype=np.float64),
        b_coef=np.array([-3.0], dtype=np.float64),
        c_coef=np.array([2.0], dtype=np.float64),
    )
    assert roots.shape == (1, 2)
    np.testing.assert_allclose(np.sort(np.real(roots[0])), np.array([1.0, 2.0]))
    np.testing.assert_allclose(np.imag(roots[0]), 0.0, atol=1.0e-12)


def test_cubic_solve_known_roots():
    """cubic_solve should recover the roots of x^3 - 6x^2 + 11x - 6."""
    roots = clubb_api.cubic_solve(
        nz=1,
        a_coef=np.array([1.0], dtype=np.float64),
        b_coef=np.array([-6.0], dtype=np.float64),
        c_coef=np.array([11.0], dtype=np.float64),
        d_coef=np.array([-6.0], dtype=np.float64),
    )
    assert roots.shape == (1, 3)
    np.testing.assert_allclose(np.sort(np.real(roots[0])), np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(np.imag(roots[0]), 0.0, atol=1.0e-12)
