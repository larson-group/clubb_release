"""Test wrappers for routines from CLUBB_core/lapack_interfaces.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_lapack_gtsv_and_gtsvx_solve_diagonal_system():
    """Raw tridiagonal LAPACK wrappers should solve a diagonal system."""
    n = 4
    nrhs = 2
    diag = np.array([2.0, 4.0, 8.0, 16.0], dtype=np.float64)
    expected = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]], dtype=np.float64, order="F")
    rhs = (diag[:, None] * expected).copy(order="F")

    soln_simple, info_simple = clubb_api.lapack_gtsv(n=n, nrhs=nrhs, dl=np.zeros(n - 1), d=diag, du=np.zeros(n - 1), b=rhs)
    soln_expert, rcond, ferr, berr, info_expert = clubb_api.lapack_gtsvx(
        n=n, nrhs=nrhs, dl=np.zeros(n - 1), d=diag, du=np.zeros(n - 1), b=rhs
    )

    np.testing.assert_allclose(soln_simple, expected)
    np.testing.assert_allclose(soln_expert, expected)
    assert int(info_simple) == 0
    assert int(info_expert) == 0
    assert float(rcond) > 0.0
    assert np.all(ferr >= 0.0)
    assert np.all(berr >= 0.0)


def test_lapack_gbsv_and_gbsvx_solve_diagonal_band_system():
    """Raw banded LAPACK wrappers should solve a diagonal system."""
    n = 4
    nrhs = 2
    diag = np.array([3.0, 6.0, 9.0, 12.0], dtype=np.float64)
    expected = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]], dtype=np.float64, order="F")
    rhs = (diag[:, None] * expected).copy(order="F")
    ab_simple = diag.reshape(1, -1).copy(order="F")

    soln_simple, ipiv, info_simple = clubb_api.lapack_gbsv(n=n, kl=0, ku=0, nrhs=nrhs, ab=ab_simple, b=rhs)
    soln_expert, rcond, ferr, berr, equed, info_expert = clubb_api.lapack_gbsvx(
        n=n, kl=0, ku=0, nrhs=nrhs, ab=ab_simple, b=rhs
    )

    np.testing.assert_allclose(soln_simple, expected)
    np.testing.assert_allclose(soln_expert, expected)
    assert ipiv.shape == (n,)
    assert int(info_simple) == 0
    assert int(info_expert) == 0
    assert float(rcond) > 0.0
    assert np.all(ferr >= 0.0)
    assert np.all(berr >= 0.0)
    assert str(equed).strip() in {"N", "Y"}


def test_lapack_factorization_and_eigen_helpers_work_on_spd_matrix():
    """potrf, poequ, laqsy, and syev wrappers should behave consistently."""
    a = np.array(
        [[4.0, 1.0, 1.0],
         [1.0, 3.0, 0.5],
         [1.0, 0.5, 2.5]],
        dtype=np.float64,
        order="F",
    )

    a_factor, info_potrf = clubb_api.lapack_potrf(n=3, a=a)
    s, scond, amax, info_poequ = clubb_api.lapack_poequ(n=3, a=a)
    a_equed, equed = clubb_api.lapack_laqsy(n=3, a=a, s=s, scond=scond, amax=amax)
    a_out, w, info_syev = clubb_api.lapack_syev(n=3, a=a)

    np.testing.assert_allclose(np.tril(a_factor) @ np.tril(a_factor).T, a, rtol=1e-10, atol=1e-10)
    assert int(info_potrf) == 0
    assert int(info_poequ) == 0
    assert np.all(s > 0.0)
    assert float(scond) > 0.0
    assert float(amax) > 0.0
    assert a_equed.shape == a.shape
    assert str(equed).strip() in {"N", "Y"}
    assert int(info_syev) == 0
    np.testing.assert_allclose(np.sort(w), np.linalg.eigvalsh(a), rtol=1e-10, atol=1e-10)
    assert a_out.shape == a.shape
