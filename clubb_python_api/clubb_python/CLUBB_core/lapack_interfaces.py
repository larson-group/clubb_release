"""User-facing wrappers for routines from CLUBB_core/lapack_interfaces.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def lapack_gtsv(n: int, nrhs: int, dl, d, du, b):
    """Solve a tridiagonal system through the raw LAPACK interface."""
    return clubb_f2py.f2py_lapack_gtsv(
        f_arr(dl),
        f_arr(d),
        f_arr(du),
        f_arr(b),
        n=int(n),
        nrhs=int(nrhs),
    )


def lapack_gtsvx(n: int, nrhs: int, dl, d, du, b):
    """Solve a tridiagonal system with reciprocal-condition diagnostics."""
    return clubb_f2py.f2py_lapack_gtsvx(
        f_arr(dl),
        f_arr(d),
        f_arr(du),
        f_arr(b),
        n=int(n),
        nrhs=int(nrhs),
    )


def lapack_gbsv(n: int, kl: int, ku: int, nrhs: int, ab, b):
    """Solve a banded system through the raw LAPACK interface."""
    return clubb_f2py.f2py_lapack_gbsv(
        int(kl),
        int(ku),
        f_arr(ab),
        f_arr(b),
        n=int(n),
        nrhs=int(nrhs),
    )


def lapack_gbsvx(n: int, kl: int, ku: int, nrhs: int, ab, b):
    """Solve a banded system with reciprocal-condition diagnostics."""
    soln, rcond, ferr, berr, equed, info = clubb_f2py.f2py_lapack_gbsvx(
        int(kl),
        int(ku),
        f_arr(ab),
        f_arr(b),
        n=int(n),
        nrhs=int(nrhs),
    )
    return soln, rcond, ferr, berr, equed.decode().strip(), info


def lapack_potrf(uplo: str = "Lower", a=None, n: int | None = None):
    """Compute a Cholesky factorization with the raw LAPACK interface."""
    if a is None:
        raise ValueError("a must be provided for lapack_potrf")
    if n is None:
        n = int(f_arr(a).shape[0])
    return clubb_f2py.f2py_lapack_potrf(str(uplo), f_arr(a), n=int(n))


def lapack_poequ(n: int, a):
    """Compute equilibration scaling factors for a symmetric positive-definite matrix."""
    return clubb_f2py.f2py_lapack_poequ(f_arr(a), n=int(n))


def lapack_laqsy(uplo: str = "Lower", a=None, s=None, scond: float = 0.0, amax: float = 0.0, n: int | None = None):
    """Equilibrate a symmetric matrix using precomputed scaling factors."""
    if a is None or s is None:
        raise ValueError("a and s must be provided for lapack_laqsy")
    if n is None:
        n = int(f_arr(a).shape[0])
    a_equed, equed = clubb_f2py.f2py_lapack_laqsy(
        str(uplo), f_arr(a), f_arr(s), float(scond), float(amax), n=int(n)
    )
    return a_equed, equed.decode().strip()


def lapack_syev(jobz: str = "No eigenvectors", uplo: str = "Lower", a=None, n: int | None = None):
    """Compute symmetric-matrix eigenvalues or eigenvectors."""
    if a is None:
        raise ValueError("a must be provided for lapack_syev")
    if n is None:
        n = int(f_arr(a).shape[0])
    return clubb_f2py.f2py_lapack_syev(str(jobz), str(uplo), f_arr(a), n=int(n))
