"""JAX port of ``src/CLUBB_core/matrix_solver_wrapper.F90``.

Porting deviations:
  * Fortran provides generic ``band_solve`` and ``tridiag_solve`` interfaces;
    JAX keeps explicit rank-dispatch wrappers.
  * Fortran mutates ``lhs``, ``rhs``, ``err_info``, and ``soln``.  JAX treats
    array inputs functionally and returns ``err_info``, ``soln``, and ``rcond``.
  * LAPACK and ``penta_bicgstab`` solve branches are intentionally left
    unsupported in JAX.  Those branches print a warning and set ``err_info``
    fatal rather than silently using a different solver.
  * The optional Fortran ``rcond`` path calls LAPACK diagnostic solvers.  JAX
    computes an equivalent reciprocal-condition estimate for active paths:
    band solves mirror LAPACK band equilibration plus condition estimation, and
    tridiagonal solves use a dense 1-norm reciprocal condition number.

Observed Fortran differences retained in this diagnostic pass:
- The rank-2/3 NaN checks retain per-column fatal flags; Fortran marks all columns fatal when any NaN is found.
"""

from __future__ import annotations

import jax
import sys
from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.model_flags import (
    lapack,
    penta_bicgstab,
    penta_lu,
    tridiag_lu,
)
from clubb_jax.src.CLUBB_core.penta_lu_solver import (
    penta_lu_solve_multiple_rhs_lhs,
    penta_lu_solve_single_rhs_multiple_lhs,
)
from clubb_jax.src.CLUBB_core.tridiag_lu_solver import (
    tridiag_lu_solve_multiple_rhs_lhs,
    tridiag_lu_solve_single_rhs_multiple_lhs,
)


# -------------------------------------------------------------------
#                     Band Solvers Procedures
# -------------------------------------------------------------------
# JAX port note: unsupported Fortran solver choices below intentionally fail
# fatal rather than silently falling back to a different linear solver.


def _optional_rcond(l_use_rcond, estimate, disabled):
    """Evaluate a requested reciprocal-condition estimate only when enabled.

    A Python false flag removes an unrequested diagnostic during tracing. A
    runtime flag retains a scalar conditional around the expensive estimate.
    ``disabled`` supplies the existing zero output with the required shape.
    """
    if isinstance(l_use_rcond, bool):
        return estimate() if l_use_rcond else disabled
    return jax.lax.cond(l_use_rcond, lambda _: estimate(), lambda _: disabled, None)


def _reciprocal_condition_number(a):
    norm_a = jnp.max(jnp.sum(jnp.abs(a), axis=-2), axis=-1)
    inv_a = jnp.linalg.inv(a)
    norm_inv_a = jnp.max(jnp.sum(jnp.abs(inv_a), axis=-2), axis=-1)
    cond = norm_a * norm_inv_a
    l_valid = (cond > 0.0) & jnp.isfinite(cond)
    rcond = jnp.where(l_valid, 1.0 / jnp.where(l_valid, cond, 1.0), 0.0)
    return rcond


def _matrix_one_norm(a):
    one_norm = jnp.max(jnp.sum(jnp.abs(a), axis=-2), axis=-1)
    return one_norm


def _solve_matrix_vector(a, x, l_transpose=False):
    lhs = jnp.swapaxes(a, -1, -2) if l_transpose else a
    return jnp.linalg.solve(lhs, x[..., None])[..., 0]


def _sign_one(x):
    sign_one = jnp.where(x >= 0.0, 1.0, -1.0)
    return sign_one


def _reciprocal_condition_number_estimate(a):
    n = a.shape[-1]
    anorm = _matrix_one_norm(a)

    if n == 1:
        l_nonzero = jnp.abs(a[..., 0, 0]) > 0.0
        ainvnm = jnp.where(l_nonzero, 1.0 / jnp.where(l_nonzero, jnp.abs(a[..., 0, 0]), 1.0), 0.0)
        l_valid = (anorm > 0.0) & (ainvnm > 0.0)
        rcond = jnp.where(
            l_valid,
            (1.0 / jnp.where(l_valid, ainvnm, 1.0)) / jnp.where(l_valid, anorm, 1.0),
            0.0,
        )
        return rcond

    x = jnp.ones(a.shape[:-1]) / float(n)
    x = _solve_matrix_vector(a, x)
    est = jnp.sum(jnp.abs(x), axis=-1)
    isgn = _sign_one(x)
    x = _solve_matrix_vector(a, isgn, l_transpose=True)
    j = jnp.argmax(jnp.abs(x), axis=-1)
    active = jnp.ones(a.shape[:-2], dtype=bool)

    batch_axis = jnp.arange(a.shape[0])
    for iteration in range(2, 6):
        x_basis = jax.nn.one_hot(j, n)
        x_solved = _solve_matrix_vector(a, x_basis)
        est_new = jnp.sum(jnp.abs(x_solved), axis=-1)
        sign_x = _sign_one(x_solved)
        same_sign = jnp.all(sign_x == isgn, axis=-1)
        improved = est_new > est
        pre_transpose = active & (~same_sign) & improved

        est = jnp.where(active, est_new, est)
        isgn = jnp.where(pre_transpose[..., None], sign_x, isgn)

        x_transpose = _solve_matrix_vector(a, sign_x, l_transpose=True)
        jlast = j
        j_new = jnp.argmax(jnp.abs(x_transpose), axis=-1)
        x_jlast = x_transpose[batch_axis, jlast]
        x_j_new = x_transpose[batch_axis, j_new]
        active = pre_transpose & (x_jlast != jnp.abs(x_j_new)) & (iteration < 5)
        j = jnp.where(active, j_new, j)

    alt_sign = jnp.where(jnp.arange(n) % 2 == 0, 1.0, -1.0)
    alt_x = alt_sign * (1.0 + jnp.arange(n, dtype=core_rknd) / float(n - 1))
    x = _solve_matrix_vector(a, jnp.broadcast_to(alt_x, a.shape[:-1]))
    temp = 2.0 * (jnp.sum(jnp.abs(x), axis=-1) / float(3 * n))
    ainvnm = jnp.where(temp > est, temp, est)
    l_valid = (anorm > 0.0) & (ainvnm > 0.0) & jnp.isfinite(ainvnm)
    rcond = jnp.where(
        l_valid,
        (1.0 / jnp.where(l_valid, ainvnm, 1.0)) / jnp.where(l_valid, anorm, 1.0),
        0.0,
    )
    return rcond


def _lapack_equilibrated_band_matrix(a):
    threshold = 0.1
    finfo = jnp.finfo(core_rknd)
    small = finfo.tiny / finfo.eps
    large = 1.0 / small

    row_max = jnp.max(jnp.abs(a), axis=-1)
    row_min = jnp.min(row_max, axis=-1)
    row_max_all = jnp.max(row_max, axis=-1)
    row_scale = 1.0 / jnp.clip(row_max, small, large)
    rowcnd = jnp.maximum(row_min, small) / jnp.minimum(row_max_all, large)
    amax = row_max_all

    col_max = jnp.max(jnp.abs(a) * row_scale[..., :, None], axis=-2)
    col_min = jnp.min(col_max, axis=-1)
    col_max_all = jnp.max(col_max, axis=-1)
    col_scale = 1.0 / jnp.clip(col_max, small, large)
    colcnd = jnp.maximum(col_min, small) / jnp.minimum(col_max_all, large)

    row_equ = (rowcnd < threshold) | (amax < small) | (amax > large)
    col_equ = colcnd < threshold
    row_factor = jnp.where(row_equ[..., None], row_scale, 1.0)
    col_factor = jnp.where(col_equ[..., None], col_scale, 1.0)
    equilibrated_a = row_factor[..., :, None] * a * col_factor[..., None, :]
    return equilibrated_a


def _band_rcond(ndim: int, lhs):
    a = jnp.zeros((lhs.shape[1], ndim, ndim))
    diag = jnp.arange(ndim)
    a = a.at[:, diag, diag].set(lhs[2, :, :])

    if ndim > 1:
        offset = jnp.arange(ndim - 1)
        a = a.at[:, offset, offset + 1].set(lhs[1, :, : ndim - 1])
        a = a.at[:, offset + 1, offset].set(lhs[3, :, 1:ndim])

    if ndim > 2:
        offset = jnp.arange(ndim - 2)
        a = a.at[:, offset, offset + 2].set(lhs[0, :, : ndim - 2])
        a = a.at[:, offset + 2, offset].set(lhs[4, :, 2:ndim])

    return _reciprocal_condition_number_estimate(_lapack_equilibrated_band_matrix(a))


def _tridiag_rcond_single_lhs(ndim: int, lhs):
    a = jnp.zeros((ndim, ndim))
    diag = jnp.arange(ndim)
    a = a.at[diag, diag].set(lhs[1, :])

    if ndim > 1:
        offset = jnp.arange(ndim - 1)
        a = a.at[offset, offset + 1].set(lhs[0, : ndim - 1])
        a = a.at[offset + 1, offset].set(lhs[2, 1:ndim])

    return _reciprocal_condition_number(a)[None]


def _tridiag_rcond_multiple_lhs(ndim: int, lhs):
    a = jnp.zeros((lhs.shape[1], ndim, ndim))
    diag = jnp.arange(ndim)
    a = a.at[:, diag, diag].set(lhs[1, :, :])

    if ndim > 1:
        offset = jnp.arange(ndim - 1)
        a = a.at[:, offset, offset + 1].set(lhs[0, :, : ndim - 1])
        a = a.at[:, offset + 1, offset].set(lhs[2, :, 1:ndim])

    return _reciprocal_condition_number(a)


def band_solve_single_rhs_multiple_lhs(
    solve_name: str,
    penta_solve_method: int,
    ngrdcol: int,
    nsup: int,
    nsub: int,
    ndim: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
    old_soln=None,
    l_use_rcond: bool = False,
):
    """Band solve written for single RHS and multiple LHS.

    Fortran declarations (band_solve_single_rhs_multiple_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          penta_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          nsup,     & ! Number of superdiagonals
          nsub,     & ! Number of subdiagonals
          ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(out) :: &
          soln

        ! ----------------------- Optional In -----------------------
        real( kind = core_rknd ), optional, dimension(ngrdcol,ndim), intent(in) :: &
          old_soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    del old_soln

    # ----------------------- Begin Code -----------------------

    # The estimate of the reciprocal condition number of matrix
    # after equilibration (if done).
    rcond = _optional_rcond(
        l_use_rcond,
        lambda: _band_rcond(ndim, lhs),
        jnp.zeros((ngrdcol,)),
    )

    if nsup != 2 or nsub != 2:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_single_rhs_multiple_lhs "
            f"only supports pentadiagonal systems; got nsup={nsup}, nsub={nsub}.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if penta_solve_method == lapack:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_single_rhs_multiple_lhs "
            "does not yet support the LAPACK band solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if penta_solve_method == penta_lu:
        # Solve the system with a penta-diagonal specific LU decomp
        soln = penta_lu_solve_single_rhs_multiple_lhs(
            ndim, ngrdcol,
            lhs[:, :, :], rhs[:, :],
        )

    elif penta_solve_method == penta_bicgstab:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_single_rhs_multiple_lhs "
            "does not yet support the penta_bicgstab solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    else:
        # The solve method should match one of the above
        if clubb_at_least_debug_level(0):
            jax.debug.callback(
                lambda info: print(info.global_error_header_host(), "Error in band_solve_single_rhs_multiple_lhs: ",
                                   "  no case for penta_solve_method = ", penta_solve_method, sep="\n", file=sys.stderr),
                err_info, ordered=True,
            )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if l_implemented:
        err_info = check_nan_soln_2d(
            solve_name, penta_solve_method, ngrdcol, ndim, soln, err_info
        )

    return err_info, soln, rcond


def band_solve_multiple_rhs_lhs(
    solve_name: str,
    penta_solve_method: int,
    ngrdcol: int,
    nsup: int,
    nsub: int,
    ndim: int,
    nrhs: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
    old_soln=None,
    l_use_rcond: bool = False,
):
    """Band solve written for multiple RHS and multiple LHS.

    Fortran declarations (band_solve_multiple_rhs_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          penta_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          nsup,     & ! Number of superdiagonals
          nsub,     & ! Number of subdiagonals
          ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
          nrhs        ! Number of RHS's to back substitute for

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
          soln

        ! ----------------------- Optional In -----------------------
        real( kind = core_rknd ), optional, dimension(ngrdcol,ndim,nrhs), intent(in) :: &
          old_soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    del old_soln

    # ----------------------- Begin Code -----------------------

    # The estimate of the reciprocal condition number of matrix
    # after equilibration (if done).
    rcond = _optional_rcond(
        l_use_rcond,
        lambda: _band_rcond(ndim, lhs),
        jnp.zeros((ngrdcol,)),
    )

    if nsup != 2 or nsub != 2:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_multiple_rhs_lhs "
            f"only supports pentadiagonal systems; got nsup={nsup}, nsub={nsub}.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if penta_solve_method == lapack:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_multiple_rhs_lhs "
            "does not yet support the LAPACK band solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if penta_solve_method == penta_lu:
        # Solve the system with a penta-diagonal specific LU decomp
        soln = penta_lu_solve_multiple_rhs_lhs(
            ndim, nrhs, ngrdcol,
            lhs[:, :, :], rhs[:, :, :],
        )

    elif penta_solve_method == penta_bicgstab:
        print(
            "WARNING: matrix_solver_wrapper.band_solve_multiple_rhs_lhs "
            "does not yet support the penta_bicgstab solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    else:
        # The solve method should match one of the above
        if clubb_at_least_debug_level(0):
            jax.debug.callback(
                lambda info: print(info.global_error_header_host(), "Error in band_solve_multiple_rhs_lhs: ",
                                   "  no case for penta_solve_method = ", penta_solve_method, sep="\n", file=sys.stderr),
                err_info, ordered=True,
            )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if l_implemented:
        err_info = check_nan_soln_3d(
            solve_name, penta_solve_method, ngrdcol, ndim, nrhs, soln, err_info
        )

    return err_info, soln, rcond


def band_solve(
    solve_name: str,
    penta_solve_method: int,
    ngrdcol: int,
    nsup: int,
    nsub: int,
    ndim: int,
    nrhs: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
    old_soln=None,
    l_use_rcond: bool = False,
):
    """Generic band solve interface.

    Fortran declarations (band_solve_multiple_rhs_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          penta_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          nsup,     & ! Number of superdiagonals
          nsub,     & ! Number of subdiagonals
          ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
          nrhs        ! Number of RHS's to back substitute for

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
          soln

        ! ----------------------- Optional In -----------------------
        real( kind = core_rknd ), optional, dimension(ngrdcol,ndim,nrhs), intent(in) :: &
          old_soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    if rhs.ndim == 2:
        return band_solve_single_rhs_multiple_lhs(
            solve_name,
            penta_solve_method,
            ngrdcol,
            nsup,
            nsub,
            ndim,
            l_implemented,
            lhs,
            rhs,
            err_info,
            old_soln=old_soln,
            l_use_rcond=l_use_rcond,
        )
    if rhs.ndim == 3:
        return band_solve_multiple_rhs_lhs(
            solve_name,
            penta_solve_method,
            ngrdcol,
            nsup,
            nsub,
            ndim,
            nrhs,
            l_implemented,
            lhs,
            rhs,
            err_info,
            old_soln=old_soln,
            l_use_rcond=l_use_rcond,
        )
    raise ValueError("rhs must be rank-2 or rank-3 for band_solve")


# -------------------------------------------------------------------
#                    Tridiag Solver Procedures
# -------------------------------------------------------------------
# JAX port note: keep LAPACK branches fatal here until those Fortran solver paths
# have true JAX equivalents.


def tridiag_solve_single_rhs_multiple_lhs(
    solve_name: str,
    tridiag_solve_method: int,
    ngrdcol: int,
    ndim: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
    l_use_rcond: bool = False,
):
    """Tridiagonal solve written for single RHS and multiple LHS.

    Fortran declarations (tridiag_solve_single_rhs_multiple_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          tridiag_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(3,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(out) :: &
          soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    # ----------------------- Begin Code -----------------------

    # The estimate of the reciprocal condition number of matrix
    # after equilibration (if done).
    rcond = _optional_rcond(
        l_use_rcond,
        lambda: _tridiag_rcond_multiple_lhs(ndim, lhs),
        jnp.zeros((ngrdcol,)),
    )

    if tridiag_solve_method == lapack:
        print(
            "WARNING: matrix_solver_wrapper.tridiag_solve_single_rhs_multiple_lhs "
            "does not yet support the LAPACK tridiagonal solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if tridiag_solve_method == tridiag_lu:
        # Solve the system with a tri-diagonal specific LU decomp
        soln = tridiag_lu_solve_single_rhs_multiple_lhs(
            ndim, ngrdcol,
            lhs, rhs,
        )
    else:
        # The solve method should match one of the above
        if clubb_at_least_debug_level(0):
            jax.debug.callback(
                lambda info: print(info.global_error_header_host(), "Error in tridiag_solve_single_rhs_multiple_lhs: ",
                                   "  no case for tridiag_solve_method = ", tridiag_solve_method, sep="\n", file=sys.stderr),
                err_info, ordered=True,
            )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if l_implemented:
        err_info = check_nan_soln_2d(
            solve_name, tridiag_solve_method, ngrdcol, ndim, soln, err_info
        )

    return err_info, soln, rcond


def tridiag_solve_multiple_rhs_lhs(
    solve_name: str,
    tridiag_solve_method: int,
    ngrdcol: int,
    ndim: int,
    nrhs: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
    l_use_rcond: bool = False,
):
    """Tridiagonal solve written for multiple RHS and multiple LHS.

    Fortran declarations (tridiag_solve_multiple_rhs_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          tridiag_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          ndim,     & ! The order of the LHS Matrix, i.e. the # of linear equations
          nrhs        ! Number of RHS's to back substitute for

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(3,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim,nrhs), intent(out) :: &
          soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    # ----------------------- Begin Code -----------------------

    # The estimate of the reciprocal condition number of matrix
    # after equilibration (if done).
    rcond = _optional_rcond(
        l_use_rcond,
        lambda: _tridiag_rcond_multiple_lhs(ndim, lhs),
        jnp.zeros((ngrdcol,)),
    )

    if tridiag_solve_method == lapack:
        print(
            "WARNING: matrix_solver_wrapper.tridiag_solve_multiple_rhs_lhs "
            "does not yet support the LAPACK tridiagonal solver path.",
            file=sys.stderr,
        )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if tridiag_solve_method == tridiag_lu:
        # Solve the system with a tri-diagonal specific LU decomp
        soln = tridiag_lu_solve_multiple_rhs_lhs(
            ndim, nrhs, ngrdcol,
            lhs, rhs,
        )
    else:
        # The solve method should match one of the above
        if clubb_at_least_debug_level(0):
            jax.debug.callback(
                lambda info: print(info.global_error_header_host(), "Error in tridiag_solve_multiple_rhs_lhs: ",
                                   "  no case for tridiag_solve_method = ", tridiag_solve_method, sep="\n", file=sys.stderr),
                err_info, ordered=True,
            )
        return err_info.set_fatal(), jnp.zeros_like(rhs), rcond

    if l_implemented:
        err_info = check_nan_soln_3d(
            solve_name, tridiag_solve_method, ngrdcol, ndim, nrhs, soln, err_info
        )

    return err_info, soln, rcond


def tridiag_solve(
    solve_name: str,
    tridiag_solve_method: int,
    ngrdcol: int,
    ndim: int,
    lhs,
    rhs,
    err_info,
    l_use_rcond: bool = False,
    l_implemented: bool = False,
):
    """Generic tridiagonal solve interface.

    Fortran declarations (tridiag_solve_single_rhs_multiple_lhs in src/CLUBB_core/matrix_solver_wrapper.F90):

        ! ----------------------- Input Variables -----------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          tridiag_solve_method

        integer, intent(in) :: &
          ngrdcol,  & ! Number of grid columns
          ndim        ! The order of the LHS Matrix, i.e. the # of linear equations

        logical, intent(in) :: &
          l_implemented ! True if CLUBB is being implemented and run in a host model

        real( kind = core_rknd ), dimension(3,ngrdcol,ndim), intent(inout) :: &
          lhs ! Left hand side

        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(inout) :: &
          rhs ! Right hand side(s)

        type(err_info_type), intent(inout) :: &
          err_info        ! err_info struct containing err_code and err_header

        ! ----------------------- Output Variables -----------------------
        real( kind = core_rknd ), dimension(ngrdcol,ndim), intent(out) :: &
          soln

        ! ----------------------- Optional Out -----------------------

        ! The estimate of the reciprocal condition number of matrix
        ! after equilibration (if done).
        real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
          rcond
    """
    if lhs.ndim == 3 and rhs.ndim == 2:
        return tridiag_solve_single_rhs_multiple_lhs(
            solve_name,
            tridiag_solve_method,
            ngrdcol,
            ndim,
            l_implemented,
            lhs,
            rhs,
            err_info,
            l_use_rcond=l_use_rcond,
        )
    if lhs.ndim == 3 and rhs.ndim == 3:
        return tridiag_solve_multiple_rhs_lhs(
            solve_name,
            tridiag_solve_method,
            ngrdcol,
            ndim,
            rhs.shape[2],
            l_implemented,
            lhs,
            rhs,
            err_info,
            l_use_rcond=l_use_rcond,
        )
    raise ValueError("Unsupported lhs/rhs ranks for tridiag_solve")


def check_nan_soln_2d(
    solve_name: str, solve_method: int, ngrdcol: int, ndim: int, soln, err_info
):
    """
    Fortran declarations (check_nan_soln_2d in src/CLUBB_core/matrix_solver_wrapper.F90):

        !--------------------------- Input Variables ---------------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          solve_method, ngrdcol, ndim

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim) :: &
          soln

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info
    """
    del ndim
    l_has_nan = jnp.any(jnp.isnan(soln), axis=1)
    # Host I/O adaptation of source rank-specific loops; fatal state unchanged.
    def print_nan_solution(values, info):
        import numpy as np
        for k, i in np.argwhere(np.isnan(values).T):
            print("ERROR: NaN detected: soln(", i+1, ",", k+1, ") = ", values[i,k], file=sys.stderr)
        print(info.global_error_header_host(), file=sys.stderr)
        print("solve_method =", solve_method, file=sys.stderr)
        print(solve_name.strip() + " NaNs in solution", file=sys.stderr)
    jax.lax.cond(
        jnp.any(l_has_nan),
        lambda _: jax.debug.callback(print_nan_solution, soln, err_info, ordered=True),
        lambda _: None, operand=None,
    )
    return err_info.set_fatal(mask=l_has_nan[:ngrdcol])


def check_nan_soln_3d(
    solve_name: str,
    solve_method: int,
    ngrdcol: int,
    ndim: int,
    nrhs: int,
    soln,
    err_info,
):
    """
    Fortran declarations (check_nan_soln_3d in src/CLUBB_core/matrix_solver_wrapper.F90):

        !--------------------------- Input Variables ---------------------------
        character(len=*), intent(in) :: &
          solve_name

        integer, intent(in) :: &
          solve_method, ngrdcol, ndim, nrhs

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim,nrhs) :: &
          soln

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info
    """
    del ndim, nrhs
    l_has_nan = jnp.any(jnp.isnan(soln), axis=(1, 2))
    # Host I/O adaptation of source rank-specific loops; fatal state unchanged.
    def print_nan_solution(values, info):
        import numpy as np
        for j, k, i in np.argwhere(np.isnan(values).transpose(2, 1, 0)):
            print("ERROR: NaN detected: soln(", i+1, ",", k+1, ",", j+1, ") = ", values[i,k,j], file=sys.stderr)
        print(info.global_error_header_host(), file=sys.stderr)
        print("solve_method =", solve_method, file=sys.stderr)
        print(solve_name.strip() + " NaNs in solution", file=sys.stderr)
    jax.lax.cond(
        jnp.any(l_has_nan),
        lambda _: jax.debug.callback(print_nan_solution, soln, err_info, ordered=True),
        lambda _: None, operand=None,
    )
    return err_info.set_fatal(mask=l_has_nan[:ngrdcol])
