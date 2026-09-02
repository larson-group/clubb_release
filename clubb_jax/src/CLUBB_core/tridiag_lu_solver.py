"""JAX implementation of ``tridiag_lu_solver.F90``.

Description:
  These routines solve lhs*soln=rhs using LU decomp.

  LHS is stored in band diagonal form.
    lhs = | lhs(0,1)  lhs(-1,1)      0          0          0          0      0
          | lhs(1,2)  lhs( 0,2)  lhs(-1,2)      0          0          0      0
          |     0     lhs( 1,3)  lhs( 0,3)  lhs(-1,3)      0          0      0
          |     0         0      lhs( 1,4)  lhs( 0,4)  lhs(-1,4)      0      0
          |     0         0          0      lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  0 ...
          |    ...

   U is stored in band diagonal form
    U = |   1    upper(1)  0       0       0       0       0
        |   0      1     upper(2)  0       0       0       0
        |   0      0       1     upper(3)  0       0       0
        |   0      0       0       1     upper(4)  0       0
        |   0      0       0       0       1     upper(5)  0  ...
        |  ...

  L is also stored in band diagonal form, but the lowest most band is equivalent to the
  lowermost band of LHS, thus we don't need to store it
    L = | l_diag(1)       0           0            0          0        0
        |  lower(2)      l_diag(2)    0            0          0        0
        |    0           lower(3)   l_diag(3)      0          0        0
        |    0            0         lower(4)    l_diag(4)     0        0
        |    0            0           0         lower(5)    l_diag(5)   0   ...
        |  ...

  To perform the LU decomposition, we go element by element.
  First we start by noting that we want lhs=LU, so the first step of calculating
  L*U, by multiplying the first row of L by the columns of U, gives us

    l_diag(1)*1        = lhs( 0,1)  =>  l_diag(1)   = lhs( 0,1)
    l_diag(1)*upper(1) = lhs(-1,1)  =>  upper(1)    = lhs(-1,1) / l_diag(1)

  Now that we're passed the k=1 step, each following step uses all the bands,
  allowing us to write the general step

    lower(k)*1 = lhs(1,k)                         =>  lower(k)    = lhs( 1,k)
    lower(k)*upper(k-1)+l_diag(k)*1  = lhs( 0,k)  =>  l_diag(k) = lhs( 0,k)
                                                                  - lower(k)*upper(k-1)
    l_diag(k)*upper(k)               = lhs(-1,k)  =>  upper(k)  = lhs(-1,k) / l_diag(k)

  This general step is done for k from 2 to ndim-1 (do k = 2, ndim-), and the last
  step is tweaked similarly to the first one, where we disclude the upper band
  since it becomes no longer relevant. Note from this general step that the lower band
  is always equivalent to first subdiagonal band of lhs, thus we do not need to
  calculate or store lower. Also note that we only ever need l_diag so that we can divide
  by it, so instead we compute lower_diag_invrs to reduce divide operations.

  After L and U are computed, normally we do forward substitution using L,
  then backward substitution using U to find the solution. This is replicated
  for every right hand side we want to solve for.

References:
  none

Porting deviations:
The Fortran generic interface is implemented as a Python rank-dispatch wrapper.
The explicit Fortran loops are represented by ``jax.lax.scan`` sweeps.  Fortran
mutates ``soln`` output arrays; the JAX routines return the solution arrays.
The band dimension uses Python indices ``0, 1, 2`` for Fortran bands
``-1, 0, 1``.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp


@partial(jax.jit, static_argnames=("ndim",))
def tridiag_lu_solve_single_rhs_lhs(ndim: int, lhs, rhs):
    """Solve one tridiagonal system with one right-hand side.

    Description:
      Written for single RHS and single LHS.
    """
    del ndim

    # Matrices to solve, stored using band diagonal vectors
    # -2 is the uppermost band, 2 is the lower most band, 0 is diagonal
    upper_band = lhs[0]
    diag_band = lhs[1]
    lower_band = lhs[2]

    # ----------------------- Begin Code -----------------------

    # Inverse of the diagonal of L
    lower_diag_invrs_0 = 1.0 / diag_band[0]

    # First U band
    upper_0 = lower_diag_invrs_0 * upper_band[0]

    def lu_step(upper_prev, x):
        diag_k, lower_k, upper_k_lhs = x
        lower_diag_invrs_k = 1.0 / (diag_k - lower_k * upper_prev)
        upper_k = lower_diag_invrs_k * upper_k_lhs
        return upper_k, (lower_diag_invrs_k, upper_k)

    upper_last, (lower_diag_invrs_mid, upper_mid) = jax.lax.scan(
        lu_step,
        upper_0,
        (diag_band[1:-1], lower_band[1:-1], upper_band[1:-1]),
    )

    lower_diag_invrs_last = 1.0 / (
        diag_band[-1] - lower_band[-1] * upper_last
    )
    lower_diag_invrs = jnp.concatenate(
        [
            lower_diag_invrs_0[None],
            lower_diag_invrs_mid,
            lower_diag_invrs_last[None],
        ],
        axis=0,
    )
    upper = jnp.concatenate([upper_0[None], upper_mid], axis=0)

    soln_0 = lower_diag_invrs[0] * rhs[0]

    def forward_substitution(soln_prev, x):
        lower_diag_invrs_k, lower_k, rhs_k = x
        soln_k = lower_diag_invrs_k * (rhs_k - lower_k * soln_prev)
        return soln_k, soln_k

    _, soln_rest = jax.lax.scan(
        forward_substitution,
        soln_0,
        (lower_diag_invrs[1:], lower_band[1:], rhs[1:]),
    )
    soln = jnp.concatenate([soln_0[None], soln_rest], axis=0)

    def backward_substitution(soln_next, x):
        soln_k_old, upper_k = x
        soln_k = soln_k_old - upper_k * soln_next
        return soln_k, soln_k

    _, soln_reversed = jax.lax.scan(
        backward_substitution,
        soln[-1],
        (soln[:-1][::-1], upper[::-1]),
    )
    return jnp.concatenate([soln_reversed[::-1], soln[-1:]], axis=0)


@partial(jax.jit, static_argnames=("ndim", "ngrdcol"))
def tridiag_lu_solve_single_rhs_multiple_lhs(ndim: int, ngrdcol: int, lhs, rhs):
    """Solve multiple tridiagonal systems with one right-hand side each.

    Description:
      Written for single RHS and multiple LHS.
    """
    del ndim, ngrdcol

    # Matrices to solve, stored using band diagonal vectors
    # -2 is the uppermost band, 2 is the lower most band, 0 is diagonal
    upper_band = lhs[0].T
    diag_band = lhs[1].T
    lower_band = lhs[2].T
    rhs_t = rhs.T

    # ----------------------- Begin Code -----------------------

    # Inverse of the diagonal of L
    lower_diag_invrs_0 = 1.0 / diag_band[0]

    # First U band
    upper_0 = lower_diag_invrs_0 * upper_band[0]

    def lu_step(upper_prev, x):
        diag_k, lower_k, upper_k_lhs = x
        lower_diag_invrs_k = 1.0 / (diag_k - lower_k * upper_prev)
        upper_k = lower_diag_invrs_k * upper_k_lhs
        return upper_k, (lower_diag_invrs_k, upper_k)

    upper_last, (lower_diag_invrs_mid, upper_mid) = jax.lax.scan(
        lu_step,
        upper_0,
        (diag_band[1:-1], lower_band[1:-1], upper_band[1:-1]),
    )

    lower_diag_invrs_last = 1.0 / (
        diag_band[-1] - lower_band[-1] * upper_last
    )
    lower_diag_invrs = jnp.concatenate(
        [
            lower_diag_invrs_0[None],
            lower_diag_invrs_mid,
            lower_diag_invrs_last[None],
        ],
        axis=0,
    )
    upper = jnp.concatenate([upper_0[None], upper_mid], axis=0)

    soln_0 = lower_diag_invrs[0] * rhs_t[0]

    def forward_substitution(soln_prev, x):
        lower_diag_invrs_k, lower_k, rhs_k = x
        soln_k = lower_diag_invrs_k * (rhs_k - lower_k * soln_prev)
        return soln_k, soln_k

    _, soln_rest = jax.lax.scan(
        forward_substitution,
        soln_0,
        (lower_diag_invrs[1:], lower_band[1:], rhs_t[1:]),
    )
    soln = jnp.concatenate([soln_0[None], soln_rest], axis=0)

    def backward_substitution(soln_next, x):
        soln_k_old, upper_k = x
        soln_k = soln_k_old - upper_k * soln_next
        return soln_k, soln_k

    _, soln_reversed = jax.lax.scan(
        backward_substitution,
        soln[-1],
        (soln[:-1][::-1], upper[::-1]),
    )
    return jnp.concatenate([soln_reversed[::-1], soln[-1:]], axis=0).T


@partial(jax.jit, static_argnames=("ndim", "nrhs", "ngrdcol"))
def tridiag_lu_solve_multiple_rhs_lhs(
    ndim: int, nrhs: int, ngrdcol: int, lhs, rhs
):
    """Solve multiple tridiagonal systems with multiple right-hand sides.

    Description:
      Written for multiple RHS and multiple LHS.
    """
    del ndim, nrhs, ngrdcol

    # Matrices to solve, stored using band diagonal vectors
    # -2 is the uppermost band, 2 is the lower most band, 0 is diagonal
    upper_band = lhs[0].T
    diag_band = lhs[1].T
    lower_band = lhs[2].T
    rhs_t = jnp.transpose(rhs, (1, 0, 2))

    # ----------------------- Begin Code -----------------------

    # Inverse of the diagonal of L
    lower_diag_invrs_0 = 1.0 / diag_band[0]

    # First U band
    upper_0 = lower_diag_invrs_0 * upper_band[0]

    def lu_step(upper_prev, x):
        diag_k, lower_k, upper_k_lhs = x
        lower_diag_invrs_k = 1.0 / (diag_k - lower_k * upper_prev)
        upper_k = lower_diag_invrs_k * upper_k_lhs
        return upper_k, (lower_diag_invrs_k, upper_k)

    upper_last, (lower_diag_invrs_mid, upper_mid) = jax.lax.scan(
        lu_step,
        upper_0,
        (diag_band[1:-1], lower_band[1:-1], upper_band[1:-1]),
    )

    lower_diag_invrs_last = 1.0 / (
        diag_band[-1] - lower_band[-1] * upper_last
    )
    lower_diag_invrs = jnp.concatenate(
        [
            lower_diag_invrs_0[None],
            lower_diag_invrs_mid,
            lower_diag_invrs_last[None],
        ],
        axis=0,
    )
    upper = jnp.concatenate([upper_0[None], upper_mid], axis=0)

    soln_0 = lower_diag_invrs[0, :, None] * rhs_t[0]

    def forward_substitution(soln_prev, x):
        lower_diag_invrs_k, lower_k, rhs_k = x
        soln_k = lower_diag_invrs_k[:, None] * (
            rhs_k - lower_k[:, None] * soln_prev
        )
        return soln_k, soln_k

    _, soln_rest = jax.lax.scan(
        forward_substitution,
        soln_0,
        (lower_diag_invrs[1:], lower_band[1:], rhs_t[1:]),
    )
    soln = jnp.concatenate([soln_0[None], soln_rest], axis=0)

    def backward_substitution(soln_next, x):
        soln_k_old, upper_k = x
        soln_k = soln_k_old - upper_k[:, None] * soln_next
        return soln_k, soln_k

    _, soln_reversed = jax.lax.scan(
        backward_substitution,
        soln[-1],
        (soln[:-1][::-1], upper[::-1]),
    )
    return jnp.transpose(
        jnp.concatenate([soln_reversed[::-1], soln[-1:]], axis=0),
        (1, 0, 2),
    )


def tridiag_lu_solve(ndim: int, lhs, rhs):
    """Generic tridiagonal LU solve interface."""
    if lhs.ndim == 2 and rhs.ndim == 1:
        return tridiag_lu_solve_single_rhs_lhs(ndim, lhs, rhs)
    if lhs.ndim == 3 and rhs.ndim == 2:
        return tridiag_lu_solve_single_rhs_multiple_lhs(
            ndim, lhs.shape[1], lhs, rhs
        )
    if lhs.ndim == 3 and rhs.ndim == 3:
        return tridiag_lu_solve_multiple_rhs_lhs(
            ndim, rhs.shape[2], lhs.shape[1], lhs, rhs
        )
    raise ValueError("Unsupported lhs/rhs ranks for tridiag_lu_solve")


__all__ = [
    "tridiag_lu_solve",
    "tridiag_lu_solve_single_rhs_lhs",
    "tridiag_lu_solve_single_rhs_multiple_lhs",
    "tridiag_lu_solve_multiple_rhs_lhs",
]
