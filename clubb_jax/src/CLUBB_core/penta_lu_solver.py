"""JAX implementation of ``penta_lu_solver.F90``.

Description:
  These routines solve lhs*soln=rhs using LU decomp.

  LHS is stored in band diagonal form.
    lhs = | lhs(0,1)  lhs(-1,1)  lhs(-2,1)     0           0          0          0
          | lhs(1,2)  lhs( 0,2)  lhs(-1,2)  lhs(-2,2)      0          0          0
          | lhs(2,3)  lhs( 1,3)  lhs( 0,3)  lhs(-1,3)  lhs(-2,3)      0          0
          |     0     lhs( 2,4)  lhs( 1,4)  lhs( 0,4)  lhs(-1,4)  lhs(-2,4)      0
          |     0         0      lhs( 2,5)  lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  lhs(-2,5) ...
          |    ...

   U is stored in band diagonal form
    U = |   1    upper_1(1)  upper_2(1)      0           0           0           0
        |   0        1       upper_1(2)  upper_2(2)      0           0           0
        |   0        0           1       upper_1(3)  upper_2(3)      0           0
        |   0        0           0           1       upper_1(4)  upper_2(4)      0
        |   0        0           0           0           1       upper_1(5)  upper_2(5)  ...
        |  ...

  L is also stored in band diagonal form, but the lowest most band is equivalent to the
  lowermost band of LHS, thus we don't need to store it
    L = | l_diag(1)        0            0           0             0       0
        | lower_1(2)    l_diag(2)       0           0             0       0
        |   l_2(3)      lower_1(3)   l_diag(3)      0             0       0
        |     0           l_2(4)     lower_1(4)  l_diag(4)        0       0
        |     0             0        l_2(5)      lower_1(5)    l_diag(5)  0   ...
        |  ...

  To perform the LU decomposition, we go element by element.
  First we start by noting that we want lhs=LU, so the first step of calculating
  L*U, by multiplying the first row of L by the columns of U, gives us

    l_diag(1)*1          = lhs( 0,1)  =>  l_diag(1)  = lhs( 0,1)
    l_diag(1)*upper_1(1) = lhs(-1,1)  =>  upper_1(1) = lhs(-1,1) / l_diag(1)
    l_diag(1)*upper_2(1) = lhs(-2,1)  =>  upper_2(1) = lhs(-2,1) / l_diag(1)

  Multiplying the second row of L by U now we get

    lower_1(2)*1                               = lhs(1,2)   =>  lower_1(2)  = lhs(1,2)
    lower_1(2)*upper_1(1)+l_diag(2)*1          = lhs(0,2)   =>  l_diag(2)
                                                                = lhs(0,2)
                                                                  - lower_1(2)*upper_1(1)
    lower_1(2)*upper_2(1)+l_diag(2)*upper_1(2) = lhs(-1,2)  =>  upper_1(2)
                                                                = ( lhs(-1,2)
                                                                  - lower_1(2)*upper_2(1) )
                                                                  / l_diag(2)
    l_diag(2)*upper_2(2)                       = lhs(-2,2)  =>  upper_2(2)  =
                                                                            lhs(-2,2)
                                                                            / l_diag(2)

  Now that we're passed the k=1 and k=2 steps, each following step uses all the bands,
  allowing us to write the general step

    l_2(k)*1                                    = lhs(2,k)   =>  l_2(k)     = lhs(2,k)
    l_2(k)*upper_1(k-2)+lower_1(k)*1            = lhs(1,k)   =>  lower_1(k)
                                                                 = lhs(1,k)
                                                                   - l_2(k)*upper_1(k-2)
    l_2(k)*upper_2(k-2)+lower_1(k)*upper_1(k-1) = lhs( 0,k)  =>  l_diag(k)  = lhs(0,k)
                   +l_diag(k)*1                                               - l_2(k)
                                                                              *upper_2(k-2)
                                                                              + lower_1(k)
                                                                              *upper_1(k-1)

    lower_1(k)*upper_2(k-1)
    + l_diag(k)*upper_1(k) = lhs(-1,k)  =>  upper_1(k) = ( lhs(-1,k) - lower_1(k)*upper_2(k-1) )
                                                         / l_diag(k)
    l_diag(k)*upper_2(k)   = lhs(-2,k)  =>  upper_2(k) = lhs(-2,k) / l_diag(k)

  This general step is done for k from 3 to ndim-2 (do k = 3, ndim-2), and the last two
  steps are tweaked similarly to the first two, where we disclude one then two bands
  since they become no longer relevant. Note from this general step that the l_2 band
  is always equivalent to second subdiagonal band of lhs, thus we do not need to
  calculate or store l_2. Also note that we only ever need l_diag so that we can divide
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
The band dimension uses Python indices ``0, 1, 2, 3, 4`` for Fortran bands
``-2, -1, 0, 1, 2``.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp


@partial(jax.jit, static_argnames=("ndim", "ngrdcol"))
def penta_lu_solve_single_rhs_multiple_lhs(ndim: int, ngrdcol: int, lhs, rhs):
    """Solve multiple pentadiagonal systems with one right-hand side each.

    Description:
      Written for single RHS and multiple LHS.
    """
    del ngrdcol

    # Matrices to solve, stored using band diagonal vectors
    # -2 is the uppermost band, 2 is the lower most band, 0 is diagonal
    upper_2_band = lhs[0].T
    upper_1_band = lhs[1].T
    diag_band = lhs[2].T
    lower_1_band = lhs[3].T
    lower_2_band = lhs[4].T
    rhs_t = rhs.T

    # ----------------------- Begin Code -----------------------

    # Inverse of the diagonal of L
    lower_diag_invrs_0 = 1.0 / diag_band[0]

    # First U band
    upper_1_0 = lower_diag_invrs_0 * upper_1_band[0]

    # Second U band
    upper_2_0 = lower_diag_invrs_0 * upper_2_band[0]

    # First L band
    lower_1_0 = jnp.zeros_like(lower_diag_invrs_0)

    # Second L band
    lower_2_0 = jnp.zeros_like(lower_diag_invrs_0)

    lower_1_1 = lower_1_band[1]
    lower_2_1 = jnp.zeros_like(lower_diag_invrs_0)
    lower_diag_invrs_1 = 1.0 / (
        diag_band[1] - lower_1_1 * upper_1_0
    )
    upper_1_1 = lower_diag_invrs_1 * (
        upper_1_band[1] - lower_1_1 * upper_2_0
    )
    upper_2_1 = lower_diag_invrs_1 * upper_2_band[1]

    def lu_step(carry, x):
        upper_1_km1, upper_1_km2, upper_2_km1, upper_2_km2 = carry
        upper_2_lhs, upper_1_lhs, diag_lhs, lower_1_lhs, lower_2_lhs = x
        lower_2_k = lower_2_lhs
        lower_1_k = lower_1_lhs - lower_2_k * upper_1_km2
        lower_diag_invrs_k = 1.0 / (
            diag_lhs
            - lower_2_k * upper_2_km2
            - lower_1_k * upper_1_km1
        )
        upper_1_k = lower_diag_invrs_k * (
            upper_1_lhs - lower_1_k * upper_2_km1
        )
        upper_2_k = lower_diag_invrs_k * upper_2_lhs
        return (
            upper_1_k,
            upper_1_km1,
            upper_2_k,
            upper_2_km1,
        ), (
            lower_diag_invrs_k,
            lower_1_k,
            lower_2_k,
            upper_1_k,
            upper_2_k,
        )

    _, (
        lower_diag_invrs_rest,
        lower_1_rest,
        lower_2_rest,
        upper_1_rest,
        upper_2_rest,
    ) = jax.lax.scan(
        lu_step,
        (upper_1_1, upper_1_0, upper_2_1, upper_2_0),
        (
            upper_2_band[2:],
            upper_1_band[2:],
            diag_band[2:],
            lower_1_band[2:],
            lower_2_band[2:],
        ),
    )

    lower_diag_invrs = jnp.concatenate(
        [
            lower_diag_invrs_0[None],
            lower_diag_invrs_1[None],
            lower_diag_invrs_rest,
        ],
        axis=0,
    )
    lower_1 = jnp.concatenate(
        [lower_1_0[None], lower_1_1[None], lower_1_rest],
        axis=0,
    )
    lower_2 = jnp.concatenate(
        [lower_2_0[None], lower_2_1[None], lower_2_rest],
        axis=0,
    )
    upper_1 = jnp.concatenate(
        [upper_1_0[None], upper_1_1[None], upper_1_rest],
        axis=0,
    )
    upper_2 = jnp.concatenate(
        [upper_2_0[None], upper_2_1[None], upper_2_rest],
        axis=0,
    )

    soln_0 = lower_diag_invrs[0] * rhs_t[0]
    soln_1 = lower_diag_invrs[1] * (rhs_t[1] - lower_1[1] * soln_0)

    def forward_substitution(carry, x):
        soln_km2, soln_km1 = carry
        rhs_k, lower_diag_invrs_k, lower_1_k, lower_2_k = x
        soln_k = lower_diag_invrs_k * (
            rhs_k - lower_2_k * soln_km2 - lower_1_k * soln_km1
        )
        return (soln_km1, soln_k), soln_k

    _, soln_rest = jax.lax.scan(
        forward_substitution,
        (soln_0, soln_1),
        (rhs_t[2:], lower_diag_invrs[2:], lower_1[2:], lower_2[2:]),
    )
    soln = jnp.concatenate([soln_0[None], soln_1[None], soln_rest], axis=0)

    soln_ndim_minus_1 = (
        soln[ndim - 2] - upper_1[ndim - 2] * soln[ndim - 1]
    )

    def backward_substitution(carry, x):
        soln_kp1, soln_kp2 = carry
        soln_k_old, upper_1_k, upper_2_k = x
        soln_k = soln_k_old - upper_1_k * soln_kp1 - upper_2_k * soln_kp2
        return (soln_k, soln_kp1), soln_k

    _, soln_reversed = jax.lax.scan(
        backward_substitution,
        (soln_ndim_minus_1, soln[ndim - 1]),
        (
            soln[: ndim - 2][::-1],
            upper_1[: ndim - 2][::-1],
            upper_2[: ndim - 2][::-1],
        ),
    )
    return jnp.concatenate(
        [
            soln_reversed[::-1],
            soln_ndim_minus_1[None],
            soln[ndim - 1 : ndim],
        ],
        axis=0,
    ).T


@partial(jax.jit, static_argnames=("ndim", "nrhs", "ngrdcol"))
def penta_lu_solve_multiple_rhs_lhs(
    ndim: int, nrhs: int, ngrdcol: int, lhs, rhs
):
    """Solve multiple pentadiagonal systems with multiple right-hand sides.

    Description:
      Written for multiple RHS and multiple LHS.
    """
    del nrhs, ngrdcol

    # Matrices to solve, stored using band diagonal vectors
    # -2 is the uppermost band, 2 is the lower most band, 0 is diagonal
    upper_2_band = lhs[0].T
    upper_1_band = lhs[1].T
    diag_band = lhs[2].T
    lower_1_band = lhs[3].T
    lower_2_band = lhs[4].T
    rhs_t = jnp.transpose(rhs, (1, 0, 2))

    # ----------------------- Begin Code -----------------------

    # Inverse of the diagonal of L
    lower_diag_invrs_0 = 1.0 / diag_band[0]

    # First U band
    upper_1_0 = lower_diag_invrs_0 * upper_1_band[0]

    # Second U band
    upper_2_0 = lower_diag_invrs_0 * upper_2_band[0]

    # First L band
    lower_1_0 = jnp.zeros_like(lower_diag_invrs_0)

    # Second L band
    lower_2_0 = jnp.zeros_like(lower_diag_invrs_0)

    lower_1_1 = lower_1_band[1]
    lower_2_1 = jnp.zeros_like(lower_diag_invrs_0)
    lower_diag_invrs_1 = 1.0 / (
        diag_band[1] - lower_1_1 * upper_1_0
    )
    upper_1_1 = lower_diag_invrs_1 * (
        upper_1_band[1] - lower_1_1 * upper_2_0
    )
    upper_2_1 = lower_diag_invrs_1 * upper_2_band[1]

    def lu_step(carry, x):
        upper_1_km1, upper_1_km2, upper_2_km1, upper_2_km2 = carry
        upper_2_lhs, upper_1_lhs, diag_lhs, lower_1_lhs, lower_2_lhs = x
        lower_2_k = lower_2_lhs
        lower_1_k = lower_1_lhs - lower_2_k * upper_1_km2
        lower_diag_invrs_k = 1.0 / (
            diag_lhs
            - lower_2_k * upper_2_km2
            - lower_1_k * upper_1_km1
        )
        upper_1_k = lower_diag_invrs_k * (
            upper_1_lhs - lower_1_k * upper_2_km1
        )
        upper_2_k = lower_diag_invrs_k * upper_2_lhs
        return (
            upper_1_k,
            upper_1_km1,
            upper_2_k,
            upper_2_km1,
        ), (
            lower_diag_invrs_k,
            lower_1_k,
            lower_2_k,
            upper_1_k,
            upper_2_k,
        )

    _, (
        lower_diag_invrs_rest,
        lower_1_rest,
        lower_2_rest,
        upper_1_rest,
        upper_2_rest,
    ) = jax.lax.scan(
        lu_step,
        (upper_1_1, upper_1_0, upper_2_1, upper_2_0),
        (
            upper_2_band[2:],
            upper_1_band[2:],
            diag_band[2:],
            lower_1_band[2:],
            lower_2_band[2:],
        ),
    )

    lower_diag_invrs = jnp.concatenate(
        [
            lower_diag_invrs_0[None],
            lower_diag_invrs_1[None],
            lower_diag_invrs_rest,
        ],
        axis=0,
    )
    lower_1 = jnp.concatenate(
        [lower_1_0[None], lower_1_1[None], lower_1_rest],
        axis=0,
    )
    lower_2 = jnp.concatenate(
        [lower_2_0[None], lower_2_1[None], lower_2_rest],
        axis=0,
    )
    upper_1 = jnp.concatenate(
        [upper_1_0[None], upper_1_1[None], upper_1_rest],
        axis=0,
    )
    upper_2 = jnp.concatenate(
        [upper_2_0[None], upper_2_1[None], upper_2_rest],
        axis=0,
    )

    soln_0 = lower_diag_invrs[0, :, None] * rhs_t[0]
    soln_1 = lower_diag_invrs[1, :, None] * (
        rhs_t[1] - lower_1[1, :, None] * soln_0
    )

    def forward_substitution(carry, x):
        soln_km2, soln_km1 = carry
        rhs_k, lower_diag_invrs_k, lower_1_k, lower_2_k = x
        soln_k = lower_diag_invrs_k[:, None] * (
            rhs_k
            - lower_2_k[:, None] * soln_km2
            - lower_1_k[:, None] * soln_km1
        )
        return (soln_km1, soln_k), soln_k

    _, soln_rest = jax.lax.scan(
        forward_substitution,
        (soln_0, soln_1),
        (rhs_t[2:], lower_diag_invrs[2:], lower_1[2:], lower_2[2:]),
    )
    soln = jnp.concatenate([soln_0[None], soln_1[None], soln_rest], axis=0)

    soln_ndim_minus_1 = (
        soln[ndim - 2] - upper_1[ndim - 2, :, None] * soln[ndim - 1]
    )

    def backward_substitution(carry, x):
        soln_kp1, soln_kp2 = carry
        soln_k_old, upper_1_k, upper_2_k = x
        soln_k = (
            soln_k_old
            - upper_1_k[:, None] * soln_kp1
            - upper_2_k[:, None] * soln_kp2
        )
        return (soln_k, soln_kp1), soln_k

    _, soln_reversed = jax.lax.scan(
        backward_substitution,
        (soln_ndim_minus_1, soln[ndim - 1]),
        (
            soln[: ndim - 2][::-1],
            upper_1[: ndim - 2][::-1],
            upper_2[: ndim - 2][::-1],
        ),
    )
    return jnp.transpose(
        jnp.concatenate(
            [
                soln_reversed[::-1],
                soln_ndim_minus_1[None],
                soln[ndim - 1 : ndim],
            ],
            axis=0,
        ),
        (1, 0, 2),
    )


def penta_lu_solve(ndim: int, ngrdcol: int, lhs, rhs):
    """Generic pentadiagonal LU solve interface."""
    if rhs.ndim == 2:
        return penta_lu_solve_single_rhs_multiple_lhs(ndim, ngrdcol, lhs, rhs)
    if rhs.ndim == 3:
        return penta_lu_solve_multiple_rhs_lhs(
            ndim, rhs.shape[2], ngrdcol, lhs, rhs
        )
    raise ValueError("rhs must be rank-2 or rank-3 for penta_lu_solve")


__all__ = [
    "penta_lu_solve",
    "penta_lu_solve_single_rhs_multiple_lhs",
    "penta_lu_solve_multiple_rhs_lhs",
]
