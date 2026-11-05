"""JAX implementations of selected routines from ``matrix_operations.F90``.

Porting deviations:
The JAX code implements only the routines currently used by the JAX PDF path:
``Cholesky_factor`` and ``mirror_lower_triangular_matrix``.  Fortran helpers
``symm_covar_matrix_2_corr_matrix``, ``row_mult_lower_tri_matrix``,
``set_lower_triangular_matrix``, ``get_lower_triangular_matrix``,
``print_lower_triangular_matrix``, and ``Symm_matrix_eigenvalues`` are omitted
because no JAX source currently calls them.  LAPACK ``poequ``/``laqsy``/``potrf``
are represented with JAX/NumPy-equivalent algebra.  The tau-on-diagonal fallback
only triggers on a concrete non-positive-definite input and is skipped under a
JAX trace; the positive-definite path remains differentiable.
With ``CLUBB_JAX_PORTABLE_CHOLESKY=1``, small PDF matrices use a native JAX
Cholesky factorization instead of the vendor's batched solver kernel. Scaling,
dtype, lower-triangle ownership, and the concrete diagonal-retry path are retained.

Observed Fortran differences retained in this diagnostic pass:
- LAPACK argument status/failing-minor numbers and eigenvalue/correlation diagnostic helpers are absent; only available matrices and retry messages are printed.
- The concrete retry path allows an initial factorization plus ten retries;
  Fortran limits its factorization loop to ten attempts total.
"""
import os
import sys
from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level

import jax
import jax.numpy as jnp

_ITERMAX = 10      # matrix_operations.F90:137
_D_COEF = 0.1      # matrix_operations.F90:139
_THRESH = 0.1      # LAPACK dlaqsy equilibration threshold


def _lower_mask(n):
    """Boolean (n,n) mask, True on and below the diagonal."""
    r = jnp.arange(n)
    l_lower = r[:, None] >= r[None, :]
    return l_lower


def mirror_lower_triangular_matrix(matrix):
    """Mirrors the elements of a lower triangular matrix to the upper
    triangle so that it is symmetric.

    Returns:
        matrix

    Fortran declarations (mirror_lower_triangular_matrix in src/CLUBB_core/matrix_operations.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nvars ! Number of variables in each dimension of square matrix

        !--------------------------- Input/Output Variables ---------------------------
        real( kind = core_rknd ), dimension(nvars,nvars), intent(inout) :: &
          matrix  ! Lower triangluar square matrix
    """
    #----- Begin Code -----
    lower = jnp.tril(matrix)
    matrix = lower + jnp.tril(matrix, -1).T
    return matrix


def _portable_cholesky(matrix):
    matrix = (matrix + matrix.T) * 0.5
    factor = jnp.zeros_like(matrix)
    for column in range(matrix.shape[0]):
        diagonal = jnp.sqrt(matrix[column, column] - jnp.sum(factor[column, :column] ** 2))
        factor = factor.at[column, column].set(diagonal)
        product = jnp.sum(factor[column + 1:, :column] * factor[column, :column], axis=-1)
        entries = (matrix[column + 1:, column] - product) / diagonal
        factor = factor.at[column + 1:, column].set(entries)
    diagonal = jnp.diagonal(factor)
    factor = jnp.where(jnp.any(jnp.isnan(diagonal) | (diagonal <= 0.0)), jnp.nan, factor)
    return factor


def Cholesky_factor(a_input):
    """Create a Cholesky factorization of a_input.
    If the factorization fails we use a modified a_input matrix and attempt
    to factorize again.

    References:

    Returns:
        (a_scaling, a_cholesky, l_scaled)

    Fortran declarations (Cholesky_factor in src/CLUBB_core/matrix_operations.F90):

           ! Coefficient applied if the decomposition doesn't work

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ndim

        real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: a_input

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ndim), intent(out) :: a_scaling

        real( kind = core_rknd ), dimension(ndim,ndim), intent(out) :: a_Cholesky

        logical, intent(out) :: l_scaled
    """
    n = a_input.shape[0]
    diag = jnp.diagonal(a_input)

    # ---- Begin code ----

    # Copy input array into output array

    # Compute scaling for a_input, using Lapack routine spoequ for single precision,
    # or dpoequ for double precision
    a_scaling = 1.0 / jnp.sqrt(diag)
    amax = jnp.max(diag)
    amin = jnp.min(diag)
    scond = jnp.sqrt(amin / amax)

    fi = jnp.finfo(core_rknd)
    small, large = fi.tiny, 1.0 / fi.tiny

    # Apply scaling to a_input, using Lapack routine slaqsy for single precision,
    # or dlaqsy for double precision
    l_scaled = ~((scond >= _THRESH) & (amax >= small) & (amax <= large))

    a_scaled = a_scaling[:, None] * a_input * a_scaling[None, :]
    a_work = jnp.where(l_scaled, a_scaled, a_input)

    mask = _lower_mask(n)

    def _factor(m):
        # Lapack Cholesky factorization, spotrf for single or dpotrf for double precision
        if os.environ.get("CLUBB_JAX_PORTABLE_CHOLESKY") == "1":
            factor = _portable_cholesky(m)
        else:
            factor = jnp.linalg.cholesky(m)
        return jnp.where(mask, factor, m)

    a_cholesky = _factor(a_work)
    if clubb_at_least_debug_level(1):
        # Host formatting only; LAPACK's failing-minor info is unavailable in JAX.
        def print_cholesky_error(original, factored):
            print("factorization failed.", file=sys.stderr)
            for label, values in (("a_input=", original), ("a_Cholesky=", factored)):
                print(label, file=sys.stderr)
                for i in range(values.shape[0]):
                    print(" ".join(f"{value:10.3g}" for value in values[i, :i+1]), file=sys.stderr)
        jax.lax.cond(
            jnp.any(jnp.isnan(jnp.diagonal(a_cholesky))),
            lambda _: jax.debug.callback(print_cholesky_error, a_input, a_cholesky, ordered=True),
            lambda _: None, operand=None,
        )
    # TODO: Skip LAPACK argument/minor info and eigenvalue/correlation dumps;
    # those diagnostic capabilities are absent and exceed print-only restoration.


    try:
        failed = bool(jnp.isnan(jnp.diagonal(a_cholesky)).any())
    except (jax.errors.TracerBoolConversionError, jax.errors.ConcretizationTypeError):
        failed = False
    if failed:
        d_smallest = float(jnp.min(jnp.diagonal(a_work)))
        for it in range(1, _ITERMAX + 1):
            # This shouldn't happen now that the s and t Mellor(chi/eta) elements have been
            # modified to never be perfectly correlated, but it's here just in case.
            # -dschanen 10 Sept 2010

            # The number used for tau here is case specific to the Sigma covariance
            # matrix in the latin hypercube code and is not at all general.
            # Tau should be a number that is small relative to the other diagonal
            # elements of the matrix to have keep the error caused by modifying 'a' low.
            # -dschanen 30 Aug 2010

            # Use the smallest element * d_coef * iteration
            if clubb_at_least_debug_level(1):
                print("Fatal error in Cholesky_factor", file=sys.stderr)
                print("Attempting to modify matrix to allow factorization.", file=sys.stderr)
            tau = d_smallest * _D_COEF * it
            m = a_work + tau * jnp.eye(n)
            cand = _factor(m)
            if not bool(jnp.isnan(jnp.diagonal(cand)).any()):
                a_cholesky = cand
                if clubb_at_least_debug_level(1):
                    print("a_factored (worked)=", file=sys.stderr)
                    for i in range(n):
                        print(" ".join(f"{float(value):10.3g}" for value in cand[i, :i+1]), file=sys.stderr)
                break
            a_cholesky = cand
            if it == _ITERMAX:
                print("iteration =", it, "itermax =", _ITERMAX, file=sys.stderr)
                print("Fatal error in Cholesky_factor", file=sys.stderr)

    return a_scaling, a_cholesky, l_scaled
