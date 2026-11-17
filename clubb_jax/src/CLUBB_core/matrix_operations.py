"""JAX counterparts of the matrix helpers used by the hydrometeor PDF.

Cholesky_factor follows Fortran's scaling, lower-triangle factorization and
at most ten attempts, including under JIT. cholesky_lower returns the same
one-based failing-pivot status. Matrix order is inferred at the existing public
Cholesky_factor boundary; the private kernel retains the Fortran argument.
Independent rows are vectorized and the small, shape-static column/reduction
loops are unrolled. Host callbacks format diagnostics; NumPy eigenvalues are
used only in debug-level-2 failure dumps. Other unused Fortran matrix helpers
remain outside this port.
"""
import sys

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import core_rknd
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level


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


def Cholesky_factor(a_input):
    """Factor the lower triangle of a_input, preserving the upper triangle.
    Return reciprocal-diagonal scaling factors and whether scaling was applied.
    If the factorization fails we use a modified a_input matrix and attempt
    to factorize again.

    References:
    <http://www.netlib.org/lapack/explore-html/a00860.html> dpoequ
    <http://www.netlib.org/lapack/explore-html/a00753.html> dlaqsy

    Returns:
        (a_scaling, a_Cholesky, l_scaled)

    Fortran declarations (Cholesky_factor in src/CLUBB_core/matrix_operations.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ndim ! Matrix order [-]

        real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: a_input ! Input covariance or correlation matrix [units vary]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ndim), intent(out) :: a_scaling ! Reciprocal diagonal square roots [units vary]

        real( kind = core_rknd ), dimension(ndim,ndim), intent(out) :: a_Cholesky ! Factor in the lower triangle [units vary]

        logical, intent(out) :: l_scaled ! Whether diagonal equilibration was applied
    """
    ndim = a_input.shape[0]
    itermax = 10
    d_coef = 0.1

    # ---- Begin code ----
    a_scaled = a_input
    l_scaled = jnp.asarray(False)
    if ndim == 0:
        a_scaling = jnp.empty((0,), dtype=a_input.dtype)
        a_Cholesky = a_input
        return a_scaling, a_Cholesky, l_scaled

    # Compute reciprocal-diagonal scaling, retaining the diagonal on invalid input.
    a_scaling = jnp.diagonal(a_input)
    amin = jnp.min(a_scaling)
    amax = jnp.max(a_scaling)
    l_positive = amin > 0.0
    # Safe operands also protect inactive branches when vmap batches a condition.
    diagonal = jnp.where(l_positive, a_scaling, 1.0)
    a_scaling = jnp.where(l_positive, 1.0 / jnp.sqrt(diagonal), a_scaling)
    scond = jnp.sqrt(jnp.where(l_positive, amin, 1.0)) / jnp.sqrt(jnp.where(l_positive, amax, 1.0))
    small = jnp.finfo(core_rknd).tiny / jnp.finfo(core_rknd).eps
    large = 1.0 / small
    l_scaled = l_positive & ~((scond >= 0.1) & (amax >= small) & (amax <= large))
    indices = jnp.arange(ndim)
    lower = indices[:, None] >= indices[None, :]
    scaling = jnp.where(l_scaled, a_scaling, 1.0)
    a_scaled = jnp.where(lower, scaling[None, :] * scaling[:, None] * a_input, a_input)
    a_Cholesky = a_scaled

    # Host formatting is an adaptation of the Fortran write statements. The
    # predicate remains explicit because vmap can evaluate both cond branches.
    def report(a_input_host, iter, info, l_active, factored, modified):
        """Print the source's success, failed-pivot and retry diagnostics."""
        if not l_active:
            return
        if info == 0:
            if clubb_at_least_debug_level(1) and iter > 1:
                print("a_factored (worked)=", file=sys.stderr)
                for i in range(ndim):
                    print(" ".join(f"{value:10.3g}" for value in factored[i, :i+1]), file=sys.stderr)
            return
        if clubb_at_least_debug_level(1):
            print("Cholesky_factor: leading minor of order", info, "is not positive definite.", file=sys.stderr)
            print("factorization failed.", file=sys.stderr)
            for label, values in (("a_input=", a_input_host), ("a_Cholesky=", factored)):
                print(label, file=sys.stderr)
                for i in range(ndim):
                    print(" ".join(f"{value:10.3g}" for value in values[i, :i+1]), file=sys.stderr)
        if clubb_at_least_debug_level(2):
            import numpy as np
            original = np.tril(a_input_host) + np.tril(a_input_host, -1).T
            print("a_eigenvalues=", np.linalg.eigvalsh(original), file=sys.stderr)
            print("a_correlations=", file=sys.stderr)
            for i in range(ndim):
                print(" ".join(f"{original[i,j] / np.sqrt(original[i,i] * original[j,j]):10.3g}"
                               for j in range(i+1)), file=sys.stderr)
        if iter == itermax:
            print("iteration =", iter, "itermax =", itermax, file=sys.stderr)
            print("Fatal error in Cholesky_factor", file=sys.stderr)
        elif clubb_at_least_debug_level(1):
            print("Fatal error in Cholesky_factor", file=sys.stderr)
            print("Attempting to modify matrix to allow factorization.", file=sys.stderr)
        if clubb_at_least_debug_level(2):
            symmetric = np.tril(modified) + np.tril(modified, -1).T
            print("a_modified eigenvalues=", np.linalg.eigvalsh(symmetric), file=sys.stderr)

    def iteration(iter, state):
        """Carry Fortran's early-exit state through a fixed differentiable loop."""
        a_Cholesky, l_done = state
        factored, info = cholesky_lower(ndim, a_Cholesky)
        l_failed = info > 0
        # The number used for tau here is case specific to the Sigma covariance
        # matrix in the latin hypercube code and is not at all general.
        # Tau should be a number that is small relative to the other diagonal
        # elements of the matrix to have keep the error caused by modifying 'a' low.
        # -dschanen 30 Aug 2010
        d_smallest = jnp.min(jnp.diagonal(a_scaled))
        # Use the smallest element * d_coef * iteration
        tau = d_smallest * d_coef * jnp.asarray(iter, dtype=a_input.dtype)
        modified = a_scaled.at[indices, indices].add(tau)

        def emit(_):
            jax.debug.callback(report, a_input, iter, info, ~l_done, factored, modified, ordered=True)

        l_report = (~l_done) & (clubb_at_least_debug_level(1) | (l_failed & (iter == itermax)))
        jax.lax.cond(l_report, emit, lambda _: None, operand=None)
        candidate = jnp.where(l_failed, modified, factored)
        a_Cholesky = jnp.where(l_done, a_Cholesky, candidate)
        l_done = l_done | ~l_failed
        return a_Cholesky, l_done

    def advance(iter, state):
        """Skip completed factorizations while retaining static loop bounds for AD."""
        return jax.lax.cond(state[1], lambda state: state,
                            lambda state: iteration(iter, state), state)

    a_Cholesky, _ = jax.lax.fori_loop(1, itermax + 1, advance, (a_Cholesky, jnp.asarray(False)))
    return a_scaling, a_Cholesky, l_scaled


def cholesky_lower(ndim, a):
    """Factor a symmetric positive-definite matrix in place, reading and updating
    only its lower triangle. Return info = 0 on success, or the one-based index
    of the first non-positive or NaN pivot; earlier columns remain factored.

    References:
    A = L L^T, with each column computed from previously factored columns.

    Returns:
        (a, info)

    Fortran declarations (cholesky_lower in src/CLUBB_core/matrix_operations.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ndim ! Matrix order [-]

        !--------------------------- Input/Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ndim,ndim), intent(inout) :: &
          a ! Matrix on entry; lower-triangular factor on success [units vary]

        !--------------------------- Output Variables ---------------------------
        integer, intent(out) :: info ! Zero on success, failing pivot, or -1 for ndim < 0
    """
    info = jnp.asarray(0, dtype=jnp.int32)
    if ndim < 0:
        info = jnp.asarray(-1, dtype=jnp.int32)
        return a, info

    for j in range(ndim):
        l_active = info == 0
        value = jnp.where(l_active, a[j, j], 1.0)
        for k in range(j):
            entry = jnp.where(l_active, a[j, k], 0.0)
            value = value - entry * entry
        a = a.at[j, j].set(jnp.where(l_active, value, a[j, j]))
        l_failed = (value <= 0.0) | jnp.isnan(value)
        info = jnp.where(l_active & l_failed, j + 1, info)
        l_active = info == 0
        diagonal = jnp.sqrt(jnp.where(l_active, value, 1.0))
        a = a.at[j, j].set(jnp.where(l_active, diagonal, a[j, j]))

        # Rows below the pivot are independent; preserve k's reduction order.
        value = jnp.where(l_active, a[j+1:, j], 0.0)
        for k in range(j):
            entry = jnp.where(l_active, a[j+1:, k], 0.0)
            pivot_entry = jnp.where(l_active, a[j, k], 0.0)
            value = value - entry * pivot_entry
        a = a.at[j+1:, j].set(jnp.where(l_active, value / diagonal, a[j+1:, j]))
    return a, info
