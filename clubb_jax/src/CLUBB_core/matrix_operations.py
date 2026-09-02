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
"""
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_ITERMAX = 10      # matrix_operations.F90:137
_D_COEF = 0.1      # matrix_operations.F90:139
_THRESH = 0.1      # LAPACK dlaqsy equilibration threshold


def _lower_mask(n):
    """Boolean (n,n) mask, True on and below the diagonal."""
    r = jnp.arange(n)
    return r[:, None] >= r[None, :]


def mirror_lower_triangular_matrix(matrix):
    """Mirror the lower triangular entries to the upper triangle.

    Description:
      Mirrors the elements of a lower triangular matrix to the upper
      triangle so that it is symmetric.

    References:
      None
    """
    m = jnp.asarray(matrix, dtype=jnp.float64)

    #----- Begin Code -----
    lower = jnp.tril(m)
    return lower + jnp.tril(m, -1).T


def Cholesky_factor(a_input):
    """Create a Cholesky factorization of a_input.

    Description:
      Create a Cholesky factorization of a_input.
      If the factorization fails we use a modified a_input matrix and attempt
      to factorize again.

    References:
      <http://www.netlib.org/lapack/explore-html/a00868.html> dpotrf
      <http://www.netlib.org/lapack/explore-html/a00860.html> dpoequ
      <http://www.netlib.org/lapack/explore-html/a00753.html> dlaqsy
    """
    a = jnp.asarray(a_input, dtype=jnp.float64)
    n = a.shape[0]
    diag = jnp.diagonal(a)

    # ---- Begin code ----

    # Copy input array into output array

    # Compute scaling for a_input, using Lapack routine spoequ for single precision,
    # or dpoequ for double precision
    a_scaling = 1.0 / jnp.sqrt(diag)
    amax = jnp.max(diag)
    amin = jnp.min(diag)
    scond = jnp.sqrt(amin / amax)

    fi = jnp.finfo(jnp.float64)
    small, large = fi.tiny, 1.0 / fi.tiny

    # Apply scaling to a_input, using Lapack routine slaqsy for single precision,
    # or dlaqsy for double precision
    l_scaled = ~((scond >= _THRESH) & (amax >= small) & (amax <= large))

    a_scaled = a_scaling[:, None] * a * a_scaling[None, :]
    a_work = jnp.where(l_scaled, a_scaled, a)

    mask = _lower_mask(n)

    def _factor(m):
        # Lapack Cholesky factorization, spotrf for single or dpotrf for double precision
        L = jnp.linalg.cholesky(m)
        return jnp.where(mask, L, m)

    a_cholesky = _factor(a_work)

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
            tau = d_smallest * _D_COEF * it
            m = a_work + tau * jnp.eye(n)
            cand = _factor(m)
            if not bool(jnp.isnan(jnp.diagonal(cand)).any()):
                a_cholesky = cand
                break
            a_cholesky = cand

    return a_scaling, a_cholesky, l_scaled
