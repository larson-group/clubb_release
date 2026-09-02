"""JAX port of CLUBB_core/diagnose_correlations_module.F90.

Porting deviations:
- ``iiPDF_w`` follows the Fortran-facing API and is 1-based.  Callers using
  JAX metadata must pass ``hm_metadata.iiPDF_w + 1``.
- The Fortran ``approx_w_corr`` and ``approx_w_covar`` routines are commented
  out in the source and are not ported.  Therefore ``l_calc_w_corr=True``
  raises ``NotImplementedError`` instead of silently taking an incomplete path.
- Fortran intent(out/inout) arguments are returned as values.
- Fortran error-stops in correlation assertion checks are exposed as concrete
  boolean helpers instead of host-side side effects inside JAX array code.
"""

import numpy as np
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import max_mag_correlation


_SQRT_FLOOR = jnp.finfo(jnp.float64).tiny


def _safe_corr_sqrt(value):
    """sqrt(max(value, 0)) with finite gradient at exactly zero."""
    value = jnp.asarray(value, dtype=jnp.float64)
    return jnp.where(value > 0.0, jnp.sqrt(jnp.maximum(value, _SQRT_FLOOR)), 0.0)


def diagnose_correlations(pdf_dim, iiPDF_w, corr_array_pre, l_calc_w_corr=False):
    """This subroutine diagnoses the correlation matrix in order to feed it
    into SILHS microphysics.

    References:
      Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
      (see CLUBB Trac ticket#514)
    """
    if l_calc_w_corr:
        raise NotImplementedError(
            "diagnose_correlations: l_calc_w_corr=True is not ported; the "
            "underlying Fortran approx_w_corr path is commented out."
        )

    # Initialize sigma2_on_mu2_ip_array

    # Swap the w-correlations to the first row for the prescribed correlations
    corr_array_pre_swapped = rearrange_corr_array(pdf_dim, iiPDF_w, corr_array_pre)
    corr_array_swapped = corr_array_pre_swapped

    # diagnose correlations
    corr_array_swapped = diagnose_corr(
        pdf_dim,
        jnp.zeros((pdf_dim,), dtype=jnp.float64),
        corr_array_pre_swapped,
        corr_array_swapped,
    )
    # Swap rows back
    return rearrange_corr_array(pdf_dim, iiPDF_w, corr_array_swapped)


def diagnose_corr(n_variables, sqrt_sigma2_on_mu2_ip,
                  corr_matrix_prescribed, corr_matrix_approx):
    """This subroutine diagnoses the correlation matrix for each timestep.

    References:
      Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
      (see CLUBB Trac ticket#514)
    """
    del sqrt_sigma2_on_mu2_ip

    # calculate all square roots
    s_1j = _safe_corr_sqrt(1.0 - corr_matrix_approx[:, 0] ** 2)
    f_ij = jnp.clip(
        corr_matrix_prescribed,
        -max_mag_correlation,
        max_mag_correlation,
    )
    diagnosed = (
        # formula (15) in the ref. paper (Larson et al. (2011))
        corr_matrix_approx[:, 0, None] * corr_matrix_approx[:, 0][None, :]
        + f_ij * s_1j[:, None] * s_1j[None, :]
    )

    rows = jnp.arange(n_variables)[:, None]
    cols = jnp.arange(n_variables)[None, :]
    mask = (rows > cols) & (cols >= 1) & (cols <= n_variables - 2)
    return jnp.where(mask, diagnosed, corr_matrix_approx)


def calc_w_corr(wpxp, stdev_w, stdev_x, w_tol, x_tol):
    """Compute the correlations of w with the hydrometeors.

    References:
      clubb:ticket:514
    """
    calc_w_corr_value = wpxp / (
        jnp.maximum(stdev_x, x_tol) * jnp.maximum(stdev_w, w_tol)
    )
    # Make sure the correlation is in [-1,1]
    return jnp.clip(
        calc_w_corr_value,
        -max_mag_correlation,
        max_mag_correlation,
    )


def calc_varnce(mixt_frac, x1, x2, xm, x1p2, x2p2):
    """Calculate the variance xp2 from the components x1, x2.

    References:
      Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
      page 3535
    """
    return (
        mixt_frac * ((x1 - xm) ** 2 + x1p2)
        + (1.0 - mixt_frac) * ((x2 - xm) ** 2 + x2p2)
    )


def calc_mean(mixt_frac, x1, x2):
    """Calculate the mean xm from the components x1, x2.

    References:
      Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
      page 3535
    """
    return mixt_frac * x1 + (1.0 - mixt_frac) * x2


def calc_cholesky_corr_mtx_approx(n_variables, iiPDF_w, corr_matrix):
    """This subroutine calculates the transposed correlation cholesky matrix
    from the correlation matrix

    References:
      1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
      2 CLUBB Trac ticket#514
    """
    corr_mtx_swap = rearrange_corr_array(n_variables, iiPDF_w, corr_matrix)
    corr_cholesky_mtx_swap = setup_corr_cholesky_mtx(
        n_variables,
        corr_mtx_swap,
    )
    corr_cholesky_mtx = rearrange_corr_array(
        n_variables,
        iiPDF_w,
        corr_cholesky_mtx_swap,
    )
    corr_mtx_approx_swap = cholesky_to_corr_mtx_approx(
        n_variables,
        corr_cholesky_mtx_swap,
    )
    corr_mtx_approx = rearrange_corr_array(
        n_variables,
        iiPDF_w,
        corr_mtx_approx_swap,
    )

    # Fortran error-stops in corr_array_assertion_checks here. The JAX port
    # keeps that as an explicit concrete helper rather than adding host-side
    # side effects to this jittable routine.
    # Set lower triangle to zero for conformity
    rows = jnp.arange(n_variables)[:, None]
    cols = jnp.arange(n_variables)[None, :]
    corr_mtx_approx = jnp.where(rows < cols, 0.0, corr_mtx_approx)
    return corr_cholesky_mtx, corr_mtx_approx


def setup_corr_cholesky_mtx(n_variables, corr_matrix):
    """This subroutine calculates the transposed correlation cholesky matrix
    from the correlation matrix

    References:
      1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
      2 CLUBB Trac ticket#514
    """
    # calculate all necessary square roots
    s = _safe_corr_sqrt(1.0 - corr_matrix ** 2)

    # calculate transposed correlation cholesky matrix; ref 1 formula 10

    # initialize matrix to zero
    rows = jnp.arange(n_variables)[:, None]
    cols = jnp.arange(n_variables)[None, :]
    # initialize upper triangle and diagonal to one
    corr_cholesky_mtx_t = jnp.where(rows >= cols, 1.0, 0.0)

    # set diagonal elements
    for j in range(1, n_variables):
        diag_value = 1.0
        for i in range(j):
            diag_value = diag_value * s[j, i]
        corr_cholesky_mtx_t = corr_cholesky_mtx_t.at[j, j].set(diag_value)

    # set first row
    for j in range(1, n_variables):
        corr_cholesky_mtx_t = corr_cholesky_mtx_t.at[j, 0].set(
            corr_matrix[j, 0]
        )

    # set upper triangle
    for i in range(1, n_variables - 1):
        for j in range(i + 1, n_variables):
            upper_value = 1.0
            for k in range(i):
                upper_value = upper_value * s[j, k]
            upper_value = upper_value * corr_matrix[j, i]
            corr_cholesky_mtx_t = corr_cholesky_mtx_t.at[j, i].set(upper_value)

    return corr_cholesky_mtx_t


def cholesky_to_corr_mtx_approx(n_variables, corr_cholesky_mtx_t=None):
    """This subroutine approximates the correlation matrix from the correlation
    cholesky matrix

    References:
      1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
      2 CLUBB Trac ticket#514
    """
    if corr_cholesky_mtx_t is None:
        corr_cholesky_mtx_t = n_variables
    else:
        del n_variables
    return corr_cholesky_mtx_t @ corr_cholesky_mtx_t.T


def corr_array_assertion_checks(n_variables, corr_array=None):
    """Check the correlation matrix elements and diagonal.

    The Fortran routine stops the program on failure; this JAX/Python helper
    returns ``False`` instead.
    """
    if corr_array is None:
        corr_array = n_variables
        n_variables = np.asarray(corr_array).shape[0]

    corr_array = np.asarray(corr_array, dtype=np.float64)
    off_diagonal = corr_array[~np.eye(n_variables, dtype=bool)]
    off_diagonal_in_range = bool(
        np.all(np.abs(off_diagonal) <= max_mag_correlation)
    )
    diagonal_is_one = bool(
        np.all(np.abs(np.diagonal(corr_array) - 1.0) <= 1.0e-6)
    )
    return off_diagonal_in_range and diagonal_is_one


def rearrange_corr_array(pdf_dim, iiPDF_w, corr_array):
    """Swap the w-correlations to the first row."""
    w_idx = int(iiPDF_w) - 1
    corr_array = jnp.asarray(corr_array, dtype=jnp.float64)
    swap_array = corr_array[:, 0]

    corr_array_swapped = corr_array
    corr_array_swapped = corr_array_swapped.at[:w_idx + 1, 0].set(
        corr_array[w_idx, w_idx::-1]
    )
    corr_array_swapped = corr_array_swapped.at[w_idx + 1:pdf_dim, 0].set(
        corr_array[w_idx + 1:pdf_dim, w_idx]
    )
    corr_array_swapped = corr_array_swapped.at[w_idx, :w_idx + 1].set(
        swap_array[w_idx::-1]
    )
    corr_array_swapped = corr_array_swapped.at[w_idx + 1:pdf_dim, w_idx].set(
        swap_array[w_idx + 1:pdf_dim]
    )
    return corr_array_swapped
