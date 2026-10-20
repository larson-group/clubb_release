"""Check the JAX mirror_lower_triangular_matrix port (matrix_operations.F90)
against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.matrix_operations import mirror_lower_triangular_matrix


def test_f2py_oracle():
    import clubb_f2py
    rng = np.random.default_rng(9)
    worst = 0.0
    for nvars in (1, 2, 4, 7):
        m = rng.uniform(-3, 3, (nvars, nvars))
        # f2py mutates in-place AND returns; pass a fresh Fortran-ordered copy so the input array is untouched.
        f = np.asarray(clubb_f2py.f2py_mirror_lower_triangular_matrix(np.asfortranarray(m.copy())))
        g = np.asarray(mirror_lower_triangular_matrix(m))
        worst = max(worst, max_abs_diff("mirror_lower_triangular", g, f))
    assert worst < 1e-15, f"mirror_lower_triangular_matrix f2py mismatch {worst:.2e}"
