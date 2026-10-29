"""Check the JAX compute_gamma_Skw port (Skx_module.F90:compute_gamma_Skw)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.skx_module as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.Skx_module import compute_gamma_Skw

NG, NZM = 2, 9
NZT = NZM - 1
NPARAMS = 102
I_GC, I_GB, I_GCF = 57, 58, 59   # 1-based igamma_coef / coefb / coefc


def _params(gc, gb, gcf):
    p = np.zeros((NG, NPARAMS))
    p[:, I_GC - 1] = gc
    p[:, I_GB - 1] = gb
    p[:, I_GCF - 1] = gcf
    return p


def test_f2py_oracle():
    """compute_gamma_Skw matches Fortran on zm and zt levels, including degenerate coefficients.

    Fortran's l_gamma_Skw is a compile-time parameter (.true. in model_flags.F90), so only the
    Skw-dependent branch exists in the reference; the JAX constant branch has no Fortran counterpart.
    """
    rng = np.random.default_rng(4)
    worst = 0.0
    for params in (_params(0.32, 0.08, 1.2), _params(0.25, 0.25, 1.0)):   # varying; degenerate coefs
        for nz in (NZM, NZT):
            skw = rng.uniform(-3, 3, (NG, nz))
            worst = max(worst, max_abs_diff("compute_gamma_Skw", compute_gamma_Skw(nz, NG, skw, params, True),
                                            fortran.compute_gamma_skw(nz, NG, skw, params)))
    assert worst < 1e-12, f"compute_gamma_Skw Fortran mismatch {worst:.2e}"
