"""Check the JAX Nc_in_cloud_to_Ncnm port (Nc_Ncn_eqns.F90) against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_rel_diff
from clubb_jax.src.CLUBB_core.Nc_Ncn_eqns import Nc_in_cloud_to_Ncnm


def test_Ncnm_bit_to_bit_vs_f2py():
    """Nc_in_cloud_to_Ncnm matches the Fortran f2py oracle over a random sweep."""
    import clubb_f2py as f
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(3000):
        mu1, mu2 = rng.uniform(-1e-3, 1e-3), rng.uniform(-1e-3, 1e-3)
        s1 = rng.choice([0.0, rng.uniform(0, 5e-4)])      # exercise sigma_chi=0 branch
        s2 = rng.choice([0.0, rng.uniform(0, 5e-4)])
        mf, ncl = rng.uniform(0.01, 0.99), rng.uniform(1e6, 3e8)
        cf1, cf2 = rng.uniform(0, 1), rng.uniform(0, 1)
        cnp2 = rng.choice([0.0, rng.uniform(0, 2.0)])     # exercise const-Ncn vs varying
        ccc = rng.choice([0.0, rng.uniform(-0.9, 0.9)])
        fo = float(f.f2py_nc_in_cloud_to_ncnm(mu1, mu2, s1, s2, mf, ncl, cf1, cf2, cnp2, ccc))
        jo = float(Nc_in_cloud_to_Ncnm(mu1, mu2, s1, s2, mf, ncl, cf1, cf2, cnp2, ccc))
        worst = max(worst, max_rel_diff("Nc_in_cloud_to_Ncnm", jo, fo, 1e-300))
    # Machine-level (residual is the JAX-vs-Fortran erfc implementation difference).
    assert worst < 1e-12, f"Ncnm vs f2py worst rel {worst:.2e}"
