"""Check advance_xm_wpxp_module.py:diagnose_upxp against the Fortran formula."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.grid_class as fortran_grid

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_xm_wpxp_module import diagnose_upxp

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_f2py_oracle():
    """diagnose_upxp matches the Fortran formula (Andre et al. 1978) evaluated with Fortran's ddzt.

    There is no F2PY wrapper for diagnose_upxp itself, so the reference applies the Fortran
    expression to vertical derivatives computed by the Fortran grid operator."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(12):
        ypwp = rng.uniform(-0.5, 0.5, (ng, nzm))
        xm   = rng.uniform(280.0, 320.0, (ng, nzt))      # e.g. thlm (zt)
        wpxp = rng.uniform(-0.5, 0.5, (ng, nzm))
        ym   = rng.uniform(-15.0, 15.0, (ng, nzt))       # e.g. um (zt)
        C6x  = rng.uniform(2.0, 6.0, (ng, nzm))
        tau  = rng.uniform(100.0, 900.0, (ng, nzm))
        C7   = rng.uniform(0.3, 0.8, (ng, nzm))
        j = diagnose_upxp(nzm, nzt, ng, jgr, jnp.asarray(ypwp), jnp.asarray(xm), jnp.asarray(wpxp),
                          jnp.asarray(ym), jnp.asarray(C6x), jnp.asarray(tau), jnp.asarray(C7))
        ddzt_xm = np.asarray(fortran_grid.ddzt(nzm, nzt, ng, fgr, xm))
        ddzt_ym = np.asarray(fortran_grid.ddzt(nzm, nzt, ng, fgr, ym))
        # Fortran evaluates interior momentum levels and sets both boundary levels to zero.
        ref = np.zeros((ng, nzm))
        ref[:, 1:-1] = ((tau / C6x) * (-ypwp * ddzt_xm - (1.0 - C7) * (wpxp * ddzt_ym)))[:, 1:-1]
        worst = max(worst, max_abs_diff("diagnose_upxp", j, ref))
    assert worst < 1e-12, f"diagnose_upxp vs Fortran-ddzt formula mismatch {worst:.2e}"
