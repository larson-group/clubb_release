"""Check Microphys/cloud_sed_module.py:cloud_drop_sed against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.microphys as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_rel_diff, no_stats, setup_matching_grids
from clubb_jax.src.Microphys.cloud_sed_module import cloud_drop_sed

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_SIGMA_G = 1.5


def test_f2py_oracle():
    """cloud_drop_sed rcm/thlm tendencies match Fortran (relative error; rcm_mc is O(1e-6)).

    The sedimentation flux Fcsed is zero at both boundaries and rcm_mc is its
    vertical divergence, so matching rcm_mc also pins Fcsed on every level."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(9)
    worst = {}
    for _ in range(10):
        # Cloudy in a mid layer and clear elsewhere, so both sides of the rcm>0 and Ncm>0 mask run.
        rcm = np.zeros((ng, nzt)); Ncm = np.zeros((ng, nzt))
        rcm[:, 8:20] = rng.uniform(1e-5, 1e-3, (ng, 12))
        Ncm[:, 8:20] = rng.uniform(1e7, 1e8, (ng, 12))
        rho_zm = rng.uniform(0.5, 1.2, (ng, nzm))
        rho = rng.uniform(0.5, 1.2, (ng, nzt))
        exner = rng.uniform(0.6, 1.0, (ng, nzt))
        # Nonzero incoming tendencies check that the routine adds to them.
        rcm_mc = rng.uniform(-1e-8, 1e-8, (ng, nzt))
        thlm_mc = rng.uniform(-1e-5, 1e-5, (ng, nzt))
        _stats, *j_out = cloud_drop_sed(
            jgr, ng, *(jnp.asarray(x) for x in (rcm, Ncm, rho_zm, rho, exner)), _SIGMA_G,
            no_stats(ng), jnp.asarray(rcm_mc), jnp.asarray(thlm_mc))
        f_out = fortran.cloud_drop_sed(nzm, nzt, ng, fgr, rcm, Ncm, rho_zm, rho, exner, _SIGMA_G,
                                       rcm_mc.copy(), thlm_mc.copy())
        assert len(j_out) == len(f_out) == 2, f"output counts: JAX {len(j_out)}, Fortran {len(f_out)}"
        for name, j, f in zip(("rcm_mc", "thlm_mc"), j_out, f_out, strict=True):
            worst[name] = max(worst.get(name, 0.0), max_rel_diff(name, j, f, 1e-300))
    bad = {k: v for k, v in worst.items() if not v < 1e-12}
    assert not bad, f"cloud_drop_sed Fortran mismatch: {bad}"
