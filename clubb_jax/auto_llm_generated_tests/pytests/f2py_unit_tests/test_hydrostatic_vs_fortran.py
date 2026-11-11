"""Check the hydrostatic-pressure ports against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.calc_pressure as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.calc_pressure import init_pressure
from clubb_jax.src.Input_fields.hydrostatic_module import hydrostatic

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _pressure_mismatch(jgr, fgr, rng, n_draws):
    """Return the worst relative JAX-vs-Fortran difference of init_pressure and hydrostatic outputs."""
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    w_ip = w_hy = 0.0
    for _ in range(n_draws):
        thvm = np.sort(rng.uniform(290.0, 330.0, (ng, nzt)), axis=1)[:, ::-1].copy()  # decreasing with height
        psfc = rng.uniform(9e4, 1.02e5, ng)
        ji = init_pressure(jnp.asarray(thvm), jnp.asarray(psfc), jgr)            # (p, exner, p_zm, exner_zm)
        fi = fortran.init_pressure(ng, fgr, thvm, psfc)
        w_ip = max(w_ip, max(max_rel_diff("hydrostatic", ji[k], fi[k], 1.0) for k in range(4)))
        jh = hydrostatic(jnp.asarray(thvm), jnp.asarray(psfc), jgr)              # (p, p_zm, exner, exner_zm, rho, rho_zm)
        fh = fortran.hydrostatic(ng, fgr, thvm, psfc)
        w_hy = max(w_hy, max(max_rel_diff("hydrostatic", jh[k], fh[k], 1.0) for k in range(6)))
    return w_ip, w_hy


def test_f2py_oracle():
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    w_ip, w_hy = _pressure_mismatch(jgr, fgr, np.random.default_rng(1), 15)
    assert w_ip < 1e-12, f"init_pressure Fortran rel mismatch {w_ip:.2e}"
    assert w_hy < 1e-12, f"hydrostatic Fortran rel mismatch {w_hy:.2e}"


def test_f2py_oracle_stretched():
    """Re-validate init_pressure + hydrostatic on a STRETCHED grid (non-uniform dz, like CGILS). The hydrostatic
    integration accumulates `p -= rho*g*dz` over the column, so a stretched grid exercises the integration over varying
    dz that the uniform test (test_f2py_oracle) doesn't. Both sides use grid_type=2, as in the stretched grid-operator
    test (a JAX grid_type=3 vs Fortran grid_type=2 would define the grid differently)."""
    nzt0 = 18
    dz = np.linspace(20.0, 120.0, nzt0)
    zt1 = np.concatenate([[10.0], 10.0 + np.cumsum(dz)])[:nzt0]
    jgr, fgr = setup_matching_grids(_NG, _DZ, float(zt1[-1] + 50.0), grid_type=2,
                                    thermodynamic_heights=np.tile(zt1, (_NG, 1)))
    assert max_abs_diff("hydrostatic", fgr.zm, jgr.zm) == 0.0, "stretched grid zm mismatch (setup)"
    w_ip, w_hy = _pressure_mismatch(jgr, fgr, np.random.default_rng(1), 12)
    assert w_ip < 1e-12 and w_hy < 1e-12, f"stretched hydrostatic mismatch: init_pressure {w_ip:.2e}, hydrostatic {w_hy:.2e}"
