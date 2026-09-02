#!/usr/bin/env python3
"""test_hydrostatic.py — validate the hydrostatic-pressure ports vs the f2py Fortran oracle.

  init_pressure  (calc_pressure.F90)       — hydrostatic p / exner on zt + zm from thvm + surface pressure.
  hydrostatic    (hydrostatic_module.F90)  — same, plus the dry-static density rho / rho_zm.
Both build the column pressure by integrating the hydrostatic relation over the grid; the JAX takes the grid object,
the f2py uses the module-global grid set to the JAX grid heights. SKIPs if clubb_f2py / clubb_python are unbuilt.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.calc_pressure import init_pressure
from clubb_jax.src.Input_fields.hydrostatic_module import hydrostatic
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _setup_f2py_grid():
    import clubb_f2py
    from clubb_python import clubb_api
    from clubb_python.derived_types.err_info import ErrInfo
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    return clubb_f2py, jgr, ng, nzm


def test_f2py_oracle():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py hydrostatic oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    rng = np.random.default_rng(1)
    w_ip = w_hy = 0.0
    for _ in range(15):
        thvm = np.sort(rng.uniform(290.0, 330.0, (ng, nzt)), axis=1)[:, ::-1].copy()  # decreasing with height
        psfc = rng.uniform(9e4, 1.02e5, ng)
        ji = init_pressure(jnp.asarray(thvm), jnp.asarray(psfc), jgr)            # (p, exner, p_zm, exner_zm)
        fi = clubb_f2py.f2py_init_pressure(nzm, thvm, psfc)
        w_ip = max(w_ip, max(float(np.max(np.abs(np.asarray(ji[k]) - np.asarray(fi[k]))
                                          / np.maximum(np.abs(np.asarray(fi[k])), 1.0))) for k in range(4)))
        jh = hydrostatic(jnp.asarray(thvm), jnp.asarray(psfc), jgr)              # (p, p_zm, exner, exner_zm, rho, rho_zm)
        fh = clubb_f2py.f2py_hydrostatic(nzm, thvm, psfc)
        w_hy = max(w_hy, max(float(np.max(np.abs(np.asarray(jh[k]) - np.asarray(fh[k]))
                                          / np.maximum(np.abs(np.asarray(fh[k])), 1.0))) for k in range(6)))
    assert w_ip < 1e-12, f"init_pressure f2py rel mismatch {w_ip:.2e}"
    assert w_hy < 1e-12, f"hydrostatic f2py rel mismatch {w_hy:.2e}"
    print(f"  f2py init_pressure (4 outputs) + hydrostatic (6 outputs, 15 cases): rel-match, worst {max(w_ip, w_hy):.2e}  PASS")


def test_f2py_oracle_stretched():
    """Re-validate init_pressure + hydrostatic on a STRETCHED grid (non-uniform dz, like CGILS). The hydrostatic
    integration accumulates `p -= rho*g*dz` over the column, so a stretched grid exercises the integration over varying
    dz that the uniform test (test_f2py_oracle) doesn't. Same consistent-grid_type=2 setup as the stretched grid-operator
    test (a JAX grid_type=3 vs f2py grid_type=2 would define the grid differently — a setup trap). SKIPs if unbuilt.
    (iter 482)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py hydrostatic (stretched) oracle: SKIP ({type(e).__name__})")
        return
    nzt0 = 18
    dz = np.linspace(20.0, 120.0, nzt0)
    zt1 = np.concatenate([[10.0], 10.0 + np.cumsum(dz)])[:nzt0]
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=float(zt1[-1] + 50.0),
                     grid_type=2, thermodynamic_heights=np.tile(zt1, (_NG, 1)))
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    g = clubb_f2py.get_grid(ng, nzm, nzt)
    assert float(np.max(np.abs(np.asarray(g[0]) - np.asarray(jgr.zm)))) == 0.0, "stretched grid zm mismatch (setup)"
    rng = np.random.default_rng(1)
    w_ip = w_hy = 0.0
    for _ in range(12):
        thvm = np.sort(rng.uniform(290.0, 330.0, (ng, nzt)), axis=1)[:, ::-1].copy()
        psfc = rng.uniform(9e4, 1.02e5, ng)
        ji = init_pressure(jnp.asarray(thvm), jnp.asarray(psfc), jgr); fi = clubb_f2py.f2py_init_pressure(nzm, thvm, psfc)
        w_ip = max(w_ip, max(float(np.max(np.abs(np.asarray(ji[k]) - np.asarray(fi[k]))
                                          / np.maximum(np.abs(np.asarray(fi[k])), 1.0))) for k in range(4)))
        jh = hydrostatic(jnp.asarray(thvm), jnp.asarray(psfc), jgr); fh = clubb_f2py.f2py_hydrostatic(nzm, thvm, psfc)
        w_hy = max(w_hy, max(float(np.max(np.abs(np.asarray(jh[k]) - np.asarray(fh[k]))
                                          / np.maximum(np.abs(np.asarray(fh[k])), 1.0))) for k in range(6)))
    assert w_ip < 1e-12 and w_hy < 1e-12, f"stretched hydrostatic mismatch: init_pressure {w_ip:.2e}, hydrostatic {w_hy:.2e}"
    print(f"  f2py init_pressure + hydrostatic on a STRETCHED grid (12 cases): rel-match, worst {max(w_ip, w_hy):.2e}  PASS")


def main():
    print("test_hydrostatic:")
    test_f2py_oracle()
    test_f2py_oracle_stretched()
    print("All hydrostatic checks PASSED")


if __name__ == "__main__":
    main()
