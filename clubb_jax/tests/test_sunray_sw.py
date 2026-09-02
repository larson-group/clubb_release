#!/usr/bin/env python3
"""test_sunray_sw.py — validate the JAX simplified-shortwave-radiation port vs the f2py Fortran oracle.

`sunray_sw` (rad_lwsw_module.F90) is the two-stream shortwave flux for the simplified radiation scheme: per-layer
optical depth from rcm·rho, then a downward/upward two-stream solve to a column SW flux Frad_SW on momentum levels.

This routine was previously case-validated-only ("divergent f2py signature"). The divergence was illusory: the f2py
wrapper `f2py_sunray_sw(fs0, amu0, rho, rcm)` simply hides the radiation constants inside `parameters_radiation`
(eff_drop_radius/alvdr/gc/omega) and the grid inside `stored_grid`, and passes `amu0` as the JAX's `xi_abs` arg. Setting
those constants deterministically via `f2py_set_simplified_radiation_params` and feeding the SAME values + grid to the JAX
yields a clean bit-match. SKIPs if clubb_f2py / clubb_python are unbuilt. (iter 444)
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

from clubb_jax.src.Radiation.rad_lwsw_module import sunray_sw
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_RADIUS, _ALVDR, _GC, _OMEGA = 1.0e-5, 0.1, 0.85, 0.9965   # eff_drop_radius / alvdr / gc / omega


def test_jit_multicol():
    """The source's sequential vertical taupath compiles for multiple columns."""
    gr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = gr.nzt
    result = sunray_sw(
        _NG, nzt,
        jnp.full((_NG, nzt), 1.0e-4), jnp.full((_NG, nzt), 1.0),
        0.5, gr.dzt, gr.zm, gr.zt,
        _RADIUS, _ALVDR, _GC, 1000.0, _OMEGA, True,
    )
    assert result.shape == (_NG, nzt + 1)
    assert np.isfinite(np.asarray(result)).all()
    print("  sunray_sw direct JIT: multicolumn flux is finite  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py sunray_sw oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    # Pin the parameters_radiation constants the f2py wrapper reads (fs_values/cos_solar_zen are fixed dim(20)).
    z20 = np.zeros(20)
    clubb_f2py.f2py_set_simplified_radiation_params(0.0, 0.0, 0.0, _RADIUS, _ALVDR, _GC, _OMEGA,
                                                    False, True, True, 1, z20.copy(), z20.copy(), z20.copy())
    dzt = np.asarray(jgr.dzt); zm = np.asarray(jgr.zm); zt = np.asarray(jgr.zt)
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(15):
        rcm = rng.uniform(0.0, 1e-3, (ng, nzt)); rho = rng.uniform(0.5, 1.2, (ng, nzt))
        fs0 = float(rng.uniform(100.0, 400.0)); amu0 = float(rng.uniform(0.3, 0.9))
        j = np.asarray(sunray_sw(ng, nzt, rcm, rho, amu0, dzt, zm, zt, _RADIUS, _ALVDR, _GC, fs0, _OMEGA, True))
        f = np.asarray(clubb_f2py.f2py_sunray_sw(fs0, amu0, rho, rcm))
        worst = max(worst, float(np.max(np.abs(j - f) / np.maximum(np.abs(f), 1e-30))))
    assert worst < 1e-12, f"sunray_sw f2py rel mismatch {worst:.2e}"
    print(f"  f2py sunray_sw (simplified SW two-stream flux, 15 cases): rel-match, worst {worst:.2e}  PASS")


def main():
    print("test_sunray_sw:")
    test_jit_multicol()
    test_f2py_oracle()
    print("All sunray_sw checks PASSED")


if __name__ == "__main__":
    main()
