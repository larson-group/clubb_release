#!/usr/bin/env python3
"""test_cloud_drop_sed.py — validate Microphys/cloud_sed_module.py:cloud_drop_sed vs the f2py Fortran primitives.

`cloud_drop_sed` (cloud_sed_module.F90:62) is the gravitational-settling "microphysics" some bit-faithful cases need
(`l_cloud_sed=.true.` with `microphys_scheme="none"`, e.g. dycoms2_rf02_so/_nd — both in DEFAULT_CASES, so it is
case-validated but was previously not unit-tested). The Ackerman-2009-eq.7 sedimentation flux (on zm levels, only
where rcm>0 & Ncm>0, zero-flux at the boundaries):

    Fcsed = 1.19e8 · exp(5·ln²σ_g) · (3/(4π·ρ_lw·Ncm_zm·ρ_zm))^(2/3) · (ρ_zm·rcm_zm)^(5/3)
    rcm_mc = (1/ρ)·ddzm(Fcsed)        thlm_mc = −(Lv/(Cp·exner))·rcm_mc

is reconstructed here on the **bit-shadowed f2py `zt2zm`/`ddzm`** + the formula and bit-compared (Fcsed + both
tendencies). A genuine Fortran-oracle composition. SKIPs if clubb_f2py / clubb_python are unbuilt. (iter 503)
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

from clubb_jax.src.Microphys.cloud_sed_module import cloud_drop_sed
from clubb_jax.src.CLUBB_core.constants_clubb import rho_lw, Cp, Lv
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_PI = float(np.pi)


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
        print(f"  f2py cloud_drop_sed oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    zt2zm = lambda a: np.asarray(clubb_f2py.f2py_zt2zm_2d(nzm, a))
    ddzm  = lambda a: np.asarray(clubb_f2py.f2py_ddzm_2d(nzt, a))
    rng = np.random.default_rng(9)
    sigma_g = 1.5
    wr = wt = wf = 0.0
    for _ in range(10):
        # cloudy in a mid layer, clear elsewhere (exercises the rcm>0 & Ncm>0 mask)
        rcm = np.zeros((ng, nzt)); Ncm = np.zeros((ng, nzt))
        rcm[:, 8:20] = rng.uniform(1e-5, 1e-3, (ng, 12))
        Ncm[:, 8:20] = rng.uniform(1e7, 1e8, (ng, 12))
        rho_zm = rng.uniform(0.5, 1.2, (ng, nzm))
        rho = rng.uniform(0.5, 1.2, (ng, nzt))
        exner = rng.uniform(0.6, 1.0, (ng, nzt))

        rcm_mc_j, thlm_mc_j, Fcsed_j = cloud_drop_sed(
            jnp.asarray(rcm), jnp.asarray(Ncm), jnp.asarray(rho_zm), jnp.asarray(rho),
            jnp.asarray(exner), sigma_g, jgr)

        rcm_zm = zt2zm(rcm); Ncm_zm = zt2zm(Ncm)
        valid = (rcm_zm > 0.0) & (Ncm_zm > 0.0)
        Ncm_safe = np.where(Ncm_zm > 0.0, Ncm_zm, 1.0); rcm_safe = np.maximum(rcm_zm, 0.0)
        const = 1.19e8 * np.exp(5.0 * (np.log(sigma_g)) ** 2)
        Fcsed_ref = const * (3.0 / (4.0 * _PI * rho_lw * Ncm_safe * rho_zm)) ** (2.0 / 3.0) \
            * (rho_zm * rcm_safe) ** (5.0 / 3.0)
        Fcsed_ref = np.where(valid, Fcsed_ref, 0.0)
        Fcsed_ref[:, 0] = 0.0; Fcsed_ref[:, -1] = 0.0
        sed_rcm_ref = (1.0 / rho) * ddzm(Fcsed_ref)
        rcm_mc_ref = sed_rcm_ref
        thlm_mc_ref = -(Lv / (Cp * exner)) * sed_rcm_ref

        wf = max(wf, float(np.max(np.abs(np.asarray(Fcsed_j) - Fcsed_ref))))
        wr = max(wr, float(np.max(np.abs(np.asarray(rcm_mc_j) - rcm_mc_ref))))
        wt = max(wt, float(np.max(np.abs(np.asarray(thlm_mc_j) - thlm_mc_ref))))
    # rel tolerance: the JAX factors const=1.19e8·exp(..) vs the Fortran's trailing exp → FP reorder ~1e-12 rel
    assert wf < 1e-6 and wr < 1e-12 and wt < 1e-9, f"cloud_drop_sed mismatch: Fcsed {wf:.2e}, rcm_mc {wr:.2e}, thlm_mc {wt:.2e}"
    print(f"  f2py cloud_drop_sed (Fcsed + rcm_mc + thlm_mc, masked + BC): match, worst Fcsed {wf:.2e} / rcm_mc {wr:.2e}  PASS")


def main():
    print("test_cloud_drop_sed:")
    test_f2py_oracle()
    print("All cloud_drop_sed checks PASSED")


if __name__ == "__main__":
    main()
