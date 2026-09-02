#!/usr/bin/env python3
"""test_pvertinterp.py — validate the JAX pvertinterp port (advance_helper_module.F90:2072).

pvertinterp interpolates a field on pressure levels to a single output pressure, clamping to the bottom/top
value outside the column. Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_pvertinterp on a stored grid matching the JAX grid — interior + both
     clamp branches. SKIPs if clubb_f2py/clubb_python are unbuilt.
  2. Pressure-weighted linear-interpolation identity at an interior target.
  3. A finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import pvertinterp
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _p_profile(nzt):
    # Descending pressure with height (index): surface ~1000 hPa down to ~870 hPa over the column.
    return np.tile(np.linspace(100000.0, 87000.0, nzt)[::-1][None, :], (_NG, 1)).copy()


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py pvertinterp oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(4)
    p_mid = np.asfortranarray(_p_profile(nzt))
    input_var = np.asfortranarray(rng.standard_normal((ng, nzt)))
    worst = 0.0
    # An interior pressure, plus targets outside both ends (exercise the clamp branches).
    for p_out in (95000.0, 88000.0, 101000.0, 80000.0):
        ref = np.asarray(clubb_f2py.f2py_pvertinterp(p_mid.copy(), float(p_out), input_var.copy()))
        got = np.asarray(pvertinterp(p_mid, p_out, input_var))
        worst = max(worst, np.max(np.abs(got - ref)))
    assert worst < 1e-11, f"pvertinterp f2py mismatch {worst:.2e}"
    print(f"  f2py pvertinterp: match (interior + both clamps), worst {worst:.2e}  PASS")


def test_pressure_weighted_identity():
    nzt = 10
    p_mid = np.tile(np.linspace(100000.0, 80000.0, nzt)[::-1][None, :], (1, 1)).copy()
    # p_mid descending in k: index 0 = 100000, index nzt-1 = 80000.
    p_mid = np.linspace(100000.0, 80000.0, nzt)[None, :]
    input_var = np.random.default_rng(2).standard_normal((1, nzt))
    p_out = 91234.0
    got = float(np.asarray(pvertinterp(p_mid, p_out, input_var))[0])
    # Find the bracket (p_mid descending) and apply the Fortran pressure-weighted formula.
    k = next(i for i in range(nzt - 1) if p_mid[0, i] > p_out >= p_mid[0, i + 1])
    dpu = p_mid[0, k] - p_out
    dpl = p_out - p_mid[0, k + 1]
    ref = (input_var[0, k] * dpl + input_var[0, k + 1] * dpu) / (dpl + dpu)
    assert abs(got - ref) < 1e-12, f"pressure-weighted mismatch {abs(got-ref):.2e}"
    print("  pressure-weighted linear interpolation identity  PASS")


def test_differentiable():
    nzt = 12
    p_mid = jnp.asarray(np.linspace(100000.0, 80000.0, nzt)[None, :])
    def loss(v):
        return jnp.sum(pvertinterp(p_mid, 90000.0, v) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(1).standard_normal((1, nzt)))))
    assert np.isfinite(g).all(), "non-finite grad through pvertinterp"
    print(f"  jax.grad through pvertinterp: finite ({g.size} entries)  PASS")


def main():
    print("test_pvertinterp:")
    for t in (test_f2py_oracle, test_pressure_weighted_identity, test_differentiable):
        t()
    print("All pvertinterp checks PASSED")


if __name__ == "__main__":
    main()
