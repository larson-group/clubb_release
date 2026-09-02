#!/usr/bin/env python3
"""test_mean_adv.py — validate the mean_adv.F90 ports (term_ma_zt_lhs / term_ma_zm_lhs) vs the f2py Fortran oracle.

term_ma_zt_lhs is the mean-advection LHS for a zt-level variable; iter 395 re-merged the JAX's previously-split
centered/upwind variants back into the single Fortran-named subroutine (internal l_upwind_xm_ma branch). This test
bit-shadows BOTH branches vs the oracle, plus term_ma_zm_lhs — confirming that merge is faithful. SKIPs if clubb_f2py
is unbuilt.
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

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs, term_ma_zm_lhs


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py mean_adv oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=2, deltaz=40.0, zm_init=0.0, zm_top=1200.0, grid_type=1)
    ng, nzm = jgr.zm.shape
    idzt = np.asarray(jgr.invrs_dzt); idzm = np.asarray(jgr.invrs_dzm)
    w_zt2zm = np.asarray(jgr.weights_zt2zm); w_zm2zt = np.asarray(jgr.weights_zm2zt); gd = float(jgr.grid_dir)
    rng = np.random.default_rng(1)
    wzt = wzm = 0.0
    for _ in range(20):
        wm_zt = rng.uniform(-2.0, 2.0, (ng, nzm - 1)); wm_zm = rng.uniform(-2.0, 2.0, (ng, nzm))
        for lup in (False, True):   # term_ma_zt_lhs centered + upwind branches (the iter-395 merge)
            j = np.asarray(term_ma_zt_lhs(
                jgr.nzm, jgr.nzt, jgr.ngrdcol, jnp.asarray(wm_zt),
                jgr.weights_zt2zm, jgr.invrs_dzt, jgr.invrs_dzm, lup, gd,
            ))
            f = np.asarray(clubb_f2py.f2py_term_ma_zt_lhs(wm_zt, w_zt2zm, idzt, idzm, int(lup), gd))
            wzt = max(wzt, float(np.max(np.abs(j - f))))
        jm = np.asarray(term_ma_zm_lhs(
            jgr.nzm, jgr.nzt, jgr.ngrdcol, jnp.asarray(wm_zm),
            jgr.invrs_dzm, jgr.weights_zm2zt,
        ))
        fm = np.asarray(clubb_f2py.f2py_term_ma_zm_lhs(wm_zm, idzm, w_zm2zt))
        wzm = max(wzm, float(np.max(np.abs(jm - fm))))
    assert wzt < 1e-12, f"term_ma_zt_lhs f2py mismatch {wzt:.2e}"
    assert wzm < 1e-12, f"term_ma_zm_lhs f2py mismatch {wzm:.2e}"
    print(f"  f2py term_ma_zt_lhs (centered+upwind) + term_ma_zm_lhs (20 cases): bit-match, worst {max(wzt, wzm):.2e}  PASS")


def main():
    print("test_mean_adv:")
    test_f2py_oracle()
    print("All mean_adv checks PASSED")


if __name__ == "__main__":
    main()
