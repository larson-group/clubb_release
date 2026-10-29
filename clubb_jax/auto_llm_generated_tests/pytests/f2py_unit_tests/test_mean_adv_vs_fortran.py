"""Check the mean_adv.F90 ports (term_ma_zt_lhs / term_ma_zm_lhs)
against Fortran."""

import jax.numpy as jnp
import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs, term_ma_zm_lhs


def test_f2py_oracle():
    import clubb_f2py
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
            wzt = max(wzt, float(max_abs_diff("mean_adv", j, f)))
        jm = np.asarray(term_ma_zm_lhs(
            jgr.nzm, jgr.nzt, jgr.ngrdcol, jnp.asarray(wm_zm),
            jgr.invrs_dzm, jgr.weights_zm2zt,
        ))
        fm = np.asarray(clubb_f2py.f2py_term_ma_zm_lhs(wm_zm, idzm, w_zm2zt))
        wzm = max(wzm, float(max_abs_diff("mean_adv", jm, fm)))
    assert wzt < 1e-12, f"term_ma_zt_lhs f2py mismatch {wzt:.2e}"
    assert wzm < 1e-12, f"term_ma_zm_lhs f2py mismatch {wzm:.2e}"
