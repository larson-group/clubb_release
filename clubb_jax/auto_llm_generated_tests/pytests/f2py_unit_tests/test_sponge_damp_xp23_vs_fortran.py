"""Check the JAX sponge_layer_damping.py ports (sponge_damp_xm/xp2/xp3)
against Fortran."""

from clubb_jax.auto_llm_generated_tests.pytests.array_inputs import call_jax

from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np
import pytest

import clubb_python.CLUBB_core.sponge_layer_damping as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.sponge_layer_damping import (
    sponge_damp_xm,
    sponge_damp_xp2,
    sponge_damp_xp3,
    initialize_tau_sponge_damp_api, sponge_damp_settings, sponge_damp_profile,
)

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_TAU, _DEPTH, _DT, _XTOLSQ = 100.0, 300.0, 60.0, 1.0e-12


def test_f2py_oracle():
    """sponge_damp_xp2 / sponge_damp_xp3 damp only the top sponge layer, exactly as Fortran does."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    zm = np.asarray(jgr.zm); zt = np.asarray(jgr.zt)
    rng = np.random.default_rng(4)
    xp2 = rng.uniform(0.01, 2.0, (ng, nzm)); xp3 = rng.uniform(-1, 1, (ng, nzt))
    depth = np.full(ng, _DEPTH)
    f2 = fortran.sponge_damp_xp2(ng, nzm, fgr, _DT, zm, xp2, _XTOLSQ,
                                 SimpleNamespace(tau_sponge_damp=np.full((ng, nzm), _TAU), sponge_layer_depth=depth))
    g2 = call_jax(sponge_damp_xp2, ng, nzm, jgr, _DT, zm, xp2, _XTOLSQ, sponge_damp_profile(jnp.full((ng, nzm), _TAU), jnp.full(ng, _DEPTH)))
    f3 = fortran.sponge_damp_xp3(ng, nzm, nzt, fgr, _DT, zt, zm, xp3,
                                 SimpleNamespace(tau_sponge_damp=np.full((ng, nzt), _TAU), sponge_layer_depth=depth))
    g3 = call_jax(sponge_damp_xp3, ng, nzm, nzt, jgr, _DT, zt, zm, xp3, sponge_damp_profile(jnp.full((ng, nzt), _TAU), jnp.full(ng, _DEPTH)))
    w = max(max_abs_diff("sponge_damp_xp2", g2, f2), max_abs_diff("sponge_damp_xp3", g3, f3))
    assert w < 1e-12, f"sponge_damp xp2/xp3 Fortran mismatch {w:.2e}"
    assert np.any(np.asarray(g2) != xp2) and np.any(np.asarray(g3) != xp3), "sponge layer not exercised"


@pytest.mark.parametrize("different_tops", [False, True])
def test_sponge_damp_xm_f2py(different_tops):
    """sponge_damp_xm nudges xm toward xm_ref in the top sponge exactly as Fortran does.

    Both implementations use column-specific time scales and layer depths."""
    if different_tops:
        heights = np.array([[0., 100., 600., 800., 1000.], [0., 200., 800., 1500., 2000.]])
        jgr, fgr = setup_matching_grids(2, 0., [1000., 2000.], grid_type=3, momentum_heights=heights)
    else:
        jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    zt = np.asarray(jgr.zt); zm = np.asarray(jgr.zm); zm_top = float(jgr.zm[0, -1])
    tau_min, tau_max, depth_frac = 120.0, 1800.0, 0.25
    settings = sponge_damp_settings(tau_min, tau_max, depth_frac)
    jprofile = initialize_tau_sponge_damp_api(ng, nzt, jgr, _DT, jgr.zt, settings)
    profile = SimpleNamespace(tau_sponge_damp=np.asarray(jprofile.tau_sponge_damp), sponge_layer_depth=jprofile.sponge_layer_depth)
    rng = np.random.default_rng(2)
    worst = 0.0; damped = False
    for _ in range(10):
        xm = rng.uniform(280.0, 320.0, (ng, nzt)); xm_ref = rng.uniform(280.0, 320.0, (ng, nzt))
        j = np.asarray(call_jax(sponge_damp_xm, ng, nzm, nzt, jgr, _DT, zt, zm, xm_ref, xm, jprofile))
        f = np.asarray(fortran.sponge_damp_xm(ng, nzm, nzt, fgr, _DT, zt, zm, xm_ref, xm, profile))
        worst = max(worst, max_rel_diff("sponge_damp_xm", j, f, 1.0))
        damped = damped or bool(np.any(j != xm))
    assert worst < 1e-12, f"sponge_damp_xm Fortran rel mismatch {worst:.2e}"
    assert damped, "sponge layer not exercised for xm"
