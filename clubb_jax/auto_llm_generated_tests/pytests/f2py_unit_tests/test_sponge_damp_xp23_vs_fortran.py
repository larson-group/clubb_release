"""Check the JAX sponge_layer_damping.py ports (sponge_damp_xm/xp2/xp3)
against Fortran."""

from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.sponge_layer_damping as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.sponge_layer_damping import (
    sponge_damp_xm,
    sponge_damp_xp2,
    sponge_damp_xp3,
    initialize_tau_sponge_damp,
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
    f2 = fortran.sponge_damp_xp2(fgr, nzm, _DT, zm, xp2, np.full(ng, _XTOLSQ),
                                 SimpleNamespace(tau_sponge_damp=np.full((ng, nzm), _TAU), sponge_layer_depth=depth))
    g2 = sponge_damp_xp2(_DT, zm, xp2, np.full((ng, 1), _XTOLSQ), np.full((ng, nzm), _TAU), depth[:, None])
    f3 = fortran.sponge_damp_xp3(fgr, nzm, nzt, _DT, zt, zm, xp3,
                                 SimpleNamespace(tau_sponge_damp=np.full((ng, nzt), _TAU), sponge_layer_depth=depth))
    g3 = sponge_damp_xp3(_DT, zt, zm, xp3, np.full((ng, nzt), _TAU), depth[:, None])
    w = max(max_abs_diff("sponge_damp_xp2", g2, f2), max_abs_diff("sponge_damp_xp3", g3, f3))
    assert w < 1e-12, f"sponge_damp xp2/xp3 Fortran mismatch {w:.2e}"
    assert np.any(np.asarray(g2) != xp2) and np.any(np.asarray(g3) != xp3), "sponge layer not exercised"


def test_sponge_damp_xm_f2py():
    """sponge_damp_xm nudges xm toward xm_ref in the top sponge exactly as Fortran does.

    The JAX relies on tau = inf below the sponge (from initialize_tau_sponge_damp) instead of a
    geometric gate, so both sides receive that tau profile."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    zt = np.asarray(jgr.zt); zm = np.asarray(jgr.zm); zm_top = float(jgr.zm[0, -1])
    tau_min, tau_max, depth_frac = 120.0, 1800.0, 0.25
    tau = np.zeros((ng, nzt))
    for i in range(ng):
        tau[i], _ = initialize_tau_sponge_damp(zt[i], _DT, zm_top, tau_min, tau_max, depth_frac)
    sld = depth_frac * zm_top
    profile = SimpleNamespace(tau_sponge_damp=tau, sponge_layer_depth=np.full(ng, sld))
    rng = np.random.default_rng(2)
    worst = 0.0; damped = False
    for _ in range(10):
        xm = rng.uniform(280.0, 320.0, (ng, nzt)); xm_ref = rng.uniform(280.0, 320.0, (ng, nzt))
        j = np.asarray(sponge_damp_xm(jnp.asarray(xm), jnp.asarray(xm_ref), zt, zm_top, jnp.asarray(tau), sld, _DT))
        f = np.asarray(fortran.sponge_damp_xm(fgr, nzm, nzt, _DT, zt, zm, xm_ref, xm, profile))
        worst = max(worst, max_rel_diff("sponge_damp_xm", j, f, 1.0))
        damped = damped or bool(np.any(j != xm))
    assert worst < 1e-12, f"sponge_damp_xm Fortran rel mismatch {worst:.2e}"
    assert damped, "sponge layer not exercised for xm"
