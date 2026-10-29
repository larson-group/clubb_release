#!/usr/bin/env python3
"""validate the JAX sponge_layer_damping.py ports (sponge_damp_xm/xp2/xp3)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.sponge_layer_damping import sponge_damp_xp2, sponge_damp_xp3

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_TAU, _DEPTH, _DT, _XTOLSQ = 100.0, 300.0, 60.0, 1.0e-12


def test_below_sponge_noop():
    nzm = 12
    zm = np.tile(np.linspace(0, 1200, nzm)[None, :], (1, 1))
    zt = 0.5 * (zm[:, :-1] + zm[:, 1:])
    xp2 = np.ones((1, nzm)) * 0.5; xp3 = np.ones((1, nzm - 1)) * 0.3
    tau2 = np.full((1, nzm), _TAU); tau3 = np.full((1, nzm - 1), _TAU)
    g2 = np.asarray(sponge_damp_xp2(_DT, zm, xp2, _XTOLSQ, tau2, _DEPTH))
    g3 = np.asarray(sponge_damp_xp3(_DT, zt, zm, xp3, tau3, _DEPTH))
    far = (zm[0, -1] - zm[0]) >= _DEPTH
    assert np.allclose(g2[0, far], xp2[0, far]), "xp2 modified below sponge"
    fart = (zm[0, -1] - zt[0]) >= _DEPTH
    assert np.allclose(g3[0, fart], xp3[0, fart]), "xp3 modified below sponge"
    print("  below-sponge no-op for both  PASS")
