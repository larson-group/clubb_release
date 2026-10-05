#!/usr/bin/env python3
"""validate the JAX simplified-shortwave-radiation port."""


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
