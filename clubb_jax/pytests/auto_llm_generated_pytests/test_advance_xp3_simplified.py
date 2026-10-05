#!/usr/bin/env python3
"""validate advance_xp3_module.py:advance_xp3_simplified + term_tp_rhs/term_ac_rhs."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_xp3_module import term_tp_rhs, term_ac_rhs

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_term_rhs_formulas():
    """term_tp_rhs / term_ac_rhs vs their literal F90 formulas (advance_xp3_module.F90:934-935 / 1004) + grad."""
    rng = np.random.default_rng(11)
    shp = (_NG, 9)
    xp2_zt = rng.uniform(1e-4, 1.0, shp); wpxpp1 = rng.standard_normal(shp); wpxp = rng.standard_normal(shp)
    rho_p1 = rng.uniform(0.5, 1.2, shp); rho = rng.uniform(0.5, 1.2, shp)
    irho = rng.uniform(0.8, 2.0, shp); idzt = rng.uniform(0.01, 0.05, shp)
    xm_p1 = rng.standard_normal(shp); xm = rng.standard_normal(shp); wpxp2 = rng.standard_normal(shp)

    tp = np.asarray(term_tp_rhs(*(jnp.asarray(a) for a in
                    (xp2_zt, wpxpp1, wpxp, rho_p1, rho, irho, idzt))))
    tp_ref = 3.0 * xp2_zt * irho * idzt * (rho_p1 * wpxpp1 - rho * wpxp)
    ac = np.asarray(term_ac_rhs(*(jnp.asarray(a) for a in (xm_p1, xm, wpxp2, idzt))))
    ac_ref = -3.0 * wpxp2 * idzt * (xm_p1 - xm)
    assert float(np.max(np.abs(tp - tp_ref))) < 1e-15, "term_tp_rhs != literal F90 formula"
    assert float(np.max(np.abs(ac - ac_ref))) < 1e-15, "term_ac_rhs != literal F90 formula"
    g = jax.grad(lambda a: jnp.sum(term_tp_rhs(a, jnp.asarray(wpxpp1), jnp.asarray(wpxp),
                jnp.asarray(rho_p1), jnp.asarray(rho), jnp.asarray(irho), jnp.asarray(idzt)) ** 2))(jnp.asarray(xp2_zt))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad through term_tp_rhs"
    print("  term_tp_rhs / term_ac_rhs == literal F90 formula (exact) + grad finite  PASS")
