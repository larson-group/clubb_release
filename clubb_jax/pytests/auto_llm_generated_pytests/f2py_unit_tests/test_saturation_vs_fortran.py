"""Check the JAX saturation ports (saturation.F90) against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.saturation as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_rel_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.saturation import (
    sat_mixrat_liq,
    sat_mixrat_ice,
    SATURATION_FLATAU,
    SATURATION_BOLTON,
)


def test_sat_mixrat_f2py():
    """Bit-shadow the saturation mixing ratios vs the f2py Fortran oracle. The single elementwise JAX
    sat_mixrat_liq / sat_mixrat_ice cover both the Fortran rank-1 (_k) and rank-2 (_2D) generic-interface
    procedures; validate the 2D form against f2py_sat_mixrat_liq_2d (both saturation_formula leaves) and
    f2py_sat_mixrat_ice_2d. The leaf SVP polynomials match coefficient-for-coefficient, so this is bit-faithful
    (the ice path's slightly larger residual is the polynomial evaluation order)."""
    import clubb_f2py
    # sat_mixrat_liq's wrapper takes a grid; a 12-level column matches the test arrays.
    _, fgr = setup_matching_grids(2, 40.0, 480.0)
    rng = np.random.default_rng(0)
    worst_liq = worst_ice = 0.0
    for formula in (SATURATION_FLATAU, SATURATION_BOLTON):
        for _ in range(15):
            ng, nz = 2, 12
            p = rng.uniform(2e4, 1.02e5, (ng, nz)); T = rng.uniform(220.0, 310.0, (ng, nz))
            j = np.asarray(sat_mixrat_liq(jnp.asarray(p), jnp.asarray(T), formula))
            f = np.asarray(fortran.sat_mixrat_liq(nz, ng, fgr, p, T, formula))
            worst_liq = max(worst_liq, max_rel_diff("saturation", j, f, 1e-30))
    for _ in range(15):
        ng, nz = 2, 12
        p = rng.uniform(2e4, 1.02e5, (ng, nz)); T = rng.uniform(220.0, 273.0, (ng, nz))
        j = np.asarray(sat_mixrat_ice(jnp.asarray(p), jnp.asarray(T)))
        f = np.asarray(clubb_f2py.f2py_sat_mixrat_ice_2d(p, T, SATURATION_FLATAU))
        worst_ice = max(worst_ice, max_rel_diff("saturation", j, f, 1e-30))
    assert worst_liq < 1e-12, f"sat_mixrat_liq f2py rel mismatch {worst_liq:.2e}"
    assert worst_ice < 1e-11, f"sat_mixrat_ice f2py rel mismatch {worst_ice:.2e}"
