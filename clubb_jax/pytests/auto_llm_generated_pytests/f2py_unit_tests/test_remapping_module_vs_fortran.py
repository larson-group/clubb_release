"""Check the JAX mass-conserving remapping port (remapping_module.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.remapping_module as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.remapping_module import remap_vals_to_target_same_grid

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def _rho_spline(ng, ztop):
    # Density spline covering [0, ztop+pad], decreasing with height; ascending levels.
    levels = np.linspace(0.0, ztop + 100.0, 24)
    vals = 1.2 * np.exp(-levels / 8000.0)
    return np.tile(vals, (ng, 1)), np.tile(levels, (ng, 1))


def test_f2py_oracle():
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rho_vals, rho_levels = _rho_spline(ng, float(jgr.zm[0, -1]))
    p_sfc = np.full(ng, 1.0e5)
    rng = np.random.default_rng(7)

    worst = 0.0
    for l_zt in (True, False):
        ncell = nzt if l_zt else nzm
        src = rng.uniform(-2, 2, (ng, ncell))
        f = np.asarray(fortran.remap_vals_to_target(
            ng, fgr, fgr, ncell, src, ncell, rho_vals.shape[1], rho_vals, rho_levels, 1, p_sfc, 1, l_zt))
        g = np.asarray(remap_vals_to_target_same_grid(
            src, np.asarray(jgr.zm), np.asarray(jgr.zt), rho_vals, rho_levels, p_sfc,
            grid_remap_method=1, l_zt_variable=l_zt))
        worst = max(worst, max_abs_diff("remapping_module", g, f))
    assert worst < 1e-12, f"remap_vals_to_target_same_grid f2py mismatch {worst:.2e}"
