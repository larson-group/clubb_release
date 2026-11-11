"""Check the JAX PPM (method 2) remapping port (remapping_module.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.remapping_module as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.remapping_module import remap_vals_to_target_same_grid

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def _rho_spline(ng, ztop):
    levels = np.linspace(0.0, ztop + 100.0, 24)
    vals = 1.2 * np.exp(-levels / 8000.0)
    return np.tile(vals, (ng, 1)), np.tile(levels, (ng, 1))


def test_f2py_same_grid_identity():
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rho_vals, rho_levels = _rho_spline(ng, float(jgr.zm[0, -1]))
    p_sfc = np.full(ng, 1.0e5)
    rng = np.random.default_rng(3)
    worst = 0.0
    for iv in (1, 0, -1):
        src = rng.uniform(0.2, 2.0, (ng, nzt))
        f = np.asarray(fortran.remap_vals_to_target(
            ng, fgr, fgr, nzt, src, nzt, rho_vals.shape[1], rho_vals, rho_levels, iv, p_sfc, 2, True))   # grid_remap_method=2 (PPM)
        g = np.asarray(remap_vals_to_target_same_grid(
            src, np.asarray(jgr.zm), np.asarray(jgr.zt), rho_vals, rho_levels, p_sfc,
            grid_remap_method=2, l_zt_variable=True, iv=iv))
        worst = max(worst, max_abs_diff("remapping_ppm", g, f))
        # Remapping onto the identical grid must reproduce the input (independent of either implementation).
        assert max_abs_diff("Fortran same-grid identity", f, src) < 1e-12, f"Fortran PPM identity fails (iv={iv})"
        assert max_abs_diff("JAX same-grid identity", g, src) < 1e-12, f"JAX PPM identity fails (iv={iv})"
    assert worst < 1e-11, f"PPM same-grid f2py mismatch {worst:.2e}"
