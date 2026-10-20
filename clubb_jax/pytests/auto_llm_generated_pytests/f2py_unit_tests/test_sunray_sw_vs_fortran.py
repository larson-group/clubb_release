"""Check the JAX simplified-shortwave-radiation port against Fortran."""

import numpy as np

import clubb_python.radiation as fortran_radiation

from clubb_jax.pytests.f2py_unit_test_helpers import max_rel_diff
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.Radiation.rad_lwsw_module import sunray_sw

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_RADIUS, _ALVDR, _GC, _OMEGA = 1.0e-5, 0.1, 0.85, 0.9965   # eff_drop_radius / alvdr / gc / omega


def test_f2py_oracle():
    """sunray_sw (simplified two-stream shortwave flux) matches Fortran for both l_center settings."""
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    dzt = np.asarray(jgr.dzt); zm = np.asarray(jgr.zm); zt = np.asarray(jgr.zt)
    rng = np.random.default_rng(0)
    worst = 0.0
    for it in range(16):
        rcm = rng.uniform(0.0, 1e-3, (ng, nzt)); rho = rng.uniform(0.5, 1.2, (ng, nzt))
        fs0 = float(rng.uniform(100.0, 400.0)); amu0 = float(rng.uniform(0.3, 0.9))
        args = (ng, nzt, rcm, rho, amu0, dzt, zm, zt, _RADIUS, _ALVDR, _GC, fs0, _OMEGA, bool(it % 2))
        j, f = np.asarray(sunray_sw(*args)), np.asarray(fortran_radiation.sunray_sw(*args))
        worst = max(worst, max_rel_diff("sunray_sw", j, f, 1e-30))
    assert worst < 1e-12, f"sunray_sw Fortran rel mismatch {worst:.2e}"
