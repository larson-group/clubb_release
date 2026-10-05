#!/usr/bin/env python3
"""validate the JAX smooth_min + calc_xpwp ports (advance_helper_module)."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os

_ROOT = str(_REPO_ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import smooth_min as _smooth_min, calc_xpwp as _calc_xpwp
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_TUNABLE_PARAMS = os.path.join(_ROOT, "input", "parameter_and_flag_configs", "default", "tunable_parameters.in")

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _smooth_dims(a, b):
    shape_a = getattr(a, "shape", ())
    shape_b = getattr(b, "shape", ())
    shape = shape_a if len(shape_a) > 0 else shape_b
    if len(shape) >= 2:
        return shape[1], shape[0]
    return shape[0] if len(shape) == 1 else 1, 1


def smooth_min(a, b, coef):
    nz, ngrdcol = _smooth_dims(a, b)
    return _smooth_min(nz, ngrdcol, a, b, coef)


def calc_xpwp(Km_zm, xm, invrs_dzm):
    Km_zm = jnp.asarray(Km_zm, dtype=jnp.float64)
    xm = jnp.asarray(xm, dtype=jnp.float64)
    invrs_dzm = jnp.asarray(invrs_dzm, dtype=jnp.float64)
    nzm = Km_zm.shape[-1]
    ngrdcol = 1 if Km_zm.ndim == 1 else Km_zm.shape[0]
    # Supply the existing kernel's grid field; the test adapter does no flux arithmetic.
    gr = setup_grid(ngrdcol, 1.0, 0.0, float(nzm - 1))
    gr = gr._replace(invrs_dzm=jnp.broadcast_to(invrs_dzm, (ngrdcol, nzm)))
    return _calc_xpwp(nzm, xm.shape[-1], ngrdcol, gr, Km_zm, xm)


def test_smooth_min_closed_form():
    a = np.array([[1.0, 5.0, -2.0]])
    coef = 1e-3
    out = np.asarray(smooth_min(a, 3.0, coef))
    ref = 0.5 * ((a + 3.0) - np.sqrt((a - 3.0) ** 2 + coef ** 2))
    assert np.max(np.abs(out - ref)) < 1e-14
    assert np.all(out <= np.minimum(a, 3.0) + 1e-9), "smooth_min must be <= min"
    print("  smooth_min: closed-form + (smooth_min <= min) bound  PASS")


def test_calc_xpwp_identity():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]; nzt = nzm - 1
    rng = np.random.default_rng(3)
    Km = np.abs(rng.standard_normal((1, nzm))) + 0.1
    xm = rng.standard_normal((1, nzt))
    invrs_dzm = np.asarray(jgr.invrs_dzm)
    out = np.asarray(calc_xpwp(Km, xm, invrs_dzm))
    for k in range(1, nzm - 1):
        expect = Km[0, k] * invrs_dzm[0, k] * (xm[0, k] - xm[0, k - 1])
        assert abs(out[0, k] - expect) < 1e-13, f"xpwp identity failed at k={k}"
    assert out[0, 0] == 0.0 and out[0, nzm - 1] == 0.0, "boundary levels must be zero"
    print("  calc_xpwp: down-gradient identity on interior + zero boundaries  PASS")
