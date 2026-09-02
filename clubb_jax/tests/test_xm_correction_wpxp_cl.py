"""Unit test for `xm_correction_wpxp_cl` (advance_xm_wpxp_module.py) — the JAX mirror
of the Fortran `xm_correction_wpxp_cl` (advance_xm_wpxp_module.F90:5766), which corrects xm
when w'x' was clipped.

The f2py oracle exposes only the whole `advance_xm_wpxp` (not this sub-subroutine), and the
routine is a no-op for every bit-faithful integration case (the covariance clip does not fire),
so it cannot be positively validated through the whole-driver gate. This test validates the math
+ the per-column `l_clipping_needed` (|wpxp_chnge| > eps) gate against an independent NumPy
transcription of the Fortran loop.
"""
import numpy as np
from jax import config
config.update("jax_enable_x64", True)   # match the production driver (float64)
import jax.numpy as jnp

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.advance_xm_wpxp_module import xm_correction_wpxp_cl
from clubb_jax.src.CLUBB_core.constants_clubb import eps


def _ref(xm, wpxp_chnge, invrs_dzt, dt):
    """Independent NumPy transcription of the Fortran xm_correction_wpxp_cl loop."""
    ng, nzt = xm.shape
    out = xm.copy()
    for i in range(ng):
        if np.any(np.abs(wpxp_chnge[i]) > eps):          # l_clipping_needed(i)
            for k in range(nzt):
                tndcy = -invrs_dzt[i, k] * (wpxp_chnge[i, k + 1] - wpxp_chnge[i, k])
                out[i, k] = xm[i, k] + tndcy * dt
    return out


def test_xm_correction_wpxp_cl_matches_reference():
    rng = np.random.default_rng(0)
    ng, nzt = 3, 12
    nzm = nzt + 1
    xm = rng.normal(size=(ng, nzt))
    invrs_dzt = rng.uniform(0.5, 2.0, size=(ng, nzt))
    dt = 60.0

    wpxp_chnge = rng.normal(size=(ng, nzm)) * 1e-3      # col 0: clipping fired
    wpxp_chnge[1] = 0.0                                  # col 1: no clipping → no-op
    wpxp_chnge[2] = 0.5 * eps                            # col 2: sub-eps → gated off (no-op)

    got = np.asarray(xm_correction_wpxp_cl(
        jnp.asarray(xm), jnp.asarray(wpxp_chnge), jnp.asarray(invrs_dzt), dt))
    ref = _ref(xm, wpxp_chnge, invrs_dzt, dt)

    assert np.allclose(got, ref, atol=1e-14), f"max |Δ| = {np.abs(got - ref).max():.3e}"
    # columns with no (or sub-eps) clipping must be left exactly untouched
    assert np.array_equal(got[1], xm[1]), "no-clip column was modified"
    assert np.array_equal(got[2], xm[2]), "sub-eps column was not gated off"
    # column 0 (clipping fired) must actually change
    assert not np.allclose(got[0], xm[0]), "clip column unchanged — correction not applied"


def test_xm_correction_wpxp_cl_all_zero_is_noop():
    rng = np.random.default_rng(1)
    ng, nzt = 2, 8
    xm = rng.normal(size=(ng, nzt))
    invrs_dzt = rng.uniform(0.5, 2.0, size=(ng, nzt))
    got = np.asarray(xm_correction_wpxp_cl(
        jnp.asarray(xm), jnp.zeros((ng, nzt + 1)), jnp.asarray(invrs_dzt), 60.0))
    assert np.array_equal(got, xm), "zero wpxp_chnge should be an exact no-op"


if __name__ == "__main__":
    test_xm_correction_wpxp_cl_matches_reference()
    test_xm_correction_wpxp_cl_all_zero_is_noop()
    print("PASS: xm_correction_wpxp_cl matches the NumPy reference + gates correctly")
