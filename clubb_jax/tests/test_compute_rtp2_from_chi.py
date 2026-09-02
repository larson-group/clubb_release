"""Validate setup_clubb_pdf_params.compute_rtp2_from_chi against a literal NumPy transcription
of the Fortran (setup_clubb_pdf_params.F90:compute_rtp2_from_chi).

The Fortran computes this only as the optional `rtp2_from_chi` stats diagnostic (so there is no f2py
trajectory oracle for it), but it is a pure algebraic combination — the per-component rt variance
implied by the chi/eta PDF, binormal-combined — so a transcription check is exact, and a jax.grad
check confirms it stays differentiable.
"""
import numpy as np
import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import compute_rtp2_from_chi


def _ref(sc1, sc2, se1, se2, rt1, rt2, c1, c2, mf, ce1, ce2):
    v1 = (ce1 * sc1 * se1 + 0.5 * sc1 ** 2 + 0.5 * se1 ** 2) / (2.0 * c1 ** 2)
    v2 = (ce2 * sc2 * se2 + 0.5 * sc2 ** 2 + 0.5 * se2 ** 2) / (2.0 * c2 ** 2)
    rtm = mf * rt1 + (1.0 - mf) * rt2
    s1, s2 = np.sqrt(v1), np.sqrt(v2)
    return mf * (s1 ** 2 + (rt1 - rtm) ** 2) + (1.0 - mf) * (s2 ** 2 + (rt2 - rtm) ** 2)


def test_compute_rtp2_from_chi_vs_reference():
    rng = np.random.default_rng(3)
    n = 50
    sc1 = rng.uniform(1e-4, 1e-3, n); sc2 = rng.uniform(1e-4, 1e-3, n)
    se1 = rng.uniform(1e-4, 1e-3, n); se2 = rng.uniform(1e-4, 1e-3, n)
    rt1 = rng.uniform(0.01, 0.02, n); rt2 = rng.uniform(0.01, 0.02, n)
    c1 = rng.uniform(0.5, 1.5, n); c2 = rng.uniform(0.5, 1.5, n)
    mf = rng.uniform(0.3, 0.7, n)
    ce1 = rng.uniform(-0.5, 0.5, n); ce2 = rng.uniform(-0.5, 0.5, n)
    got = np.asarray(compute_rtp2_from_chi(
        jnp.asarray(sc1), jnp.asarray(sc2), jnp.asarray(se1), jnp.asarray(se2),
        jnp.asarray(rt1), jnp.asarray(rt2), jnp.asarray(c1), jnp.asarray(c2),
        jnp.asarray(mf), jnp.asarray(ce1), jnp.asarray(ce2)))
    exp = _ref(sc1, sc2, se1, se2, rt1, rt2, c1, c2, mf, ce1, ce2)
    d = float(np.max(np.abs(got - exp)))
    assert d == 0.0, f"compute_rtp2_from_chi vs reference: {d:.2e}"
    assert np.all(got >= 0.0), "rtp2 must be non-negative"
    print(f"  compute_rtp2_from_chi vs NumPy transcription: bit-exact ({d:.1e}), all >= 0 ({n} pts)  PASS")


def test_compute_rtp2_from_chi_grad():
    def f(sc1):
        return jnp.sum(compute_rtp2_from_chi(
            jnp.atleast_1d(sc1), jnp.atleast_1d(5e-4), jnp.atleast_1d(5e-4), jnp.atleast_1d(5e-4),
            jnp.atleast_1d(0.015), jnp.atleast_1d(0.012), jnp.atleast_1d(1.0), jnp.atleast_1d(1.0),
            jnp.atleast_1d(0.5), jnp.atleast_1d(0.2), jnp.atleast_1d(0.2)))
    g = jax.grad(f)(5e-4)
    assert bool(jnp.isfinite(g)) and float(g) != 0.0, "grad must be finite + nonzero"
    print(f"  compute_rtp2_from_chi grad wrt sigma_chi_1: finite, nonzero ({float(g):.3e})  PASS")


if __name__ == "__main__":
    print("setup_clubb_pdf_params.compute_rtp2_from_chi vs Fortran transcription:")
    test_compute_rtp2_from_chi_vs_reference()
    test_compute_rtp2_from_chi_grad()
    print("Done.")
