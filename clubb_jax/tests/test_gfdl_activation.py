#!/usr/bin/env python3
"""test_gfdl_activation.py — validate the CLUBB-side pieces of GFDL droplet activation.

erff is validated against math.erf (the true independent oracle); updraft_weights against a literal NumPy
transcription of the Fortran (including its sequential-assignment normalization quirk). Plus a finite jax.grad.
"""
import os
import sys
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.Microphys.gfdl_activation import (
    erff, updraft_weights, aer_act_clubb_ndrop, _WP2_EPS, _P_UPDRAFT_EPS)


def test_erff_vs_math():
    xs = np.linspace(-4.0, 4.0, 2001)
    got = np.asarray(erff(jnp.asarray(xs)))
    ref = np.array([math.erf(x) for x in xs])
    err = np.max(np.abs(got - ref))
    assert err < 1e-6, f"erff vs math.erf max abs err {err:.2e}"
    print(f"  erff vs math.erf (2001 points, [-4,4]): max abs err {err:.1e}  PASS")


def _ref_updraft(w1, w2, vw1, vw2, mixt, cf1, cf2):
    """Literal transcription of gfdl_activation.F90:109-143 (one level), quirk included. Uses the ported `erff`
    (not math.erf) so this isolates the TRANSCRIPTION logic; erff's own accuracy is checked separately."""
    _erf = lambda x: float(erff(jnp.asarray(x)))
    if vw1 > _WP2_EPS:
        P1 = (0.5 + 0.5 * _erf(w1 / math.sqrt(2.0 * vw1))) * mixt * cf1
    elif w1 > 0.0:
        P1 = mixt * cf1
    else:
        P1 = 0.0
    if vw2 > _WP2_EPS:
        P2 = (0.5 + 0.5 * _erf(w2 / math.sqrt(2.0 * vw2))) * (1.0 - mixt) * cf2
    elif w2 > 0.0:
        P2 = (1.0 - mixt) * cf2
    else:
        P2 = 0.0
    if P1 + P2 > _P_UPDRAFT_EPS:
        P1 = P1 / (P1 + P2)          # reassigned first
        P2 = P2 / (P1 + P2)          # denominator uses the NEW P1 (the Fortran quirk)
    else:
        P1 = 0.0; P2 = 0.0
    return P1, P2


def test_updraft_weights_vs_literal():
    rng = np.random.default_rng(5)
    n = 500
    w1 = rng.standard_normal(n); w2 = rng.standard_normal(n)
    vw1 = np.abs(rng.standard_normal(n)) * 0.5      # mix of >wp2_eps and (rarely) small
    vw2 = np.abs(rng.standard_normal(n)) * 0.5
    vw1[::5] = 1e-6; vw2[::7] = 1e-6                # force the degenerate (varnce<=eps) branch
    mixt = rng.uniform(0.1, 0.9, n)
    cf1 = rng.uniform(0.0, 1.0, n); cf2 = rng.uniform(0.0, 1.0, n)
    P1, P2 = updraft_weights(*(jnp.asarray(a) for a in (w1, w2, vw1, vw2, mixt, cf1, cf2)))
    P1 = np.asarray(P1); P2 = np.asarray(P2)
    worst = 0.0
    for i in range(n):
        r1, r2 = _ref_updraft(w1[i], w2[i], vw1[i], vw2[i], mixt[i], cf1[i], cf2[i])
        worst = max(worst, abs(P1[i] - r1), abs(P2[i] - r2))
    assert worst < 1e-13, f"updraft_weights vs literal {worst:.2e}"
    print(f"  updraft_weights vs literal NumPy (500 levels, incl. degenerate branch + quirk): max diff {worst:.1e}  PASS")


def test_differentiable():
    gE = jax.grad(lambda x: jnp.sum(erff(x)))(jnp.linspace(-2.0, 2.0, 11))
    # d/dx erf(x) = 2/sqrt(pi) exp(-x^2); spot-check at 0 -> 2/sqrt(pi)
    g0 = float(jax.grad(lambda x: erff(x))(jnp.asarray(0.0)))
    assert abs(g0 - 2.0 / math.sqrt(math.pi)) < 1e-3, f"d erf/dx at 0 = {g0}"
    gw = jax.grad(lambda v: jnp.sum(updraft_weights(
        jnp.asarray([0.3, -0.2]), jnp.asarray([0.1, 0.4]), v, jnp.asarray([0.2, 0.3]),
        jnp.asarray([0.4, 0.6]), jnp.asarray([0.7, 0.5]), jnp.asarray([0.6, 0.8]))[0]))(jnp.asarray([0.3, 0.5]))
    assert np.isfinite(np.asarray(gE)).all() and np.isfinite(np.asarray(gw)).all(), "non-finite grad"
    print(f"  jax.grad: d(erf)/dx|0 = {g0:.5f} (=2/sqrt(pi)); updraft_weights grad finite  PASS")


def test_ndrop_assembly():
    rng = np.random.default_rng(8)
    n = 100
    d1 = rng.uniform(0, 3e8, n); d2 = rng.uniform(0, 3e8, n)
    P1 = rng.uniform(0, 1, n); P2 = rng.uniform(0, 1, n)
    mixt = rng.uniform(0.2, 0.8, n); cf1 = rng.uniform(0, 1, n); cf2 = rng.uniform(0, 1, n)
    got = np.asarray(aer_act_clubb_ndrop(*(jnp.asarray(a) for a in (d1, d2, P1, P2, mixt, cf1, cf2))))
    ref = (d1 * P1 + d2 * P2) * (mixt * cf1 + (1 - mixt) * cf2)
    assert np.max(np.abs(got - ref)) < 1e-6, "ndrop vs literal"
    # Ndrop = 0 when no cloud (cf1=cf2=0); non-negative for non-negative inputs.
    assert float(aer_act_clubb_ndrop(1e8, 1e8, 0.5, 0.5, 0.5, 0.0, 0.0)) == 0.0, "no cloud -> 0"
    assert np.all(got >= 0.0), "Ndrop should be non-negative"
    g = jax.grad(lambda d: aer_act_clubb_ndrop(d, 1e8, 0.6, 0.4, 0.5, 0.8, 0.6))(jnp.asarray(2e8))
    assert np.isfinite(float(g)), "non-finite grad"
    print("  aer_act_clubb_ndrop vs literal + no-cloud→0 + non-negative + grad  PASS")


def main():
    print("test_gfdl_activation:")
    for t in (test_erff_vs_math, test_updraft_weights_vs_literal, test_ndrop_assembly, test_differentiable):
        t()
    print("All gfdl_activation checks PASSED")


if __name__ == "__main__":
    main()
