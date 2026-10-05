#!/usr/bin/env python3
"""validate the JAX compute_mean_binormal / compute_variance_binormal ports."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, compute_variance_binormal


def test_monte_carlo():
    rng = np.random.default_rng(11)
    mu1, mu2, s1, s2, a = 1.0, -2.0, 0.5, 1.3, 0.35
    n = 4_000_000
    pick1 = rng.random(n) < a
    samples = np.where(pick1, rng.normal(mu1, s1, n), rng.normal(mu2, s2, n))
    xm = float(compute_mean_binormal(mu1, mu2, a))
    xp2 = float(compute_variance_binormal(xm, mu1, mu2, s1, s2, a))
    assert abs(xm - samples.mean()) < 5e-3, f"mean {xm} vs MC {samples.mean()}"
    assert abs(xp2 - samples.var()) < 1e-2, f"variance {xp2} vs MC {samples.var()}"
    print(f"  Monte-Carlo: formula mean/variance match samples (xm={xm:.3f}, xp2={xp2:.3f})  PASS")
