#!/usr/bin/env python3
"""validate the JAX calc_comp_corrs_binormal + smooth_corr_quotient ports."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_utilities import calc_comp_corrs_binormal, smooth_corr_quotient

NG, NZ = 2, 6


def test_round_trip():
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(100):
        a = rng.uniform(0.3, 0.7)
        mu_x_1, mu_x_2, mu_y_1, mu_y_2 = rng.uniform(-2, 2, 4)
        sx1, sx2, sy1, sy2 = rng.uniform(0.3, 1.5, 4)   # std devs
        xm = a * mu_x_1 + (1 - a) * mu_x_2
        ym = a * mu_y_1 + (1 - a) * mu_y_2
        corr = rng.uniform(-0.9, 0.9)
        # Forward covariance assembly.
        xpyp = (a * ((mu_x_1 - xm) * (mu_y_1 - ym) + corr * sx1 * sy1)
                + (1 - a) * ((mu_x_2 - xm) * (mu_y_2 - ym) + corr * sx2 * sy2))
        g1, g2 = calc_comp_corrs_binormal(
            np.array([[xpyp]]), np.array([[xm]]), np.array([[ym]]),
            np.array([[mu_x_1]]), np.array([[mu_x_2]]), np.array([[mu_y_1]]), np.array([[mu_y_2]]),
            np.array([[sx1 ** 2]]), np.array([[sx2 ** 2]]), np.array([[sy1 ** 2]]), np.array([[sy2 ** 2]]),
            np.array([[a]]))
        v1 = float(np.asarray(g1).item()); v2 = float(np.asarray(g2).item())
        worst = max(worst, abs(v1 - corr), abs(v1 - v2))
    assert worst < 1e-9, f"round-trip recovery {worst:.2e}"
    print(f"  round-trip: recover corr from assembled <x'y'>, worst {worst:.2e}  PASS")


def test_bound():
    # Huge covariance with tiny variances -> quotient would blow up; smoothing must cap |corr| <= 0.99.
    q = float(np.asarray(smooth_corr_quotient(np.array([[1.0e6]]), np.array([[1.0e-3]]), 1.0e-10)).item())
    assert abs(q) <= 0.99 + 1e-9, f"smooth_corr_quotient exceeded max_mag_correlation: {q}"
    qn = float(np.asarray(smooth_corr_quotient(np.array([[-1.0e6]]), np.array([[1.0e-3]]), 1.0e-10)).item())
    assert abs(qn) <= 0.99 + 1e-9
    print(f"  smooth_corr_quotient bound: |corr| <= 0.99 for a huge covariance ({q:.4f})  PASS")
