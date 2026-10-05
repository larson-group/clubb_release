#!/usr/bin/env python3
"""validate the JAX diagnose_correlations_module port."""

# _ROOT first so `import clubb_jax` resolves to this checkout's package.

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.diagnose_correlations_module import calc_mean, calc_varnce, calc_w_corr, corr_array_assertion_checks


def _corr_matrix(n, seed):
    """A symmetric matrix with unit diagonal and off-diagonals in (-0.9, 0.9) — a valid 'prescribed' input."""
    rng = np.random.default_rng(seed)
    a = rng.uniform(-0.9, 0.9, size=(n, n))
    a = 0.5 * (a + a.T)
    np.fill_diagonal(a, 1.0)
    return np.asfortranarray(a)


def test_helpers():
    # calc_mean(a,x1,x2) = a*x1 + (1-a)*x2
    assert abs(float(calc_mean(0.3, 5.0, 2.0)) - (0.3 * 5 + 0.7 * 2)) < 1e-14
    # calc_varnce = a*((x1-xm)^2+x1p2) + (1-a)*((x2-xm)^2+x2p2)
    a, x1, x2, xm, x1p2, x2p2 = 0.4, 3.0, -1.0, 0.6, 0.5, 0.2
    ref = a * ((x1 - xm) ** 2 + x1p2) + (1 - a) * ((x2 - xm) ** 2 + x2p2)
    assert abs(float(calc_varnce(a, x1, x2, xm, x1p2, x2p2)) - ref) < 1e-13
    # calc_w_corr = clip(wpxp/(max(sx,xt)*max(sw,wt)), ±0.99); the clip must fire for a large cov
    assert abs(float(calc_w_corr(0.5, 1.0, 1.0, 1e-2, 1e-2)) - 0.5) < 1e-14
    assert abs(float(calc_w_corr(100.0, 1.0, 1.0, 1e-2, 1e-2)) - 0.99) < 1e-14   # clipped
    assert abs(float(calc_w_corr(-100.0, 1.0, 1.0, 1e-2, 1e-2)) + 0.99) < 1e-14  # clipped
    print("  helpers calc_mean/calc_varnce/calc_w_corr (incl. ±0.99 clip): closed-form  PASS")


def test_corr_array_assertion_checks():
    c = _corr_matrix(5, 9)                       # valid: symmetric, unit diag, off-diag in (-0.9,0.9)
    assert corr_array_assertion_checks(c) is True, "valid correlation matrix rejected"
    bad = c.copy(); bad[0, 1] = bad[1, 0] = 1.5  # off-diagonal out of [-0.99, 0.99]
    assert corr_array_assertion_checks(bad) is False, "out-of-range off-diagonal accepted"
    nd = c.copy(); nd[2, 2] = 0.9                # diagonal != 1
    assert corr_array_assertion_checks(nd) is False, "non-unit-diagonal accepted"
    print("  corr_array_assertion_checks: valid PASS / out-of-range + non-unit-diag FAIL  PASS")
