"""Verification of pdf_utilities.py — lognormal<->normal moment/correlation conversions."""

import numpy as np
import jax

jax.config.update("jax_enable_x64", True)


from clubb_jax.src.CLUBB_core.pdf_utilities import stdev_L2N, corr_NL2NN, corr_LL2NN, MAX_MAG_CORRELATION


def test_corr_NL2NN_vs_montecarlo():
    """corr_NL2NN(corr(x,y)) reproduces the measured corr(x, ln y) for x normal, y lognormal."""
    rng = np.random.default_rng(0)
    n = 4_000_000
    worst = 0.0
    for sigma_y_n, rho_target in [(0.3, 0.6), (0.5, -0.4), (0.8, 0.3)]:
        # build x ~ N(0,1), ln y = mu_yn + sigma_yn*(rho_n*x + sqrt(1-rho_n^2)*z) so that
        # corr(x, ln y) = rho_n; then MEASURE the linear corr(x, y) and check the JAX
        # conversion maps that measured linear corr back to rho_n.
        rho_n = rho_target
        x = rng.standard_normal(n)
        z = rng.standard_normal(n)
        mu_yn = 2.0
        lny = mu_yn + sigma_y_n * (rho_n * x + np.sqrt(1 - rho_n**2) * z)
        y = np.exp(lny)
        corr_xy = np.corrcoef(x, y)[0, 1]            # linear correlation (the input)
        mu_y, sig_y = y.mean(), y.std()
        y_s2m2 = (sig_y / mu_y) ** 2
        got = float(corr_NL2NN(corr_xy, sigma_y_n, y_s2m2))
        worst = max(worst, abs(got - rho_n))
    assert worst < 5e-3, f"corr_NL2NN vs Monte-Carlo worst |Δ| {worst:.2e}"
    print(f"  corr_NL2NN vs Monte-Carlo: worst |Δ| {worst:.1e}  PASS")


def test_corr_LL2NN_vs_montecarlo():
    """corr_LL2NN(corr(x,y)) reproduces the measured corr(ln x, ln y) for x,y lognormal."""
    rng = np.random.default_rng(1)
    n = 4_000_000
    worst = 0.0
    for sx_n, sy_n, rho_target in [(0.3, 0.4, 0.5), (0.6, 0.5, -0.3), (0.7, 0.2, 0.2)]:
        a = rng.standard_normal(n)
        b = rng.standard_normal(n)
        lnx = 1.0 + sx_n * a
        lny = 2.0 + sy_n * (rho_target * a + np.sqrt(1 - rho_target**2) * b)
        x, y = np.exp(lnx), np.exp(lny)
        corr_xy = np.corrcoef(x, y)[0, 1]
        x_s2m2 = (x.std() / x.mean()) ** 2
        y_s2m2 = (y.std() / y.mean()) ** 2
        got = float(corr_LL2NN(corr_xy, sx_n, sy_n, x_s2m2, y_s2m2))
        worst = max(worst, abs(got - rho_target))
    assert worst < 5e-3, f"corr_LL2NN vs Monte-Carlo worst |Δ| {worst:.2e}"
    print(f"  corr_LL2NN vs Monte-Carlo: worst |Δ| {worst:.1e}  PASS")


def test_corr_clip_and_zero_sigma():
    """Correlation conversions clip to +/-0.99 and fall back to corr_x_y at sigma_n=0."""
    # inconsistent inputs that would overshoot |corr|>1 -> clipped
    big = float(corr_NL2NN(0.95, 0.01, 1.0))   # sqrt(1)/0.01 * 0.95 huge -> clip +0.99
    assert abs(big - MAX_MAG_CORRELATION) < 1e-15
    # sigma_y_n == 0 -> returns corr_x_y unchanged
    assert float(corr_NL2NN(0.42, 0.0, 0.0)) == 0.42
    assert float(corr_LL2NN(0.37, 0.0, 0.5, 0.0, 0.3)) == 0.37
    # differentiable
    g = float(jax.grad(lambda s: stdev_L2N(s))(0.5))
    assert np.isfinite(g) and g > 0
    print(f"  clip / zero-sigma fallback / differentiable: PASS")
