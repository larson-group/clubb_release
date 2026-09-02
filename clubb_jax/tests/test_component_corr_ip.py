#!/usr/bin/env python3
"""test_component_corr_ip.py — validate the JAX component_corr_*_ip ports.

These complete the in-precip PDF-component hydrometeor-correlation family in setup_clubb_pdf_params.F90:
  * component_corr_w_hm_n_ip   (:2669) — w & ln(hm); pass-through if l_calc_w_corr else cloud/below by rc_i
  * component_corr_x_hm_n_ip   (:2754) — x(=chi/eta) & ln(hm); cloud/below by rc_i
  * component_corr_hmx_hmy_n_ip(:2823) — ln(hmx) & ln(hmy); cloud/below by rc_i
  * component_corr_eta_hm_n_ip (:2892) — eta & ln(hm) estimated as corr_chi_eta * corr_chi_hm_n

Oracle: literal NumPy transcription of each (rc_tol cloud/below selection, the l_calc_w_corr branch, and the
product identity), over a grid that straddles rc_tol; plus a finite jax.grad through the eta-product estimate.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import (
    component_corr_w_hm_n_ip, component_corr_x_hm_n_ip,
    component_corr_hmx_hmy_n_ip, component_corr_eta_hm_n_ip,
    component_corr_w_x, component_corr_chi_eta)
from clubb_jax.src.CLUBB_core.constants_clubb import rc_tol, max_mag_correlation


def _rc_fields(seed):
    rng = np.random.default_rng(seed)
    # Values straddling rc_tol (1e-6): some below, some above.
    rc_1 = rng.choice([0.0, 0.5 * rc_tol, 2.0 * rc_tol, 1e-3], size=(3, 7))
    rc_2 = rng.choice([0.0, 0.5 * rc_tol, 2.0 * rc_tol, 1e-3], size=(3, 7))
    return rc_1, rc_2


def _ref_cloud_below(rc_1, rc_2, c_cloud, c_below):
    return (np.where(rc_1 > rc_tol, c_cloud, c_below),
            np.where(rc_2 > rc_tol, c_cloud, c_below))


def test_x_hm():
    rc_1, rc_2 = _rc_fields(1)
    c_cloud, c_below = 0.7, -0.3
    g1, g2 = component_corr_x_hm_n_ip(rc_1, rc_2, c_cloud, c_below)
    r1, r2 = _ref_cloud_below(rc_1, rc_2, c_cloud, c_below)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2)
    print("  component_corr_x_hm_n_ip: cloud/below by rc_tol (exact vs NumPy)  PASS")


def test_hmx_hmy():
    rc_1, rc_2 = _rc_fields(2)
    c_cloud, c_below = 0.9, 0.1
    g1, g2 = component_corr_hmx_hmy_n_ip(rc_1, rc_2, c_cloud, c_below)
    r1, r2 = _ref_cloud_below(rc_1, rc_2, c_cloud, c_below)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2)
    print("  component_corr_hmx_hmy_n_ip: cloud/below by rc_tol (exact vs NumPy)  PASS")


def test_w_hm_both_branches():
    rc_1, rc_2 = _rc_fields(3)
    c_cloud, c_below = 0.6, -0.2
    rng = np.random.default_rng(4)
    corr_in_1 = rng.uniform(-0.9, 0.9, size=rc_1.shape)
    corr_in_2 = rng.uniform(-0.9, 0.9, size=rc_2.shape)
    # l_calc_w_corr=False -> cloud/below selection (corr_in ignored).
    g1, g2 = component_corr_w_hm_n_ip(corr_in_1, rc_1, corr_in_2, rc_2, c_cloud, c_below, False)
    r1, r2 = _ref_cloud_below(rc_1, rc_2, c_cloud, c_below)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2), "l_calc=False branch"
    # l_calc_w_corr=True -> pass-through of the diagnosed correlations.
    p1, p2 = component_corr_w_hm_n_ip(corr_in_1, rc_1, corr_in_2, rc_2, c_cloud, c_below, True)
    assert np.array_equal(np.asarray(p1), corr_in_1) and np.array_equal(np.asarray(p2), corr_in_2), "passthrough"
    print("  component_corr_w_hm_n_ip: both l_calc_w_corr branches (cloud/below + passthrough)  PASS")


def test_eta_product_and_grad():
    rng = np.random.default_rng(5)
    ce1 = rng.uniform(-0.9, 0.9, size=(3, 7))
    ce2 = rng.uniform(-0.9, 0.9, size=(3, 7))
    ch1 = rng.uniform(-0.9, 0.9, size=(3, 7))
    ch2 = rng.uniform(-0.9, 0.9, size=(3, 7))
    g1, g2 = component_corr_eta_hm_n_ip(ce1, ch1, ce2, ch2)
    assert np.max(np.abs(np.asarray(g1) - ce1 * ch1)) < 1e-15
    assert np.max(np.abs(np.asarray(g2) - ce2 * ch2)) < 1e-15
    # Differentiable: d(corr_eta_1)/d(corr_chi_eta_1) = corr_chi_hm_n_1.
    import jax.numpy as jnp
    def loss(ce):
        a, _ = component_corr_eta_hm_n_ip(ce, jnp.asarray(ch1), jnp.asarray(ce2), jnp.asarray(ch2))
        return jnp.sum(a)
    grad = np.asarray(jax.grad(loss)(jnp.asarray(ce1)))
    assert np.isfinite(grad).all() and np.max(np.abs(grad - ch1)) < 1e-12, "grad wrong"
    print("  component_corr_eta_hm_n_ip: product identity + grad == corr_chi_hm_n (exact)  PASS")


def test_w_x():
    rc_1, rc_2 = _rc_fields(6)
    c_cloud, c_below = 0.4, -0.1
    # ADG standards on + an ADG PDF type -> zeros regardless of rc.
    for iipdf in (1, 2, 7):
        g1, g2 = component_corr_w_x(rc_1, rc_2, c_cloud, c_below, iipdf, True)
        assert np.all(np.asarray(g1) == 0.0) and np.all(np.asarray(g2) == 0.0), f"ADG type {iipdf} not zeroed"
    # ADG standards off -> cloud/below selection even for an ADG type.
    g1, g2 = component_corr_w_x(rc_1, rc_2, c_cloud, c_below, 1, False)
    r1, r2 = _ref_cloud_below(rc_1, rc_2, c_cloud, c_below)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2), "ADG-off should select"
    # Non-ADG type (e.g. 6 = LY93) -> cloud/below regardless of the flag.
    g1, g2 = component_corr_w_x(rc_1, rc_2, c_cloud, c_below, 6, True)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2), "non-ADG should select"
    print("  component_corr_w_x: ADG-zero branch + prescribed cloud/below branch  PASS")


def test_chi_eta():
    rc_1, rc_2 = _rc_fields(7)
    # Cloud value exceeds max_mag_correlation -> with l_limit it must clamp; without, pass through.
    c_cloud, c_below = 1.3, -1.5
    g1, g2 = component_corr_chi_eta(rc_1, rc_2, c_cloud, c_below, False)
    r1, r2 = _ref_cloud_below(rc_1, rc_2, c_cloud, c_below)
    assert np.array_equal(np.asarray(g1), r1) and np.array_equal(np.asarray(g2), r2), "no-limit branch"
    l1, l2 = component_corr_chi_eta(rc_1, rc_2, c_cloud, c_below, True)
    lr1 = np.clip(r1, -max_mag_correlation, max_mag_correlation)
    lr2 = np.clip(r2, -max_mag_correlation, max_mag_correlation)
    assert np.array_equal(np.asarray(l1), lr1) and np.array_equal(np.asarray(l2), lr2), "limit/clamp branch"
    assert np.max(np.abs(np.asarray(l1))) <= max_mag_correlation + 1e-15, "clamp failed"
    print("  component_corr_chi_eta: cloud/below + optional ±max_mag_correlation clamp  PASS")


def main():
    print("test_component_corr_ip:")
    for t in (test_x_hm, test_hmx_hmy, test_w_hm_both_branches, test_eta_product_and_grad,
              test_w_x, test_chi_eta):
        t()
    print("All component_corr_*_ip checks PASSED")


if __name__ == "__main__":
    main()
