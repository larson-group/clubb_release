#!/usr/bin/env python3
"""validate the new-hybrid PDF leaf routines (new_hybrid_pdf.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.new_hybrid_pdf import calculate_responder_params

NG, NZ = 2, 8


def test_responder_zero_cov_gaussian():
    sh = (NG, NZ)
    rng = np.random.default_rng(2)
    xm = rng.uniform(285, 305, sh); xp2 = rng.uniform(0.05, 1, sh)
    Skx = rng.uniform(-2, 2, sh); wp2 = rng.uniform(0.05, 2, sh)
    F_w = rng.uniform(0.05, 0.95, sh); mf = rng.uniform(0.2, 0.8, sh)
    wpxp = np.zeros(sh)                                  # |<w'x'>| = 0 → single Gaussian
    mu1, mu2, s1, s2, c1, c2 = (np.asarray(o) for o in
                                calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mf))
    assert np.allclose(mu1, xm) and np.allclose(mu2, xm), "zero-cov means not collapsed to xm"
    assert np.allclose(s1, xp2) and np.allclose(s2, xp2), "zero-cov variances not xp2"
    assert np.allclose(c1, 1.0) and np.allclose(c2, 1.0), "zero-cov coefs not 1"
    print("  responder |<w'x'>|=0 → single Gaussian  PASS")
