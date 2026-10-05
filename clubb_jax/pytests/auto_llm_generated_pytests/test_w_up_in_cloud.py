#!/usr/bin/env python3
"""validate the JAX calc_w_up_in_cloud port."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_w_up_in_cloud

NG, NZ = 2, 8


def _fields(seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0.2, 0.8, (NG, NZ))
    cf1 = rng.uniform(0.0, 1.0, (NG, NZ)); cf2 = rng.uniform(0.0, 1.0, (NG, NZ))
    w1 = rng.uniform(-3, 3, (NG, NZ)); w2 = rng.uniform(-3, 3, (NG, NZ))
    v1 = rng.uniform(0.01, 1.0, (NG, NZ)); v2 = rng.uniform(0.01, 1.0, (NG, NZ))
    return a, cf1, cf2, w1, w2, v1, v2


def test_invariants():
    a, cf1, cf2, w1, w2, v1, v2 = _fields(5)
    w_up, w_down, uf, df = (np.asarray(x) for x in calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2))
    assert np.all(uf >= -1e-12) and np.all(uf <= 1.0 + 1e-12), "updraft frac out of [0,1]"
    assert np.all(df >= -1e-12) and np.all(df <= 1.0 + 1e-12), "downdraft frac out of [0,1]"
    assert np.all(w_up >= -1e-9), "mean cloudy updraft should be >= 0"
    assert np.all(w_down <= 1e-9), "mean cloudy downdraft should be <= 0"
    # All-updraft shortcut: a strongly positive w with tiny variance -> updraft_frac ~ 1, w_up ~ w.
    a1 = np.ones((1, 1)); cf = np.ones((1, 1))
    w_up1, _, uf1, _ = (np.asarray(x) for x in
                        calc_w_up_in_cloud(a1, cf, cf, 5.0 * np.ones((1, 1)), np.zeros((1, 1)),
                                           1e-4 * np.ones((1, 1)), 1e-4 * np.ones((1, 1))))
    assert abs(uf1[0, 0] - 1.0) < 1e-6 and abs(w_up1[0, 0] - 5.0) < 1e-3, "all-updraft shortcut wrong"
    print("  invariants: fractions in [0,1], w_up>=0>=w_down, all-updraft shortcut  PASS")
