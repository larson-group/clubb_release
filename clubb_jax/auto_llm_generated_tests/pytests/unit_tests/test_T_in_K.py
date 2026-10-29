#!/usr/bin/env python3
"""validate the T_in_K_module port (thlm <-> absolute-T conversions)."""


import numpy as np

from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K, T_in_K2thlm
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv


def test_closed_form():
    thlm, exner, rcm = 300.0, 0.9, 1e-3
    assert np.isclose(thlm2T_in_K(thlm, exner, rcm), thlm * exner + Lv * rcm / Cp, rtol=0, atol=0)
    T = 285.0
    assert np.isclose(T_in_K2thlm(T, exner, rcm), (T - Lv / Cp * rcm) / exner, rtol=0, atol=0)
    print("  thlm2T_in_K / T_in_K2thlm closed-form  PASS")


def test_exact_inverse():
    rng = np.random.default_rng(0)
    thlm = rng.uniform(270, 330, 5000)
    exner = rng.uniform(0.5, 1.0, 5000)
    rcm = rng.uniform(0.0, 3e-3, 5000)
    rt = T_in_K2thlm(thlm2T_in_K(thlm, exner, rcm), exner, rcm)
    worst = np.max(np.abs(rt - thlm) / np.abs(thlm))
    assert worst < 1e-13, f"round-trip worst rel {worst:.2e}"
    print(f"  T_in_K2thlm is the exact inverse of thlm2T_in_K: {len(thlm)} cases, worst rel {worst:.1e}  PASS")
