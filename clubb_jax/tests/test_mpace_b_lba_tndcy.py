#!/usr/bin/env python3
"""test_mpace_b_lba_tndcy.py — validate the M-PACE B large-scale forcing + the (zero) LBA forcing.

mpace_b_tndcy is checked vs a literal NumPy transcription + physical invariants (subsidence sign, capped
cooling, wm_zm boundary conditions) + a finite jax.grad. lba_tndcy is identically zero.
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
import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.mpace_b import mpace_b_tndcy
from clubb_jax.src.Benchmark_cases.lba import lba_tndcy
from clubb_jax.src.CLUBB_core.constants_clubb import Rd, Cp, grav, g_per_kg, sec_per_day
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.CLUBB_core.grid_class import setup_grid


def _grid(nzt=40, dz=100.0):
    return setup_grid(1, dz, dz, dz * (nzt + 1))


def _ref(p, thvm, gr):
    D, p_sfc, pinv = 5.8e-6, 101000.0, 85000.0
    omega = np.minimum(D * (p_sfc - p), D * (p_sfc - pinv))
    wm_zt = -omega * Rd * thvm / p / grav
    wm_zm = np.array(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, jnp.asarray(wm_zt)))   # writable copy
    wm_zm[:, 0] = 0.0; wm_zm[:, -1] = 0.0
    tt = np.minimum(-4.0, -15.0 * (1.0 - (p_sfc - p) / 21818.0))
    thlm_f = (tt * (p_sfc / p) ** (Rd / Cp)) / sec_per_day
    rtm_f = np.minimum(0.164, -3.0 * (1.0 - (p_sfc - p) / 15171.0)) / g_per_kg / sec_per_day
    return wm_zt, wm_zm, thlm_f, rtm_f


_SCLR_IDX = SclrIdx(-1, -1, -1, -1, -1, -1)


def test_mpace_b_vs_literal():
    gr = _grid()
    nzt = gr.zt.shape[1]
    rng = np.random.default_rng(4)
    p = np.linspace(101000.0, 40000.0, nzt)[None, :]      # descending pressure with height
    thvm = (285.0 + 0.004 * gr.zt[0] + rng.standard_normal(nzt) * 0.3)[None, :]
    out = mpace_b_tndcy(1, 0, 0, _SCLR_IDX, gr, p, thvm)
    rwzt, rwzm, rth, rrt = _ref(p, thvm, gr)
    for name, got, ref in (('wm_zt', out[0], rwzt), ('wm_zm', out[1], rwzm),
                           ('thlm_forcing', out[2], rth), ('rtm_forcing', out[3], rrt)):
        rel = np.max(np.abs(np.asarray(got) - ref) / (np.abs(ref) + 1e-30))
        assert rel < 1e-12, f"{name} vs literal rel {rel:.2e}"
    print("  mpace_b_tndcy vs literal NumPy (all 4 outputs): rel <1e-12  PASS")


def test_mpace_b_physical():
    gr = _grid()
    nzt = gr.zt.shape[1]
    p = np.linspace(101000.0, 40000.0, nzt)[None, :]
    thvm = (288.0 + 0.004 * gr.zt[0])[None, :]
    out = mpace_b_tndcy(1, 0, 0, _SCLR_IDX, gr, p, thvm)
    # Subsidence: omega>0 (descent) → wm_zt < 0 wherever p < p_sfc.
    assert np.all(np.asarray(out[0])[:, 1:] <= 0.0), "subsidence wm_zt should be <= 0"
    # wm_zm boundary conditions.
    assert out[1][0, 0] == 0.0 and out[1][0, -1] == 0.0, "wm_zm BCs not zero"
    # Cooling capped at -4 K/day → thlm_forcing <= -4/86400 * exner (negative).
    assert np.all(np.asarray(out[2]) < 0.0), "thlm forcing should be cooling (<0)"
    print("  mpace_b_tndcy physical: subsidence<=0, wm_zm BCs=0, cooling<0  PASS")


def test_lba_tndcy_zero():
    gr = _grid(nzt=30)
    th, rt, sclr, edsclr = lba_tndcy(1, 0, 0, _SCLR_IDX, gr)
    assert not np.any(np.asarray(th)) and not np.any(np.asarray(rt)), "lba forcing should be identically zero"
    assert sclr.shape[-1] == 0 and edsclr.shape[-1] == 0
    print("  lba_tndcy: identically zero (surface-driven LBA, no LS forcing)  PASS")


def test_differentiable():
    gr = _grid()
    nzt = gr.zt.shape[1]
    p = jnp.asarray(np.linspace(101000.0, 40000.0, nzt)[None, :])
    g = jax.grad(lambda thvm: jnp.sum(mpace_b_tndcy(1, 0, 0, _SCLR_IDX, gr, p, thvm)[0] ** 2))(
        jnp.asarray((288.0 + 0.004 * np.asarray(gr.zt[0]))[None, :]))
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad"
    print("  jax.grad(mpace_b_tndcy) wrt thvm: finite  PASS")


def main():
    print("test_mpace_b_lba_tndcy:")
    for t in (test_mpace_b_vs_literal, test_mpace_b_physical, test_lba_tndcy_zero, test_differentiable):
        t()
    print("All mpace_b/lba tndcy checks PASSED")


if __name__ == "__main__":
    main()
