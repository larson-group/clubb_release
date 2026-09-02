#!/usr/bin/env python3
"""test_fill_holes_sliding_window.py — validate the windowed mass-conserving hole-fill (fill_holes_type=2).

`fill_holes_sliding_window(field, rho_dz, threshold, lower_k, upper_k, num_draw=2)` (fill_holes.py ↔
fill_holes.F90:fill_holes_sliding_window) fills sub-threshold holes by sweeping a width-(2·num_draw+1) window over
the levels k ∈ [lower_k+num_draw, upper_k−num_draw], doing a per-window mass-conserving fill (same clip+rescale as
fill_holes_global), then a GLOBAL fallback (fill_holes_global over [lower_k, upper_k]) if any hole remains. Each
window operation is mass-neutral and the windows stay within [lower_k, upper_k], so the ρ·dz-weighted mass is
conserved overall and (via the fallback) no holes remain. Only fill_holes_vertical was tested. This pins those two
defining properties + that out-of-range values are untouched + a finite grad. Companion to iter-571's
fill_holes_global. Oracle-independent; never SKIPs. (iter 572)
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

from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_sliding_window

_NZ = 12


def _setup(rng):
    rho_dz = rng.uniform(0.5, 2.0, (1, _NZ))
    field = rng.uniform(0.8, 3.0, (1, _NZ))
    for k in (2, 4, 7, 9):
        field[0, k] = rng.uniform(-0.3, 0.05)   # holes scattered through the range
    return field, rho_dz


def test_conservation_and_no_holes():
    rng = np.random.default_rng(572)
    field, rho_dz = _setup(rng)
    thr = 0.1
    lo, hi = 1, 10           # window k sweeps 3..8, fallback over [1,10]
    out = np.asarray(fill_holes_sliding_window(jnp.asarray(field), jnp.asarray(rho_dz), thr, lo, hi))
    m_in = np.sum(field[0, lo:hi + 1] * rho_dz[0, lo:hi + 1])
    m_out = np.sum(out[0, lo:hi + 1] * rho_dz[0, lo:hi + 1])
    assert abs(m_in - m_out) / (abs(m_in) + 1e-300) < 1e-12, f"mass not conserved: {m_in} vs {m_out}"
    assert np.all(out[0, lo:hi + 1] >= thr - 1e-12), f"holes remain: {out[0, lo:hi+1].min()}"
    assert out[0, 0] == field[0, 0] and out[0, 11] == field[0, 11], "out-of-range modified"
    print(f"  windowed fill: ρ·dz mass conserved (Δ {abs(m_in-m_out):.1e}); no holes ≥thr; out-of-range untouched  PASS")


def test_no_holes_is_noop():
    # field already everywhere ≥ threshold -> unchanged (mass trivially conserved).
    rng = np.random.default_rng(3)
    rho_dz = rng.uniform(0.5, 2.0, (1, _NZ))
    field = rng.uniform(0.5, 3.0, (1, _NZ))   # all > thr=0.1
    out = np.asarray(fill_holes_sliding_window(jnp.asarray(field), jnp.asarray(rho_dz), 0.1, 1, 10))
    assert np.max(np.abs(out - field)) < 1e-13, "no-hole field should be unchanged"
    print("  no-hole field: unchanged (noop)  PASS")


def test_grad_finite():
    rng = np.random.default_rng(9)
    field, rho_dz = _setup(rng)
    g = jax.grad(lambda f: jnp.sum(fill_holes_sliding_window(f, jnp.asarray(rho_dz), 0.1, 1, 10) ** 2))(jnp.asarray(field))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad"
    print("  jax.grad through fill_holes_sliding_window finite  PASS")


def main():
    print("test_fill_holes_sliding_window:")
    test_conservation_and_no_holes()
    test_no_holes_is_noop()
    test_grad_finite()
    print("All fill_holes_sliding_window checks PASSED")


if __name__ == "__main__":
    main()
