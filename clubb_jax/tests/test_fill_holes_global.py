#!/usr/bin/env python3
"""test_fill_holes_global.py — validate the mass-conserving global hole-fill.

`fill_holes_global(field, rho_dz, threshold, lower_k, upper_k)` (fill_holes.py ↔ fill_holes.F90:fill_holes_global)
raises sub-threshold "holes" up to the threshold over [lower_k, upper_k] while CONSERVING the ρ·dz-weighted mass:
it clips to max(threshold, field), then rescales (threshold + mass_frac·(clipped−threshold)) so the weighted mean is
unchanged. It is the conservative-fill primitive underlying the hole-filling family but was untested (only
fill_holes_vertical is). This pins its DEFINING properties: (1) ρ·dz-weighted MASS CONSERVATION over the fill range
(exact); (2) no holes remain (all in-range values ≥ threshold when there is enough mass); (3) values OUTSIDE the
range are untouched; (4) a finite jax.grad. Oracle-independent; never SKIPs. (iter 571)
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

from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_global

_NZ = 10


def _setup(rng):
    rho_dz = rng.uniform(0.5, 2.0, (1, _NZ))
    field = rng.uniform(0.5, 3.0, (1, _NZ))     # mean well above threshold so there is mass to fill
    field[0, 2] = -0.4                          # holes (sub-threshold)
    field[0, 5] = 0.05
    field[0, 8] = -0.1
    return field, rho_dz


def test_mass_conservation_and_no_holes():
    rng = np.random.default_rng(571)
    field, rho_dz = _setup(rng)
    thr = 0.1
    lo, hi = 1, 8
    out = np.asarray(fill_holes_global(jnp.asarray(field), jnp.asarray(rho_dz), thr, lo, hi))
    # (1) ρ·dz-weighted mass conserved over [lo, hi]
    m_in = np.sum(field[0, lo:hi + 1] * rho_dz[0, lo:hi + 1])
    m_out = np.sum(out[0, lo:hi + 1] * rho_dz[0, lo:hi + 1])
    assert abs(m_in - m_out) / (abs(m_in) + 1e-300) < 1e-13, f"mass not conserved: {m_in} vs {m_out}"
    # (2) no holes remain in range (weighted mean is above threshold -> all raised to >= threshold)
    assert np.all(out[0, lo:hi + 1] >= thr - 1e-12), f"holes remain: {out[0, lo:hi+1].min()}"
    # (3) outside [lo, hi] untouched
    assert out[0, 0] == field[0, 0] and out[0, 9] == field[0, 9], "out-of-range field modified"
    print(f"  mass-conserving (Δ {abs(m_in-m_out):.1e}); no holes ≥thr in [{lo},{hi}]; out-of-range untouched  PASS")


def test_grad_finite():
    rng = np.random.default_rng(7)
    field, rho_dz = _setup(rng)
    g = jax.grad(lambda f: jnp.sum(fill_holes_global(f, jnp.asarray(rho_dz), 0.1, 1, 8) ** 2))(jnp.asarray(field))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad"
    print("  jax.grad through fill_holes_global finite  PASS")


def main():
    print("test_fill_holes_global:")
    test_mass_conservation_and_no_holes()
    test_grad_finite()
    print("All fill_holes_global checks PASSED")


if __name__ == "__main__":
    main()
