#!/usr/bin/env python3
"""test_component_precip_frac_specify.py — property-pin the per-component precip-fraction split.

`component_precip_frac_specify` (precipitation_fraction.py ↔ precipitation_fraction.F90, precip_frac_calc_type=2,
the fixed default) splits the overall precip fraction `pf` into the two PDF-component fractions (pf1, pf2) via the
specified-upsilon method. Its GENERAL (0<upsilon<1) branch is bit-exact-validated end-to-end vs the rico oracle at
upsilon=0.55 (test_precip_fraction.py), but the `upsilon==1` branch (rico uses 0.55) was never exercised, and the
branches carry intricate clip + recompute-and-double-check sub-cases (too conditional to faithfully transcribe).

Instead of transcribing, this pins the DEFINING property — component-fraction CONSERVATION: a weighted average of the
component fractions returns the overall fraction, `mixt_frac·pf1 + (1−mixt_frac)·pf2 == pf`, which holds EXACTLY in the
non-clipping regime for BOTH branches (it is how pf2 is defined from pf1). Tested for the general branch AND the
upsilon==1 branch (so that branch is now exercised — a crash/NaN there is caught). Plus validity bounds (pf1,pf2 ∈
[0,1] across a broad sweep including clip regimes) and a finite jax.grad. Oracle-independent; never SKIPs. (iter 548)
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

from clubb_jax.src.CLUBB_core.precipitation_fraction import component_precip_frac_specify

_NG, _NZT = 2, 10
_TOL = 1.0e-3


def _call(pf, mf, upsilon):
    tol = np.full((_NG, 1), _TOL)
    return component_precip_frac_specify(jnp.asarray(pf), jnp.asarray(mf), jnp.asarray(tol), float(upsilon))


def test_conservation_general_branch_no_clip():
    """General branch (upsilon=0.55), inputs chosen so no clip fires: pf<mf ⇒ pf1=upsilon·pf/mf<1 and >tol."""
    rng = np.random.default_rng(548)
    mf = rng.uniform(0.55, 0.85, (_NG, _NZT))
    pf = mf * rng.uniform(0.2, 0.9, (_NG, _NZT))        # pf < mf
    pf1, pf2 = _call(pf, mf, 0.55)
    pf1, pf2 = np.asarray(pf1), np.asarray(pf2)
    recon = mf * pf1 + (1.0 - mf) * pf2
    worst = float(np.max(np.abs(recon - pf)))
    assert worst < 1e-13, f"general branch not conserving: mf·pf1+(1−mf)·pf2 vs pf {worst:.2e}"
    print(f"  general (upsilon=0.55), no clip: mf·pf1+(1−mf)·pf2 == pf (worst {worst:.1e})  PASS")


def test_conservation_upsilon_one_branch():
    """upsilon==1 branch — exercised here for the first time. pf<=mf: pf1=pf/mf,pf2=0; pf>mf (small): pf1=1,pf2=(pf−mf)/omf."""
    rng = np.random.default_rng(11)
    mf = rng.uniform(0.4, 0.7, (_NG, _NZT))
    # half the points pf<=mf (comp-1 only), half pf>mf but pf<mf+(1-mf) so pf2 unclipped
    pf = np.where(rng.random((_NG, _NZT)) < 0.5,
                  mf * rng.uniform(0.2, 0.95, (_NG, _NZT)),
                  mf + (1.0 - mf) * rng.uniform(0.05, 0.8, (_NG, _NZT)))
    pf1, pf2 = _call(pf, mf, 1.0)
    pf1, pf2 = np.asarray(pf1), np.asarray(pf2)
    assert np.all(np.isfinite(pf1)) and np.all(np.isfinite(pf2)), "upsilon==1 branch produced non-finite output"
    recon = mf * pf1 + (1.0 - mf) * pf2
    worst = float(np.max(np.abs(recon - pf)))
    assert worst < 1e-13, f"upsilon==1 branch not conserving {worst:.2e}"
    print(f"  upsilon==1 branch (now exercised): conserves, finite (worst {worst:.1e})  PASS")


def test_validity_bounds_broad_sweep():
    """Across a broad random sweep (incl. clip regimes), both component fractions stay valid: pf1,pf2 ∈ [0,1]."""
    rng = np.random.default_rng(3)
    worst_lo = 0.0; worst_hi = 0.0
    for ups in (0.3, 0.55, 0.9, 1.0):
        mf = rng.uniform(0.1, 0.9, (_NG, _NZT))
        pf = rng.uniform(0.0, 1.0, (_NG, _NZT))
        pf1, pf2 = _call(pf, mf, ups)
        for a in (np.asarray(pf1), np.asarray(pf2)):
            assert np.all(np.isfinite(a)), f"non-finite component fraction at upsilon={ups}"
            worst_lo = min(worst_lo, float(np.min(a)))
            worst_hi = max(worst_hi, float(np.max(a)))
    assert worst_lo >= -1e-12 and worst_hi <= 1.0 + 1e-12, f"component fractions escaped [0,1]: [{worst_lo},{worst_hi}]"
    print(f"  validity: pf1,pf2 ∈ [0,1] across upsilon∈{{0.3,0.55,0.9,1}} (range [{worst_lo:.2f},{worst_hi:.2f}])  PASS")


def test_grad_finite():
    rng = np.random.default_rng(99)
    mf = jnp.asarray(rng.uniform(0.4, 0.8, (_NG, _NZT)))
    pf0 = jnp.asarray(mf * rng.uniform(0.2, 0.9, (_NG, _NZT)))
    tol = jnp.asarray(np.full((_NG, 1), _TOL))
    g = jax.grad(lambda p: jnp.sum(sum(x ** 2 for x in component_precip_frac_specify(p, mf, tol, 0.55))))(pf0)
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad wrt pf"
    print("  jax.grad(component_precip_frac_specify) wrt pf finite  PASS")


def main():
    print("test_component_precip_frac_specify:")
    test_conservation_general_branch_no_clip()
    test_conservation_upsilon_one_branch()
    test_validity_bounds_broad_sweep()
    test_grad_finite()
    print("All component_precip_frac_specify checks PASSED")


if __name__ == "__main__":
    main()
