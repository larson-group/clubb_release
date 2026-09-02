#!/usr/bin/env python3
"""test_hydrometp2_zt.py — pin the precipitating-hydrometeor overall-variance formula.

`hydrometp2_zt` (setup_clubb_pdf_params.py ↔ setup_clubb_pdf_params.F90:449) converts a hydrometeor's
PRESCRIBED in-precip variance ratio (hmp2_ip_on_hmm2_ip) into the OVERALL grid-box variance <hm'²> used to seed the
hydrometeor PDF for the KK / Morrison microphysics. It is on the live microphysics PDF-param path (callers:
kk_microphys_step, morrison_microphys_step, kk_microphys_driver) but has no f2py wrapper and no isolation test —
only end-to-end coverage. The Fortran formula (the non-zero branch; the caller zeroes sub-tolerance levels):

    hydrometp2_zt = ( (hmp2_ip_on_hmm2_ip + 1) / precip_frac − 1 ) · hmm²

This pins it against an INDEPENDENT transcription over positive precip fractions (**exact**), checks the safe-division
guard (precip_frac=0 must not NaN/Inf — the JAX clamps the denominator, the caller overwrites those levels), the
in-cloud limit (precip_frac→1, ratio→0 ⇒ <hm'²>→0), and a finite jax.grad. Oracle-independent; never SKIPs. (iter 547)
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

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import hydrometp2_zt

_NG, _NZT = 2, 8


def test_matches_f90_formula():
    rng = np.random.default_rng(547)
    hmm = rng.uniform(1e-6, 1e-3, (_NG, _NZT))         # hydrometeor mean (e.g. rain mixing ratio)
    precip_frac = rng.uniform(0.05, 1.0, (_NG, _NZT))  # positive precip fraction (the meaningful branch)
    ratio = rng.uniform(0.0, 5.0, (_NG, _NZT))         # hmp2_ip_on_hmm2_ip
    got = np.asarray(hydrometp2_zt(jnp.asarray(hmm), jnp.asarray(precip_frac), jnp.asarray(ratio)))
    ref = ((ratio + 1.0) / precip_frac - 1.0) * hmm ** 2
    worst = float(np.max(np.abs(got - ref)))
    assert worst < 1e-18, f"hydrometp2_zt mismatch vs F90 formula {worst:.2e}"
    print(f"  hydrometp2_zt = ((ratio+1)/precip_frac − 1)·hmm² matches F90 (worst {worst:.1e})  PASS")


def test_safe_division_and_in_cloud_limit():
    # precip_frac=0 must not produce NaN/Inf (the safe-division guard; caller overwrites these levels).
    hmm = np.full((_NG, _NZT), 1e-4)
    pf0 = np.zeros((_NG, _NZT))
    ratio = np.full((_NG, _NZT), 2.0)
    out0 = np.asarray(hydrometp2_zt(jnp.asarray(hmm), jnp.asarray(pf0), jnp.asarray(ratio)))
    assert np.all(np.isfinite(out0)), "precip_frac=0 produced non-finite output (safe-division guard broken)"
    # In-cloud limit: precip_frac→1 and ratio→0 ⇒ <hm'²> → (1/1 − 1)·hmm² = 0.
    one = np.ones((_NG, _NZT)); zero = np.zeros((_NG, _NZT))
    out1 = np.asarray(hydrometp2_zt(jnp.asarray(hmm), jnp.asarray(one), jnp.asarray(zero)))
    assert np.max(np.abs(out1)) < 1e-20, "precip_frac=1, ratio=0 must give zero variance (fully in-cloud, no spread)"
    print("  safe-division at precip_frac=0 (finite) + in-cloud limit (pf=1,ratio=0 ⇒ 0)  PASS")


def test_grad_finite():
    rng = np.random.default_rng(7)
    hmm = jnp.asarray(rng.uniform(1e-6, 1e-3, (_NG, _NZT)))
    pf = jnp.asarray(rng.uniform(0.05, 1.0, (_NG, _NZT)))
    ratio = jnp.asarray(rng.uniform(0.0, 5.0, (_NG, _NZT)))
    g = jax.grad(lambda h: jnp.sum(hydrometp2_zt(h, pf, ratio) ** 2))(hmm)
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad wrt hmm"
    print("  jax.grad(hydrometp2_zt) wrt hmm finite  PASS")


def main():
    print("test_hydrometp2_zt:")
    test_matches_f90_formula()
    test_safe_division_and_in_cloud_limit()
    test_grad_finite()
    print("All hydrometp2_zt checks PASSED")


if __name__ == "__main__":
    main()
