#!/usr/bin/env python3
"""pin the precipitating-hydrometeor overall-variance formula."""


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
