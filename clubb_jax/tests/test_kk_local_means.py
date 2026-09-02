#!/usr/bin/env python3
"""test_kk_local_means.py — validate the JAX KK_local_means port (grid-mean KK rates).

No f2py oracle exists for the local-mean variants (the gated KK path uses the upscaled means), so these
power laws are validated against an INDEPENDENT NumPy transcription of the documented `parameters_KK.F90`
exponents, the supersaturation/cloud-presence branch (`s<=0` evap, `rc>0` auto/accr), and a finite gradient.
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

from clubb_jax.src.Microphys.KK_microphys.KK_local_means import (
    KK_evap_local_mean, KK_auto_local_mean, KK_accr_local_mean, KK_mvr_local_mean)


def test_formulas_and_branches():
    rng = np.random.default_rng(20260604)
    n = 200
    s = rng.standard_normal(n) * 1e-3          # mix of >0 and <=0
    rrm = np.abs(rng.standard_normal(n)) * 1e-4 + 1e-8
    Nrm = np.abs(rng.standard_normal(n)) * 1e6 + 1.0
    Ncm = np.abs(rng.standard_normal(n)) * 1e8 + 1.0
    rc = np.abs(rng.standard_normal(n)) * 1e-4   # >=0, but inject some zeros/negatives
    rc[::7] = 0.0
    evap_c, auto_c, accr_c, mvr_c, Nr_tol = 7.0e-4, 1350.0, 67.0, 3.6e-3, 1.0e-3

    # Independent NumPy references (exponents from parameters_KK.F90).
    ref_evap = np.where(s <= 0.0, evap_c * s ** 1.0 * rrm ** (1 / 3) * Nrm ** (2 / 3), 0.0)
    ref_auto = np.where(rc > 0.0, auto_c * np.maximum(rc, 0.0) ** 2.47 * Ncm ** (-1.79), 0.0)
    ref_accr = np.where(rc > 0.0, accr_c * np.maximum(rc, 0.0) ** 1.15 * rrm ** 1.15, 0.0)
    ref_mvr = mvr_c * rrm ** (1 / 3) * np.maximum(Nrm, Nr_tol) ** (-1 / 3)

    got_evap = np.asarray(KK_evap_local_mean(jnp.asarray(s), jnp.asarray(rrm), jnp.asarray(Nrm), evap_c))
    got_auto = np.asarray(KK_auto_local_mean(jnp.asarray(rc), jnp.asarray(Ncm), auto_c))
    got_accr = np.asarray(KK_accr_local_mean(jnp.asarray(rc), jnp.asarray(rrm), accr_c))
    got_mvr = np.asarray(KK_mvr_local_mean(jnp.asarray(rrm), jnp.asarray(Nrm), mvr_c, Nr_tol))

    for name, got, ref in (("evap", got_evap, ref_evap), ("auto", got_auto, ref_auto),
                           ("accr", got_accr, ref_accr), ("mvr", got_mvr, ref_mvr)):
        rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
        assert rel < 1e-13, f"KK_{name}_local_mean rel {rel:.2e}"
    # Branch coverage: evap is 0 where s>0; auto/accr are 0 where rc<=0.
    assert np.all(got_evap[s > 0.0] == 0.0) and np.all(got_auto[rc <= 0.0] == 0.0) \
        and np.all(got_accr[rc <= 0.0] == 0.0), "rate branch not respected"
    print("  KK evap/auto/accr/mvr local means vs independent NumPy + branches: rel <1e-13  PASS")


def test_differentiable():
    rrm = jnp.asarray(np.abs(np.random.default_rng(1).standard_normal(50)) * 1e-4 + 1e-6)
    Nrm = jnp.asarray(np.abs(np.random.default_rng(2).standard_normal(50)) * 1e6 + 1.0)
    g = jax.grad(lambda r: jnp.sum(KK_mvr_local_mean(r, Nrm, 3.6e-3, 1e-3)))(rrm)
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad through KK_mvr_local_mean"
    # auto rate grad w.r.t. positive rc
    rc = jnp.asarray(np.abs(np.random.default_rng(3).standard_normal(50)) * 1e-4 + 1e-6)
    Ncm = jnp.asarray(np.full(50, 1e8))
    ga = jax.grad(lambda c: jnp.sum(KK_auto_local_mean(c, Ncm, 1350.0)))(rc)
    assert np.isfinite(np.asarray(ga)).all(), "non-finite grad through KK_auto_local_mean"
    print(f"  jax.grad through KK_mvr/KK_auto local means: finite  PASS")


def main():
    print("test_kk_local_means:")
    for t in (test_formulas_and_branches, test_differentiable):
        t()
    print("All KK_local_means checks PASSED")


if __name__ == "__main__":
    main()
