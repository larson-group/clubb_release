#!/usr/bin/env python3
"""test_KK_Nrm_tendencies.py — pin the local KK rain-number (N_r) tendency formulas.

`KK_Nrm_auto_mean` and `KK_Nrm_evap_local_mean` (KK_Nrm_tendencies.py ↔ KK_Nrm_tendencies.F90:296/205) convert the
KK rain-MASS tendencies into rain-NUMBER tendencies. They are on the live KK path (rico) but had no isolation test:
  KK_Nrm_auto_mean(rate)            = rate / ( (4/3)·π·ρ_lw·r₀³ )          — one new drop of radius r₀ per unit mass
  KK_Nrm_evap_local_mean(ev,Nr,rr)  = dt^(ν−1) · (Nr / rr^ν) · ev^ν,  ν=1  ⇒  (Nr/rr)·ev
This pins both against an INDEPENDENT transcription of the F90 (so the ρ_lw=1000 / r₀=2.5e-5 drop-mass constant and the
ν=1 reduction are checked directly), the safe-division guard (rr=0 → finite), the auto rate being exactly linear (one
drop per unit mass), and a finite jax.grad. Oracle-independent; never SKIPs. (iter 557)
"""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = str(_REPO_ROOT)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import (
    KK_Nrm_auto_mean, KK_Nrm_evap_local_mean)

_RHO_LW, _R0 = 1000.0, 2.5e-5
_PI = 3.141592654  # constants_clubb.F90 deliberately uses this rounded literal.
_DROP_MASS = (4.0 / 3.0) * _PI * _RHO_LW * _R0 ** 3   # mass of one radius-r0 drop


def test_auto_mean_drop_mass_constant():
    worst = 0.0
    for rate in (1.0e-8, 5.0e-7, -2.0e-7, 0.0):
        got = float(KK_Nrm_auto_mean(rate))
        ref = rate / _DROP_MASS
        worst = max(worst, abs(got - ref))
        np.testing.assert_allclose(got, ref, rtol=2e-15, atol=1e-20)
    # linearity: doubling the mass rate doubles the number rate (one drop per unit mass).
    assert abs(float(KK_Nrm_auto_mean(2.0e-7)) - 2.0 * float(KK_Nrm_auto_mean(1.0e-7))) < 1e-12
    print(f"  KK_Nrm_auto_mean = rate/((4/3)πρ_lw r0³) (drop-mass {_DROP_MASS:.3e}); linear  PASS")


def test_evap_local_mean_formula_and_safe_div():
    # ν=1 ⇒ (Nr/rr)·evap, dt^(ν−1)=dt^0=1 (so dt-independent at ν=1).
    worst = 0.0
    rng = np.random.default_rng(557)
    for _ in range(20):
        ev = rng.uniform(-1e-6, 0.0)          # evaporation is a sink (negative)
        Nr = rng.uniform(1e3, 1e7)
        rr = rng.uniform(1e-6, 1e-3)
        for dt in (60.0, 300.0):
            got = float(KK_Nrm_evap_local_mean(ev, Nr, rr, dt))
            ref = (Nr / rr) * ev          # ν=1 reduction
            worst = max(worst, abs(got - ref) / (abs(ref) + 1e-300))
    assert worst < 1e-13, f"KK_Nrm_evap_local_mean vs (Nr/rr)·ev {worst:.2e}"
    # dt-independence at ν=1
    assert abs(float(KK_Nrm_evap_local_mean(-1e-7, 1e5, 1e-4, 60.0))
               - float(KK_Nrm_evap_local_mean(-1e-7, 1e5, 1e-4, 999.0))) < 1e-18, "should be dt-independent at ν=1"
    # safe division: rr=0 must not NaN/Inf.
    assert np.isfinite(float(KK_Nrm_evap_local_mean(-1e-7, 1e5, 0.0, 60.0))), "rr=0 not finite (safe-div broken)"
    print(f"  KK_Nrm_evap_local_mean = (Nr/rr)·ev (ν=1, dt-independent); safe at rr=0 (worst {worst:.1e})  PASS")


def test_grad_finite():
    g1 = jax.grad(lambda r: KK_Nrm_auto_mean(r) ** 2)(1.0e-7)
    g2 = jax.grad(lambda rr: KK_Nrm_evap_local_mean(-1e-7, 1e5, rr, 60.0) ** 2)(1.0e-4)
    assert np.isfinite(float(g1)) and np.isfinite(float(g2)), "non-finite grad"
    # Differentiate the independent drop-mass and nu=1 formulas at nonzero smooth inputs.
    np.testing.assert_allclose(g1, 2.0 * 1.0e-7 / _DROP_MASS**2, rtol=1e-13)
    np.testing.assert_allclose(g2, -2.0 * (-1e-7 * 1e5)**2 / 1.0e-4**3, rtol=1e-13)
    print("  N_r gradients match independent analytic derivatives  PASS")


def main():
    print("test_KK_Nrm_tendencies:")
    test_auto_mean_drop_mass_constant()
    test_evap_local_mean_formula_and_safe_div()
    test_grad_finite()
    print("All KK_Nrm tendency checks PASSED")


if __name__ == "__main__":
    main()
