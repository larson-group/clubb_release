#!/usr/bin/env python3
"""test_lba_sfclyr.py — validate the LBA diurnal surface scheme.

Oracle: a literal NumPy transcription of the Fortran (the analytic diurnal factor + bulk-flux formulas;
diag_ustar is the already-validated MOST routine), plus the physical diurnal structure (peak at 5.25 h, zero
outside the cosine window, fluxes >= 0, ustar > 0) and a finite jax.grad.
"""
import os
import sys
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.lba import lba_sfclyr
from clubb_jax.src.Benchmark_cases.diag_ustar_module import diag_ustar
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, grav, sec_per_hr


def _ref(t, z, rho, thlm, ubar):
    ft = max(0.0, math.cos(0.5 * math.pi * (5.25 - t / sec_per_hr) / 5.25))
    wpthlp = (270.0 * ft ** 1.5) / (rho * Cp)
    wprtp = (554.0 * ft ** 1.3) / (rho * Lv)
    bflx = grav / thlm * wpthlp
    ustar = diag_ustar(z, bflx, ubar, 0.035)
    return wpthlp, wprtp, ustar


def test_vs_literal():
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(200):
        t = rng.uniform(0.0, 12.0 * 3600.0)
        z, rho, thlm, ubar = rng.uniform(10, 40), rng.uniform(1.0, 1.2), rng.uniform(295, 305), rng.uniform(1, 8)
        wth, wrt, ust = lba_sfclyr(1, t, 0.0, z, rho, thlm, ubar)
        rth, rrt, rust = _ref(t, z, rho, thlm, ubar)
        worst = max(worst, abs(float(wth) - rth), abs(float(wrt) - rrt), abs(float(ust) - rust))
    assert worst < 1e-12, f"lba_sfclyr vs literal {worst:.2e}"
    print(f"  lba_sfclyr vs literal NumPy (200 cases): max diff {worst:.1e}  PASS")


def test_diurnal_structure():
    # Peak at t = 5.25 h (ft = cos(0) = 1).
    diurnal_factor = lambda t: max(0.0, math.cos(0.5 * math.pi * (5.25 - t / sec_per_hr) / 5.25))
    assert abs(diurnal_factor(5.25 * 3600.0) - 1.0) < 1e-12, "peak not at 5.25 h"
    # ft = 0 before t such that 0.5pi(5.25-t_hr)/5.25 >= pi/2 -> t_hr <= 0 ; and after the symmetric cutoff.
    assert diurnal_factor(0.0) >= 0.0 and diurnal_factor(0.0) < 1e-12, "ft(0) should be ~0"
    assert diurnal_factor(12.0 * 3600.0) == 0.0, "ft should be 0 well after the peak window"
    # Fluxes non-negative; ustar positive at the peak.
    wth, wrt, ust = lba_sfclyr(1, 5.25 * 3600.0, 0.0, 25.0, 1.1, 300.0, 4.0)
    assert float(wth) > 0 and float(wrt) > 0 and float(ust) > 0, "peak fluxes/ustar should be positive"
    print(f"  lba_sfclyr diurnal: peak@5.25h, ft(0)≈0, zero outside window, peak fluxes/ustar>0  PASS")


def test_differentiable():
    # Grad w.r.t. the physical inputs (rho, thlm, ubar) at the diurnal peak (ft constant → fractional power safe).
    g = jax.grad(lambda r: lba_sfclyr(1, 5.25 * 3600.0, 0.0, 25.0, r, 300.0, 4.0)[0])(jnp.asarray(1.1))
    gt = jax.grad(lambda T: lba_sfclyr(1, 5.25 * 3600.0, 0.0, 25.0, 1.1, T, 4.0)[2])(jnp.asarray(300.0))
    assert np.isfinite(float(g)) and np.isfinite(float(gt)), "non-finite grad"
    print(f"  jax.grad(lba_sfclyr): d wpthlp/d rho = {float(g):+.3e}, d ustar/d thlm = {float(gt):+.3e}: finite  PASS")


def main():
    print("test_lba_sfclyr:")
    for t in (test_vs_literal, test_diurnal_structure, test_differentiable):
        t()
    print("All lba_sfclyr checks PASSED")


if __name__ == "__main__":
    main()
