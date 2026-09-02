#!/usr/bin/env python3
"""test_calc_corr_w_hm_n.py — validate the JAX calc_corr_w_hm_n port.

calc_corr_w_hm_n (setup_clubb_pdf_params.F90:3428) diagnoses the PDF-component correlation of w and
ln(hm) in-precip from the overall w-hm flux <w'hm'>. It is the *inverse* of the flux assembly, which gives
a strong independent oracle:

  1. Round-trip (both-vary branch): pick a correlation corr in (-0.99, 0.99), build the overall flux
     wphydrometp from the forward formula
         <w'hm'> = Σ_i mixt_i precip_frac_i [ (mu_w_i - <w>) mu_hm_i + corr σ_w_i σ_hm_i_n mu_hm_i ],
     then calc_corr_w_hm_n must recover corr to machine precision.
  2. Literal transcription of the 4-way branch (both vary / comp1 only / comp2 only / neither) — checked
     against a direct NumPy re-implementation over randomized degenerate-sigma configs.
  3. The ±max_mag_correlation clamp fires for an out-of-range flux.
  4. A finite jax.grad through the diagnosed correlation.
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

from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import calc_corr_w_hm_n
from clubb_jax.src.CLUBB_core.constants_clubb import max_mag_correlation, w_tol

_HM_TOL = 1.0e-12


def _forward_flux(corr, wm, mu_w_1, mu_w_2, mu_hm_1, mu_hm_2,
                  sigma_w_1, sigma_w_2, sigma_hm_1_n, sigma_hm_2_n, mixt_frac, pf1, pf2):
    """Assemble the overall <w'hm'> flux from a single shared component correlation (both-vary case)."""
    w1, w2 = mixt_frac * pf1, (1.0 - mixt_frac) * pf2
    return (w1 * ((mu_w_1 - wm) * mu_hm_1 + corr * sigma_w_1 * sigma_hm_1_n * mu_hm_1)
            + w2 * ((mu_w_2 - wm) * mu_hm_2 + corr * sigma_w_2 * sigma_hm_2_n * mu_hm_2))


def test_round_trip():
    rng = np.random.default_rng(2026_0603)
    worst = 0.0
    for _ in range(200):
        corr = rng.uniform(-0.95, 0.95)
        wm = rng.uniform(-1.0, 1.0)
        mu_w_1, mu_w_2 = rng.uniform(-2.0, 2.0, 2)
        mu_hm_1, mu_hm_2 = rng.uniform(0.5, 5.0, 2)          # in-precip means are positive
        sigma_w_1, sigma_w_2 = rng.uniform(0.1, 1.0, 2)
        sigma_hm_1_n, sigma_hm_2_n = rng.uniform(0.1, 1.0, 2)
        mixt_frac = rng.uniform(0.2, 0.8)
        pf1, pf2 = rng.uniform(0.3, 1.0, 2)
        # sigma_hm (linear) only gates the branch; keep > tol so the both-vary branch is taken.
        sigma_hm_1, sigma_hm_2 = 1.0, 1.0
        wphm = _forward_flux(corr, wm, mu_w_1, mu_w_2, mu_hm_1, mu_hm_2,
                             sigma_w_1, sigma_w_2, sigma_hm_1_n, sigma_hm_2_n, mixt_frac, pf1, pf2)
        c1, c2 = calc_corr_w_hm_n(wm, wphm, mu_w_1, mu_w_2, mu_hm_1, mu_hm_2,
                                  sigma_w_1, sigma_w_2, sigma_hm_1, sigma_hm_2,
                                  sigma_hm_1_n, sigma_hm_2_n, mixt_frac, pf1, pf2, _HM_TOL)
        c1, c2 = float(c1), float(c2)
        assert abs(c1 - c2) < 1e-14, "both-vary: components must be equal"
        worst = max(worst, abs(c1 - corr))
    assert worst < 1e-12, f"round-trip recovery rel {worst:.2e}"
    print(f"  round-trip: recover corr from assembled <w'hm'> over 200 configs, worst {worst:.2e}  PASS")


def _np_ref(wm, wphm, mu_w_1, mu_w_2, mu_hm_1, mu_hm_2, sw1, sw2, shm1, shm2,
            shm1n, shm2n, a, pf1, pf2, hm_tol):
    """Literal NumPy transcription of the Fortran branch ladder."""
    w1, w2 = a * pf1, (1.0 - a) * pf2
    num = wphm - w1 * (mu_w_1 - wm) * mu_hm_1 - w2 * (mu_w_2 - wm) * mu_hm_2
    clamp = lambda x: max(-max_mag_correlation, min(max_mag_correlation, x))
    c1v = (sw1 > w_tol) and (shm1 > hm_tol)
    c2v = (sw2 > w_tol) and (shm2 > hm_tol)
    if c1v and c2v:
        cn = clamp(num / (w1 * sw1 * shm1n * mu_hm_1 + w2 * sw2 * shm2n * mu_hm_2))
        return cn, cn
    if c1v:
        return clamp(num / (w1 * sw1 * shm1n * mu_hm_1)), 0.0
    if c2v:
        return 0.0, clamp(num / (w2 * sw2 * shm2n * mu_hm_2))
    return 0.0, 0.0


def test_branches_vs_numpy():
    rng = np.random.default_rng(7)
    worst = 0.0
    seen = set()
    for _ in range(400):
        # Randomly switch each sigma between "varies" (>tol) and "constant" (0) to exercise all 4 branches.
        sw1 = rng.uniform(0.1, 1.0) if rng.random() < 0.6 else 0.0
        sw2 = rng.uniform(0.1, 1.0) if rng.random() < 0.6 else 0.0
        shm1 = 1.0 if rng.random() < 0.6 else 0.0
        shm2 = 1.0 if rng.random() < 0.6 else 0.0
        args = dict(wm=rng.uniform(-1, 1), wphm=rng.uniform(-2, 2),
                    mu_w_1=rng.uniform(-2, 2), mu_w_2=rng.uniform(-2, 2),
                    mu_hm_1=rng.uniform(0.5, 5), mu_hm_2=rng.uniform(0.5, 5),
                    sw1=sw1, sw2=sw2, shm1=shm1, shm2=shm2,
                    shm1n=rng.uniform(0.1, 1), shm2n=rng.uniform(0.1, 1),
                    a=rng.uniform(0.2, 0.8), pf1=rng.uniform(0.3, 1), pf2=rng.uniform(0.3, 1),
                    hm_tol=_HM_TOL)
        seen.add(((sw1 > w_tol) and (shm1 > _HM_TOL), (sw2 > w_tol) and (shm2 > _HM_TOL)))
        r1, r2 = _np_ref(**args)
        g1, g2 = calc_corr_w_hm_n(args['wm'], args['wphm'], args['mu_w_1'], args['mu_w_2'],
                                  args['mu_hm_1'], args['mu_hm_2'], args['sw1'], args['sw2'],
                                  args['shm1'], args['shm2'], args['shm1n'], args['shm2n'],
                                  args['a'], args['pf1'], args['pf2'], args['hm_tol'])
        worst = max(worst, abs(float(g1) - r1), abs(float(g2) - r2))
    assert worst < 1e-13, f"branch transcription mismatch {worst:.2e}"
    assert len(seen) == 4, f"did not exercise all 4 branches, only {seen}"
    print(f"  4-way branch ladder vs literal NumPy: all 4 branches hit, worst {worst:.2e}  PASS")


def test_clamp():
    # A huge flux -> correlation would exceed 1; must clamp to +max_mag_correlation (and -ve for -ve flux).
    base = dict(mu_w_1=0.0, mu_w_2=0.0, mu_hm_1=1.0, mu_hm_2=1.0,
                sigma_w_1=0.5, sigma_w_2=0.5, sigma_hm_1=1.0, sigma_hm_2=1.0,
                sigma_hm_1_n=0.5, sigma_hm_2_n=0.5, mixt_frac=0.5, precip_frac_1=1.0,
                precip_frac_2=1.0, hm_tol=_HM_TOL)
    cp1, cp2 = calc_corr_w_hm_n(0.0, 1.0e6, **base)
    cm1, cm2 = calc_corr_w_hm_n(0.0, -1.0e6, **base)
    assert abs(float(cp1) - max_mag_correlation) < 1e-14 and abs(float(cp2) - max_mag_correlation) < 1e-14
    assert abs(float(cm1) + max_mag_correlation) < 1e-14 and abs(float(cm2) + max_mag_correlation) < 1e-14
    print("  clamp: |corr| limited to max_mag_correlation for out-of-range flux  PASS")


def test_differentiable():
    base = dict(mu_w_1=0.3, mu_w_2=-0.2, mu_hm_1=2.0, mu_hm_2=1.5,
                sigma_w_1=0.5, sigma_w_2=0.4, sigma_hm_1=1.0, sigma_hm_2=1.0,
                sigma_hm_1_n=0.6, sigma_hm_2_n=0.5, mixt_frac=0.4, precip_frac_1=0.8,
                precip_frac_2=0.9, hm_tol=_HM_TOL)
    def loss(wphm):
        c1, c2 = calc_corr_w_hm_n(0.1, wphm, **base)
        return c1 ** 2 + c2 ** 2
    g = float(jax.grad(loss)(0.05))
    assert np.isfinite(g), "non-finite grad through calc_corr_w_hm_n"
    print(f"  jax.grad through calc_corr_w_hm_n: finite (d/dwphm = {g:.3e})  PASS")


def main():
    print("test_calc_corr_w_hm_n:")
    for t in (test_round_trip, test_branches_vs_numpy, test_clamp, test_differentiable):
        t()
    print("All calc_corr_w_hm_n checks PASSED")


if __name__ == "__main__":
    main()
