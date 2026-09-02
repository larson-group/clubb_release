#!/usr/bin/env python3
"""test_calc_liquid_cloud_frac_component.py — validate the per-component Gaussian liquid-cloud closure.

`calc_liquid_cloud_frac_component(mean_chi, stdev_chi)` (pdf_closure_module.py ↔ pdf_closure_module.F90:2453) returns
the liquid cloud fraction and liquid water of ONE PDF component from the Gaussian distribution of χ (extended liquid
water): for χ ~ N(mean_chi, stdev²),
    cloud_frac = P(χ > 0) = ½(1 + erf(mean_chi/(stdev·√2)))      (the Gaussian CDF at 0)
    rc         = E[max(χ, 0)] = mean_chi·cloud_frac + stdev·φ(mean_chi/stdev)   (φ = std-normal pdf)
with ±max_num_stdevs(=5)·stdev truncation to the clear (cf=0, rc=0) / fully-cloudy (cf=1, rc=mean_chi) limits. It
feeds the ADG1 cloud-fraction closure but was validated only end-to-end. This pins both outputs against (1) the
closed-form probability/conditional-expectation, and (2) an INDEPENDENT Monte-Carlo (draw χ, measure the cloudy
fraction and mean positive part) — so the formulas are validated, not just transcribed — plus the truncation limits and
a finite grad. (iter 562)
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_liquid_cloud_frac_component


def test_closed_form_and_monte_carlo():
    rng = np.random.default_rng(562)
    N = 12_000_000
    worst_cf = 0.0; worst_rc = 0.0
    # interior cases (within ~±3σ so the ±5σ truncation doesn't fire)
    for mean_chi, stdev in [(2.0e-4, 1.0e-3), (-5.0e-4, 8.0e-4), (0.0, 1.0e-3), (1.5e-3, 1.0e-3)]:
        cf, rc = calc_liquid_cloud_frac_component(mean_chi, stdev)
        cf, rc = float(cf), float(rc)
        z = mean_chi / stdev
        cf_cf = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        rc_cf = mean_chi * cf_cf + stdev * math.exp(-0.5 * z ** 2) / math.sqrt(2.0 * math.pi)
        assert abs(cf - cf_cf) < 1e-13 and abs(rc - rc_cf) < 1e-15, "closed-form mismatch"
        chi = rng.normal(mean_chi, stdev, N)
        cf_mc = np.mean(chi > 0.0)
        rc_mc = np.mean(np.maximum(chi, 0.0))
        worst_cf = max(worst_cf, abs(cf - cf_mc))
        worst_rc = max(worst_rc, abs(rc - rc_mc) / (abs(rc_mc) + 1e-300))
    assert worst_cf < 5e-4 and worst_rc < 5e-3, f"MC: cf {worst_cf:.2e}, rc {worst_rc:.2e}"
    print(f"  cloud_frac=P(χ>0), rc=E[max(χ,0)]: closed form + {N//1_000_000}M MC (cf Δ{worst_cf:.1e}, rc rel {worst_rc:.1e})  PASS")


def test_truncation_limits():
    s = 1.0e-3
    # mean_chi > 5σ -> fully cloudy: cf=1, rc=mean_chi
    cf, rc = calc_liquid_cloud_frac_component(6.0 * s, s)
    assert float(cf) == 1.0 and float(rc) == 6.0 * s, "fully-cloudy limit"
    # mean_chi < -5σ -> clear: cf=0, rc=0
    cf, rc = calc_liquid_cloud_frac_component(-6.0 * s, s)
    assert float(cf) == 0.0 and float(rc) == 0.0, "clear limit"
    # zero mean, zero stdev -> clear (the |chi|<=eps & stdev<=chi_tol branch)
    cf, rc = calc_liquid_cloud_frac_component(0.0, 0.0)
    assert float(cf) == 0.0 and float(rc) == 0.0, "deterministic clear"
    print("  truncation: >5σ -> (1, mean_chi); <−5σ -> (0,0); deterministic dry -> (0,0)  PASS")


def test_grad_finite():
    g = jax.grad(lambda mc: calc_liquid_cloud_frac_component(mc, 1.0e-3)[1])(2.0e-4)
    assert np.isfinite(float(g)), "non-finite grad of rc wrt mean_chi"
    print("  jax.grad(rc) wrt mean_chi finite  PASS")


def main():
    print("test_calc_liquid_cloud_frac_component:")
    test_closed_form_and_monte_carlo()
    test_truncation_limits()
    test_grad_finite()
    print("All calc_liquid_cloud_frac_component checks PASSED")


if __name__ == "__main__":
    main()
