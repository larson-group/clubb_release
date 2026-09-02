#!/usr/bin/env python3
"""test_calc_ice_cloud_frac_component.py — validate the per-component ice-supersaturation fraction.

`calc_ice_cloud_frac_component(mean_chi, stdev_chi, crt, rsatl, tl, cf_liq, p_in_Pa)` (pdf_closure_module.F90:2490)
returns one PDF component's ice-supersaturation fraction. Below freezing it is the fraction of the Gaussian χ PDF
SUPERSATURATED w.r.t. ice — P(χ > χ_at_ice_sat), χ_at_ice_sat = crt·(rsat_ice − rsatl), rsat_ice = sat_mixrat_ice(p,tl)
— i.e. ½(1+erf((mean_chi − χ_at_ice_sat)/(σ√2))) with ±5σ truncation; ABOVE freezing it equals the liquid cloud
fraction `cf_liq` exactly. It feeds the ADG1 ice-supersat closure but was validated only end-to-end. This pins it
against (1) the closed-form shifted Gaussian CDF and (2) an INDEPENDENT Monte-Carlo (draw χ, measure the fraction
above χ_at_ice_sat), plus the above-freezing passthrough, the truncation limits, and a finite grad. (iter 563)
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

from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_ice_cloud_frac_component
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_ice
from clubb_jax.src.CLUBB_core.constants_clubb import T_freeze_K

_CRT, _RSATL, _P = 1.0, 2.0e-3, 7.0e4


def _thresh(tl):
    rsat_ice = float(sat_mixrat_ice(_P, tl))   # Flatau (the formula the JAX hardcodes)
    return _CRT * (rsat_ice - _RSATL)


def test_below_freezing_closed_form_and_mc():
    tl = 260.0   # below freezing
    chi_ice = _thresh(tl)
    rng = np.random.default_rng(563)
    N = 12_000_000
    worst_cf = 0.0
    for mean_chi, stdev in [(chi_ice + 3e-4, 1.0e-3), (chi_ice - 2e-4, 8e-4), (chi_ice, 1.0e-3)]:
        ssf = float(calc_ice_cloud_frac_component(mean_chi, stdev, _CRT, _RSATL, tl, 0.123, _P))
        z = (mean_chi - chi_ice) / stdev
        cf = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        assert abs(ssf - cf) < 1e-12, f"closed form: {ssf} vs {cf}"
        chi = rng.normal(mean_chi, stdev, N)
        mc = np.mean(chi > chi_ice)
        worst_cf = max(worst_cf, abs(ssf - mc))
    assert worst_cf < 5e-4, f"MC mismatch {worst_cf:.2e}"
    print(f"  below freezing: ssf = P(χ > χ_ice_sat) closed form + {N//1_000_000}M MC (Δ {worst_cf:.1e})  PASS")


def test_above_freezing_returns_cf_liq():
    # tl > T_freeze_K -> returns cf_liq EXACTLY, regardless of the chi inputs.
    for cf_liq in (0.0, 0.37, 1.0):
        got = float(calc_ice_cloud_frac_component(5.0e-4, 1.0e-3, _CRT, _RSATL, 290.0, cf_liq, _P))
        assert got == cf_liq, f"above freezing should passthrough cf_liq={cf_liq}, got {got}"
    print("  above freezing (tl>T_freeze): returns cf_liq exactly  PASS")


def test_truncation_limits():
    tl = 255.0; chi_ice = _thresh(tl); s = 1.0e-3
    assert float(calc_ice_cloud_frac_component(chi_ice + 6 * s, s, _CRT, _RSATL, tl, 0.0, _P)) == 1.0
    assert float(calc_ice_cloud_frac_component(chi_ice - 6 * s, s, _CRT, _RSATL, tl, 0.0, _P)) == 0.0
    print("  truncation: Δ>5σ -> 1, Δ<−5σ -> 0  PASS")


def test_grad_finite():
    g = jax.grad(lambda mc: calc_ice_cloud_frac_component(mc, 1.0e-3, _CRT, _RSATL, 260.0, 0.1, _P))(_thresh(260.0) + 1e-4)
    assert np.isfinite(float(g)), "non-finite grad wrt mean_chi"
    print("  jax.grad(ssf) wrt mean_chi finite  PASS")


def main():
    print("test_calc_ice_cloud_frac_component:")
    test_below_freezing_closed_form_and_mc()
    test_above_freezing_returns_cf_liq()
    test_truncation_limits()
    test_grad_finite()
    print("All calc_ice_cloud_frac_component checks PASSED")


if __name__ == "__main__":
    main()
