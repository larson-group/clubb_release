#!/usr/bin/env python3
"""test_spec_hum_to_mixing_ratio.py — validate Benchmark_cases/spec_hum_to_mixing_ratio.py.

Mirrors spec_hum_to_mixing_ratio.F90 — converts a surface flux / large-scale forcing expressed in specific
humidity q_t into total-water mixing ratio r_t via the linearised Jacobian dr_t/dq_t = (1 + r_t)^2. Used by the
q_t-forced cases (bomex, gabls2, …) in prescribe_forcings.py (`use spec_hum_to_mixing_ratio`); case-active but
previously untested. Pinned vs the literal F90 formulas:

  flux_spec_hum_to_mixing_ratio  = (1 + rtm_sfc)^2 · w'q_t'     (F90:64)
  force_spec_hum_to_mixing_ratio = (1 + rtm)^2     · qtm_forcing (F90:109)

Pure-Python (the F90 formulas are the oracle), so it never SKIPs. (iter 509)
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

from clubb_jax.src.Benchmark_cases.spec_hum_to_mixing_ratio import (
    flux_spec_hum_to_mixing_ratio, force_spec_hum_to_mixing_ratio,
)


def test_conversions():
    rng = np.random.default_rng(59)
    rtm = rng.uniform(1e-4, 2.5e-2, (3, 20))     # total water mixing ratio
    wpqtp = rng.uniform(-1e-4, 1e-3, (3, 20))
    qtm_forcing = rng.uniform(-1e-7, 1e-7, (3, 20))

    jf = np.asarray(flux_spec_hum_to_mixing_ratio(3, rtm, wpqtp))
    jg = np.asarray(force_spec_hum_to_mixing_ratio(3, 20, rtm, qtm_forcing))
    assert float(np.max(np.abs(jf - (1.0 + rtm) ** 2 * wpqtp))) == 0.0, "flux conversion != F90"
    assert float(np.max(np.abs(jg - (1.0 + rtm) ** 2 * qtm_forcing))) == 0.0, "force conversion != F90"
    # physical: the factor reduces to identity as r_t -> 0 (factor (1+0)^2 = 1)
    identity = np.asarray(flux_spec_hum_to_mixing_ratio(3, np.zeros_like(wpqtp), wpqtp))
    assert float(np.max(np.abs(identity - wpqtp))) == 0.0, "r_t=0 must give identity"
    print("  flux_/force_spec_hum_to_mixing_ratio == F90 ((1+r_t)^2 Jacobian, exact)  PASS")


def main():
    print("test_spec_hum_to_mixing_ratio:")
    test_conversions()
    print("All spec_hum_to_mixing_ratio checks PASSED")


if __name__ == "__main__":
    main()
