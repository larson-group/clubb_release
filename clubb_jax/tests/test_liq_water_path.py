#!/usr/bin/env python3
"""test_liq_water_path.py — validate Radiation/simple_rad_module.py:liq_water_path.

`liq_water_path` (simple_rad_module.F90:liq_water_path) is the cumulative-from-the-top liquid water path that the
simplified cloud-top-cooling radiation (`simple_rad`, used by the bit-faithful dycoms2_rf01/rf02 cases) integrates
the longwave flux against (Frad_LW = F0·exp(−κ·LWP) + F1·exp(−κ·(LWP_bottom − LWP))). It was validated only
implicitly by those cases. Pinned here vs an independent top-down transcription of the F90 loop:

    LWP(nzm) = 0;  LWP(k) = LWP(k+1) + rcm(k)·rho(k)/invrs_dzt(k)   for k = nzt..1   (F90:38-39)

Pure-Python (the F90 loop is the oracle), so it never SKIPs. (iter 504)
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

from clubb_jax.src.Radiation.simple_rad_module import liq_water_path


def test_liq_water_path():
    ng, nzt = 2, 24
    nzm = nzt + 1
    rng = np.random.default_rng(54)
    for _ in range(8):
        rho = rng.uniform(0.5, 1.2, (ng, nzt))
        # cloudy mid layer, clear elsewhere
        rcm = np.zeros((ng, nzt)); rcm[:, 6:16] = rng.uniform(1e-5, 1e-3, (ng, 10))
        invrs_dzt = rng.uniform(0.01, 0.05, (ng, nzt))
        j = np.asarray(liq_water_path(ng, nzm, nzt, rho, rcm, invrs_dzt))
        # independent top-down transcription
        contrib = rcm * rho / invrs_dzt
        ref = np.zeros((ng, nzm))
        for k in range(nzt - 1, -1, -1):
            ref[:, k] = ref[:, k + 1] + contrib[:, k]
        assert j.shape == ref.shape, f"shape {j.shape} != {ref.shape}"
        assert float(np.max(np.abs(j - ref))) < 1e-15, f"liq_water_path != F90 top-down loop ({np.max(np.abs(j-ref)):.2e})"
        # physical: monotonic non-increasing with height, top = 0
        assert float(np.max(np.abs(j[:, -1]))) == 0.0, "top boundary LWP must be 0"
        assert np.all(np.diff(j, axis=1) <= 1e-30), "LWP must be non-increasing upward"
    print("  liq_water_path == F90 top-down transcription (rev-cumsum, top=0, monotone)  PASS")


def main():
    print("test_liq_water_path:")
    test_liq_water_path()
    print("All liq_water_path checks PASSED")


if __name__ == "__main__":
    main()
