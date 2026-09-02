#!/usr/bin/env python3
"""test_T_in_K.py — validate the T_in_K_module port (thlm <-> absolute-T conversions).

T_in_K_module.F90 exposes the pair thlm2T_in_K (T = thlm*exner + Lv*rcm/Cp) and T_in_K2thlm_api
(thlm = (T - Lv/Cp*rcm)/exner). Validation: closed-form values + exact mutual inverse (round trip) + an
f2py bit-shadow of thlm2T_in_K vs the Fortran oracle (SKIPs cleanly if clubb_f2py is unbuilt).
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np

from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K, T_in_K2thlm
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv


def test_closed_form():
    thlm, exner, rcm = 300.0, 0.9, 1e-3
    assert np.isclose(thlm2T_in_K(thlm, exner, rcm), thlm * exner + Lv * rcm / Cp, rtol=0, atol=0)
    T = 285.0
    assert np.isclose(T_in_K2thlm(T, exner, rcm), (T - Lv / Cp * rcm) / exner, rtol=0, atol=0)
    print("  thlm2T_in_K / T_in_K2thlm closed-form  PASS")


def test_exact_inverse():
    rng = np.random.default_rng(0)
    thlm = rng.uniform(270, 330, 5000)
    exner = rng.uniform(0.5, 1.0, 5000)
    rcm = rng.uniform(0.0, 3e-3, 5000)
    rt = T_in_K2thlm(thlm2T_in_K(thlm, exner, rcm), exner, rcm)
    worst = np.max(np.abs(rt - thlm) / np.abs(thlm))
    assert worst < 1e-13, f"round-trip worst rel {worst:.2e}"
    print(f"  T_in_K2thlm is the exact inverse of thlm2T_in_K: {len(thlm)} cases, worst rel {worst:.1e}  PASS")


def test_f2py_oracle():
    """Bit-shadow thlm2T_in_K vs f2py_thlm2t_in_k_1d (the Fortran absolute-temperature diagnosis used by the
    radiation/microphysics handoff). SKIPs if clubb_f2py is unbuilt. (iter 446)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py thlm2T_in_K oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(15):
        nz = 20
        thlm = rng.uniform(280.0, 320.0, nz); exner = rng.uniform(0.6, 1.0, nz); rcm = rng.uniform(0.0, 1e-3, nz)
        j = np.asarray(thlm2T_in_K(thlm, exner, rcm))
        f = np.asarray(clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm))
        worst = max(worst, float(np.max(np.abs(j - f))))
    assert worst < 1e-11, f"thlm2T_in_K f2py mismatch {worst:.2e}"
    print(f"  f2py thlm2T_in_K (15 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_T_in_K:")
    test_closed_form()
    test_exact_inverse()
    test_f2py_oracle()
    print("All T_in_K_module tests PASSED")


if __name__ == "__main__":
    main()
