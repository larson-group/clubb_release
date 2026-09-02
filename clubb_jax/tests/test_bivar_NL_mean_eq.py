#!/usr/bin/env python3
"""test_bivar_NL_mean_eq.py — pin the 4-way variance DISPATCH of the KK bivariate normal-lognormal mean integral.

`bivar_NL_mean_eq` (KK_upscaled_means.py ↔ KK_upscaled_means.F90:864) computes the per-component moment
<chi^alpha · y^beta>_i (y = N_cn for autoconv / r_r for accretion) by selecting among four PDF-integral primitives
according to whether chi and/or y have (near-)zero variance — mirroring the Fortran if/elseif/elseif/else ladder:

    x1_const = sigma_chi ≤ CHI_TOL  OR  |s_c| > 49 (parabolic-cylinder overflow guard)
    x2_const = sigma_y   ≤ y_tol
      x1_const & x2_const -> bivar_NL_mean_const_all
      x1_const            -> bivar_NL_mean_const_x1
      x2_const            -> bivar_NL_mean_const_x2
      else                -> bivar_NL_mean   (both vary)

The four primitives are themselves validated by test_pdf_integrals_all_mm.py, and the whole KK rate path is bit-exact
vs the rico oracle (test_kk_rico_oracle.py) — but the DISPATCH wiring (which primitive + with which arg subset for
each regime) was never pinned in isolation. This drives each of the 4 regimes and asserts `bivar_NL_mean_eq` returns
EXACTLY the corresponding primitive called with the documented args (so a mis-wired branch / wrong-arg bug is caught
directly). Oracle-independent; never SKIPs. (iter 549)
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

from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import bivar_NL_mean_eq, CHI_TOL
from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_means import (
    bivar_NL_mean, bivar_NL_mean_const_x1, bivar_NL_mean_const_x2, bivar_NL_mean_const_all)

_Y_TOL = 1.0e-7
_ALPHA, _BETA = 2.39, 0.68     # KK-autoconversion-like exponents, both ≥ 0 (so mu_x2 = mu_y)


def _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho):
    return float(bivar_NL_mean_eq(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho, _Y_TOL, _ALPHA, _BETA))


def test_both_vary_selects_general():
    # sigma_chi > CHI_TOL, sigma_y > y_tol, and |s_c| ≤ 49 (moderate mu_chi/sigma_chi keeps s_c small).
    mu_chi, mu_y, mu_y_n = 0.5, 1.0e-3, -7.0
    sig_chi, sig_y, sig_y_n, rho = 0.3, 2.0e-4, 0.5, 0.2
    got = _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho)
    sig1_safe = max(sig_chi, CHI_TOL)
    ref = float(bivar_NL_mean(mu_chi, mu_y_n, sig1_safe, sig_y_n, rho, _ALPHA, _BETA))
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"both-vary != bivar_NL_mean: {got} vs {ref}"
    print("  both vary -> bivar_NL_mean (general)  PASS")


def test_x1_const_selects_const_x1():
    # sigma_chi ≤ CHI_TOL → x1_const; sigma_y stays > y_tol.
    mu_chi, mu_y, mu_y_n = 0.5, 1.0e-3, -7.0
    sig_chi, sig_y, sig_y_n, rho = 0.0, 2.0e-4, 0.5, 0.2
    got = _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho)
    ref = float(bivar_NL_mean_const_x1(mu_chi, mu_y_n, sig_y_n, _ALPHA, _BETA))
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"x1_const != const_x1: {got} vs {ref}"
    print("  sigma_chi≤CHI_TOL -> bivar_NL_mean_const_x1  PASS")


def test_x2_const_selects_const_x2():
    # sigma_y ≤ y_tol → x2_const; sigma_chi stays > CHI_TOL (and s_c small).
    mu_chi, mu_y, mu_y_n = 0.5, 1.0e-3, -7.0
    sig_chi, sig_y, sig_y_n, rho = 0.3, 0.0, 0.5, 0.2
    got = _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho)
    sig1_safe = max(sig_chi, CHI_TOL)
    ref = float(bivar_NL_mean_const_x2(mu_chi, mu_y, sig1_safe, _ALPHA, _BETA))
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"x2_const != const_x2: {got} vs {ref}"
    print("  sigma_y≤y_tol -> bivar_NL_mean_const_x2  PASS")


def test_both_const_selects_const_all():
    mu_chi, mu_y, mu_y_n = 0.5, 1.0e-3, -7.0
    sig_chi, sig_y, sig_y_n, rho = 0.0, 0.0, 0.5, 0.2
    got = _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho)
    ref = float(bivar_NL_mean_const_all(mu_chi, mu_y, _ALPHA, _BETA))
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"both-const != const_all: {got} vs {ref}"
    print("  both σ≤tol -> bivar_NL_mean_const_all  PASS")


def test_parab_cyl_overflow_forces_const_x1():
    # |s_c| > 49 (huge mu_chi/sigma_chi) must force the const_x1 branch even though sigma_chi > CHI_TOL.
    mu_chi, mu_y, mu_y_n = 100.0, 1.0e-3, -7.0
    sig_chi, sig_y, sig_y_n, rho = 1.0e-3, 2.0e-4, 0.5, 0.2   # mu_chi/sig_chi = 1e5 ≫ 49
    got = _f(mu_chi, mu_y, mu_y_n, sig_chi, sig_y, sig_y_n, rho)
    ref = float(bivar_NL_mean_const_x1(mu_chi, mu_y_n, sig_y_n, _ALPHA, _BETA))
    assert abs(got - ref) <= 1e-12 * (abs(ref) + 1e-300), f"|s_c|>49 didn't force const_x1: {got} vs {ref}"
    print("  |s_c|>49 (parab-cyl overflow guard) -> bivar_NL_mean_const_x1  PASS")


def main():
    print("test_bivar_NL_mean_eq:")
    test_both_vary_selects_general()
    test_x1_const_selects_const_x1()
    test_x2_const_selects_const_x2()
    test_both_const_selects_const_all()
    test_parab_cyl_overflow_forces_const_x1()
    print("All bivar_NL_mean_eq dispatch checks PASSED")


if __name__ == "__main__":
    main()
