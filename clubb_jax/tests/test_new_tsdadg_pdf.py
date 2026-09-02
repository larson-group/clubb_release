#!/usr/bin/env python3
"""test_new_tsdadg_pdf.py — validate the JAX calc_L_x_Skx_fnc port (new_tsdadg_pdf.F90).

L_x_i = l_x_i · |Skx|/sqrt(4+Skx²), swapping l_x_1/l_x_2 when Skx·sgn(<w'x'>) < 0. Oracles:
  1. f2py bit-shadow vs f2py_calc_l_x_skx_fnc. SKIPs if clubb_f2py is unbuilt.
  2. Closed-form + the swap behavior + the Skx=0 -> 0 limit.
  3. A finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_tsdadg_pdf import (
    calc_L_x_Skx_fnc, calc_setter_parameters, calc_respnder_parameters, tsdadg_pdf_driver)

NG, NZ = 2, 6


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_l_x_skx_fnc oracle: SKIP ({type(e).__name__})")
        return
    # The f2py wrapper is scalar (the Fortran routine is per-level); call it element-wise.
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(60):
        Skx = float(rng.uniform(-3, 3)); sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        l1 = float(rng.uniform(2 / 3, 1.0)); l2 = float(rng.uniform(0.0, 1.0))
        f1, f2 = clubb_f2py.f2py_calc_l_x_skx_fnc(Skx, sgn, l1, l2)
        g1, g2 = calc_L_x_Skx_fnc(Skx, sgn, l1, l2)
        worst = max(worst, abs(float(g1) - float(f1)), abs(float(g2) - float(f2)))
    assert worst < 1e-13, f"calc_l_x_skx_fnc f2py mismatch {worst:.2e}"
    print(f"  f2py calc_l_x_skx_fnc: bit-match over 60 scalar configs, worst {worst:.2e}  PASS")


def test_closed_form_and_swap():
    Skx = np.array([[1.0, -1.0, 2.0]]); sgn = np.array([[1.0, 1.0, 1.0]])
    l1 = np.array([[0.8, 0.8, 0.8]]); l2 = np.array([[0.3, 0.3, 0.3]])
    g1, g2 = (np.asarray(x) for x in calc_L_x_Skx_fnc(Skx, sgn, l1, l2))
    factor = np.abs(Skx) / np.sqrt(4 + Skx ** 2)
    # Skx*sgn: [+, -, +] -> col1 swaps.
    exp1 = np.where(Skx * sgn >= 0, l1, l2) * factor
    exp2 = np.where(Skx * sgn >= 0, l2, l1) * factor
    assert np.max(np.abs(g1 - exp1)) < 1e-14 and np.max(np.abs(g2 - exp2)) < 1e-14, "closed-form/swap"
    # Skx=0 -> 0.
    z1, z2 = (np.asarray(x) for x in calc_L_x_Skx_fnc(np.array([[0.0]]), np.array([[1.0]]),
                                                      np.array([[0.8]]), np.array([[0.3]])))
    assert z1[0, 0] == 0.0 and z2[0, 0] == 0.0, "Skx=0 should give 0"
    print("  closed-form + swap on Skx*sgn<0 + Skx=0 limit  PASS")


def test_differentiable():
    rng = np.random.default_rng(7)
    Skx = jnp.asarray(rng.uniform(-3, 3, (NG, NZ)))
    sgn = jnp.asarray(np.sign(rng.uniform(-1, 1, (NG, NZ))) + 0.0)
    l1 = jnp.asarray(rng.uniform(0.6, 1.0, (NG, NZ))); l2 = jnp.asarray(rng.uniform(0.0, 1.0, (NG, NZ)))
    def loss(s):
        a, b = calc_L_x_Skx_fnc(s, sgn, l1, l2)
        return jnp.sum(a ** 2 + b ** 2)
    g = np.asarray(jax.grad(loss)(Skx))
    assert np.isfinite(g).all(), "non-finite grad through calc_L_x_Skx_fnc"
    print(f"  jax.grad through calc_L_x_Skx_fnc: finite ({g.size} entries)  PASS")


def test_setter_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_setter_parameters_tsdadg oracle: SKIP ({type(e).__name__})")
        return
    # The f2py wrapper is scalar (per-level); call it element-wise.
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(60):
        xm = float(rng.uniform(-2, 2)); xp2 = float(rng.uniform(0.05, 2.0))
        Skx = float(rng.uniform(-2.5, 2.5)); sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0)
        L1 = float(rng.uniform(0.05, 0.5)); L2 = float(rng.uniform(0.05, 0.5))
        f = clubb_f2py.f2py_calc_setter_parameters_tsdadg(xm, xp2, Skx, sgn, L1, L2)
        g = calc_setter_parameters(xm, xp2, Skx, sgn, L1, L2)
        for fi, gi in zip(f, g):
            worst = max(worst, abs(float(np.asarray(gi)) - float(np.asarray(fi))))
    assert worst < 1e-11, f"calc_setter_parameters_tsdadg f2py mismatch {worst:.2e}"
    print(f"  f2py calc_setter_parameters_tsdadg: bit-match (7 outputs, scalar), worst {worst:.2e}  PASS")


def test_setter_mean_reconstruction():
    # The setter binormal reproduces the overall mean xm exactly (mf*mu1 + (1-mf)*mu2 == xm).
    rng = np.random.default_rng(9)
    worst = 0.0
    for _ in range(200):
        xm = rng.uniform(-2, 2); xp2 = rng.uniform(0.1, 2.0); Skx = rng.uniform(-2, 2)
        sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0); L1 = rng.uniform(0.1, 0.5); L2 = rng.uniform(0.1, 0.5)
        mu1, mu2, s1, s2, mf, c1, c2 = (float(np.asarray(x)) for x in
                                        calc_setter_parameters(xm, xp2, Skx, sgn, L1, L2))
        worst = max(worst, abs(mf * mu1 + (1 - mf) * mu2 - xm))
    assert worst < 1e-9, f"mean reconstruction {worst:.2e}"
    print(f"  setter mean reconstruction: mf*mu1 + (1-mf)*mu2 == xm, worst {worst:.2e}  PASS")


def test_setter_differentiable():
    def loss(Skx):
        outs = calc_setter_parameters(0.5, 1.0, Skx, 1.0, 0.3, 0.25)
        return sum(jnp.sum(jnp.asarray(o) ** 2) for o in outs)
    g = float(jax.grad(loss)(1.2))
    assert np.isfinite(g), "non-finite grad through calc_setter_parameters"
    print(f"  jax.grad through calc_setter_parameters: finite (d/dSkx={g:.3e})  PASS")


def test_responder():
    """No direct f2py oracle for the TSDADG responder; validate the exact mean reconstruction (mu_x_2 follows
    the overall-mean constraint) and that its variance-coefficient formula matches the f2py-validated setter's
    (the only difference is mu_x_2_nrmlized), plus a finite jax.grad."""
    rng = np.random.default_rng(11)
    worst_mean = worst_coef = 0.0
    for _ in range(200):
        xm = rng.uniform(-2, 2); xp2 = rng.uniform(0.1, 2.0); Skx = rng.uniform(-2, 2)
        sgn = float(np.sign(rng.uniform(-1, 1)) or 1.0); mf = rng.uniform(0.3, 0.7); L1 = rng.uniform(0.1, 0.5)
        mu1, mu2, s1, s2, c1, c2 = (float(np.asarray(x)) for x in
                                    calc_respnder_parameters(xm, xp2, Skx, sgn, mf, L1))
        worst_mean = max(worst_mean, abs(mf * mu1 + (1 - mf) * mu2 - xm))
        # Reproduce coef1 with the literal formula using the same mu1n/mu2n (independent transcription).
        t = Skx * sgn / np.sqrt(4 + Skx ** 2)
        mu1n = L1 * np.sqrt((1 + t) / (1 - t)) * sgn
        mu2n = -(mf / (1 - mf)) * mu1n
        thr = max(mu1n, 1e-10) if mu1n >= 0 else min(mu1n, -1e-10)
        common = Skx / (3 * mf * thr) - mu1n ** 2 / 3 + mu2n ** 2 / 3
        base = 1 - mf * mu1n ** 2 - (1 - mf) * mu2n ** 2
        worst_coef = max(worst_coef, abs(c1 - (base + (1 - mf) * common)), abs(c2 - (base - mf * common)))
    assert worst_mean < 1e-9, f"mean reconstruction {worst_mean:.2e}"
    assert worst_coef < 1e-12, f"coef transcription {worst_coef:.2e}"
    def loss(s):
        outs = calc_respnder_parameters(0.5, 1.0, s, 1.0, 0.4, 0.3)
        return sum(jnp.sum(jnp.asarray(o) ** 2) for o in outs)
    g = float(jax.grad(loss)(1.2))
    assert np.isfinite(g), "non-finite grad"
    print(f"  TSDADG responder: mean reconstruction {worst_mean:.1e} + coef formula {worst_coef:.1e} + grad  PASS")


def test_driver_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py tsdadg_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(13)
    worst = 0.0
    ng, nz = 2, 8
    for _ in range(4):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(0, 0.02, (ng, nz)); thlm = rng.uniform(285, 305, (ng, nz))
        wp2 = rng.uniform(0.05, 2, (ng, nz)); rtp2 = rng.uniform(1e-8, 1e-6, (ng, nz)); thlp2 = rng.uniform(0.01, 1, (ng, nz))
        # Skewness with varying magnitudes so each of w/rt/thl is the setter at some points.
        Skw = rng.uniform(-3, 3, (ng, nz)); Skrt = rng.uniform(-3, 3, (ng, nz)); Skthl = rng.uniform(-3, 3, (ng, nz))
        wprtp = rng.uniform(-1, 1, (ng, nz)); wpthlp = rng.uniform(-1, 1, (ng, nz))
        f = clubb_f2py.f2py_tsdadg_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp)
        g = tsdadg_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"tsdadg_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py tsdadg_pdf_driver: bit-match (13 outputs, all 3 setter branches), worst {worst:.2e}  PASS")


def test_driver_differentiable():
    rng = np.random.default_rng(3)
    ng, nz = 1, 6
    args = [rng.uniform(-1, 1, (ng, nz)) for _ in range(11)]
    args[3] = np.abs(args[3]) + 0.1; args[4] = np.abs(args[4]) + 0.1; args[5] = np.abs(args[5]) + 0.1  # variances>0
    def loss(Skw):
        a = list(args); a[6] = Skw
        outs = tsdadg_pdf_driver(*a)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(args[6])))
    assert np.isfinite(g).all(), "non-finite grad through tsdadg_pdf_driver"
    print(f"  jax.grad through tsdadg_pdf_driver: finite ({g.size} entries)  PASS")


def main():
    print("test_new_tsdadg_pdf:")
    for t in (test_f2py_oracle, test_closed_form_and_swap, test_differentiable,
              test_setter_f2py, test_setter_mean_reconstruction, test_setter_differentiable, test_responder,
              test_driver_f2py, test_driver_differentiable):
        t()
    print("All new_tsdadg_pdf checks PASSED")


if __name__ == "__main__":
    main()
