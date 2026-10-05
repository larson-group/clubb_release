#!/usr/bin/env python3
"""validate the JAX calc_L_x_Skx_fnc port (new_tsdadg_pdf.F90)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_tsdadg_pdf import calc_L_x_Skx_fnc, calc_setter_parameters, calc_respnder_parameters

NG, NZ = 2, 6


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


def test_responder():
    """Validate exact mean reconstruction (mu_x_2 follows the overall-mean constraint)
    and compare the variance coefficient with the setter's
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
