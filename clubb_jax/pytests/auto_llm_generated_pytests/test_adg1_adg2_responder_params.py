#!/usr/bin/env python3
"""validate the JAX ADG1_ADG2_responder_params port."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import ADG1_ADG2_responder_params, zero_threshold

NG, NZ = 2, 8


def _ref(xm, xp2, wp2, sqrt_wp2, wpxp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta):
    """Independent numpy per-(i,k) transcription of ADG1_ADG2_responder_params (F90:1184-1225)."""
    x_1 = np.empty((NG, NZ)); x_2 = np.empty((NG, NZ))
    vx_1 = np.empty((NG, NZ)); vx_2 = np.empty((NG, NZ)); alpha = np.empty((NG, NZ))
    for k in range(NZ):
        for i in range(NG):
            x_1[i, k] = xm[i, k] - wpxp[i, k] / (sqrt_wp2[i, k] * w_2_n[i, k])
            x_2[i, k] = xm[i, k] - wpxp[i, k] / (sqrt_wp2[i, k] * w_1_n[i, k])
            a = 0.5 * (1.0 - wpxp[i, k] * wpxp[i, k]
                       / ((1.0 - sigma_sqd_w[i, k]) * wp2[i, k] * xp2[i, k]))
            a = max(min(a, 1.0), zero_threshold)
            alpha[i, k] = a
            wf1 = (2.0 / 3.0) * beta[i] + 2.0 * mixt_frac[i, k] * (1.0 - (2.0 / 3.0) * beta[i])
            vx_1[i, k] = wf1 * xp2[i, k] * a / mixt_frac[i, k]
            vx_2[i, k] = (2.0 - wf1) * xp2[i, k] * a / (1.0 - mixt_frac[i, k])
    return x_1, x_2, vx_1, vx_2, alpha


def _inputs(rng):
    xm = rng.uniform(280.0, 300.0, (NG, NZ))           # thl-like responder mean
    xp2 = rng.uniform(0.1, 4.0, (NG, NZ))              # variance > x_tol^2
    wp2 = rng.uniform(0.05, 1.5, (NG, NZ))
    sqrt_wp2 = np.sqrt(wp2)
    wpxp = rng.uniform(-0.4, 0.4, (NG, NZ))
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))        # strictly interior (0,1)
    sigma_sqd_w = rng.uniform(0.1, 0.6, (NG, NZ))      # < 1
    # normalized w-component means: w_1_n>0, w_2_n<0 (the two ADG plumes), away from 0
    w_1_n = rng.uniform(0.4, 1.4, (NG, NZ))
    w_2_n = -rng.uniform(0.4, 1.4, (NG, NZ))
    beta = rng.uniform(0.5, 2.5, (NG,))
    return xm, xp2, wp2, sqrt_wp2, wpxp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta


def test_responder_params_matches_reference():
    rng = np.random.default_rng(20240524)
    args = _inputs(rng)
    # Source signature: (xm, xp2, wp2, sqrt_wp2, wpxp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta) — pass straight.
    x_1, x_2, vx_1, vx_2, alpha = ADG1_ADG2_responder_params(
        jnp.asarray(args[0]), jnp.asarray(args[1]), jnp.asarray(args[2]), jnp.asarray(args[3]),
        jnp.asarray(args[4]), jnp.asarray(args[5]), jnp.asarray(args[6]),
        jnp.asarray(args[7]), jnp.asarray(args[8]), jnp.asarray(args[9]))
    rx_1, rx_2, rvx_1, rvx_2, ralpha = _ref(*args)
    worst = 0.0
    for got, ref, nm in ((x_1, rx_1, "x_1"), (x_2, rx_2, "x_2"), (vx_1, rvx_1, "varnce_x_1"),
                         (vx_2, rvx_2, "varnce_x_2"), (alpha, ralpha, "alpha_x")):
        rel = float(np.max(np.abs(np.asarray(got) - ref) / (np.abs(ref) + 1e-30)))
        worst = max(worst, rel)
        assert rel < 1e-12, f"{nm} rel-mismatch {rel:.2e} vs the F90 transcription"
    print(f"  ADG1_ADG2_responder_params: 5 outputs match F90 transcription (worst rel {worst:.2e})  PASS")


def test_alpha_x_clip_to_unit_interval():
    """alpha_x = max(min(., 1), zero_threshold) — force the clip both ways and pin the saturated values."""
    rng = np.random.default_rng(7)
    args = list(_inputs(rng))
    # Drive alpha well above 1 (tiny |wpxp|) on column 0, and below 0 (large |wpxp|) on column 1.
    args[4] = args[4].copy()
    args[4][0, :] = 0.0                               # wpxp=0 -> alpha = 0.5 (interior, sanity)
    args[4][1, :] = np.sqrt((1.0 - args[8][1, :]) * args[2][1, :] * args[1][1, :]) * 1.5  # |wpxp|^2 huge -> alpha<0 -> clip 0
    _, _, _, _, alpha = ADG1_ADG2_responder_params(*[jnp.asarray(a) for a in args])
    alpha = np.asarray(alpha)
    assert np.all(alpha >= zero_threshold - 1e-15) and np.all(alpha <= 1.0 + 1e-15), "alpha_x escaped [0,1]"
    assert np.allclose(alpha[0, :], 0.5), "wpxp=0 must give alpha_x=0.5"
    assert np.allclose(alpha[1, :], zero_threshold), "large wpxp must clip alpha_x to zero_threshold"
    print("  alpha_x clipped to [zero_threshold, 1] (0.5 at wpxp=0, floor at large wpxp)  PASS")
