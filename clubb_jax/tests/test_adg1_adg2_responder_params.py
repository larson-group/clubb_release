#!/usr/bin/env python3
"""test_adg1_adg2_responder_params.py — validate the JAX ADG1_ADG2_responder_params port.

`ADG1_ADG2_responder_params` (adg1_adg2_3d_luhar_pdf.py:73 ↔ adg1_adg2_3d_luhar_pdf.F90:1069) computes the
PDF component means (x_1, x_2), component variances (varnce_x_1, varnce_x_2), and the normalized-variance
factor alpha_x for a "responder" variable (rt, thl, sclr) given the ADG1/ADG2 w-closure. It is the only
routine in that module WITHOUT a direct behavioral test — its siblings (ADG1/ADG2/Luhar drivers, calc/close
Luhar, max_cubic_root, backsolve) each have one. It is exercised indirectly through the wired ADG1 path
(full-case gates), but never unit-pinned.

There is no f2py wrapper for it (the f2py `*responder_params*` wrappers are the NEW-PDF responders, a different
signature), so this is an independent per-level transcription of the Fortran do-loop (F90:1204-1220):

    x_1        = xm - wpxp / (sqrt_wp2 * w_2_n)
    x_2        = xm - wpxp / (sqrt_wp2 * w_1_n)
    alpha_x    = max( min( 0.5*(1 - wpxp^2 / ((1-sigma_sqd_w)*wp2*xp2)), 1 ), zero_threshold )
    width_factor_1 = (2/3)*beta + 2*mixt_frac*(1 - (2/3)*beta)
    varnce_x_1 = width_factor_1     * xp2 * alpha_x / mixt_frac
    varnce_x_2 = (2 - width_factor_1)* xp2 * alpha_x / (1 - mixt_frac)

Oracle-independent (pure numpy reference), never SKIPs; plus a finite jax.grad check (mirror the sibling tests).
(iter 524)
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


def test_responder_params_grad_finite():
    rng = np.random.default_rng(99)
    args = _inputs(rng)
    jargs = [jnp.asarray(a) for a in args]

    def loss(wpxp):
        x_1, x_2, vx_1, vx_2, alpha = ADG1_ADG2_responder_params(
            jargs[0], jargs[1], jargs[2], jargs[3], wpxp, jargs[5], jargs[6],
            jargs[7], jargs[8], jargs[9])
        return jnp.sum(x_1 ** 2 + x_2 ** 2 + vx_1 + vx_2 + alpha)

    g = jax.grad(loss)(jargs[4])
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad wrt wpxp"
    print("  jax.grad(ADG1_ADG2_responder_params) wrt wpxp is finite  PASS")


def main():
    print("test_adg1_adg2_responder_params:")
    test_responder_params_matches_reference()
    test_alpha_x_clip_to_unit_interval()
    test_responder_params_grad_finite()
    print("All ADG1_ADG2_responder_params checks PASSED")


if __name__ == "__main__":
    main()
