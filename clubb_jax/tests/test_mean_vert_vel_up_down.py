#!/usr/bin/env python3
"""test_mean_vert_vel_up_down.py — validate the PDF mean up/down vertical velocity (mono flux limiter).

`calc_mean_w_up_down_component(w_i, varnce_w_i, wm)` (mono_flux_limiter.F90:calc_mean_w_up_down_component) returns the
mean DOWNWARD and UPWARD vertical velocity of ONE Gaussian PDF component on zm — the Gaussian truncated means
    mwu = E[max(w,0)] = w_i·Φ(w_i/σ) + σ·φ(w_i/σ),   mwd = E[min(w,0)] = w_i − mwu   (w ~ N(w_i, σ²))
with three short-circuit branches (too_weak: |w_i|+3σ ≤ wm → (0,0); all_dn: w_i+3σ ≤ 0 → (w_i,0); all_up:
w_i−3σ ≥ 0 → (0,w_i)) and the domain boundaries (k=0, k=−1) zeroed. `mean_vert_vel_up_down` is the mixt_frac-weighted
two-component combine. These feed the monotonic flux limiter (live xm_wpxp path) but were validated only end-to-end.
This pins (1) the interior branch's truncated means vs an INDEPENDENT Monte-Carlo (mean max/min of Gaussian draws),
(2) mwu+mwd == w_i, (3) the three branches + boundary zeroing, (4) the mixt_frac combine. (`calc_mean_w_up_down_component`
returns numpy — a flux-limiter RANGE diagnostic, off the grad path — so no grad test.) (iter 569)
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
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.mono_flux_limiter import (
    calc_mean_w_up_down_component, mean_vert_vel_up_down)

_NZM = 6   # interior levels 1..4 (0 and 5 are zeroed boundaries)


def _col(val):
    a = np.full((1, _NZM), val)
    return a


def test_interior_truncated_means_vs_mc():
    wi, sig = 0.2, 0.5            # |wi| < 3σ (else branch), wm=0 so not too_weak
    var = sig ** 2
    mwd, mwu = calc_mean_w_up_down_component(_col(wi), _col(var), _col(0.0))
    mwd, mwu = float(mwd[0, 2]), float(mwu[0, 2])    # an interior level
    # closed form
    Phi = 0.5 * (1 + math.erf(wi / (sig * math.sqrt(2))))
    phi = math.exp(-0.5 * (wi / sig) ** 2) / math.sqrt(2 * math.pi)
    cf_up = wi * Phi + sig * phi
    cf_dn = wi - cf_up
    assert abs(mwu - cf_up) < 1e-12 and abs(mwd - cf_dn) < 1e-12, "closed-form truncated means"
    assert abs((mwu + mwd) - wi) < 1e-12, "mwu + mwd must equal w_i"
    # Monte-Carlo
    rng = np.random.default_rng(569)
    w = rng.normal(wi, sig, 14_000_000)
    mc_up, mc_dn = np.mean(np.maximum(w, 0.0)), np.mean(np.minimum(w, 0.0))
    assert abs(mwu - mc_up) < 5e-4 and abs(mwd - mc_dn) < 5e-4, f"MC: up {abs(mwu-mc_up):.2e}, dn {abs(mwd-mc_dn):.2e}"
    print(f"  interior: mwu=E[max(w,0)], mwd=E[min(w,0)] closed form + MC (Δ {max(abs(mwu-mc_up),abs(mwd-mc_dn)):.1e}); mwu+mwd=w_i  PASS")


def test_branches_and_boundary():
    sig = 0.3; var = sig ** 2
    # too_weak: |wi|+3σ <= wm -> (0,0)
    mwd, mwu = calc_mean_w_up_down_component(_col(0.1), _col(var), _col(5.0))
    assert mwd[0, 2] == 0.0 and mwu[0, 2] == 0.0, "too_weak -> (0,0)"
    # all_dn: wi+3σ <= 0 -> (wi, 0)
    mwd, mwu = calc_mean_w_up_down_component(_col(-2.0), _col(var), _col(0.0))
    assert mwd[0, 2] == -2.0 and mwu[0, 2] == 0.0, "all_dn -> (wi, 0)"
    # all_up: wi-3σ >= 0 -> (0, wi)
    mwd, mwu = calc_mean_w_up_down_component(_col(2.0), _col(var), _col(0.0))
    assert mwd[0, 2] == 0.0 and mwu[0, 2] == 2.0, "all_up -> (0, wi)"
    # boundaries zeroed
    mwd, mwu = calc_mean_w_up_down_component(_col(0.2), _col(var), _col(0.0))
    assert mwd[0, 0] == 0.0 and mwu[0, 0] == 0.0 and mwd[0, -1] == 0.0 and mwu[0, -1] == 0.0, "boundary zero"
    print("  branches (too_weak/all_dn/all_up) + boundary (k=0,−1) zeroing  PASS")


def test_combine():
    # NB calc_mean_w_up_down_component returns numpy (a flux-limiter RANGE diagnostic, off the grad path) — no grad test.
    var1, var2 = 0.25, 0.16
    mwd1, mwu1 = calc_mean_w_up_down_component(_col(0.2), _col(var1), _col(0.0))
    mwd2, mwu2 = calc_mean_w_up_down_component(_col(-0.3), _col(var2), _col(0.0))
    mf = 0.4
    cwd, cwu = mean_vert_vel_up_down(_col(0.2), _col(-0.3), _col(var1), _col(var2), _col(mf), _col(0.0))
    assert np.allclose(np.asarray(cwd), mf * mwd1 + (1 - mf) * mwd2) and \
           np.allclose(np.asarray(cwu), mf * mwu1 + (1 - mf) * mwu2), "mixt_frac combine"
    print("  mean_vert_vel_up_down = mixt_frac-weighted combine of the two components  PASS")


def main():
    print("test_mean_vert_vel_up_down:")
    test_interior_truncated_means_vs_mc()
    test_branches_and_boundary()
    test_combine()
    print("All mean_vert_vel_up_down checks PASSED")


if __name__ == "__main__":
    main()
