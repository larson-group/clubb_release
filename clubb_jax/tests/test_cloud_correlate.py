#!/usr/bin/env python3
"""test_cloud_correlate.py — validate the bugs_ctot generalized cloud-overlap total cloud amount.

Oracles (no f2py for BUGSrad): a LITERAL NumPy transcription of the Fortran nested-loop algorithm (the right
check for the cumulative-product refactor of the cld_below recurrence), the analytic random-overlap limit
(l_c<1m → c_tot = 1-prod(1-cloud), which is also the maximum possible c_tot), the physical bounds
max(cloud) <= c_tot <= random and monotonic decrease in l_c, the edge branches, and an FD-correct jax.grad of
the differentiable overlap core (the discrete cloudy-layer selection is not a differentiable black box).
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

from clubb_jax.src.Radiation.BUGSrad.cloud_correlate import (
    bugs_ctot, ctot_from_cloudy_layers, bugs_cloudfit, _midlayer_heights, MIN_CF, MIN_LC)

_R_DZ = 29.286


def _ref_ctot(pl2, tl, acld, l_c):
    """Literal NumPy transcription of cloud_correlate.F90:bugs_ctot (one column)."""
    nlm = acld.shape[0]
    Ncloud = int((acld > MIN_CF).sum())
    if Ncloud == 0:
        return 0.0
    if Ncloud == 1:
        return acld.max()
    if acld.max() >= 1.0:
        return 1.0
    z = np.zeros(nlm)
    z0 = 0.0
    for i in range(nlm - 1, -1, -1):
        dz = _R_DZ * np.log(pl2[i + 1] / pl2[i]) * tl[i]
        z[i] = z0 + dz / 2.0
        z0 += dz
    i_cld = np.nonzero(acld > MIN_CF)[0]
    cloud = acld[i_cld]
    olap = np.zeros(Ncloud - 1)
    for k in range(1, Ncloud):
        if l_c < MIN_LC:
            olap[k - 1] = cloud[k - 1] * cloud[k]
        else:
            dz = z[i_cld[k - 1]] - z[i_cld[k]]
            alpha = np.exp(-dz / l_c)
            olap[k - 1] = alpha * min(cloud[k - 1], cloud[k]) + (1 - alpha) * cloud[k - 1] * cloud[k]
    cld_below = np.zeros(Ncloud)
    cld_below[0] = cloud[0]
    cld_below[1] = cloud[1] - olap[0]
    for k in range(2, Ncloud):
        for j in range(k - 1, -1, -1):
            if k - j == 1:
                cld_below[k] = cloud[k] - olap[k - 1]
            else:
                cld_below[k] = cld_below[k] * (1.0 - (cloud[j] - olap[j]) / (1.0 - cloud[j + 1]))
    return cld_below.sum()


def _profile(seed, nlm=20):
    rng = np.random.default_rng(seed)
    pl2 = np.sort(rng.uniform(100.0, 1000.0, nlm + 1))[::-1] * 100.0   # decreasing-with-index? level pressures
    pl2 = np.sort(pl2)                                                  # ascending top->bottom (Pa)
    tl = rng.uniform(220.0, 290.0, nlm)
    acld = np.zeros(nlm)
    n_cld = rng.integers(3, nlm)
    cld_layers = rng.choice(nlm, size=n_cld, replace=False)
    acld[cld_layers] = rng.uniform(0.05, 0.95, n_cld)
    return pl2, tl, acld


def test_vs_literal_loop():
    worst = 0.0
    for seed in range(40):
        pl2, tl, acld = _profile(seed)
        for l_c in (0.5, 500.0, 5000.0):
            got = float(bugs_ctot(pl2, tl, acld, l_c))
            ref = _ref_ctot(pl2, tl, acld, l_c)
            worst = max(worst, abs(got - ref) / (abs(ref) + 1e-30))
    assert worst < 1e-12, f"bugs_ctot vs literal loop rel {worst:.2e}"
    print(f"  bugs_ctot vs literal NumPy loop (40 profiles × 3 l_c): rel {worst:.1e}  PASS")


def test_random_limit_bounds_monotonic():
    pl2, tl, acld = _profile(7)
    cl = acld[acld > MIN_CF]
    # Random overlap (l_c < MIN_LC) is the MAXIMUM possible c_tot = 1 - prod(1-cloud) (Fortran's own comment).
    rnd = float(bugs_ctot(pl2, tl, acld, 0.5))
    assert abs(rnd - (1.0 - np.prod(1.0 - cl))) < 1e-12, "random-overlap limit"
    # For any l_c, max(cloud) <= c_tot <= random, and c_tot decreases monotonically as l_c grows (more overlap).
    vals = [float(bugs_ctot(pl2, tl, acld, l)) for l in (0.5, 500.0, 2000.0, 5000.0, 1e9)]
    for v in vals:
        assert cl.max() - 1e-9 <= v <= rnd + 1e-9, f"c_tot {v} out of [max(cloud), random]"
    assert all(vals[i] >= vals[i + 1] - 1e-12 for i in range(len(vals) - 1)), f"not monotone in l_c: {vals}"
    print(f"  bugs_ctot: random=1-prod(1-cloud) ({rnd:.4f}) exact; c_tot in [max,random] & monotone-decreasing  PASS")


def test_edge_branches():
    pl2, tl, _ = _profile(3)
    nlm = tl.shape[0]
    assert float(bugs_ctot(pl2, tl, np.zeros(nlm), 1000.0)) == 0.0, "no cloud -> 0"
    one = np.zeros(nlm); one[5] = 0.6
    assert abs(float(bugs_ctot(pl2, tl, one, 1000.0)) - 0.6) < 1e-12, "one cloud -> its frac"
    full = np.zeros(nlm); full[3] = 1.0; full[8] = 0.5
    assert float(bugs_ctot(pl2, tl, full, 1000.0)) == 1.0, "max>=1 -> 1"
    print("  bugs_ctot edge branches (Ncloud=0/1, max>=1)  PASS")


def _ref_cloudfit(acld, c_tot, tol=1.0e-5):
    """Literal NumPy transcription of cloud_correlate.F90:bugs_cloudfit (one column, do-while)."""
    nlm = acld.shape[0]
    cf_max = np.zeros(nlm); cf_random = np.zeros(nlm)
    Nc = 0; min_frac = 1.0; max_frac = 0.0; frac_set = 0
    for j in range(nlm):
        if acld[j] > MIN_CF:
            frac_set = 1; Nc += 1
            min_frac = min(min_frac, acld[j]); max_frac = max(max_frac, acld[j])
    if frac_set == 0:
        return 0.0, cf_max, cf_random
    prod = 1.0
    for j in range(nlm):
        prod *= (1.0 - acld[j])
    c_tot_max = 1.0 - prod
    if c_tot_max - c_tot > tol:
        cf_stacked = 0.0; c_tot_calcd = c_tot_max; delta = max_frac / 1000.0
        while c_tot_calcd > c_tot and cf_stacked <= max_frac:
            prod = 1.0
            for j in range(nlm):
                if acld[j] > cf_stacked:
                    prod *= (1.0 - acld[j]) / (1.0 - cf_stacked)
            c_tot_calcd = cf_stacked + (1 - cf_stacked) * (1.0 - prod)
            cf_stacked += delta
        if cf_stacked > max_frac:
            cf_stacked = max_frac
        else:
            cf_stacked -= delta / 2.0
    else:
        cf_stacked = min_frac if Nc == 1 else 0.0
    c_maximal = cf_stacked
    for j in range(nlm):
        if acld[j] >= cf_stacked:
            cf_max[j] = cf_stacked; cf_random[j] = acld[j] - cf_stacked
        else:
            cf_max[j] = acld[j]; cf_random[j] = 0.0
    if 0.0 < c_maximal < 1.0:
        cf_max = cf_max / c_maximal
        cf_random = cf_random / (1.0 - c_maximal)
    return c_maximal, cf_max, cf_random


def test_cloudfit_vs_literal_loop():
    worst = 0.0
    for seed in range(40):
        pl2, tl, acld = _profile(seed)
        c_tot = float(bugs_ctot(pl2, tl, acld, 2000.0))
        cm, cfx, cfr = bugs_cloudfit(acld, c_tot)
        rcm, rcfx, rcfr = _ref_cloudfit(acld, c_tot)
        worst = max(worst, abs(cm - rcm), np.max(np.abs(cfx - rcfx)), np.max(np.abs(cfr - rcfr)))
    assert worst < 1e-12, f"bugs_cloudfit vs literal loop {worst:.2e}"
    print(f"  bugs_cloudfit vs literal NumPy do-while (40 profiles): max abs diff {worst:.1e}  PASS")


def test_cloudfit_reconstruction_invariant():
    for seed in range(40):
        pl2, tl, acld = _profile(seed + 100)
        c_tot = float(bugs_ctot(pl2, tl, acld, 2000.0))
        cm, cfx, cfr = bugs_cloudfit(acld, c_tot)
        recon = cm * cfx + (1.0 - cm) * cfr
        assert np.max(np.abs(recon - acld)) < 1e-12, f"reconstruction seed={seed}"
        assert 0.0 <= cm <= acld[acld > MIN_CF].max() + 1e-12, "c_maximal out of [0, max_frac]"
    # clear column -> zeros
    cm, cfx, cfr = bugs_cloudfit(np.zeros(20), 0.0)
    assert cm == 0.0 and not cfx.any() and not cfr.any(), "clear column"
    print("  bugs_cloudfit: c_maximal*cf_max + (1-c_maximal)*cf_random = acld (40 profiles) + clear column  PASS")


def test_differentiable_core():
    # bugs_ctot has discrete control flow (cloud-layer count/condensation) so it is not a differentiable
    # black box; the OVERLAP FORMULA is differentiable w.r.t. cloud fractions once the cloudy layers are
    # selected. Grad the differentiable core (the part that carries physical sensitivity).
    pl2, tl, acld = _profile(11)
    z = np.asarray(_midlayer_heights(jnp.asarray(pl2), jnp.asarray(tl)))
    idx = np.nonzero(acld > MIN_CF)[0]
    cloud = jnp.asarray(acld[idx]); zc = jnp.asarray(z[idx])
    g = jax.grad(lambda c: ctot_from_cloudy_layers(c, zc, 2000.0))(cloud)
    assert np.isfinite(np.asarray(g)).all() and np.any(np.asarray(g) != 0.0), "grad bad"
    # finite-difference check on one component
    eps = 1e-6
    base = float(ctot_from_cloudy_layers(cloud, zc, 2000.0))
    cloud2 = cloud.at[0].add(eps)
    fd = (float(ctot_from_cloudy_layers(cloud2, zc, 2000.0)) - base) / eps
    assert abs(float(g[0]) - fd) < 1e-4, f"grad {float(g[0])} vs FD {fd}"
    print(f"  jax.grad(ctot core) wrt cloud fractions: finite, nonzero, FD-correct (||g||={float(jnp.linalg.norm(g)):.3e})  PASS")


def main():
    print("test_cloud_correlate:")
    for t in (test_vs_literal_loop, test_random_limit_bounds_monotonic, test_edge_branches,
              test_cloudfit_vs_literal_loop, test_cloudfit_reconstruction_invariant, test_differentiable_core):
        t()
    print("All cloud_correlate checks PASSED")


if __name__ == "__main__":
    main()
