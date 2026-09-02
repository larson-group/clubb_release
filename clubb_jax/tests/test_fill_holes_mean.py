"""Unit tests for the MEAN-FIELD vertical hole-filler (the rtm_cl / thlm_cl mechanism).

These lock in the Iter186 fix: advance_xm_wpxp applies `fill_holes_vertical` to the mean fields
rtm/thlm after the solve (advance_xm_wpxp_module.F90:4974-5018, the `rtm_cl`/`thlm_cl` budget). For
a case whose mean field dips below its tolerance at a stretched dry top (rico, rt→0 above ~9 km), this
fill — NOT an IC floor — is what keeps the field >= tol each step, mass-conservingly pulling a tiny
amount from the topmost above-tol level (the rico step-1 k51 seed of 1.37e-8 = that pull × dt).

A regression that drops the rtm/thlm fill (or re-adds an IC floor) would silently reintroduce the rico
seed; these tests are the fast unit-level guard (no Fortran oracle needed).
"""
import os
import sys

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical, fill_holes_wp2_from_horz_tke

RT_TOL = 1.0e-8


def _col_mass(field, rho_ds, dz):
    """Column mass with the same rho_ds*dz weighting fill_holes conserves."""
    return float(np.sum(np.asarray(field) * np.asarray(rho_ds) * np.asarray(dz)))


def test_mean_fill_raises_dry_top_and_conserves_mass():
    """A rico-like rtm column (moist below, a sharp drop to 0 at the top) — fill must raise the
    sub-tol dry top to >= rt_tol, conserve rho_ds-weighted column mass to machine precision, and pull
    the deficit from the topmost above-tol (donor) level."""
    nzt = 24
    # Stretched-ish density (decreasing with height) and non-uniform dz.
    rho = np.linspace(1.1, 0.35, nzt)[None, :]
    dz = np.linspace(30.0, 120.0, nzt)[None, :]
    # Moist below, then a hard cliff to 0 over the top 5 levels (the dry top).
    rtm = np.full((1, nzt), 5.0e-3)
    rtm[0, -8:] = np.array([3.0e-3, 1.2e-3, 4.0e-4, 9.0e-5, 0.0, 0.0, 0.0, 0.0])
    rtm_j = jnp.asarray(rtm)

    m_before = _col_mass(rtm, rho, dz)
    out = np.asarray(fill_holes_vertical(
        field=rtm_j, rho_ds=jnp.asarray(rho), dz=jnp.asarray(dz),
        threshold=RT_TOL, lower_k=0, upper_k=nzt - 1, fill_holes_type=2))
    m_after = _col_mass(out, rho, dz)

    assert np.all(out >= RT_TOL - 1e-300), "dry top not raised to threshold"
    # Mass conservation to MACHINE PRECISION. The global redistribution conserves rho_ds-weighted column mass
    # exactly in exact arithmetic; in x64 the residual is ~1e-16 (the earlier "~1e-8 precision floor" note was a
    # float32 artifact of this file before it enabled jax_enable_x64 — see iters 432/437). The JAX bit-matches the
    # Fortran fill_holes_vertical (test_fill_holes_vertical_f2py, 3.3e-16).
    rel_mass = abs(m_after - m_before) / abs(m_before)
    assert rel_mass < 1e-12, f"column mass not conserved to machine precision: rel {rel_mass:.2e}"
    # The donor (topmost moist level, index where rtm was 9e-5) must have given up mass.
    donor = int(np.where(rtm[0] >= RT_TOL)[0][-1])
    assert out[0, donor] < rtm[0, donor], "donor level did not give up mass to fill the dry top"
    print(f"  mean-field fill: dry top raised, mass conserved (rel {rel_mass:.1e}), "
          f"donor k{donor} {rtm[0,donor]:.2e}->{out[0,donor]:.2e}  PASS")


def test_mean_fill_noop_when_all_above_threshold():
    """When every level is >= threshold (all 15 uniform-grid cases), the fill is a BITWISE no-op —
    this is why adding the rtm/thlm fill did not perturb any faithful case."""
    nzt = 20
    rho = np.linspace(1.1, 0.5, nzt)[None, :]
    dz = np.full((1, nzt), 40.0)
    rtm = jnp.asarray(np.linspace(8.0e-3, 1.0e-6, nzt)[None, :])  # all >> rt_tol
    out = fill_holes_vertical(
        field=rtm, rho_ds=jnp.asarray(rho), dz=jnp.asarray(dz),
        threshold=RT_TOL, lower_k=0, upper_k=nzt - 1, fill_holes_type=2)
    assert np.array_equal(np.asarray(out), np.asarray(rtm)), "fill changed an all-above-threshold field"
    print("  mean-field fill: bitwise no-op when all >= threshold  PASS")


def test_thlm_fill_is_noop():
    """thlm ~300 K >> thl_tol, so the thlm_cl fill is always a guaranteed no-op (mirrors the Fortran
    but never changes thlm) — guards against an accidental threshold that would corrupt thlm."""
    nzt = 16
    rho = np.linspace(1.1, 0.6, nzt)[None, :]
    dz = np.full((1, nzt), 40.0)
    thl_tol = 1.0e-2
    thlm = jnp.asarray(np.linspace(298.0, 320.0, nzt)[None, :])
    out = fill_holes_vertical(
        field=thlm, rho_ds=jnp.asarray(rho), dz=jnp.asarray(dz),
        threshold=thl_tol, lower_k=0, upper_k=nzt - 1, fill_holes_type=2)
    assert np.array_equal(np.asarray(out), np.asarray(thlm)), "thlm fill changed thlm"
    print("  thlm fill: no-op (thlm >> thl_tol)  PASS")


def test_fill_holes_vertical_f2py():
    """fill_holes_vertical (fill_holes.F90) — the mass-conserving vertical hole-filler, validated vs the f2py oracle
    for BOTH fill_holes_type=1 (single global pass) and 2 (sliding 3-pt window + global fallback). The JAX lower_k/
    upper_k are 0-based (Fortran lower/upper_hf_level = +1). This is the validation wrongly reverted at iter 411 as an
    "FP-floor": that 5.91e-08 was the SAME float32 artifact later diagnosed at iter 432 (the file lacked
    jax_enable_x64) — with x64 it is bit-faithful (~3.3e-16), NOT an FP-floor. SKIPs if clubb_f2py is unbuilt.
    (iter 437)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py fill_holes_vertical oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(0)
    worst = 0.0
    for fht in (1, 2):
        for _ in range(20):
            ng, nz = 2, 16
            field = rng.uniform(-0.05, 1.0, (ng, nz))   # negatives = holes to fill
            rho = rng.uniform(0.3, 1.2, (ng, nz)); dz = rng.uniform(20.0, 120.0, (ng, nz)); thr = 1e-3
            j = np.asarray(fill_holes_vertical(jnp.asarray(field), jnp.asarray(rho), jnp.asarray(dz),
                                               thr, 0, nz - 1, fht, 1))
            f = np.asarray(clubb_f2py.f2py_fill_holes_vertical(thr, 1, nz, dz.copy(), rho.copy(), 1, fht, field.copy()))
            worst = max(worst, float(np.max(np.abs(j - f))))
    assert worst < 1e-12, f"fill_holes_vertical f2py mismatch {worst:.2e}"
    print(f"  f2py fill_holes_vertical (fill_holes_type 1 + 2, 40 cases): bit-match, worst {worst:.2e}  PASS")


def test_fill_holes_wp2_from_horz_tke_f2py():
    """fill_holes_wp2_from_horz_tke (fill_holes.F90) — the TKE-conserving wp2 hole-fill: where wp2 < threshold and
    up2/vp2 have spare TKE, it borrows from up2/vp2 proportionally to raise wp2. Unlike fill_holes_vertical this is
    a direct proportional borrow (NOT the FP-floor mass_frac redistribution), so it bit-shadows the f2py oracle
    cleanly. The JAX lower_k/upper_k are 0-based (Fortran lower/upper_hf_level = +1). SKIPs if clubb_f2py is unbuilt.
    (iter 432)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py fill_holes_wp2_from_horz_tke oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(20):
        ng, nzm = 2, 12
        wp2 = rng.uniform(-0.1, 1.0, (ng, nzm))   # negatives / sub-threshold -> holes to fill
        up2 = rng.uniform(0.05, 1.0, (ng, nzm)); vp2 = rng.uniform(0.05, 1.0, (ng, nzm)); thr = 1e-3
        j = fill_holes_wp2_from_horz_tke(jnp.asarray(wp2), jnp.asarray(up2), jnp.asarray(vp2), thr, 0, nzm - 1)
        f = clubb_f2py.f2py_fill_holes_wp2_from_horz_tke(thr, 1, nzm, np.asfortranarray(wp2.copy()),
                                                         np.asfortranarray(up2.copy()), np.asfortranarray(vp2.copy()))
        for k in range(3):   # wp2, up2, vp2
            worst = max(worst, float(np.max(np.abs(np.asarray(j[k]) - np.asarray(f[k])))))
    assert worst < 1e-12, f"fill_holes_wp2_from_horz_tke f2py mismatch {worst:.2e}"
    print(f"  f2py fill_holes_wp2_from_horz_tke (wp2/up2/vp2 TKE-conserving fill, 20 cases): bit-match, worst {worst:.2e}  PASS")


if __name__ == "__main__":
    print("Mean-field fill_holes (rtm_cl / thlm_cl) tests:")
    test_mean_fill_raises_dry_top_and_conserves_mass()
    test_mean_fill_noop_when_all_above_threshold()
    test_thlm_fill_is_noop()
    test_fill_holes_vertical_f2py()
    test_fill_holes_wp2_from_horz_tke_f2py()
    print("All mean-field fill_holes tests PASSED.")
