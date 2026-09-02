#!/usr/bin/env python3
"""test_ice_dfsn.py — validate the JAX ice_dfsn port (cloud-water depletion by diffusional growth of ice).

Oracles:
  1. `thlm2T_in_K` is bit-validated vs the f2py oracle `f2py_thlm2t_in_k_1d` (SKIPs if unbuilt).
  2. `ice_dfsn` is validated vs a LITERAL NumPy transcription of the Fortran top-to-bottom mass-integration
     loop (the algorithm is sequential, so a direct Python loop is a faithful, independent reference) — rel ~1e-14.
  3. Branch coverage: tendencies are 0 outside in-cloud/below-freezing; the over-depletion cap is exercised.
  4. The thlm tendency equals -(Lv/(Cp*exner)) * rcm tendency exactly.
  5. A finite `jax.grad` through the whole routine.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
# Append (not prepend) the f2py directory so it cannot shadow this checkout's package.
_F2PY = os.path.join(_ROOT, "clubb_python_api")
if _F2PY not in sys.path:
    sys.path.append(_F2PY)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp, Lv, ep, Rv, Lf, Ls, T_freeze_K, cm_per_m)
from clubb_jax.src.Microphys.ice_dfsn_module import ice_dfsn
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_SAT_FORMULA = 3   # SATURATION_FLATAU

# ice_dfsn_module.F90 parameters (duplicated here as an independent reference).
_N_I, _MASS_INIT, _RCM_THR = 2000.0, 1.0e-11, 1.0e-5
_A, _B, _KU, _Q, _N = 2.05e-3, 1.8, 55.0, 0.17, 0.70


def _make_grid_and_profile(nzt=40, deltaz=100.0):
    gr = setup_grid(1, deltaz, deltaz, deltaz * (nzt + 1))
    rng = np.random.default_rng(2026)
    # A mixed-phase-ish profile: cold aloft, some cloud water in a mid layer.
    z = np.arange(nzt) * deltaz
    T = 285.0 - 6.5e-3 * z + rng.standard_normal(nzt) * 0.5   # crosses freezing partway up
    exner = (1.0 - 2.0e-5 * z)
    thlm = T / exner
    p = 1.0e5 * exner ** (1004.67 / 287.04)
    rho = p / (287.04 * T)
    rcm = np.zeros(nzt)
    rcm[8:24] = 5.0e-5 + np.abs(rng.standard_normal(16)) * 5.0e-5   # cloud layer, > threshold
    rcm[12] = 2.0e-8   # one sub-threshold point inside the layer
    return gr, thlm, rcm, exner, p, rho


def _numpy_ice_dfsn(gr, dt, thlm, rcm, exner, p, rho):
    """Literal NumPy transcription of the Fortran ice_dfsn loop (independent reference)."""
    nzt = thlm.shape[0]
    T = thlm * exner + Lv * rcm / Cp
    mass = np.full(nzt, _MASS_INIT)
    rcm_ice = np.zeros(nzt)
    inv_dzm = np.asarray(gr.invrs_dzm[0])
    for k in range(nzt - 1, -1, -1):   # 0-based for Fortran k=nzt..1
        if rcm[k] >= _RCM_THR and T[k] < T_freeze_K:
            r_s = float(sat_mixrat_liq(p[k], T[k], _SAT_FORMULA))
            e_s = (r_s * p[k]) / (ep + r_s)
            e_i = e_s / np.exp((Lf / (Rv * T_freeze_K)) * (T_freeze_K / T[k] - 1.0))
            S_i = e_s / e_i
            Cel = T[k] - T_freeze_K
            Ka = 4.1868 * 100.0 * ((5.69 + 0.017 * Cel) * 0.00001)
            Dv = (0.221 * (T[k] / T_freeze_K) ** 1.94 * (101325.0 / p[k])) / 10000.0
            Denom = (Ls / (Rv * T[k]) - 1.0) * Ls / (Ka * T[k]) + (Rv * T[k]) / (Dv * e_i)
            fac = 4.0 * (S_i - 1.0) / Denom
            base = mass[k] / _A
            rcm_ice[k] = -(_N_I / rho[k]) * fac * base ** (1.0 / _B)
            if rcm[k] + rcm_ice[k] * dt < 0.0:
                rcm_ice[k] = -rcm[k] / dt
            inv_dzm_lag = inv_dzm[k - 1] if k - 1 >= 0 else inv_dzm[0]
            dmass = fac * (1.0 / _KU) * rho[k] ** _Q * base ** ((1.0 - _N) / _B) * (1.0 / inv_dzm_lag)
            if k > 0:
                mass[k - 1] = mass[k] + dmass
        else:
            if k > 0:
                mass[k - 1] = mass[k]
            rcm_ice[k] = 0.0
    thlm_ice = -(Lv / (Cp * exner)) * rcm_ice
    return rcm_ice, thlm_ice


def test_thlm2T_in_K_vs_f2py():
    rng = np.random.default_rng(1)
    thlm = 290.0 + rng.standard_normal(30) * 5
    exner = 0.9 + rng.standard_normal(30) * 0.02
    rcm = np.abs(rng.standard_normal(30)) * 1e-4
    got = np.asarray(thlm2T_in_K(jnp.asarray(thlm), jnp.asarray(exner), jnp.asarray(rcm)))
    try:
        import clubb_f2py
    except Exception:
        print("  thlm2T_in_K vs f2py: SKIP (oracle unbuilt)")
        return
    ref = np.asarray(clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm))
    rel = np.max(np.abs(got - ref) / (np.abs(ref) + 1e-30))
    assert rel == 0.0, f"thlm2T_in_K vs f2py rel {rel:.2e}"
    print(f"  thlm2T_in_K vs f2py_thlm2t_in_k_1d: BIT-EXACT (rel {rel:.0e})  PASS")


def test_ice_dfsn_vs_numpy():
    gr, thlm, rcm, exner, p, rho = _make_grid_and_profile()
    dt = 60.0
    g_rcm, g_thl = ice_dfsn(gr, dt, thlm, rcm, exner, p, rho, _SAT_FORMULA)
    g_rcm, g_thl = np.asarray(g_rcm), np.asarray(g_thl)
    r_rcm, r_thl = _numpy_ice_dfsn(gr, dt, thlm, rcm, exner, p, rho)
    rel_rcm = np.max(np.abs(g_rcm - r_rcm) / (np.abs(r_rcm) + 1e-30))
    rel_thl = np.max(np.abs(g_thl - r_thl) / (np.abs(r_thl) + 1e-30))
    assert rel_rcm < 1e-12, f"rcm_icedfsn rel {rel_rcm:.2e}"
    assert rel_thl < 1e-12, f"thlm_icedfsn rel {rel_thl:.2e}"
    print(f"  ice_dfsn vs literal NumPy loop: rcm rel {rel_rcm:.2e}, thlm rel {rel_thl:.2e}  PASS")


def test_branches_and_coupling():
    gr, thlm, rcm, exner, p, rho = _make_grid_and_profile()
    T = np.asarray(thlm2T_in_K(jnp.asarray(thlm), jnp.asarray(exner), jnp.asarray(rcm)))
    g_rcm, g_thl = ice_dfsn(gr, 60.0, thlm, rcm, exner, p, rho, _SAT_FORMULA)
    g_rcm = np.asarray(g_rcm)
    outside = ~((rcm >= _RCM_THR) & (T < T_freeze_K))
    assert np.all(g_rcm[outside] == 0.0), "tendency nonzero outside in-cloud/below-freezing"
    assert np.any(g_rcm < 0.0), "expected some negative (depleting) rcm tendency in cloud"
    # thlm tendency must be exactly -(Lv/(Cp*exner)) * rcm tendency.
    assert np.allclose(np.asarray(g_thl), -(Lv / (Cp * exner)) * g_rcm, rtol=0, atol=0), "thlm/rcm coupling broken"
    print("  branch (in-cloud/freezing) + exact thlm/rcm coupling  PASS")


def test_over_depletion_cap():
    # Tiny dt-relative rcm so the unclamped tendency would over-deplete -> cap to -rcm/dt.
    gr, thlm, _, exner, p, rho = _make_grid_and_profile()
    nzt = thlm.shape[0]
    rcm = np.zeros(nzt)
    rcm[8:24] = 1.01e-5   # just above the in-cloud threshold
    dt = 1.0e5            # large dt to force the unclamped tendency to over-deplete (exercises the cap)
    g_rcm, _ = ice_dfsn(gr, dt, thlm, rcm, exner, p, rho, _SAT_FORMULA)
    g_rcm = np.asarray(g_rcm)
    # No level may remove more liquid than exists: rcm + tend*dt >= 0 (within fp).
    assert np.all(rcm + g_rcm * dt >= -1e-18), "over-depletion cap failed"
    # And at least one level should hit the cap exactly (-rcm/dt).
    capped = np.isclose(g_rcm, -rcm / dt, rtol=1e-12, atol=0) & (rcm > 0)
    assert np.any(capped), "expected at least one capped level"
    print(f"  over-depletion cap (no negative liquid; {capped.sum()} levels capped)  PASS")


def test_differentiable():
    gr, thlm, rcm, exner, p, rho = _make_grid_and_profile()
    def loss(rcm_in):
        r, _ = ice_dfsn(gr, 60.0, jnp.asarray(thlm), rcm_in, jnp.asarray(exner),
                        jnp.asarray(p), jnp.asarray(rho), _SAT_FORMULA)
        return jnp.sum(r ** 2)
    g = jax.grad(loss)(jnp.asarray(rcm))
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad through ice_dfsn"
    print(f"  jax.grad(ice_dfsn) wrt rcm: finite (||g||={float(jnp.linalg.norm(g)):.3e})  PASS")


def main():
    print("test_ice_dfsn:")
    for t in (test_thlm2T_in_K_vs_f2py, test_ice_dfsn_vs_numpy, test_branches_and_coupling,
              test_over_depletion_cap, test_differentiable):
        t()
    print("All ice_dfsn checks PASSED")


if __name__ == "__main__":
    main()
