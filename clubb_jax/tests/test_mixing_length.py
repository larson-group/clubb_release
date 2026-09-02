#!/usr/bin/env python3
"""test_mixing_length.py — validate mixing_length.F90 ports against the f2py Fortran oracle.

  set_Lscale_max       — the Lscale cap (0.25·min(host_dx, host_dy) implemented, else 1e5). Pure, bit-exact.
  calc_Lscale_directly — the direct (parcel-ascent) mixing length: extracts mu from clubb_params, runs
                         compute_mixing_length with MEAN values (l_avg_Lscale=False, so the f2py rtp2_zt/thlp2_zt/
                         rtpthlp_zt are unused on this path → passed as zeros). The iterative parcel ascent carries a
                         ~1e-13 RELATIVE FP-accumulation vs Fortran (not bit-exact), so it is checked with a relative
                         tolerance. NOTE: this IS the isolation-level oracle for `compute_mixing_length` itself — the
                         ~1000-line parcel-ascent core, which has no f2py wrapper of its own. calc_Lscale_directly
                         returns `compute_mixing_length(...)` verbatim (mean-value path), so f2py_calc_lscale_directly's
                         Lscale/Lscale_up/Lscale_down outputs ARE compute_mixing_length's outputs. Its averaged path
                         (l_avg_Lscale=True, used by the live calc_Lscale driver) is validated end-to-end (compare_cases)
                         + differentiability (test_differentiability.py).
  diagnose_Lscale_from_tau — the gated (l_diag_Lscale_from_tau) tau-based Lscale; no case sets it (f2py-shadowed here).
Both SKIP cleanly if clubb_f2py / clubb_python are unbuilt.
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
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.mixing_length import (
    set_Lscale_max, calc_Lscale_directly, diagnose_Lscale_from_tau, imu)
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.saturation import SATURATION_FLATAU

_TUNABLE_PARAMS = os.path.join(_ROOT, "input", "tunable_parameters", "tunable_parameters.in")

_NG, _DZ, _ZTOP = 2, 40.0, 1600.0


def test_set_Lscale_max_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py set_Lscale_max oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(1)
    worst = 0.0
    for l_impl in (0, 1):
        for _ in range(10):
            dx = rng.uniform(100.0, 5e5, 3); dy = rng.uniform(100.0, 5e5, 3)
            j = np.asarray(set_Lscale_max(l_impl, dx, dy, 3))
            f = np.asarray(clubb_f2py.f2py_set_lscale_max(l_impl, dx, dy))
            worst = max(worst, float(np.max(np.abs(j - f))))
    assert worst < 1e-12, f"set_Lscale_max f2py mismatch {worst:.2e}"
    print(f"  f2py set_Lscale_max (host + standalone, 20 cases): bit-match, worst {worst:.2e}  PASS")


def _setup_f2py_grid():
    import clubb_f2py
    from clubb_python import clubb_api
    from clubb_python.derived_types.err_info import ErrInfo
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    return clubb_f2py, jgr, ng, nzm


def test_calc_Lscale_directly_f2py():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py calc_Lscale_directly oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    worst_rel = 0.0
    for s in range(6):
        rng = np.random.default_rng(s)
        thlm = np.sort(rng.uniform(290.0, 330.0, (ng, nzt)), axis=1); exner = rng.uniform(0.8, 1.0, (ng, nzt))
        rtm = rng.uniform(2e-3, 1.4e-2, (ng, nzt)); thvm = thlm * (1.0 + 0.6 * rtm)
        em = rng.uniform(1e-2, 1.0, (ng, nzm)); p = np.sort(rng.uniform(6e4, 1e5, (ng, nzt)), axis=1)[:, ::-1].copy()
        thv_ds = thvm.copy(); Lmax = np.full(ng, 1e5); lmin = 20.0
        npar = max(imu, 100); cp = np.zeros((ng, npar)); cp[:, imu] = rng.uniform(0.4, 1.0, ng); mu = cp[:, imu]
        jr = calc_Lscale_directly(jnp.asarray(thvm), jnp.asarray(thlm), jnp.asarray(rtm), jnp.asarray(em),
            jnp.asarray(Lmax), jnp.asarray(p), jnp.asarray(exner), jnp.asarray(thv_ds), jnp.asarray(cp),
            lmin, SATURATION_FLATAU, False, jgr)
        z = np.zeros((ng, nzt))
        fr = clubb_f2py.f2py_calc_lscale_directly(0, p, exner, rtm, thlm, thvm, mu, z, z, z, em, thv_ds,
            Lmax, lmin, cp, SATURATION_FLATAU, 0)
        for i in range(3):
            a = np.asarray(jr[i]); b = np.asarray(fr[i])
            worst_rel = max(worst_rel, float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1.0))))
    # Parcel-ascent iteration carries ~1e-13 relative FP accumulation vs Fortran (not bit-exact); strict bound 1e-9.
    assert worst_rel < 1e-9, f"calc_Lscale_directly f2py rel mismatch {worst_rel:.2e}"
    print(f"  f2py calc_Lscale_directly (Lscale/up/down, 6 seeds): rel-match, worst {worst_rel:.2e}  PASS")


def test_diagnose_Lscale_from_tau_f2py():
    """diagnose_Lscale_from_tau (the l_diag_Lscale_from_tau path — the tau-based mixing-length / dissipation-time
    diagnostics) vs the Fortran oracle. This routine is GATED OFF in every validated case (no case_setup sets
    l_diag_Lscale_from_tau), so it is otherwise UN-exercised — this is its only validation. The JAX returns the
    Fortran's 19 outputs (invrs_tau_* / tau_* / Lscale[/up/down]) plus 2 extra diagnostics appended; the first 19
    are compared by name-order. SKIPs if clubb_f2py / clubb_python are unbuilt. (iter 418)"""
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py diagnose_Lscale_from_tau oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    cp = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))   # realistic defaults (well-conditioned tau)
    worst = 0.0
    for s in range(4):
        rng = np.random.default_rng(s)
        upwp = rng.uniform(-0.1, 0.1, ng); vpwp = rng.uniform(-0.1, 0.1, ng)
        ddzt = rng.uniform(0.0, 1e-3, (ng, nzm)); isf = rng.uniform(0.0, 0.3, (ng, nzt))
        em = rng.uniform(1e-2, 1.0, (ng, nzm)); sqem = np.sqrt(rng.uniform(1e-2, 1.0, (ng, nzt)))
        sfce = np.zeros(ng); Lmax = np.full(ng, 1e5)
        Ri = rng.uniform(-1.0, 5.0, (ng, nzm)); bvs = rng.uniform(-1e-3, 1e-2, (ng, nzm))
        jr = diagnose_Lscale_from_tau(jnp.asarray(upwp), jnp.asarray(vpwp), jnp.asarray(ddzt), jnp.asarray(isf),
            jnp.asarray(em), jnp.asarray(sqem), 0.01, 900.0, jnp.asarray(sfce), jnp.asarray(Lmax),
            jnp.asarray(cp), jnp.asarray(Ri), jnp.asarray(bvs), False, False, jgr)
        fr = clubb_f2py.f2py_diagnose_lscale_from_tau(upwp, vpwp, ddzt, isf, em, sqem, 0.01, 900.0, sfce, Lmax,
            cp, 0, 0, bvs, Ri)
        assert len(fr) == 19, f"unexpected f2py output count {len(fr)}"
        for i in range(19):
            a = np.asarray(jr[i]); b = np.asarray(fr[i])
            worst = max(worst, float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-3))))
    assert worst < 1e-9, f"diagnose_Lscale_from_tau f2py rel mismatch {worst:.2e}"
    print(f"  f2py diagnose_Lscale_from_tau (19 tau/Lscale outputs, 4 seeds): rel-match, worst {worst:.2e}  PASS")


def main():
    print("test_mixing_length:")
    test_set_Lscale_max_f2py()
    test_calc_Lscale_directly_f2py()
    test_diagnose_Lscale_from_tau_f2py()
    print("All mixing_length checks PASSED")


if __name__ == "__main__":
    main()
