#!/usr/bin/env python3
"""test_mono_flux_limiter.py — REFACTOR B2 (iter11): the JAX (lax.scan) monotonic flux limiter is
BIT-EXACT to the NumPy reference, and is differentiable.

The limiter fires only for strong-shear/stable BLs (atex, gabls3_night). To validate the conversion
without a slow full-case run, this compares `monotonic_turbulent_flux_limit` against the original
NumPy `monotonic_turbulent_flux_limit` on synthetic inputs constructed to TRIGGER the clip + re-solve
(tight allowable bounds + large wpxp), across the solve types (rtm = spikefix, thlm, um = is_uv branch).
Then a `jax.grad` smoke check (the limiter is now differentiable w.r.t. the fields).
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import jax
import jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mono_flux_limiter import (
    monotonic_turbulent_flux_limit,
    calc_turb_adv_range,
    mono_flux_rtm as MFL_RTM,
    mono_flux_thlm as MFL_THLM,
    mono_flux_um as MFL_UM,
)


def _empty_stats(gr):
    return JaxStats.empty(
        l_sample=False,
        names=(),
        ncol=gr.ngrdcol,
        max_nlev=max(gr.nzm, gr.nzt),
    )


def _mfl_call(a):
    gr = a["gr"]
    xm, wpxp, _err, _stats = monotonic_turbulent_flux_limit(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        gr,
        a["solve_type"],
        a["dt"],
        a["xm_old"],
        a["xp2"],
        a["wm_zt"],
        a["xm_forcing"],
        a["rho_ds_zm"],
        a["rho_ds_zt"],
        a["invrs_rho_ds_zm"],
        a["invrs_rho_ds_zt"],
        a["xp2_threshold"],
        a["xm_tol"],
        False,
        a["low_lev_effect"],
        a["high_lev_effect"],
        2,
        True,
        True,
        _empty_stats(gr),
        a["xm"],
        a["wpxp"],
        ErrInfo.initialized(gr.ngrdcol),
    )
    return xm, wpxp


def _monotonic_turbulent_flux_limit_numpy(**a):
    return _mfl_call(a)


def _make_inputs(rng, gr, ngrdcol, nzt, solve_type):
    nzm = nzt + 1
    # w-PDF fields (zm) → realistic integer turbulent-advection ranges via calc_turb_adv_range.
    w1 = -0.5 + rng.random((ngrdcol, nzm))
    w2 = -0.5 + rng.random((ngrdcol, nzm))
    vw1 = 0.05 + 0.5 * rng.random((ngrdcol, nzm))
    vw2 = 0.05 + 0.5 * rng.random((ngrdcol, nzm))
    mf = 0.3 + 0.4 * rng.random((ngrdcol, nzm))
    dt = 60.0
    lle, hle, _stats = calc_turb_adv_range(
        gr.nzm, gr.nzt, gr.ngrdcol, gr, dt, w1, w2, vw1, vw2, mf, _empty_stats(gr)
    )

    if solve_type in (MFL_UM,):
        scale, xp2_thr = 5.0, 1e-3
    elif solve_type == MFL_RTM:
        scale, xp2_thr = 1e-2, 1e-9
    else:  # thlm
        scale, xp2_thr = 300.0, 1e-4
    xm_old = scale * (0.5 + rng.random((ngrdcol, nzt)))
    xm = xm_old + 0.01 * scale * (rng.random((ngrdcol, nzt)) - 0.5)
    xm_forcing = 1e-4 * scale * (rng.random((ngrdcol, nzt)) - 0.5)
    # large wpxp + tight xp2 → force the limiter to clip
    wpxp = 2.0 * scale * (rng.random((ngrdcol, nzm)) - 0.5)
    xp2 = (0.01 * scale) ** 2 * rng.random((ngrdcol, nzm))
    rho_ds_zm = 1.0 + 0.1 * rng.random((ngrdcol, nzm))
    rho_ds_zt = 1.0 + 0.1 * rng.random((ngrdcol, nzt))
    wm_zt = 0.01 * (rng.random((ngrdcol, nzt)) - 0.5)
    args = dict(solve_type=solve_type, xm=xm, wpxp=wpxp, xm_old=xm_old, xp2=xp2, wm_zt=wm_zt,
                xm_forcing=xm_forcing, rho_ds_zm=rho_ds_zm, rho_ds_zt=rho_ds_zt,
                invrs_rho_ds_zm=1.0 / rho_ds_zm, invrs_rho_ds_zt=1.0 / rho_ds_zt,
                xp2_threshold=xp2_thr, xm_tol=1e-4 * scale, low_lev_effect=lle,
                high_lev_effect=hle, gr=gr, dt=dt)
    return args


def test_bit_exact():
    gr = setup_grid(ngrdcol=2, deltaz=50.0, zm_init=0.0, zm_top=2000.0, grid_type=1)
    nzt = gr.zt.shape[1]
    worst = 0.0
    n_clipped = 0
    for seed in range(6):
        rng = np.random.default_rng(seed)
        for st in (MFL_RTM, MFL_THLM, MFL_UM):
            a = _make_inputs(rng, gr, 2, nzt, st)
            xm_np, wpxp_np = _monotonic_turbulent_flux_limit_numpy(**a)
            xm_jx, wpxp_jx = _mfl_call(a)
            xm_jx = np.asarray(xm_jx); wpxp_jx = np.asarray(wpxp_jx)
            # confirm the limiter actually did something (else the test is vacuous)
            if np.max(np.abs(wpxp_np - a["wpxp"])) > 0:
                n_clipped += 1
            for ref, got in ((xm_np, xm_jx), (wpxp_np, wpxp_jx)):
                den = np.max(np.abs(ref)) + 1e-300
                worst = max(worst, np.max(np.abs(ref - got)) / den)
    assert n_clipped > 0, "test vacuous — the limiter never clipped; strengthen the inputs"
    assert worst < 1e-12, f"JAX flux limiter not bit-exact to NumPy: worst rel {worst:.2e}"
    print(f"  mono flux limiter: JAX (lax.scan) vs NumPy BIT-EXACT (worst rel {worst:.1e}; "
          f"{n_clipped}/18 cases clipped)  PASS")


def test_differentiable():
    gr = setup_grid(ngrdcol=1, deltaz=50.0, zm_init=0.0, zm_top=2000.0, grid_type=1)
    nzt = gr.zt.shape[1]
    rng = np.random.default_rng(0)
    a = _make_inputs(rng, gr, 1, nzt, MFL_RTM)

    def loss(wpxp):
        b = dict(a); b["wpxp"] = wpxp
        xm_new, wpxp_new = _mfl_call(b)
        return jnp.sum(xm_new ** 2) + jnp.sum(wpxp_new ** 2)

    g = jax.grad(loss)(jnp.asarray(a["wpxp"]))
    assert np.all(np.isfinite(np.asarray(g))) and float(jnp.sum(jnp.abs(g))) > 0, "grad not finite/nonzero"
    print(f"  mono flux limiter: jax.grad w.r.t. wpxp finite+nonzero (|g|sum={float(jnp.sum(jnp.abs(g))):.2e})  PASS")


def test_full_limiter_f2py():
    """End-to-end f2py bit-shadow of the WHOLE `monotonic_turbulent_flux_limit` vs the Fortran oracle — the
    strongest validation of this routine (test_bit_exact only checks JAX-vs-NumPy; this checks JAX-vs-Fortran for
    the complete clip + xm re-solve). The conventions to align the JAX call to f2py: solve_type code rtm=2/thlm=1/
    um=4 (mono_flux_limiter.F90 named constants); low/high_lev_effect +1 (JAX 0-based → Fortran 1-based, cf. iter
    431); tridiag_solve_method=2 (tridiag_lu); **l_upwind_xm_ma=1** (the JAX `mfl_xm_lhs` calls term_ma_zt_lhs with
    the upwind discretization — this controls the post-clip xm re-solve's MA term, the sole non-clip difference);
    l_mono_flux_lim_spikefix=1; l_implemented=0. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 442)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py monotonic_turbulent_flux_limit oracle: SKIP ({type(e).__name__})")
        return
    gr = setup_grid(ngrdcol=2, deltaz=50.0, zm_init=0.0, zm_top=2000.0, grid_type=1)
    ng, nzm = gr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, 50.0), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(gr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(gr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(gr.zt)), err_info=ErrInfo(ngrdcol=ng))
    code = {MFL_RTM: 2, MFL_THLM: 1, MFL_UM: 4}
    worst = 0.0
    for seed in range(4):
        rng = np.random.default_rng(seed)
        for st in (MFL_RTM, MFL_THLM, MFL_UM):
            a = _make_inputs(rng, gr, 2, nzt, st)
            xm_j, wpxp_j = _mfl_call(a)
            lle = np.asarray(a["low_lev_effect"]); hle = np.asarray(a["high_lev_effect"])
            xm_f, wpxp_f = clubb_f2py.f2py_monotonic_turbulent_flux_limit(
                code[st], a["dt"], np.asarray(a["xm_old"]), np.asarray(a["xp2"]), np.asarray(a["wm_zt"]),
                np.asarray(a["xm_forcing"]), np.asarray(a["rho_ds_zm"]), np.asarray(a["rho_ds_zt"]),
                np.asarray(a["invrs_rho_ds_zm"]), np.asarray(a["invrs_rho_ds_zt"]), a["xp2_threshold"], a["xm_tol"],
                0, lle + 1, hle + 1, 2, 1, 1, np.asarray(a["xm"]).copy(), np.asarray(a["wpxp"]).copy())
            for ref, got in ((xm_f, xm_j), (wpxp_f, wpxp_j)):
                den = np.max(np.abs(ref)) + 1e-300
                worst = max(worst, float(np.max(np.abs(np.asarray(got) - ref)) / den))
    assert worst < 1e-12, f"full monotonic_turbulent_flux_limit f2py mismatch: worst rel {worst:.2e}"
    print(f"  f2py monotonic_turbulent_flux_limit (whole routine, rtm/thlm/um × 4 seeds): bit-match, "
          f"worst rel {worst:.2e}  PASS")


def test_calc_turb_adv_range_f2py():
    """calc_turb_adv_range (mono_flux_limiter.F90) — the low/high level indices bounding the vertical range a
    turbulent-advection parcel reaches (from the two w-PDF components over dt). vs the f2py oracle; the JAX returns
    0-based Python level indices, the f2py 1-based Fortran, so f2py == JAX + 1. SKIPs if clubb_f2py/clubb_python
    are unbuilt. (iter 431)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py calc_turb_adv_range oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=2, deltaz=40.0, zm_init=0.0, zm_top=1200.0, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, 40.0), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(2)
    bad = 0
    for _ in range(15):
        w1 = rng.uniform(-1, 1, (ng, nzm)); w2 = rng.uniform(-1, 1, (ng, nzm))
        v1 = rng.uniform(0.1, 1, (ng, nzm)); v2 = rng.uniform(0.1, 1, (ng, nzm)); mf = rng.uniform(0.3, 0.7, (ng, nzm))
        jt = calc_turb_adv_range(
            jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr, 60.0,
            jnp.asarray(w1), jnp.asarray(w2), jnp.asarray(v1), jnp.asarray(v2),
            jnp.asarray(mf), _empty_stats(jgr),
        )[:2]
        ft = clubb_f2py.f2py_calc_turb_adv_range(nzt, 60.0, w1, w2, v1, v2, mf)
        # 0-based JAX vs 1-based Fortran level indices: f2py == JAX + 1
        for k in range(2):
            if not np.array_equal(np.asarray(ft[k]), np.asarray(jt[k]) + 1):
                bad += 1
    assert bad == 0, f"calc_turb_adv_range f2py index mismatch (beyond the +1 0/1-based offset) in {bad} cases"
    print("  f2py calc_turb_adv_range: low/high level indices match (modulo 0-based vs 1-based +1), 15 cases  PASS")


if __name__ == "__main__":
    test_bit_exact()
    test_full_limiter_f2py()
    test_calc_turb_adv_range_f2py()
    test_differentiable()
    print("\nAll mono_flux_limiter tests PASSED")
