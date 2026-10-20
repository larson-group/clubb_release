"""Check the JAX monotonic flux limiter ports (mono_flux_limiter.F90) against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.mono_flux_limiter as fortran
from clubb_python.derived_types.err_info import ErrInfo as FortranErrInfo

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
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


def test_full_limiter_f2py():
    """The whole `monotonic_turbulent_flux_limit` (clip + xm re-solve) matches Fortran.

    Conventions that align the Fortran call with the JAX one: solve_type codes rtm=2/thlm=1/um=4
    (mono_flux_limiter.F90 named constants); tridiag_solve_method=2 (tridiag_lu); l_upwind_xm_ma=True,
    because the JAX `mfl_xm_lhs` uses the upwind term_ma_zt_lhs in the post-clip xm re-solve;
    l_mono_flux_lim_spikefix=True; l_implemented=False. The clubb_python wrapper converts the 0-based
    low/high_lev_effect indices to Fortran's 1-based ones."""
    gr, fgr = setup_matching_grids(2, 50.0, 2000.0)
    ng, nzm = gr.zm.shape; nzt = nzm - 1
    code = {MFL_RTM: 2, MFL_THLM: 1, MFL_UM: 4}
    worst = 0.0
    for seed in range(4):
        rng = np.random.default_rng(seed)
        for st in (MFL_RTM, MFL_THLM, MFL_UM):
            a = _make_inputs(rng, gr, 2, nzt, st)
            xm_j, wpxp_j = _mfl_call(a)
            # The wrapper converts the 0-based level indices to Fortran's 1-based ones.
            err_f, xm_f, wpxp_f = fortran.monotonic_turbulent_flux_limit(
                nzm, nzt, ng, fgr, code[st], a["dt"], *(np.asarray(a[k]) for k in (
                    "xm_old", "xp2", "wm_zt", "xm_forcing", "rho_ds_zm", "rho_ds_zt",
                    "invrs_rho_ds_zm", "invrs_rho_ds_zt")), a["xp2_threshold"], a["xm_tol"],
                False, np.asarray(a["low_lev_effect"]), np.asarray(a["high_lev_effect"]), 2, True, True,
                np.asarray(a["xm"]).copy(), np.asarray(a["wpxp"]).copy(), FortranErrInfo(ngrdcol=ng))
            assert not err_f.is_fatal(), "Fortran monotonic_turbulent_flux_limit reported a fatal error"
            for ref, got in ((xm_f, xm_j), (wpxp_f, wpxp_j)):
                den = np.max(np.abs(ref)) + 1e-300
                worst = max(worst, float(max_abs_diff("mono_flux_limiter", np.asarray(got), ref) / den))
    assert worst < 1e-12, f"full monotonic_turbulent_flux_limit f2py mismatch: worst rel {worst:.2e}"


def test_calc_turb_adv_range_f2py():
    """calc_turb_adv_range (mono_flux_limiter.F90) — the low/high level indices bounding the vertical range a
    turbulent-advection parcel reaches (from the two w-PDF components over dt). vs the f2py oracle; the JAX returns
    0-based Python level indices, and so does the clubb_python wrapper (it converts from Fortran's 1-based ones)."""
    jgr, fgr = setup_matching_grids(2, 40.0, 1200.0)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
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
        ft = fortran.calc_turb_adv_range(nzm, nzt, ng, fgr, 60.0, w1, w2, v1, v2, mf)
        for k in range(2):
            if not np.array_equal(np.asarray(ft[k]), np.asarray(jt[k])):
                bad += 1
    assert bad == 0, f"calc_turb_adv_range Fortran index mismatch in {bad} cases"
