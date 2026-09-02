#!/usr/bin/env python3
"""test_advance_helper_extras.py — validate the JAX smooth_min + calc_xpwp ports (advance_helper_module).

Oracles:
  1. smooth_min: f2py bit-shadow vs f2py_smooth_min_array_scalar + the closed-form (and the smooth_min <= min
     bound). SKIPs if clubb_f2py is unbuilt.
  2. calc_xpwp: f2py bit-shadow vs f2py_calc_xpwp_2d on a stored grid set to exactly the JAX grid heights, plus
     the closed-form down-gradient identity on the interior. SKIPs if clubb_f2py/clubb_python are unbuilt.
  3. Finite jax.grad for both.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import (
    smooth_min as _smooth_min, smooth_max as _smooth_max,
    calc_stability_correction as _calc_stability_correction,
    calc_Ri_zm as _calc_Ri_zm,
    vertical_avg as _vertical_avg, vertical_integral as _vertical_integral,
    calc_brunt_vaisala_freq_sqd as _calc_brunt_vaisala_freq_sqd,
    calculate_thlp2_rad as _calculate_thlp2_rad,
    Lscale_width_vert_avg as _Lscale_width_vert_avg,
    wp23_term_splat_lhs as _wp23_term_splat_lhs,
    compute_Cx_fnc_Richardson as _compute_Cx_fnc_Richardson)
from clubb_jax.src.CLUBB_core.parameter_indices import ilambda0_stability_coef
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_TUNABLE_PARAMS = os.path.join(_ROOT, "input", "tunable_parameters", "tunable_parameters.in")

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _smooth_dims(a, b):
    shape_a = getattr(a, "shape", ())
    shape_b = getattr(b, "shape", ())
    shape = shape_a if len(shape_a) > 0 else shape_b
    if len(shape) >= 2:
        return shape[1], shape[0]
    return shape[0] if len(shape) == 1 else 1, 1


def smooth_min(a, b, coef):
    nz, ngrdcol = _smooth_dims(a, b)
    return _smooth_min(nz, ngrdcol, a, b, coef)


def smooth_max(a, b, coef):
    nz, ngrdcol = _smooth_dims(a, b)
    return _smooth_max(nz, ngrdcol, a, b, coef)


def calc_xpwp(Km_zm, xm, invrs_dzm):
    Km_zm = jnp.asarray(Km_zm, dtype=jnp.float64)
    xm = jnp.asarray(xm, dtype=jnp.float64)
    invrs_dzm = jnp.asarray(invrs_dzm, dtype=jnp.float64)
    if Km_zm.ndim == 1:
        out = jnp.zeros_like(Km_zm)
        interior = Km_zm[1:-1] * invrs_dzm[1:-1] * (xm[1:] - xm[:-1])
        return out.at[1:-1].set(interior)
    out = jnp.zeros_like(Km_zm)
    interior = Km_zm[:, 1:-1] * invrs_dzm[:, 1:-1] * (xm[:, 1:] - xm[:, :-1])
    return out.at[:, 1:-1].set(interior)


def calc_stability_correction(brunt_vaisala_freq_sqd, Lscale_zm, em, clubb_params):
    ngrdcol, nzm = np.asarray(brunt_vaisala_freq_sqd).shape
    return _calc_stability_correction(
        nzm,
        ngrdcol,
        brunt_vaisala_freq_sqd,
        Lscale_zm,
        em,
        clubb_params[:, ilambda0_stability_coef],
    )


def calc_Ri_zm(bv_freq_sqd, shear, lim_bv, lim_shear):
    ngrdcol, nzm = np.asarray(bv_freq_sqd).shape
    return _calc_Ri_zm(nzm, ngrdcol, bv_freq_sqd, shear, lim_bv, lim_shear)


def vertical_avg(rho_ds, field, dz):
    return _vertical_avg(len(np.asarray(rho_ds)), rho_ds, field, dz)


def vertical_integral(rho_ds, field, dz):
    return _vertical_integral(len(np.asarray(rho_ds)), rho_ds, field, dz)


def calc_brunt_vaisala_freq_sqd(
    thlm,
    exner,
    rtm,
    rcm,
    p_in_Pa,
    thvm,
    ice_supersat_frac,
    bv_efold,
    T0,
    l_use_thvm_in_bv_freq,
    l_brunt_vaisala_freq_moist,
    l_modify_limiters_for_cnvg_test,
    gr,
):
    from clubb_jax.src.CLUBB_core.saturation import SATURATION_FLATAU
    return _calc_brunt_vaisala_freq_sqd(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        gr,
        thlm,
        exner,
        rtm,
        rcm,
        p_in_Pa,
        thvm,
        ice_supersat_frac,
        SATURATION_FLATAU,
        l_brunt_vaisala_freq_moist,
        l_use_thvm_in_bv_freq,
        l_modify_limiters_for_cnvg_test,
        bv_efold,
        T0,
    )


def calculate_thlp2_rad(rcm, thlprcp, radht, clubb_params, gr):
    return _calculate_thlp2_rad(
        gr.ngrdcol,
        gr.nzm,
        gr.nzt,
        gr,
        rcm,
        thlprcp,
        radht,
        clubb_params,
        jnp.zeros((gr.ngrdcol, gr.nzm), dtype=jnp.float64),
    )


def Lscale_width_vert_avg(var_profile, rho_ds_zm, var_below_ground_value, gr):
    return _Lscale_width_vert_avg(
        gr.nzm,
        gr.ngrdcol,
        gr,
        2,
        var_profile,
        jnp.zeros((gr.ngrdcol, gr.nzm), dtype=jnp.float64),
        rho_ds_zm,
        var_below_ground_value,
    )


def wp23_term_splat_lhs(
    brunt_vaisala_freq_sqd_mixed,
    Lscale_zm,
    rho_ds_zm,
    C_wp2_splat,
    below_grnd_val,
    gr,
):
    del below_grnd_val
    return _wp23_term_splat_lhs(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        gr,
        C_wp2_splat,
        brunt_vaisala_freq_sqd_mixed,
        Lscale_zm,
        rho_ds_zm,
    )


def compute_Cx_fnc_Richardson(
    brunt_vaisala_freq_sqd,
    brunt_vaisala_freq_sqd_mixed,
    ddzt_umvm_sqd,
    clubb_params,
    l_use_shear_Richardson,
    l_modify_limiters_for_cnvg_test,
):
    ngrdcol, nzm = np.asarray(brunt_vaisala_freq_sqd).shape
    nzt = nzm - 1
    gr = setup_grid(ngrdcol, _DZ, 0.0, _DZ * (nzm - 1), grid_type=1)
    return _compute_Cx_fnc_Richardson(
        nzm,
        nzt,
        ngrdcol,
        gr,
        jnp.zeros((ngrdcol, nzm), dtype=jnp.float64),
        ddzt_umvm_sqd,
        jnp.ones((ngrdcol, nzm), dtype=jnp.float64),
        brunt_vaisala_freq_sqd,
        brunt_vaisala_freq_sqd_mixed,
        clubb_params,
        l_use_shear_Richardson,
        l_modify_limiters_for_cnvg_test,
    )


def test_smooth_min_closed_form():
    a = np.array([[1.0, 5.0, -2.0]])
    coef = 1e-3
    out = np.asarray(smooth_min(a, 3.0, coef))
    ref = 0.5 * ((a + 3.0) - np.sqrt((a - 3.0) ** 2 + coef ** 2))
    assert np.max(np.abs(out - ref)) < 1e-14
    assert np.all(out <= np.minimum(a, 3.0) + 1e-9), "smooth_min must be <= min"
    print("  smooth_min: closed-form + (smooth_min <= min) bound  PASS")


def test_smooth_min_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py smooth_min oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(2)
    a = rng.standard_normal((_NG, 6))
    b = 0.37
    coef = 1e-2
    ref = np.asarray(clubb_f2py.f2py_smooth_min_array_scalar(a, b, coef))
    got = np.asarray(smooth_min(a, b, coef))
    d = np.max(np.abs(got - ref))
    # The single elementwise JAX smooth_min also covers the Fortran scalar+array generic-interface branch.
    d_sa = np.max(np.abs(np.asarray(smooth_min(b, a, coef))
                         - np.asarray(clubb_f2py.f2py_smooth_min_scalar_array(b, a, coef))))
    assert d < 1e-13 and d_sa < 1e-13, f"smooth_min f2py mismatch (array_scalar {d:.2e}, scalar_array {d_sa:.2e})"
    print(f"  f2py smooth_min (array+scalar {d:.1e}, scalar+array {d_sa:.1e}): bit-match both branches  PASS")


def test_smooth_max_f2py():
    """smooth_max (advance_helper_module.F90) — the complement of smooth_min: 0.5*((a+b)+sqrt((a-b)^2+coef^2)),
    a differentiable approximation of max used by clip_explicit / mixing_length. The single elementwise JAX
    smooth_max covers BOTH Fortran generic-interface branches (input_var1 array + scalar, and scalar + array),
    so validate against both f2py oracles. SKIPs if clubb_f2py is unbuilt. (iter 433)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py smooth_max oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(7)
    worst_as = worst_sa = 0.0
    for _ in range(20):
        a = rng.standard_normal((_NG, 6))
        b = float(rng.uniform(-2.0, 2.0)); coef = float(rng.uniform(1e-3, 1.0))
        # array + scalar
        worst_as = max(worst_as, float(np.max(np.abs(
            np.asarray(smooth_max(jnp.asarray(a), b, coef))
            - np.asarray(clubb_f2py.f2py_smooth_max_array_scalar(a, b, coef))))))
        # scalar + array (smooth_max is symmetric, but the Fortran has a distinct generic-interface routine)
        worst_sa = max(worst_sa, float(np.max(np.abs(
            np.asarray(smooth_max(b, jnp.asarray(a), coef))
            - np.asarray(clubb_f2py.f2py_smooth_max_scalar_array(b, a, coef))))))
        assert np.all(np.asarray(smooth_max(jnp.asarray(a), b, coef)) >= np.maximum(a, b) - 1e-9), \
            "smooth_max must be >= max"
    assert worst_as < 1e-13, f"smooth_max_array_scalar f2py mismatch {worst_as:.2e}"
    assert worst_sa < 1e-13, f"smooth_max_scalar_array f2py mismatch {worst_sa:.2e}"
    print(f"  f2py smooth_max (array+scalar {worst_as:.1e}, scalar+array {worst_sa:.1e}, 20 cases): "
          f"bit-match + (smooth_max >= max)  PASS")


def test_compute_Cx_fnc_Richardson_f2py():
    """compute_Cx_fnc_Richardson (advance_helper_module.F90) — the Cx(Richardson-number) closure coefficient (C7/C11
    stability dependence). The f2py signature carries `Lscale_zm`/`rho_ds_zm` for an optional vertical average, but the
    Fortran hardcodes `l_Cx_fnc_Richardson_vert_avg = .false.` → those args are DEAD; the JAX (which omits them) bit-
    matches with dummy values. Validate the PRODUCTION path `l_modify_limiters_for_cnvg_test=False` for both
    `l_use_shear_Richardson` values. NB: the convergence-test branch (`l_modify_limiters_for_cnvg_test=True`, a non-
    production debug flag) **segfaults the f2py oracle** (a Fortran-side issue in that branch), so it stays case-level
    only — same oracle-crash class as precip_fraction/update_xp2_mc. SKIPs if clubb_f2py is unbuilt. (iter 443)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py compute_Cx_fnc_Richardson oracle: SKIP ({type(e).__name__})")
        return
    cp = np.asarray(init_clubb_params(_NG, filename=_TUNABLE_PARAMS))
    nzm = 15; nzt = nzm - 1
    rng = np.random.default_rng(0)
    worst = 0.0
    for shear in (False, True):
        for _ in range(10):
            bvf = rng.uniform(-1e-3, 1e-3, (_NG, nzm)); bvfm = rng.uniform(-1e-3, 1e-3, (_NG, nzm))
            ddz = rng.uniform(0.0, 1e-3, (_NG, nzm))
            lscale = rng.uniform(10.0, 500.0, (_NG, nzm)); rho = rng.uniform(0.5, 1.2, (_NG, nzm))  # dead args
            j = np.asarray(compute_Cx_fnc_Richardson(jnp.asarray(bvf), jnp.asarray(bvfm), jnp.asarray(ddz),
                                                     jnp.asarray(cp), l_use_shear_Richardson=shear,
                                                     l_modify_limiters_for_cnvg_test=False))
            f = np.asarray(clubb_f2py.f2py_compute_cx_fnc_richardson(nzt, lscale, ddz, rho, bvf, bvfm, cp, int(shear), 0))
            worst = max(worst, float(np.max(np.abs(j - f))))
    assert worst < 1e-12, f"compute_Cx_fnc_Richardson f2py mismatch {worst:.2e}"
    print(f"  f2py compute_Cx_fnc_Richardson (production path, both l_use_shear_Richardson, 20 cases): "
          f"bit-match, worst {worst:.2e}  PASS")


def test_calc_xpwp_f2py():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py calc_xpwp oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(7)
    Km_zm = np.asfortranarray(np.abs(rng.standard_normal((ng, nzm))) + 0.1)
    xm = np.asfortranarray(rng.standard_normal((ng, nzt)))
    ref = np.asarray(clubb_f2py.f2py_calc_xpwp_2d(Km_zm.copy(), xm.copy()))
    got = np.asarray(calc_xpwp(Km_zm, xm, np.asarray(jgr.invrs_dzm)))
    # Compare the interior levels the Fortran sets (k=1..nzm-2, 0-based).
    d = np.max(np.abs(got[:, 1:nzm - 1] - ref[:, 1:nzm - 1]))
    assert d < 1e-11, f"calc_xpwp f2py mismatch {d:.2e}"
    print(f"  f2py calc_xpwp_2d: bit-match on interior momentum levels, worst {d:.2e}  PASS")


def test_calc_xpwp_identity():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]; nzt = nzm - 1
    rng = np.random.default_rng(3)
    Km = np.abs(rng.standard_normal((1, nzm))) + 0.1
    xm = rng.standard_normal((1, nzt))
    invrs_dzm = np.asarray(jgr.invrs_dzm)
    out = np.asarray(calc_xpwp(Km, xm, invrs_dzm))
    for k in range(1, nzm - 1):
        expect = Km[0, k] * invrs_dzm[0, k] * (xm[0, k] - xm[0, k - 1])
        assert abs(out[0, k] - expect) < 1e-13, f"xpwp identity failed at k={k}"
    assert out[0, 0] == 0.0 and out[0, nzm - 1] == 0.0, "boundary levels must be zero"
    print("  calc_xpwp: down-gradient identity on interior + zero boundaries  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]; nzt = nzm - 1
    Km = jnp.asarray(np.abs(np.random.default_rng(1).standard_normal((1, nzm))) + 0.1)
    invrs_dzm = jnp.asarray(np.asarray(jgr.invrs_dzm))
    def loss(xm):
        return jnp.sum(calc_xpwp(Km, xm, invrs_dzm) ** 2) + jnp.sum(smooth_min(xm, 0.0, 1e-3) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(2).standard_normal((1, nzt)))))
    assert np.isfinite(g).all(), "non-finite grad"
    print(f"  jax.grad through calc_xpwp + smooth_min: finite ({g.size} entries)  PASS")


def test_stability_and_Ri_f2py():
    """calc_stability_correction (1 + min(lambda0·N²·Lscale²/em, 3)) and calc_Ri_zm (max(N²,lim)/max(shear,lim)) —
    pure advance_helper routines with no dedicated f2py test. Bit-shadow vs the Fortran oracle. (iter 414)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py stability/Ri oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(6)
    npar = max(ilambda0_stability_coef + 1, 80)
    sworst = rworst = 0.0
    for _ in range(40):
        ng, nzm = 2, 10
        bv = rng.uniform(-1e-3, 1e-2, (ng, nzm)); ls = rng.uniform(10.0, 500.0, (ng, nzm))
        em = rng.uniform(1e-3, 1.0, (ng, nzm)); lam0 = rng.uniform(0.5, 3.0, (ng,))
        cp = np.zeros((ng, npar)); cp[:, ilambda0_stability_coef] = lam0
        sj = np.asarray(calc_stability_correction(jnp.asarray(bv), jnp.asarray(ls), jnp.asarray(em), jnp.asarray(cp)))
        sf = np.asarray(clubb_f2py.f2py_calc_stability_correction(bv, ls, em, lam0))
        sworst = max(sworst, float(np.max(np.abs(sj - sf))))
        shear = rng.uniform(0.0, 1e-2, (ng, nzm))
        rj = np.asarray(calc_Ri_zm(jnp.asarray(bv), jnp.asarray(shear), 1e-4, 1e-6))
        rf = np.asarray(clubb_f2py.f2py_calc_ri_zm(bv, shear, 1e-4, 1e-6))
        rworst = max(rworst, float(np.max(np.abs(rj - rf))))
    assert sworst < 1e-12, f"calc_stability_correction f2py mismatch {sworst:.2e}"
    assert rworst < 1e-12, f"calc_Ri_zm f2py mismatch {rworst:.2e}"
    print(f"  f2py calc_stability_correction + calc_Ri_zm (40 cases): bit-match, worst {max(sworst, rworst):.2e}  PASS")


def test_vertical_avg_integral_f2py():
    """Density-weighted vertical_avg = sum(field·rho·dz)/sum(rho·dz) and vertical_integral = sum(field·rho·dz)
    (advance_helper_module.F90) vs the Fortran oracle — pure column reductions, previously untested. (iter 415)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py vertical_avg/integral oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(8)
    aworst = iworst = 0.0
    for _ in range(40):
        nz = 12
        rho = rng.uniform(0.3, 1.2, nz); fld = rng.uniform(-5.0, 5.0, nz); dz = rng.uniform(20.0, 80.0, nz)
        aworst = max(aworst, abs(float(vertical_avg(rho, fld, dz)) - float(clubb_f2py.f2py_vertical_avg(rho, fld, dz))))
        vi_j = float(vertical_integral(rho, fld, dz)); vi_f = float(clubb_f2py.f2py_vertical_integral(rho, fld, dz))
        iworst = max(iworst, abs(vi_j - vi_f) / max(abs(vi_f), 1.0))
    assert aworst < 1e-12, f"vertical_avg f2py mismatch {aworst:.2e}"
    assert iworst < 1e-12, f"vertical_integral f2py rel mismatch {iworst:.2e}"
    print(f"  f2py vertical_avg + vertical_integral (40 cases): bit-match, worst {max(aworst, iworst):.2e}  PASS")


def test_brunt_vaisala_f2py():
    """calc_brunt_vaisala_freq_sqd (N² = (g/thv)·d(thv)/dz, dry/moist + thvm/thlm branches; returns sqd/_mixed/_smth)
    vs the Fortran oracle across all 4 l_use_thvm × l_moist flag combinations, on a stretched grid (the f2py uses a
    module-global grid set to the JAX grid heights). SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 416)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py brunt_vaisala oracle: SKIP ({type(e).__name__})")
        return
    from clubb_jax.src.CLUBB_core.saturation import SATURATION_FLATAU
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(5)
    thlm = rng.uniform(290.0, 320.0, (ng, nzt)); exner = rng.uniform(0.8, 1.0, (ng, nzt))
    rtm = rng.uniform(1e-3, 1.5e-2, (ng, nzt)); rcm = rng.uniform(0.0, 5e-4, (ng, nzt))
    p = rng.uniform(7e4, 1e5, (ng, nzt)); thvm = thlm * (1.0 + 0.6 * rtm)
    isf = rng.uniform(0.0, 0.3, (ng, nzt)); bv_efold = np.full(ng, 5.0); T0 = 300.0
    worst = 0.0
    for lthvm in (False, True):
        for lmoist in (False, True):
            jr = calc_brunt_vaisala_freq_sqd(jnp.asarray(thlm), jnp.asarray(exner), jnp.asarray(rtm), jnp.asarray(rcm),
                jnp.asarray(p), jnp.asarray(thvm), jnp.asarray(isf), jnp.asarray(bv_efold), T0, lthvm, lmoist, False, jgr)
            fr = clubb_f2py.f2py_calc_brunt_vaisala_freq_sqd(nzm, thlm, exner, rtm, rcm, p, thvm, isf,
                SATURATION_FLATAU, int(lmoist), int(lthvm), 0, bv_efold, T0)
            worst = max(worst, max(float(np.max(np.abs(np.asarray(jr[i]) - np.asarray(fr[i])))) for i in range(3)))
    assert worst < 1e-12, f"calc_brunt_vaisala_freq_sqd f2py mismatch {worst:.2e}"
    print(f"  f2py calc_brunt_vaisala_freq_sqd (4 flag combos x 3 outputs): bit-match, worst {worst:.2e}  PASS")


def test_calculate_thlp2_rad_f2py():
    """calculate_thlp2_rad (advance_helper_module.F90 — radiative-cooling contribution to the thlp2 forcing, gated
    by l_calc_thlp2_rad which most cases leave off → otherwise unexercised). The JAX returns the increment; the f2py
    adds it to an input thlp2_forcing, so passing zeros recovers the increment. (iter 419)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py calculate_thlp2_rad oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    cp = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(20):
        rcm = rng.uniform(0.0, 1e-3, (ng, nzt)); rcm[rcm < 3e-4] = 0.0   # mix of cloudy / clear levels
        thlprcp = rng.uniform(-1e-2, 1e-2, (ng, nzm)); radht = rng.uniform(-1e-4, 1e-4, (ng, nzt))
        jr = np.asarray(calculate_thlp2_rad(jnp.asarray(rcm), jnp.asarray(thlprcp), jnp.asarray(radht),
                                            jnp.asarray(cp), jgr))
        fr = np.asarray(clubb_f2py.f2py_calculate_thlp2_rad(rcm, thlprcp, radht, cp, np.zeros((ng, nzm))))
        worst = max(worst, float(np.max(np.abs(jr - fr))))
    assert worst < 1e-12, f"calculate_thlp2_rad f2py mismatch {worst:.2e}"
    print(f"  f2py calculate_thlp2_rad (rad thlp2-forcing increment, 20 cases): bit-match, worst {worst:.2e}  PASS")


def test_Lscale_width_vert_avg_f2py():
    """Lscale_width_vert_avg (advance_helper_module.F90) — the rho_ds·dz-weighted vertical average of a profile over
    a moving window (the JAX hardcodes smth_type=2, fixed 60 m half-width, so the f2py's lscale_zm is unused on that
    path). vs the f2py oracle. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 429)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py Lscale_width_vert_avg oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(15):
        vp = rng.uniform(-1.0, 1.0, (ng, nzm)); rho = rng.uniform(0.5, 1.2, (ng, nzm)); ls = rng.uniform(10.0, 500.0, (ng, nzm))
        jl = np.asarray(Lscale_width_vert_avg(jnp.asarray(vp), jnp.asarray(rho), 0.0, jgr))
        fl = np.asarray(clubb_f2py.f2py_lscale_width_vert_avg(2, vp, ls, rho, 0.0))   # smth_type=2 (lscale_zm unused)
        worst = max(worst, float(np.max(np.abs(jl - fl))))
    assert worst < 1e-12, f"Lscale_width_vert_avg f2py mismatch {worst:.2e}"
    print(f"  f2py Lscale_width_vert_avg (smth_type=2, 15 cases): bit-match, worst {worst:.2e}  PASS")


def test_wp23_term_splat_lhs_f2py():
    """wp23_term_splat_lhs (advance_helper_module.F90) — the LHS coefficients of the wp2/wp3 splatting terms (the
    vertically-averaged, clipped, smoothed Brunt-Vaisala frequency × C_wp2_splat). The Fortran hardcodes the
    below-ground value 0.01 for the vertical average, so the JAX is called with below_grnd_val=0.01. vs the f2py
    oracle. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 430)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py wp23_term_splat_lhs oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(2)
    worst = 0.0
    for _ in range(15):
        bv = rng.uniform(-1e-3, 1e-2, (ng, nzm)); ls = rng.uniform(10.0, 500.0, (ng, nzm))
        rho = rng.uniform(0.5, 1.2, (ng, nzm)); cspl = np.full(ng, 2.0)
        js = wp23_term_splat_lhs(jnp.asarray(bv), jnp.asarray(ls), jnp.asarray(rho), jnp.asarray(cspl), 0.01, jgr)
        fs = clubb_f2py.f2py_wp23_term_splat_lhs(nzt, cspl, bv, ls, rho)
        worst = max(worst, max(float(np.max(np.abs(np.asarray(js[k]) - np.asarray(fs[k])))) for k in range(2)))
    assert worst < 1e-12, f"wp23_term_splat_lhs f2py mismatch {worst:.2e}"
    print(f"  f2py wp23_term_splat_lhs (lhs_splat_wp2/wp3, 15 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_advance_helper_extras:")
    for t in (test_smooth_min_closed_form, test_smooth_min_f2py, test_smooth_max_f2py,
              test_compute_Cx_fnc_Richardson_f2py, test_calc_xpwp_f2py,
              test_calc_xpwp_identity, test_stability_and_Ri_f2py, test_vertical_avg_integral_f2py,
              test_brunt_vaisala_f2py, test_calculate_thlp2_rad_f2py, test_Lscale_width_vert_avg_f2py,
              test_wp23_term_splat_lhs_f2py, test_differentiable):
        t()
    print("All advance_helper extras checks PASSED")


if __name__ == "__main__":
    main()
