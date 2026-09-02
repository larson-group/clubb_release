#!/usr/bin/env python3
"""test_grid_operators.py — validate the JAX vertical grid operators vs the f2py Fortran oracle.

These are the model's most fundamental routines (every vertical derivative and zt<->zm regrid in the closure
goes through them), in CLUBB_core/grid_class.py ↔ grid_class.F90:

  ddzt   (ddzt_2d)      — d/dz of a zt-field onto zm levels (invrs_dzm-weighted difference)
  ddzm   (ddzm_2d)      — d/dz of a zm-field onto zt levels (invrs_dzt-weighted difference)
  zt2zm  (zt2zm_2d)     — linear-interpolate a zt-field onto zm levels
  zm2zt  (zm2zt_2d)     — linear-interpolate a zm-field onto zt levels
  zt2zm2zt (zt2zm2zt_2d) — the zt→zm→zt smoothing round-trip (optional lower clip zt_min)
  zm2zt2zm (zm2zt2zm_2d) — the zm→zt→zm smoothing round-trip (optional lower clip zm_min)

The JAX operators take the grid object (interp weights precomputed at setup_grid); the f2py uses the module-global
grid set to the same JAX heights (the test_hydrostatic pattern). SKIPs cleanly if clubb_f2py / clubb_python are unbuilt.
(iter 436)
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

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.grid_class import (
    ddzt, ddzm, zt2zm, zm2zt, zt2zm2zt, zm2zt2zm,
)

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


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


def test_f2py_oracle():
    try:
        clubb_f2py, jgr, ng, nzm = _setup_f2py_grid()
    except Exception as e:
        print(f"  f2py grid-operator oracle: SKIP ({type(e).__name__})")
        return
    nzt = nzm - 1
    rng = np.random.default_rng(0)
    worst = {}
    # zt-input operators: the leading f2py size arg is nzm
    for nm, jf, ff in (
        ("ddzt", lambda a: ddzt(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_ddzt_2d),
        ("zt2zm", lambda a: zt2zm(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_zt2zm_2d),
    ):
        w = 0.0
        for _ in range(15):
            a = rng.standard_normal((ng, nzt))
            w = max(w, float(np.max(np.abs(np.asarray(jf(jnp.asarray(a))) - np.asarray(ff(nzm, a))))))
        worst[nm] = w
    # zm-input operators: the leading f2py size arg is nzt
    for nm, jf, ff in (
        ("ddzm", lambda a: ddzm(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_ddzm_2d),
        ("zm2zt", lambda a: zm2zt(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_zm2zt_2d),
    ):
        w = 0.0
        for _ in range(15):
            a = rng.standard_normal((ng, nzm))
            w = max(w, float(np.max(np.abs(np.asarray(jf(jnp.asarray(a))) - np.asarray(ff(nzt, a))))))
        worst[nm] = w
    # round-trip smoothers, with no-clip (-1e30) and a finite lower clip (0.0)
    for zmin in (-1e30, 0.0):
        wt = wm = 0.0
        for _ in range(15):
            at = rng.standard_normal((ng, nzt)); am = rng.standard_normal((ng, nzm))
            wt = max(wt, float(np.max(np.abs(np.asarray(zt2zm2zt(nzm, nzt, ng, jgr, jnp.asarray(at), zt_min=zmin))
                                             - np.asarray(clubb_f2py.f2py_zt2zm2zt_2d(nzm, at, zmin))))))
            wm = max(wm, float(np.max(np.abs(np.asarray(zm2zt2zm(nzm, nzt, ng, jgr, jnp.asarray(am), zm_min=zmin))
                                             - np.asarray(clubb_f2py.f2py_zm2zt2zm_2d(nzt, am, zmin))))))
        worst[f"zt2zm2zt(min={zmin:g})"] = wt
        worst[f"zm2zt2zm(min={zmin:g})"] = wm
    overall = max(worst.values())
    assert overall < 1e-12, f"grid-operator f2py mismatch: {worst}"
    print(f"  f2py ddzt/ddzm/zt2zm/zm2zt + zt2zm2zt/zm2zt2zm round-trips (no-clip + clip): "
          f"bit-match, worst {overall:.2e}  PASS")


def test_stretched_grid_operators():
    """The iter-436 operator check used an EVENLY-SPACED grid (constant invrs_dz → trivial 0.5 interp weights). A
    STRETCHED grid (non-uniform dz, like the CGILS cases) exercises the invrs_dz-weighted derivative + the
    relative-spacing interpolation weights for real. Build the SAME stretched grid_type=2 (zt-input) on JAX + f2py and
    re-validate ddzt/ddzm/zt2zm/zm2zt. NB: both sides MUST use the same grid_type — JAX grid_type=3 (zm-input) vs f2py
    grid_type=2 (zt-input) defines the grid differently (zt-from-zm ≠ inverse of zm-from-zt on a stretched grid), a
    setup trap that produces a spurious ~0.1 mismatch. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 481)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py stretched-grid operators: SKIP ({type(e).__name__})")
        return
    nzt0 = 18
    dz = np.linspace(20.0, 120.0, nzt0)
    zt1 = np.concatenate([[10.0], 10.0 + np.cumsum(dz)])[:nzt0]          # non-uniform thermodynamic heights
    zt_in = np.tile(zt1, (_NG, 1))
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=float(zt1[-1] + 50.0),
                     grid_type=2, thermodynamic_heights=zt_in)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    # guard the setup: the grids must actually coincide before comparing operators
    g = clubb_f2py.get_grid(ng, nzm, nzt)
    assert float(np.max(np.abs(np.asarray(g[0]) - np.asarray(jgr.zm)))) == 0.0, "stretched grid zm mismatch (setup)"
    rng = np.random.default_rng(0)
    worst = 0.0
    for nm, jf, ff, on_zt in (
            ("ddzt", lambda a: ddzt(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_ddzt_2d, True),
            ("zt2zm", lambda a: zt2zm(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_zt2zm_2d, True),
            ("ddzm", lambda a: ddzm(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_ddzm_2d, False),
            ("zm2zt", lambda a: zm2zt(nzm, nzt, ng, jgr, a), clubb_f2py.f2py_zm2zt_2d, False)):
        for _ in range(8):
            a = rng.standard_normal((ng, nzt if on_zt else nzm))
            worst = max(worst, float(np.max(np.abs(
                np.asarray(jf(jnp.asarray(a))) - np.asarray(ff(nzm if on_zt else nzt, a))))))
    assert worst < 1e-12, f"stretched-grid operators differ from f2py: worst {worst:.2e}"
    print(f"  f2py ddzt/ddzm/zt2zm/zm2zt on a STRETCHED grid (non-uniform dz): bit-match, worst {worst:.2e}  PASS")


def test_grid_type3_construction_matches_fortran():
    """`setup_grid(grid_type=3)` is the stretched zm-input path: given the momentum heights it COMPUTES the
    thermodynamic heights (zt-from-zm) — a distinct construction from grid_type=1 (even) and grid_type=2 (zt-input → zm,
    iters 481/482). It's case-validated (dycoms2_rf01_fixed_sst/dycoms2_rf02_so/neutral are grid_type=3 + bit-faithful)
    but had no direct unit test. Here the Fortran independently runs grid_type=3 on the same (stretched) zm input and the
    constructed zt/zm match exactly. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 484)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py grid_type=3 construction: SKIP ({type(e).__name__})")
        return
    nzm0 = 20
    dz = np.linspace(20.0, 120.0, nzm0 - 1)
    zm1 = np.concatenate([[0.0], np.cumsum(dz)])                 # non-uniform momentum heights
    zm_in = np.tile(zm1, (_NG, 1))
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=float(zm1[-1]),
                     grid_type=3, momentum_heights=zm_in)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=3, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(zm1[-1])), momentum_heights=np.asfortranarray(zm_in),
                         thermodynamic_heights=np.asfortranarray(np.zeros((ng, nzt))), err_info=ErrInfo(ngrdcol=ng))
    g = clubb_f2py.get_grid(ng, nzm, nzt)
    dzm = float(np.max(np.abs(np.asarray(g[0]) - np.asarray(jgr.zm))))
    dzt = float(np.max(np.abs(np.asarray(g[1]) - np.asarray(jgr.zt))))
    assert dzm == 0.0 and dzt == 0.0, f"grid_type=3 construction differs: zm {dzm:.2e}, zt {dzt:.2e}"
    print(f"  f2py grid_type=3 construction (zt-from-zm): JAX zt/zm match Fortran exactly (diff 0.0)  PASS")


def test_descending_grid_rejected():
    """`setup_grid` fail-loud rejects `l_ascending_grid=False`: the JAX grid operators (zt2zm boundary handling) are
    ascending-only, so a descending grid would silently mis-compute the two boundary levels (verified iter 483 — interior
    interp is correct, both boundaries diverge ~0.5·Δfield from the Fortran). No case uses a descending grid; the guard
    makes the unsupported case explicit instead of a silent footgun. (iter 483)"""
    try:
        setup_grid(ngrdcol=2, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1, l_ascending_grid=False)
    except ValueError as e:
        assert "ascending" in str(e).lower(), f"unexpected rejection message: {e}"
        print("  setup_grid fail-loud rejects descending grid (l_ascending_grid=False)  PASS")
        return
    raise AssertionError("setup_grid did NOT reject l_ascending_grid=False (descending grid)")


def test_setup_grid_is_sole_grid_constructor():
    """Source-grounded guardrail for the iter-483 descending-grid guard's COVERAGE. That guard lives in `setup_grid`
    (CLUBB_core/grid_class.py) — it only protects the downstream ascending-only assumptions (zt2zm boundaries,
    pos_definite_adj, sponge_layer_damping, precipitation_fraction, fill_holes, advance_xm_wpxp/xp2_xpyp, diffusion,
    remapping) IF `setup_grid` is the SINGLE place a `Grid` is constructed. Verified iter 485 by grepping src: the only
    `Grid(...)` construction is the `return Grid(...)` inside `setup_grid`. This test pins that invariant so a future
    refactor adding a direct `Grid(...)` (bypassing the guard, re-opening the descending-grid footgun) fails loudly here.
    Pure source scan — no f2py, never SKIPs. (iter 485)"""
    import re
    src_root = os.path.join(_ROOT, "clubb_jax", "src")
    # Grid construction is allowed only in setup_grid.
    allowed = {
        os.path.normpath(os.path.join(src_root, "CLUBB_core", "grid_class.py")),
    }
    offenders = []
    for dirpath, _dirs, files in os.walk(src_root):
        for fn in files:
            if not fn.endswith(".py"):
                continue
            path = os.path.join(dirpath, fn)
            with open(path) as fh:
                for lineno, line in enumerate(fh, 1):
                    code = line.split("#", 1)[0]  # strip trailing comments
                    # a Grid(...) construction: `Grid(` preceded by `=`, `return `, `(`, etc. — NOT `class Grid(`
                    if re.search(r"(?<![A-Za-z_.])Grid\s*\(", code) and "class Grid(" not in code:
                        if os.path.normpath(path) not in allowed:
                            offenders.append(f"{os.path.relpath(path, _ROOT)}:{lineno}: {line.strip()}")
    assert not offenders, (
        "Grid(...) constructed outside setup_grid (CLUBB_core/grid_class.py) — the iter-483 descending-grid guard "
        "no longer covers all construction sites:\n  " + "\n  ".join(offenders))
    print("  setup_grid (grid_class.py) is the sole Grid constructor in src — iter-483 guard coverage intact  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zt.shape[1]
    rng = np.random.default_rng(2)
    a = jnp.asarray(rng.standard_normal((1, nzt)))
    # a chain through both a derivative and a regrid → finite grad
    g = jax.grad(
        lambda v: (
            jnp.sum(ddzt(jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr, v) ** 2)
            + jnp.sum(zt2zm(jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr, v) ** 2)
        )
    )(a)
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad through ddzt/zt2zm"
    print(f"  jax.grad through ddzt + zt2zm: finite ({g.size} entries)  PASS")


def test_grid_construction_matches_fortran():
    """The JAX `setup_grid(grid_type=1)` builds the evenly-spaced grid HEIGHTS (zt/zm from deltaz/zm_top) — distinct from
    test_f2py_oracle, which validates the *operators* given heights the f2py wrapper takes from the JAX (grid_type=2). Here
    the Fortran independently **constructs** the grid_type=1 grid from the same deltaz/zm_top (heights args ignored), and
    we compare the resulting zt/zm — so the JAX height computation itself is pinned to the oracle, not just the weights.
    SKIPs if clubb_f2py / clubb_python are unbuilt. (iter 480)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py grid-construction oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape; nzt = nzm - 1
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=1, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, _ZTOP),  # grid_type=1 → Fortran constructs; heights args ignored
                         momentum_heights=np.asfortranarray(np.zeros((ng, nzm))),
                         thermodynamic_heights=np.asfortranarray(np.zeros((ng, nzt))), err_info=ErrInfo(ngrdcol=ng))
    g = clubb_f2py.get_grid(ng, nzm, nzt)
    fzm, fzt = np.asarray(g[0]), np.asarray(g[1])
    dzm = float(np.max(np.abs(fzm - np.asarray(jgr.zm))))
    dzt = float(np.max(np.abs(fzt - np.asarray(jgr.zt))))
    assert dzm == 0.0 and dzt == 0.0, f"grid_type=1 construction differs from Fortran: zm {dzm:.2e}, zt {dzt:.2e}"
    print(f"  f2py grid_type=1 construction: JAX zt/zm heights match Fortran exactly (diff 0.0)  PASS")


def main():
    print("test_grid_operators:")
    test_f2py_oracle()
    test_grid_construction_matches_fortran()
    test_grid_type3_construction_matches_fortran()
    test_stretched_grid_operators()
    test_descending_grid_rejected()
    test_setup_grid_is_sole_grid_constructor()
    test_differentiable()
    print("All grid-operator checks PASSED")


if __name__ == "__main__":
    main()
