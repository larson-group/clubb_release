"""Check the JAX vertical grid operators against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.grid_class as fortran
from clubb_python.derived_types.err_info import ErrInfo

from clubb_jax.pytests.f2py_unit_test_helpers import init_fortran_state, max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.grid_class import setup_grid, ddzt, ddzm, zt2zm, zm2zt, zt2zm2zt, zm2zt2zm

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


_OPERATORS = (("ddzt", ddzt, True), ("zt2zm", zt2zm, True), ("ddzm", ddzm, False), ("zm2zt", zm2zt, False))


def _operator_mismatch(jgr, fgr, rng, n_draws):
    """Return the worst JAX-vs-Fortran difference of each single-step operator on random inputs."""
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    worst = {}
    for nm, jax_op, on_zt in _OPERATORS:
        fortran_op = getattr(fortran, nm)
        w = 0.0
        for _ in range(n_draws):
            a = rng.standard_normal((ng, nzt if on_zt else nzm))
            w = max(w, float(max_abs_diff("grid_operators", np.asarray(jax_op(nzm, nzt, ng, jgr, jnp.asarray(a))),
                                          np.asarray(fortran_op(nzm, nzt, ng, fgr, a)))))
        worst[nm] = w
    return worst


def test_f2py_oracle():
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(0)
    worst = _operator_mismatch(jgr, fgr, rng, 15)
    # round-trip smoothers, with no-clip (-1e30) and a finite lower clip (0.0)
    for zmin in (-1e30, 0.0):
        wt = wm = 0.0
        for _ in range(15):
            at = rng.standard_normal((ng, nzt)); am = rng.standard_normal((ng, nzm))
            wt = max(wt, float(max_abs_diff("grid_operators", np.asarray(zt2zm2zt(nzm, nzt, ng, jgr, jnp.asarray(at), zt_min=zmin)), np.asarray(fortran.zt2zm2zt(nzm, nzt, ng, fgr, at, zmin)))))
            wm = max(wm, float(max_abs_diff("grid_operators", np.asarray(zm2zt2zm(nzm, nzt, ng, jgr, jnp.asarray(am), zm_min=zmin)), np.asarray(fortran.zm2zt2zm(nzm, nzt, ng, fgr, am, zmin)))))
        worst[f"zt2zm2zt(min={zmin:g})"] = wt
        worst[f"zm2zt2zm(min={zmin:g})"] = wm
    overall = max(worst.values())
    assert overall < 1e-12, f"grid-operator f2py mismatch: {worst}"


def test_stretched_grid_operators():
    """The operators on a STRETCHED grid (non-uniform dz, like the CGILS cases), which exercises the invrs_dz-weighted
    derivative and the relative-spacing interpolation weights. Both sides use grid_type=2 (thermodynamic heights in):
    a JAX grid_type=3 against a Fortran grid_type=2 defines the grid differently (zt-from-zm is not the inverse of
    zm-from-zt on a stretched grid) and produces a spurious ~0.1 mismatch."""
    nzt0 = 18
    dz = np.linspace(20.0, 120.0, nzt0)
    zt1 = np.concatenate([[10.0], 10.0 + np.cumsum(dz)])[:nzt0]          # non-uniform thermodynamic heights
    jgr, fgr = setup_matching_grids(_NG, _DZ, float(zt1[-1] + 50.0), grid_type=2,
                                    thermodynamic_heights=np.tile(zt1, (_NG, 1)))
    # guard the setup: the grids must actually coincide before comparing operators
    assert max_abs_diff("grid_operators", fgr.zm, jgr.zm) == 0.0, "stretched grid zm mismatch (setup)"
    worst = max(_operator_mismatch(jgr, fgr, np.random.default_rng(0), 8).values())
    assert worst < 1e-12, f"stretched-grid operators differ from Fortran: worst {worst:.2e}"


def _fortran_constructed_grid(ng, grid_type, zm_top, momentum_heights, thermodynamic_heights):
    """Let Fortran setup_grid construct the grid itself and return the resulting Grid."""
    nzm = momentum_heights.shape[1]
    init_fortran_state(ng)
    fgr, err = fortran.setup_grid(
        nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False, l_ascending_grid=True,
        grid_type=grid_type, deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, zm_top),
        momentum_heights=np.asfortranarray(momentum_heights),
        thermodynamic_heights=np.asfortranarray(thermodynamic_heights), err_info=ErrInfo(ngrdcol=ng))
    assert not err.is_fatal(), "Fortran setup_grid reported a fatal error"
    return fgr


def test_grid_type3_construction_matches_fortran():
    """`setup_grid(grid_type=3)` is the stretched zm-input path: given the momentum heights it COMPUTES the
    thermodynamic heights (zt-from-zm), a distinct construction from grid_type=1 (even) and grid_type=2 (zt-input).
    Fortran independently runs grid_type=3 on the same stretched zm input and the constructed zt/zm must match."""
    nzm0 = 20
    dz = np.linspace(20.0, 120.0, nzm0 - 1)
    zm1 = np.concatenate([[0.0], np.cumsum(dz)])                 # non-uniform momentum heights
    zm_in = np.tile(zm1, (_NG, 1))
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=float(zm1[-1]),
                     grid_type=3, momentum_heights=zm_in)
    fgr = _fortran_constructed_grid(_NG, 3, float(zm1[-1]), zm_in, np.zeros((_NG, nzm0 - 1)))
    dzm = float(max_abs_diff("grid_operators", fgr.zm, jgr.zm))
    dzt = float(max_abs_diff("grid_operators", fgr.zt, jgr.zt))
    assert dzm == 0.0 and dzt == 0.0, f"grid_type=3 construction differs: zm {dzm:.2e}, zt {dzt:.2e}"


def test_grid_construction_matches_fortran():
    """`setup_grid(grid_type=1)` builds the evenly-spaced grid HEIGHTS (zt/zm from deltaz/zm_top), unlike
    test_f2py_oracle, whose Fortran grid takes its heights from the JAX grid (grid_type=2). Here Fortran
    constructs the grid_type=1 grid from the same deltaz/zm_top (height arguments ignored) and the resulting zt/zm
    must match, so the JAX height computation itself is pinned to Fortran, not just the weights."""
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]
    fgr = _fortran_constructed_grid(_NG, 1, _ZTOP, np.zeros((_NG, nzm)), np.zeros((_NG, nzm - 1)))
    dzm = float(max_abs_diff("grid_operators", fgr.zm, jgr.zm))
    dzt = float(max_abs_diff("grid_operators", fgr.zt, jgr.zt))
    assert dzm == 0.0 and dzt == 0.0, f"grid_type=1 construction differs from Fortran: zm {dzm:.2e}, zt {dzt:.2e}"
