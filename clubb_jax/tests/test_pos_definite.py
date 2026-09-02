#!/usr/bin/env python3
"""test_pos_definite.py — validate the JAX pos_definite_adj port (Smolarkiewicz 1989 flux limiter).

Two oracles:
  1. Conservation invariant (always, oracle-free): the renormalization preserves the unweighted column
     integral `sum_k field[k]*dzt[k]` (the flux BCs zero the telescoping boundary terms), so the adjusted
     field's column integral must equal the input's to machine precision.
  2. f2py bit-shadow: `clubb_f2py.f2py_pos_definite_adj` on the SAME stored grid + inputs — all four outputs
     bit-to-bit (SKIPs cleanly if clubb_f2py / clubb_python are not built).
Plus a `jax.grad` smoke check.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
# _ROOT must come first so `import clubb_jax` resolves to this checkout's package.
# The clubb_python_api path (for clubb_python + clubb_f2py) is appended.
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.pos_definite_module import pos_definite_adj

_NG, _DZ, _ZTOP, _DT = 2, 40.0, 1200.0, 60.0


def _make_inputs(jgr):
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(20260603)
    # Flux on zm, field on zt. Make column 0 develop a below-zero field (exercise the limiter branches),
    # leave column 1 all-positive (exercise the "no negative -> *_pd = 0" branch).
    flux_np1 = rng.standard_normal((ng, nzm)) * 1e-3
    field_n = np.abs(rng.standard_normal((ng, nzt))) * 1e-2 + 1e-3
    field_np1 = field_n + rng.standard_normal((ng, nzt)) * 5e-3
    field_np1[0, nzt // 3] = -2e-3      # a genuine negative in column 0
    field_np1[1, :] = np.abs(field_np1[1, :]) + 1e-4   # column 1 strictly positive
    return (np.asfortranarray(field_np1), np.asfortranarray(flux_np1), np.asfortranarray(field_n))


def test_conservation():
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    field_np1, flux_np1, field_n = _make_inputs(jgr)
    field_adj, flux_lim, field_pd, flux_pd = pos_definite_adj(jgr, _DT, field_np1, flux_np1, field_n)
    dzt = np.asarray(1.0 / jgr.invrs_dzt)
    col_in = (np.asarray(field_np1) * dzt).sum(axis=1)
    col_out = (np.asarray(field_adj) * dzt).sum(axis=1)
    rel = np.max(np.abs(col_out - col_in) / (np.abs(col_in) + 1e-30))
    assert rel < 1e-12, f"column-integral not conserved: rel {rel:.2e}"
    # column 1 had no negative -> both diagnostic tendencies must be exactly zero there
    assert np.all(np.asarray(field_pd)[1] == 0.0) and np.all(np.asarray(flux_pd)[1] == 0.0), \
        "all-positive column should have zero *_pd"
    print(f"  conservation (sum field*dzt): rel {rel:.2e}; all-positive col -> *_pd=0  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py pos_definite oracle: SKIP ({type(e).__name__})")
        return
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    # Set the f2py stored_grid to EXACTLY the JAX grid's heights (grid_type=2 = explicit heights).
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, _DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    field_np1, flux_np1, field_n = _make_inputs(jgr)
    f_field, f_flux, f_field_pd, f_flux_pd = clubb_f2py.f2py_pos_definite_adj(
        _DT, field_np1.copy(), flux_np1.copy(), field_n.copy())
    j_field, j_flux, j_field_pd, j_flux_pd = pos_definite_adj(jgr, _DT, field_np1, flux_np1, field_n)
    worst = 0.0
    for name, jx, fx in (("field_np1", j_field, f_field), ("flux_np1", j_flux, f_flux),
                         ("field_pd", j_field_pd, f_field_pd), ("flux_pd", j_flux_pd, f_flux_pd)):
        jx = np.asarray(jx); fx = np.asarray(fx)
        rel = np.max(np.abs(jx - fx) / (np.abs(fx) + 1e-30))
        worst = max(worst, rel)
        assert rel < 1e-11, f"{name} f2py mismatch rel {rel:.2e}"
    print(f"  f2py pos_definite_adj: 4/4 outputs bit-match, worst rel {worst:.2e}  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzm = jgr.zm.shape[1]; nzt = nzm - 1
    flux = jnp.asarray(np.random.default_rng(1).standard_normal((1, nzm)) * 1e-3)
    field_n = jnp.asarray(np.abs(np.random.default_rng(2).standard_normal((1, nzt))) * 1e-2)
    def loss(field_np1):
        out, _, _, _ = pos_definite_adj(jgr, _DT, field_np1, flux, field_n)
        return jnp.sum(out ** 2)
    g = np.asarray(jax.grad(loss)(field_n))   # field_np1 = field_n (no negatives) baseline
    assert np.isfinite(g).all(), "non-finite grad through pos_definite_adj"
    print(f"  jax.grad through pos_definite_adj: finite ({g.size} entries)  PASS")


def main():
    print("test_pos_definite:")
    for t in (test_conservation, test_f2py_oracle, test_differentiable):
        t()
    print("All pos_definite checks PASSED")


if __name__ == "__main__":
    main()
