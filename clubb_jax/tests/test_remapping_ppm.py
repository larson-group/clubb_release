#!/usr/bin/env python3
"""test_remapping_ppm.py — validate the JAX PPM (method 2) remapping port (remapping_module.F90).

E3SM Piecewise-Parabolic-Method conservative vertical remap (map1_ppm/ppm2m/steepz/kmppm). Oracles:
  1. f2py bit-shadow vs f2py_remap_vals_to_target_same_grid with grid_remap_method=2 — same grid → identity
     (the integral of a source cell's PPM parabola over the same cell = its mean), validating the map1_ppm
     integration + k0 search end-to-end against the Fortran. SKIPs if unbuilt.
  2. MASS CONSERVATION on a non-trivial (refined) target grid — PPM is conservative by construction, so the
     pressure-weighted column integral must be preserved (a genuine, oracle-free check of the reconstruction).
  3. A finite jax.grad.
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

from clubb_jax.src.CLUBB_core.remapping_module import remap_vals_ppm, remap_vals_to_target_same_grid
from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def _rho_spline(ng, ztop):
    levels = np.linspace(0.0, ztop + 100.0, 24)
    vals = 1.2 * np.exp(-levels / 8000.0)
    return np.tile(vals, (ng, 1)), np.tile(levels, (ng, 1))


def test_f2py_same_grid_identity():
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py PPM same-grid oracle: SKIP ({type(e).__name__})")
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
    rho_vals, rho_levels = _rho_spline(ng, float(jgr.zm[0, -1]))
    p_sfc = np.full(ng, 1.0e5)
    rng = np.random.default_rng(3)
    worst = 0.0
    for iv in (1, 0, -1):
        src = rng.uniform(0.2, 2.0, (ng, nzt))
        f = np.asarray(clubb_f2py.f2py_remap_vals_to_target_same_grid(
            nzm, nzt, src, nzt, rho_vals, rho_levels, iv, p_sfc, 2, True))   # grid_remap_method=2 (PPM)
        g = np.asarray(remap_vals_to_target_same_grid(
            src, np.asarray(jgr.zm), np.asarray(jgr.zt), rho_vals, rho_levels, p_sfc,
            grid_remap_method=2, l_zt_variable=True, iv=iv))
        worst = max(worst, np.max(np.abs(g - f)))
    assert worst < 1e-11, f"PPM same-grid f2py mismatch {worst:.2e}"
    print(f"  f2py PPM same-grid (iv=1,0,-1): bit-match, worst {worst:.2e}  PASS")


def test_mass_conservation_refined():
    # Build a top→surface pressure grid; refine each cell into 2; PPM must conserve the pressure-weighted integral.
    rng = np.random.default_rng(8)
    ncol = 2
    # source edges surface→top (decreasing pressure)
    p_src = np.cumsum(np.concatenate([[1.0e5], -rng.uniform(2000, 4000, 20)]))[None, :].repeat(ncol, 0)
    # refined target: midpoints inserted
    mids = 0.5 * (p_src[:, :-1] + p_src[:, 1:])
    p_tgt = np.zeros((ncol, p_src.shape[1] + (p_src.shape[1] - 1)))
    p_tgt[:, ::2] = p_src
    p_tgt[:, 1::2] = mids
    src = rng.uniform(0.5, 3.0, (ncol, p_src.shape[1] - 1))
    tgt = np.asarray(remap_vals_ppm(p_src, p_tgt, src, iv=1))
    # pressure-weighted integral (mass) preserved: sum(src*dp_src) == sum(tgt*dp_tgt)
    dp_src = np.abs(np.diff(p_src, axis=1)); dp_tgt = np.abs(np.diff(p_tgt, axis=1))
    m_src = np.sum(src * dp_src, axis=1); m_tgt = np.sum(tgt * dp_tgt, axis=1)
    rel = np.max(np.abs(m_tgt - m_src) / np.abs(m_src))
    assert rel < 1e-12, f"PPM not mass-conservative on refined grid: rel {rel:.2e}"
    # finite + sane: kord=4 iv=1 PPM allows bounded edge overshoot (not strictly monotone), but the remapped
    # field must stay finite and within a generous band of the source range (no blow-up).
    assert np.isfinite(tgt).all(), "PPM produced non-finite values"
    span = src.max() - src.min()
    assert tgt.min() > src.min() - span and tgt.max() < src.max() + span, "PPM produced a runaway overshoot"
    print(f"  PPM mass conservation (refined grid): rel {rel:.2e}; finite + bounded  PASS")


def test_differentiable():
    rng = np.random.default_rng(1)
    ncol = 1
    p_src = np.cumsum(np.concatenate([[1.0e5], -rng.uniform(2000, 4000, 12)]))[None, :].repeat(ncol, 0)
    mids = 0.5 * (p_src[:, :-1] + p_src[:, 1:])
    p_tgt = np.zeros((ncol, p_src.shape[1] + p_src.shape[1] - 1))
    p_tgt[:, ::2] = p_src; p_tgt[:, 1::2] = mids
    src0 = jnp.asarray(rng.uniform(0.5, 3.0, (ncol, p_src.shape[1] - 1)))
    g = np.asarray(jax.grad(lambda s: jnp.sum(remap_vals_ppm(p_src, p_tgt, s, iv=1) ** 2))(src0))
    assert np.isfinite(g).all(), "non-finite grad through PPM remap"
    print(f"  jax.grad through PPM remap: finite ({g.size} entries)  PASS")


def main():
    print("test_remapping_ppm:")
    for t in (test_f2py_same_grid_identity, test_mass_conservation_refined, test_differentiable):
        t()
    print("All PPM remapping checks PASSED")


if __name__ == "__main__":
    main()
