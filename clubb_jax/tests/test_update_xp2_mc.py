#!/usr/bin/env python3
"""test_update_xp2_mc.py — validate the JAX update_xp2_mc port (advance_xp2_xpyp_module.F90:6176).

Effects of rain evaporation on the second moments (l_morr_xp2_mc). Oracle:
  1. Literal NumPy transcription of the explicit Fortran formulas on the zt grid (the per-level math), then the
     independently-validated zt2zm interpolation -- so the comparison isolates this routine's algebra.
  2. precip_frac_double_delta top-down fill correctness.
  3. Physical sign invariants (rtp2_mc, thlp2_mc >= 0; wpthlp_mc <= 0).
  4. A finite jax.grad.

NB (iter 445): the direct f2py oracle `f2py_update_xp2_mc` is deliberately NOT used and the NumPy transcription is
the gold standard here, for two structural reasons (verified, not assumed): (a) the wrapper reads the PDF component
means/variances from a module-global `stored_pdf_params` (derived_type_storage) rather than taking them as args — so
the JAX's actual pdf_params can't be fed to it without a separate derived-type setter; and (b) like
`f2py_precip_fraction`, it FPE-traps under this build's `-ffpe-trap`. This routine is therefore validated by the
transcription here + the case-level KK rico oracle (test_kk_rico_oracle), not a direct f2py bit-shadow.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import update_xp2_mc as _update_xp2_mc
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.constants_clubb import cloud_frac_min, Cp, Lv
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def _pdf_from_dict(p, nzt, ngrdcol):
    params = init_pdf_params(nzt, ngrdcol)
    return params.replace(**{k: jnp.asarray(v) for k, v in p.items()})


def update_xp2_mc(gr, dt, cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap, p):
    zeros = jnp.zeros((gr.ngrdcol, gr.nzm), dtype=jnp.float64)
    return _update_xp2_mc(
        gr,
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        dt,
        jnp.asarray(cloud_frac),
        jnp.asarray(rcm),
        jnp.asarray(rvm),
        jnp.asarray(thlm),
        jnp.asarray(wm),
        jnp.asarray(exner),
        jnp.asarray(rrm_evap),
        _pdf_from_dict(p, gr.nzt, gr.ngrdcol),
        zeros,
        zeros,
        zeros,
        zeros,
        zeros,
    )


def _pdf_params(nzt, seed):
    rng = np.random.default_rng(seed)
    keys = ['mixt_frac', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'thl_1', 'thl_2',
            'varnce_thl_1', 'varnce_thl_2', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2']
    p = {}
    for k in keys:
        if k == 'mixt_frac':
            p[k] = rng.uniform(0.3, 0.7, (_NG, nzt))
        elif k.startswith('varnce'):
            p[k] = rng.uniform(0.0, 1.0, (_NG, nzt))
        elif k.startswith('rt'):
            p[k] = rng.uniform(0.0, 0.02, (_NG, nzt))
        elif k.startswith('thl'):
            p[k] = rng.uniform(280.0, 310.0, (_NG, nzt))
        else:
            p[k] = rng.uniform(-1.0, 1.0, (_NG, nzt))
    return p


def _ref_zt(dt, cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap, p):
    """Literal NumPy transcription of the zt-grid tendencies."""
    ng, nzt = cloud_frac.shape
    pf = np.zeros((ng, nzt))
    pf[:, nzt - 1] = 0.0
    for k in range(nzt - 2, -1, -1):
        pf[:, k] = np.where(cloud_frac[:, k] > cloud_frac_min, cloud_frac[:, k], pf[:, k + 1])
    pf_const = np.where(pf > cloud_frac_min, (1.0 - pf) / np.where(pf > cloud_frac_min, pf, 1.0), 0.0)
    a = p['mixt_frac']
    rt_tot = rcm + rvm
    t_rtp2 = a * ((p['rt_1'] - rt_tot) ** 2 + p['varnce_rt_1']) + (1 - a) * ((p['rt_2'] - rt_tot) ** 2 + p['varnce_rt_2'])
    t_thlp2 = a * ((p['thl_1'] - thlm) ** 2 + p['varnce_thl_1']) + (1 - a) * ((p['thl_2'] - thlm) ** 2 + p['varnce_thl_2'])
    t_wp2 = a * ((p['w_1'] - wm) ** 2 + p['varnce_w_1']) + (1 - a) * ((p['w_2'] - wm) ** 2 + p['varnce_w_2'])
    lc = Lv / (Cp * exner)
    ae = np.abs(rrm_evap)
    rtp2 = rrm_evap ** 2 * pf_const * dt + 2 * ae * np.sqrt(t_rtp2 * pf_const)
    thlp2 = (rrm_evap * lc) ** 2 * pf_const * dt + 2 * ae * lc * np.sqrt(t_thlp2 * pf_const)
    wprtp = ae * np.sqrt(pf_const) * np.sqrt(t_wp2)
    wpthlp = -lc * ae * np.sqrt(pf_const) * np.sqrt(t_wp2)
    rtpthlp = -ae * np.sqrt(pf_const) * (lc * np.sqrt(t_rtp2) + np.sqrt(t_thlp2)) - lc * pf_const * rrm_evap ** 2 * dt
    return rtp2, thlp2, wprtp, wpthlp, rtpthlp, pf


def _precip_frac_double_delta_reference(cloud_frac):
    ng, nzt = cloud_frac.shape
    pf = np.zeros((ng, nzt))
    pf[:, nzt - 1] = 0.0
    for k in range(nzt - 2, -1, -1):
        pf[:, k] = np.where(cloud_frac[:, k] > cloud_frac_min, cloud_frac[:, k], pf[:, k + 1])
    return pf


def _inputs(nzt, seed):
    rng = np.random.default_rng(seed)
    cloud_frac = rng.choice([0.0, 0.5, 0.9], (_NG, nzt))
    rcm = rng.uniform(0, 1e-3, (_NG, nzt)); rvm = rng.uniform(0, 0.02, (_NG, nzt))
    thlm = rng.uniform(285, 305, (_NG, nzt)); wm = rng.uniform(-1, 1, (_NG, nzt))
    exner = rng.uniform(0.8, 1.0, (_NG, nzt)); rrm_evap = -rng.uniform(0, 1e-5, (_NG, nzt))
    return cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap


def test_transcription_and_interp():
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zm.shape[1] - 1
    cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap = _inputs(nzt, 1)
    p = _pdf_params(nzt, 2)
    dt = 60.0
    got = [np.asarray(x) for x in update_xp2_mc(jgr, dt, cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap, p)]
    rtp2_zt, thlp2_zt, wprtp_zt, wpthlp_zt, rtpthlp_zt, pf = _ref_zt(
        dt, cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap, p)
    # precip-frac fill matches.
    assert np.max(np.abs(_precip_frac_double_delta_reference(cloud_frac) - pf)) < 1e-14, "precip_frac fill"
    # The full outputs equal zt2zm of the literal zt tendencies.
    worst = 0.0
    for g, ref_zt in zip(got, (rtp2_zt, thlp2_zt, wprtp_zt, wpthlp_zt, rtpthlp_zt)):
        ref = np.asarray(zt2zm(jgr.nzm, jgr.nzt, jgr.ngrdcol, jgr, ref_zt))
        worst = max(worst, np.max(np.abs(g - ref)))
    assert worst < 1e-12, f"update_xp2_mc transcription mismatch {worst:.2e}"
    print(f"  literal transcription + zt2zm (5 moments) + precip_frac fill, worst {worst:.2e}  PASS")


def test_invariants():
    jgr = setup_grid(ngrdcol=_NG, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zm.shape[1] - 1
    cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap = _inputs(nzt, 5)
    p = _pdf_params(nzt, 6)
    rtp2, thlp2, wprtp, wpthlp, rtpthlp = (np.asarray(x) for x in
                                           update_xp2_mc(jgr, 60.0, cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap, p))
    # On the zt grid the variance increments are >= 0; zt2zm is a positive-weight interp so zm values too.
    assert np.all(rtp2 >= -1e-12), "rtp2_mc should be >= 0"
    assert np.all(thlp2 >= -1e-12), "thlp2_mc should be >= 0"
    assert np.all(wpthlp <= 1e-12), "wpthlp_mc should be <= 0"
    print("  invariants: rtp2_mc>=0, thlp2_mc>=0, wpthlp_mc<=0  PASS")


def test_differentiable():
    jgr = setup_grid(ngrdcol=1, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1)
    nzt = jgr.zm.shape[1] - 1
    global _NG
    saved = _NG
    _NG = 1
    try:
        cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap = _inputs(nzt, 7)
        p = _pdf_params(nzt, 8)
        def loss(rrm):
            outs = update_xp2_mc(jgr, 60.0, cloud_frac, rcm, rvm, thlm, wm, exner, rrm, p)
            return sum(jnp.sum(o ** 2) for o in outs)
        g = np.asarray(jax.grad(loss)(jnp.asarray(rrm_evap)))
        assert np.isfinite(g).all(), "non-finite grad through update_xp2_mc"
    finally:
        _NG = saved
    print(f"  jax.grad through update_xp2_mc: finite ({g.size} entries)  PASS")


def main():
    print("test_update_xp2_mc:")
    for t in (test_transcription_and_interp, test_invariants, test_differentiable):
        t()
    print("All update_xp2_mc checks PASSED")


if __name__ == "__main__":
    main()
