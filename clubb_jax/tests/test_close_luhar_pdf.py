#!/usr/bin/env python3
"""test_close_luhar_pdf.py — validate the JAX close_Luhar_pdf port (adg1_adg2_3d_luhar_pdf.F90).

PDF component widths/means/variances for the Luhar closure. Oracles:
  1. f2py bit-shadow vs f2py_close_luhar_pdf (8 outputs), in the well-defined xp2 > x_tol_sqd regime (the
     Fortran's degenerate branch reads an uninitialized sgn, so it is not compared). SKIPs if unbuilt.
  2. Moment reconstruction: the binormal reproduces the overall mean xm and variance xp2.
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

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import close_Luhar_pdf, ADG1_w_closure

NG, NZ = 2, 6
_X_TOL_SQD = 1.0e-8


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))      # all > x_tol_sqd
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    small_m = rng.uniform(0.05, 1.0, (NG, NZ))
    wpxp = rng.uniform(-1, 1, (NG, NZ))
    return xm, xp2, mixt_frac, small_m, wpxp


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py close_luhar_pdf oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, mf, m, wpxp = _inputs(seed)
        f = clubb_f2py.f2py_close_luhar_pdf(xm, xp2, mf, m, wpxp, _X_TOL_SQD)
        g = close_Luhar_pdf(xm, xp2, mf, m, wpxp, _X_TOL_SQD)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"close_luhar_pdf f2py mismatch {worst:.2e}"
    print(f"  f2py close_luhar_pdf: bit-match (8 outputs, xp2>tol), worst {worst:.2e}  PASS")


def test_moment_reconstruction():
    xm, xp2, mf, m, wpxp = _inputs(5)
    ss1, ss2, v1, v2, x1n, x2n, x1, x2 = (np.asarray(x) for x in close_Luhar_pdf(xm, xp2, mf, m, wpxp, _X_TOL_SQD))
    xm_rec = mf * x1 + (1 - mf) * x2
    assert np.max(np.abs(xm_rec - xm)) < 1e-12, "overall mean not reproduced"
    xp2_rec = mf * ((x1 - xm) ** 2 + v1) + (1 - mf) * ((x2 - xm) ** 2 + v2)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-10, "overall variance not reproduced"
    print("  moment reconstruction: binormal reproduces overall mean & variance  PASS")


def test_differentiable():
    xm, xp2, mf, m, wpxp = _inputs(7)
    def loss(v):
        outs = close_Luhar_pdf(xm, v, mf, m, wpxp, _X_TOL_SQD)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(xp2)))
    assert np.isfinite(g).all(), "non-finite grad through close_Luhar_pdf"
    print(f"  jax.grad through close_Luhar_pdf: finite ({g.size} entries)  PASS")


def test_ADG1_w_closure_f2py():
    """ADG1_w_closure (adg1_adg2_3d_luhar_pdf.F90) — the ADG1 two-component w-PDF closure (w_1/w_2, normal-space
    w_1_n/w_2_n, varnce_w_1/2, mixt_frac from wm/wp2/Skw/sigma_sqd_w). A core every-step routine that had only a
    grad test; bit-shadow its 7 outputs vs the f2py oracle. SKIPs if clubb_f2py is unbuilt. (iter 421)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py ADG1_w_closure oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(40):
        ng, nz = 2, 8
        wp2 = rng.uniform(1e-3, 2.0, (ng, nz)); wm = rng.uniform(-1.0, 1.0, (ng, nz))
        Skw = rng.uniform(-3.0, 3.0, (ng, nz)); ssw = rng.uniform(0.0, 0.6, (ng, nz)); mfmm = 0.4
        jr = ADG1_w_closure(jnp.asarray(wm), jnp.asarray(wp2), jnp.asarray(Skw), jnp.asarray(ssw),
                            jnp.asarray(np.sqrt(wp2)), mfmm)
        fr = clubb_f2py.f2py_adg1_w_closure(wm, wp2, Skw, ssw, np.sqrt(wp2), mfmm)
        for i in range(7):
            worst = max(worst, float(np.max(np.abs(np.asarray(jr[i]) - np.asarray(fr[i])))))
    assert worst < 1e-12, f"ADG1_w_closure f2py mismatch {worst:.2e}"
    print(f"  f2py ADG1_w_closure (7 PDF outputs, 40 cases): bit-match, worst {worst:.2e}  PASS")


def test_ADG1_pdf_driver_f2py():
    """ADG1_pdf_driver (adg1_adg2_3d_luhar_pdf.F90) — THE core every-step ADG1 PDF-parameter driver: from the mean
    fields + 2nd moments it produces the two-component (w, rt, thl, u, v) means/variances + mixt_frac + alpha
    factors (25 non-scalar outputs). Case-validated end-to-end but had no unit test. Bit-shadow all 25 outputs vs
    the f2py oracle (sclr_dim=1 + l_scalar_calc=False to bypass the wrapper's size-0 scalar-array error; the JAX is
    sclr_dim=0). The JAX returns a dict; f2py outputs 0..24 are these names in order. SKIPs if unbuilt. (iter 423)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py ADG1_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import ADG1_pdf_driver
    fn = ['w_1', 'w_2', 'rt_1', 'rt_2', 'thl_1', 'thl_2', 'u_1', 'u_2', 'v_1', 'v_2',
          'varnce_w_1', 'varnce_w_2', 'varnce_rt_1', 'varnce_rt_2', 'varnce_thl_1', 'varnce_thl_2',
          'varnce_u_1', 'varnce_u_2', 'varnce_v_1', 'varnce_v_2', 'mixt_frac', 'alpha_rt', 'alpha_thl',
          'alpha_u', 'alpha_v']
    rng = np.random.default_rng(1); ng, nz = 2, 8; worst = 0.0
    for _ in range(20):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(1e-3, 1.5e-2, (ng, nz)); thlm = rng.uniform(290, 320, (ng, nz))
        um = rng.uniform(-5, 5, (ng, nz)); vm = rng.uniform(-5, 5, (ng, nz)); wp2 = rng.uniform(1e-2, 2, (ng, nz))
        rtp2 = rng.uniform(1e-7, 1e-5, (ng, nz)); thlp2 = rng.uniform(1e-3, 1, (ng, nz))
        up2 = rng.uniform(1e-2, 1, (ng, nz)); vp2 = rng.uniform(1e-2, 1, (ng, nz))
        Skw = rng.uniform(-3, 3, (ng, nz)); wprtp = rng.uniform(-1e-3, 1e-3, (ng, nz)); wpthlp = rng.uniform(-1, 1, (ng, nz))
        upwp = rng.uniform(-0.1, 0.1, (ng, nz)); vpwp = rng.uniform(-0.1, 0.1, (ng, nz))
        sqw = np.sqrt(wp2); ssw = rng.uniform(0, 0.5, (ng, nz)); beta = np.full(ng, 1.0); mfmm = 0.4
        jr = ADG1_pdf_driver(jnp.asarray(wm), jnp.asarray(rtm), jnp.asarray(thlm), jnp.asarray(um), jnp.asarray(vm),
            jnp.asarray(wp2), jnp.asarray(rtp2), jnp.asarray(thlp2), jnp.asarray(up2), jnp.asarray(vp2),
            jnp.asarray(Skw), jnp.asarray(wprtp), jnp.asarray(wpthlp), jnp.asarray(upwp), jnp.asarray(vpwp),
            jnp.asarray(sqw), jnp.asarray(ssw), jnp.asarray(beta)[:, None], mfmm)
        sm = np.zeros((ng, nz, 1))
        fr = clubb_f2py.f2py_adg1_pdf_driver(1, 1e-8, wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, Skw,
            wprtp, wpthlp, upwp, vpwp, sqw, ssw, beta, mfmm, sm, sm, sm, False)
        for i, nm in enumerate(fn):
            worst = max(worst, float(np.max(np.abs(np.asarray(jr[nm]) - np.asarray(fr[i])))))
    assert worst < 1e-12, f"ADG1_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py ADG1_pdf_driver (25 PDF-param outputs, 20 cases): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_close_luhar_pdf:")
    for t in (test_f2py_oracle, test_moment_reconstruction, test_ADG1_w_closure_f2py,
              test_ADG1_pdf_driver_f2py, test_differentiable):
        t()
    print("All close_Luhar_pdf checks PASSED")


if __name__ == "__main__":
    main()
