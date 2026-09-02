#!/usr/bin/env python3
"""test_adg2_driver.py — validate the JAX ADG2_pdf_driver port (adg1_adg2_3d_luhar_pdf.F90).

ADG2 PDF: w from the Luhar closure, rt/thl as ADG responders. Oracles:
  1. f2py bit-shadow vs f2py_adg2_pdf_driver (the 16 w/rt/thl PDF-param outputs; sclr_dim=0). SKIPs if unbuilt.
  2. A finite jax.grad.
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

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import ADG2_pdf_driver as _ADG2_pdf_driver

NG, NZ = 2, 8
_BETA = 1.75


def ADG2_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, sqrt_wp2, beta):
    ng, nz = np.asarray(wm).shape
    sclr_dim = 1
    sclr_tol = jnp.full((sclr_dim,), 1e-8)
    sclrm = jnp.zeros((ng, nz, sclr_dim), dtype=jnp.float64)
    sclrp2 = jnp.ones((ng, nz, sclr_dim), dtype=jnp.float64) * 1e-6
    wpsclrp = jnp.zeros((ng, nz, sclr_dim), dtype=jnp.float64)
    return _ADG2_pdf_driver(
        nz,
        ng,
        sclr_dim,
        sclr_tol,
        wm,
        rtm,
        thlm,
        wp2,
        rtp2,
        thlp2,
        Skw,
        wprtp,
        wpthlp,
        sqrt_wp2,
        beta,
        sclrm,
        sclrp2,
        wpsclrp,
        False,
    )


def _inputs(rng):
    sh = (NG, NZ)
    wp2 = rng.uniform(0.05, 2, sh)
    return dict(wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0, 0.02, sh), thlm=rng.uniform(285, 305, sh),
                wp2=wp2, rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
                Skw=rng.uniform(-2, 2, sh), wprtp=rng.uniform(-1e-3, 1e-3, sh), wpthlp=rng.uniform(-0.5, 0.5, sh),
                sqrt_wp2=np.sqrt(wp2))


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py adg2_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(3):
        a = _inputs(rng)
        beta_col = np.full(NG, _BETA)
        # f2py rejects 0-sized dims, so use sclr_dim=1 with dummy scalars and l_scalar_calc=False
        # (the scalar branch is then skipped and only the 16 w/rt/thl outputs are compared).
        sd = 1
        sclr_tol = np.full(sd, 1e-8)
        sm = np.zeros((NG, NZ, sd)); sp = np.ones((NG, NZ, sd)) * 1e-6; ws = np.zeros((NG, NZ, sd))
        f = clubb_f2py.f2py_adg2_pdf_driver(
            sd, sclr_tol, a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['Skw'], a['wprtp'], a['wpthlp'], a['sqrt_wp2'], beta_col, sm, sp, ws, False)
        g = ADG2_pdf_driver(a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'],
                            a['Skw'], a['wprtp'], a['wpthlp'], a['sqrt_wp2'], beta_col[:, None])
        for k in range(16):   # the 16 non-scalar PDF-param outputs
            worst = max(worst, np.max(np.abs(np.asarray(g[k]) - np.asarray(f[k]))))
    assert worst < 1e-10, f"adg2_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py adg2_pdf_driver: bit-match (16 w/rt/thl PDF-param outputs), worst {worst:.2e}  PASS")


def test_differentiable():
    rng = np.random.default_rng(3)
    a = _inputs(rng)
    beta = jnp.full((NG, 1), _BETA)
    def loss(Skw):
        outs = ADG2_pdf_driver(a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'],
                               Skw, a['wprtp'], a['wpthlp'], a['sqrt_wp2'], beta)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(a['Skw'])))
    assert np.isfinite(g).all(), "non-finite grad through ADG2_pdf_driver"
    print(f"  jax.grad through ADG2_pdf_driver: finite ({g.size} entries)  PASS")


def main():
    print("test_adg2_driver:")
    for t in (test_f2py_oracle, test_differentiable):
        t()
    print("All ADG2_pdf_driver checks PASSED")


if __name__ == "__main__":
    main()
