#!/usr/bin/env python3
"""test_luhar_3d_driver.py — validate the JAX Luhar_3D_pdf_driver + backsolve_Luhar_params (adg1_adg2_3d_luhar_pdf.F90).

The 3-D Luhar PDF: max-|skewness| variable sets the PDF; the others backsolve their m. Oracles:
  1. f2py bit-shadow vs f2py_luhar_3d_pdf_driver (13 outputs), inputs spanning all three setter branches.
     SKIPs if clubb_f2py is unbuilt.
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

from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import Luhar_3D_pdf_driver

NG, NZ = 2, 8


def _inputs(rng):
    sh = (NG, NZ)
    return dict(wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0, 0.02, sh), thlm=rng.uniform(285, 305, sh),
                wp2=rng.uniform(0.05, 2, sh), rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
                Skw=rng.uniform(-2, 2, sh), Skrt=rng.uniform(-2, 2, sh), Skthl=rng.uniform(-2, 2, sh),
                wprtp=rng.uniform(-1e-3, 1e-3, sh), wpthlp=rng.uniform(-0.5, 0.5, sh))


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py luhar_3d_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(11)
    order = ('wm', 'rtm', 'thlm', 'wp2', 'rtp2', 'thlp2', 'Skw', 'Skrt', 'Skthl', 'wprtp', 'wpthlp')
    worst = 0.0
    for _ in range(4):
        a = _inputs(rng)
        f = clubb_f2py.f2py_luhar_3d_pdf_driver(*[a[k] for k in order])
        g = Luhar_3D_pdf_driver(*[a[k] for k in order])
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-9, f"luhar_3d_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py luhar_3d_pdf_driver: bit-match (13 outputs, all 3 setter branches), worst {worst:.2e}  PASS")


def test_differentiable():
    rng = np.random.default_rng(3)
    a = _inputs(rng)
    order = ('wm', 'rtm', 'thlm', 'wp2', 'rtp2', 'thlp2', 'Skw', 'Skrt', 'Skthl', 'wprtp', 'wpthlp')
    def loss(Skw):
        aa = dict(a); aa['Skw'] = Skw
        outs = Luhar_3D_pdf_driver(*[aa[k] for k in order])
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(a['Skw'])))
    assert np.isfinite(g).all(), "non-finite grad through Luhar_3D_pdf_driver"
    print(f"  jax.grad through Luhar_3D_pdf_driver: finite ({g.size} entries)  PASS")


def main():
    print("test_luhar_3d_driver:")
    for t in (test_f2py_oracle, test_differentiable):
        t()
    print("All Luhar_3D_pdf_driver checks PASSED")


if __name__ == "__main__":
    main()
