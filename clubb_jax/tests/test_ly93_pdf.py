#!/usr/bin/env python3
"""test_ly93_pdf.py — validate the JAX calc_params_LY93 port (LY93_pdf.F90, Lewellen & Yoh 1993).

Gives the two binormal PDF-component means/variances of x from (xm, xp2, Skx, mixt_frac). Oracles:
  1. f2py bit-shadow vs f2py_calc_params_ly93. SKIPs if clubb_f2py is unbuilt.
  2. Moment reconstruction (oracle-free): the resulting binormal reproduces the overall mean xm and variance
     xp2 exactly (compute_mean_binormal / compute_variance_binormal).
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

from clubb_jax.src.CLUBB_core.LY93_pdf import calc_params_LY93, calc_mixt_frac_LY93, LY93_driver
from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, compute_variance_binormal

NG, NZ = 2, 6


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-2.5, 2.5, (NG, NZ))
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    return xm, xp2, Skx, mixt_frac


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calc_params_ly93 oracle: SKIP ({type(e).__name__})")
        return
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, Skx, mf = _inputs(seed)
        f = clubb_f2py.f2py_calc_params_ly93(xm, xp2, Skx, mf)
        g = calc_params_LY93(xm, xp2, Skx, mf)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"calc_params_ly93 f2py mismatch {worst:.2e}"
    print(f"  f2py calc_params_ly93: bit-match (4 outputs) over 3 seeds, worst {worst:.2e}  PASS")


def test_moment_reconstruction():
    # Keep |Skx| modest so the component variances stay non-negative (LY93 can yield negative component
    # variances for large skewness — a known feature clipped downstream; the moment identities still hold
    # algebraically, which we verify directly with the signed component variances).
    rng = np.random.default_rng(5)
    xm = rng.uniform(-2, 2, (NG, NZ)); xp2 = rng.uniform(0.05, 2.0, (NG, NZ))
    Skx = rng.uniform(-1.0, 1.0, (NG, NZ)); mf = rng.uniform(0.3, 0.7, (NG, NZ))
    mu1, mu2, s1sq, s2sq = (np.asarray(x) for x in calc_params_LY93(xm, xp2, Skx, mf))
    a = mf
    # Overall mean = a*mu1 + (1-a)*mu2 == xm.
    xm_rec = a * mu1 + (1 - a) * mu2
    assert np.max(np.abs(xm_rec - xm)) < 1e-12, "overall mean not reproduced"
    # Overall variance = a((mu1-xm)^2+s1sq) + (1-a)((mu2-xm)^2+s2sq) == xp2 (signed component variances).
    xp2_rec = a * ((mu1 - xm) ** 2 + s1sq) + (1 - a) * ((mu2 - xm) ** 2 + s2sq)
    assert np.max(np.abs(xp2_rec - xp2)) < 1e-12, "overall variance not reproduced"
    # Overall (unnormalized) third moment = a((mu1-xm)^3 + 3(mu1-xm)s1sq) + ... == Skx*xp2^(3/2).
    m3 = (a * ((mu1 - xm) ** 3 + 3 * (mu1 - xm) * s1sq)
          + (1 - a) * ((mu2 - xm) ** 3 + 3 * (mu2 - xm) * s2sq))
    assert np.max(np.abs(m3 - Skx * xp2 ** 1.5)) < 1e-11, "skewness not reproduced"
    print("  moment reconstruction: binormal reproduces overall mean / variance / skewness  PASS")


def test_mixt_frac_root():
    # For Sk_max > 0.84 the bisection root satisfies mf^6 = Sk_max^2 (1-mf) to ~tolerance; for <=0.84, mf=0.75.
    Sk = np.array([[0.5, 0.84, 1.0, 2.0, 5.0, 0.9]])
    mf = np.asarray(calc_mixt_frac_LY93(Sk))
    assert abs(mf[0, 0] - 0.75) < 1e-12 and abs(mf[0, 1] - 0.75) < 1e-12, "Sk_max<=0.84 -> 0.75"
    for j in (2, 3, 4, 5):
        expr = mf[0, j] ** 6 - Sk[0, j] ** 2 * (1.0 - mf[0, j])
        assert abs(expr) < 1e-4, f"bisection residual at j={j}: {expr:.2e}"
        assert 0.5 <= mf[0, j] <= 1.0, "mixt_frac out of [0.5,1]"
    print("  calc_mixt_frac_LY93: root mf^6=Sk_max^2(1-mf) within tol; 0.75 below threshold  PASS")


def test_driver_f2py():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py ly93_driver oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(8)
    worst = 0.0
    for _ in range(5):
        wm = rng.uniform(-1, 1, (NG, NZ)); rtm = rng.uniform(0, 0.02, (NG, NZ)); thlm = rng.uniform(285, 305, (NG, NZ))
        wp2 = rng.uniform(0.05, 2.0, (NG, NZ)); rtp2 = rng.uniform(1e-8, 1e-6, (NG, NZ)); thlp2 = rng.uniform(0.01, 1.0, (NG, NZ))
        # Skewness spanning both the 0.75 branch (|Sk|<=0.84) and the iterative branch.
        Skw = rng.uniform(-3, 3, (NG, NZ)); Skrt = rng.uniform(-1, 1, (NG, NZ)); Skthl = rng.uniform(-2, 2, (NG, NZ))
        f = clubb_f2py.f2py_ly93_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl)
        g = LY93_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-9, f"ly93_driver f2py mismatch {worst:.2e}"
    print(f"  f2py ly93_driver: match (13 outputs incl. iterative mixt_frac), worst {worst:.2e}  PASS")


def test_differentiable():
    xm, xp2, Skx, mf = _inputs(7)
    def loss(Skx_v):
        outs = calc_params_LY93(xm, xp2, Skx_v, mf)
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(Skx)))
    assert np.isfinite(g).all(), "non-finite grad through calc_params_LY93"
    print(f"  jax.grad through calc_params_LY93: finite ({g.size} entries)  PASS")


def main():
    print("test_ly93_pdf:")
    for t in (test_f2py_oracle, test_moment_reconstruction, test_mixt_frac_root,
              test_driver_f2py, test_differentiable):
        t()
    print("All LY93 checks PASSED")


if __name__ == "__main__":
    main()
