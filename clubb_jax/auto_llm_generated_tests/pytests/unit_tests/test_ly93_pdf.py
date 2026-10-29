#!/usr/bin/env python3
"""validate the JAX calc_params_LY93 port (LY93_pdf.F90, Lewellen & Yoh 1993)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.LY93_pdf import calc_params_LY93, calc_mixt_frac_LY93

NG, NZ = 2, 6


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
