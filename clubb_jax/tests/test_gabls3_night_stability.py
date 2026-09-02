#!/usr/bin/env python3
"""test_gabls3_night_stability.py — pin the GABLS3-night Businger-Dyer surface-layer stability functions.

`gm1/gh1/fm1/fh1` are the Monin-Obukhov stability functions
that the `landflx` surface scheme composes; they are bit-validated end-to-end (gabls3_night IS in the bit-faithful
regression set) but never unit-pinned in ISOLATION, so a coefficient typo would only surface as a full-case
divergence. This pins each against an INDEPENDENT transcription of the Fortran formulas (so the 15 / 9 / 0.74 / π/2
coefficients are checked directly):
    gm1(x) = (1 − 15x)^0.25
    gh1(x) = sqrt(1 − 9x) / 0.74
    fm1(x) = 2·log((1+x)/2) + log((1+x²)/2) − 2·atan(x) + π/2
    fh1(x) = 2·log((1+0.74x)/2)
gm1/gh1 are evaluated over the unstable regime x<0 where `landflx` calls them;
fm1/fh1 preserve the source's explicit default-real `alog` evaluation before
conversion back to CLUBB core precision. Oracle-independent; never SKIPs.
"""
import os
import sys
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import ep1

configure_jax_precision()

from clubb_jax.src.Benchmark_cases.gabls3_night import (
    fh1,
    fm1,
    gh1,
    gm1,
    landflx,
    psi_h,
)


def _ref_gm1(x): return (1.0 - 15.0 * x) ** 0.25
def _ref_gh1(x): return math.sqrt(1.0 - 9.0 * x) / 0.74
def _alog(x): return float(np.log(np.float32(x)))
def _ref_fm1(x): return 2.0 * _alog((np.float32(1.0) + np.float32(x)) / np.float32(2.0)) \
    + _alog((np.float32(1.0) + np.float32(x * x)) / np.float32(2.0)) - 2.0 * math.atan(x) + math.pi / 2.0
def _ref_fh1(x): return 2.0 * _alog( \
    (np.float32(1.0) + np.float32(0.74) * np.float32(x)) / np.float32(2.0))


def test_gm1_gh1_unstable_regime():
    worst = 0.0
    for x in np.linspace(-2.0, -1.0e-3, 50):     # unstable: z/L < 0 (where landflx calls gm1/gh1)
        worst = max(worst, abs(float(gm1(float(x))) - _ref_gm1(x)),
                    abs(float(gh1(float(x))) - _ref_gh1(x)))
    assert worst < 1e-13, f"gm1/gh1 mismatch vs F90 formula {worst:.2e}"
    print(f"  gm1=(1−15x)^¼, gh1=√(1−9x)/0.74 over x∈[−2,0): match F90 formula (worst {worst:.1e})  PASS")


def test_fm1_fh1_formula():
    worst = 0.0
    for x in np.linspace(1.0, 6.0, 50):          # gm1/gh1 outputs are > 1 in the unstable regime
        worst = max(worst, abs(float(fm1(float(x))) - _ref_fm1(x)),
                    abs(float(fh1(float(x))) - _ref_fh1(x)))
    # JAX and the host libm may differ by a few float32 ulps for ALOG.
    assert worst < 5e-7, f"fm1/fh1 mismatch vs F90 formula {worst:.2e}"
    print(f"  fm1 (log+atan+π/2), fh1=2·log((1+0.74x)/2): match F90 formula (worst {worst:.1e})  PASS")


def test_landflx_stable_default_real_boundaries():
    """Pin the stable GABLS3-night path, including its default-real ALOG calls."""
    th, ts = 290.0, 289.0
    qh, qs = 0.005, 0.006
    uh, vh = 5.0, 1.0
    h, z0 = 3.125, 0.15

    zody = _alog(np.float32(h) / np.float32(z0))
    vel_floor = math.sqrt(max(0.5, uh**2 + vh**2))
    r = 9.81 / ts * (th * (1.0 + ep1 * qh) - ts * (1.0 + ep1 * qs)) * h / vel_floor**2
    a = 4.8 * 4.8 * r - 6.35
    b = (2.0 * r * 4.8 - 1.0) * zody
    c = r * zody**2
    d = math.sqrt(b * b - 4.0 * a * c)
    xsi1 = (-b + d) / a / 2.0
    xsi2 = (-b - d) / a / 2.0
    xsi = float(max(np.float32(xsi1), np.float32(xsi2)))
    fm = zody + 4.8 * xsi
    vel = math.sqrt(uh**2 + vh**2)
    ustar = 0.4 / fm * vel
    xsi = max(1.0e-5, xsi)
    xlmo = h / xsi
    denominator = _alog(np.float32(h) / np.float32(0.25)) + 5.0 * h / xlmo - 5.0 * 0.25 / xlmo
    expected = (
        0.4 * ustar * (ts - th) / denominator,
        0.4 * ustar * (qs - qh) / denominator,
        vel,
        ustar,
    )
    actual = tuple(float(value) for value in landflx(th, ts, qh, qs, uh, vh, h, z0))
    np.testing.assert_allclose(actual, expected, rtol=1.0e-14, atol=1.0e-15)


def test_gh1_called_domain():
    """The source square-root argument is positive in the called domain."""
    for x in (-0.5, -0.01, -1.5, 0.0):
        assert abs(float(gh1(x)) - math.sqrt(1.0 - 9.0 * x) / 0.74) < 1e-13
    print("  gh1 matches the Fortran sqrt(1-9x)/0.74 in the x<=0 call domain  PASS")


def test_psi_h_stable():
    """psi_h(x, xlmo) = −5·x/xlmo — the stable-case integrated heat stability function (gabls3_night.F90:150).
    Linear in x; pin the −5 coefficient and the 1/xlmo dependence (added iter 574, completing the gabls3_night
    surface stability functions gm1/gh1/fm1/fh1/psi_h)."""
    for x, xlmo in ((0.25, 50.0), (10.0, 200.0), (2.0, 30.0)):
        assert abs(float(psi_h(x, xlmo)) - (-5.0 * x) / xlmo) < 1e-14, f"psi_h({x},{xlmo})"
    # the landflx combination psi_h(0.25,xlmo) − psi_h(h,xlmo) telescopes linearly in (0.25−h)
    assert abs((psi_h(0.25, 100.0) - psi_h(5.0, 100.0)) - (-5.0 * (0.25 - 5.0) / 100.0)) < 1e-14
    print("  psi_h = −5·x/xlmo (stable-case heat stability function)  PASS")


def main():
    print("test_gabls3_night_stability:")
    test_gm1_gh1_unstable_regime()
    test_fm1_fh1_formula()
    test_landflx_stable_default_real_boundaries()
    test_gh1_called_domain()
    test_psi_h_stable()
    print("All gabls3_night stability-function checks PASSED")


if __name__ == "__main__":
    main()
