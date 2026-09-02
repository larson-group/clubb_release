"""Validation of the Morrison special functions (module_mp_graupel.py).

The KK-playbook foundational layer: validate the self-contained transcendental/utility
functions vs known physical values + an independent reference BEFORE building the process
rates on top. POLYSVP is validated here against (a) the exact triple-point intercept and
(b) the Goff-Gratch saturation formula (the physical reference the Flatau polynomial fits).
"""
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)  # CLUBB-JAX runs in float64
import jax.numpy as jnp
import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Microphys.Morrison_microphys.module_mp_graupel import (
    polysvp, derf1, gamma)

try:
    from scipy.special import erf as _scipy_erf, gamma as _scipy_gamma
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False


def _goff_gratch_liquid(T):
    """Goff-Gratch 1946 SVP over liquid water [Pa] — the physical reference Flatau fits."""
    Ts = 373.16
    e = (-7.90298 * (Ts / T - 1.0) + 5.02808 * np.log10(Ts / T)
         - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - T / Ts)) - 1.0)
         + 8.1328e-3 * (10.0 ** (-3.49149 * (Ts / T - 1.0)) - 1.0)
         + np.log10(1013.246))
    return 10.0 ** e * 100.0  # hPa -> Pa


def _goff_gratch_ice(T):
    """Goff-Gratch 1946 SVP over ice [Pa]."""
    Ti = 273.16
    e = (-9.09718 * (Ti / T - 1.0) - 3.56654 * np.log10(Ti / T)
         + 0.876793 * (1.0 - T / Ti) + np.log10(6.1071))
    return 10.0 ** e * 100.0


def test_polysvp_triple_point():
    """At T=273.16 (dt=0) POLYSVP returns exactly a0·100 — the Flatau intercept."""
    assert abs(float(polysvp(jnp.array(273.16), 0)) - 611.239921) < 1e-4, "liquid intercept"
    assert abs(float(polysvp(jnp.array(273.16), 1)) - 611.147274) < 1e-4, "ice intercept"
    print("  POLYSVP triple-point intercepts: liquid 611.24 / ice 611.15 Pa  PASS")


def _polysvp_fortran_replica(T, itype):
    """Independent pure-Python transcription of the Fortran POLYSVP Horner nesting — bit-matching
    the JAX confirms the JAX port has no structural/coefficient typo (the exact-transcription oracle,
    since the f2py API exposes no microphysics)."""
    a = ((6.11147274, 0.503160820, 0.188439774e-1, 0.420895665e-3, 0.615021634e-5,
          0.602588177e-7, 0.385852041e-9, 0.146898966e-11, 0.252751365e-14) if itype == 1
         else (6.11239921, 0.443987641, 0.142986287e-1, 0.264847430e-3, 0.302950461e-5,
               0.206739458e-7, 0.640689451e-10, -0.952447341e-13, -0.976195544e-15))
    dt = max(-80.0, T - 273.16)
    p = a[0] + dt * (a[1] + dt * (a[2] + dt * (a[3] + dt * (a[4] + dt * (a[5] + dt * (
        a[6] + dt * (a[7] + a[8] * dt)))))))
    return p * 100.0


def test_polysvp_transcription_exact():
    """JAX POLYSVP == the pure-Python Fortran-Horner replica to machine precision, over the full
    valid range and both phases — the definitive transcription check."""
    worst = 0.0
    for itype in (0, 1):
        for T in np.linspace(180.0, 320.0, 80):
            j = float(polysvp(jnp.array(T), itype))
            f = _polysvp_fortran_replica(float(T), itype)
            worst = max(worst, abs(j - f) / (abs(f) + 1e-30))
    assert worst < 1e-13, f"POLYSVP transcription rel {worst:.2e}"
    print(f"  POLYSVP == Fortran-Horner replica (rel {worst:.1e})  PASS")


def test_polysvp_physical_plausibility():
    """POLYSVP (Flatau) tracks the Goff-Gratch reference within the formulas' mutual agreement
    (~1% at the cold extremes, <0.3% near 0 C) — a gross coefficient error would show as ≫1%."""
    Tl = np.linspace(248.16, 313.16, 40)  # -25 to +40 C: Flatau & GG agree tightly
    rel = np.abs(np.array([float(polysvp(jnp.array(t), 0)) for t in Tl])
                 - _goff_gratch_liquid(Tl)) / _goff_gratch_liquid(Tl)
    assert rel.max() < 4e-3, f"liquid plausibility rel {rel.max():.2e}"
    Ti = np.linspace(243.16, 273.16, 30)  # -30 to 0 C
    reli = np.abs(np.array([float(polysvp(jnp.array(t), 1)) for t in Ti])
                  - _goff_gratch_ice(Ti)) / _goff_gratch_ice(Ti)
    assert reli.max() < 5e-3, f"ice plausibility rel {reli.max():.2e}"
    print(f"  POLYSVP physical plausibility vs Goff-Gratch: liquid<{rel.max():.1e} ice<{reli.max():.1e}  PASS")


def test_polysvp_floor_and_monotonic():
    """The dt=max(-80,...) floor caps T below 193.16 K; SVP is monotonic increasing in T."""
    # Floor: T=150 and T=193.16 (both dt clipped to -80) give the same value
    v_cold = float(polysvp(jnp.array(150.0), 0))
    v_floor = float(polysvp(jnp.array(193.16), 0))
    assert abs(v_cold - v_floor) < 1e-10, "dt floor not applied at extreme cold"
    # Monotonic over the valid range
    T = np.linspace(193.16, 313.16, 60)
    es = np.array([float(polysvp(jnp.array(t), 0)) for t in T])
    assert np.all(np.diff(es) > 0), "POLYSVP not monotonic increasing"
    print("  POLYSVP dt-floor + monotonicity  PASS")


def test_derf1_vs_scipy():
    """DERF1 (Ooura approximation) == scipy.special.erf to ~double precision over the full range,
    including across both branch boundaries (w=2.2, 6.9) and the saturated tail. The Ooura table is
    a near-double-precision fit, so the match should be tight (≪ the gate), confirming the 130
    coefficients + branch logic were transcribed correctly."""
    if not _HAVE_SCIPY:
        print("  DERF1 vs scipy: SKIPPED (no scipy)"); return
    x = np.concatenate([np.linspace(-7.0, 7.0, 600),
                        np.array([0.0, 2.2, 2.2 - 1e-9, 6.9, 6.9 - 1e-9, -2.2, -6.9, 1.0, -1.0])])
    j = np.array(derf1(jnp.asarray(x)))
    s = _scipy_erf(x)
    err = np.abs(j - s)
    assert err.max() < 1e-12, f"DERF1 vs scipy max abs err {err.max():.2e}"
    print(f"  DERF1 == scipy.special.erf (max abs err {err.max():.1e}, incl. branch boundaries)  PASS")


def test_derf1_identities():
    """erf(0)=0, erf is odd, erf(±large)=±1, monotonic increasing."""
    assert abs(float(derf1(jnp.array(0.0)))) < 1e-15, "erf(0)≠0"
    xs = jnp.array([0.3, 1.7, 3.5, 5.0])
    odd = np.abs(np.array(derf1(xs)) + np.array(derf1(-xs)))
    assert odd.max() < 1e-14, "erf not odd"
    assert abs(float(derf1(jnp.array(8.0))) - 1.0) < 1e-15, "erf(8)≠1"
    assert abs(float(derf1(jnp.array(-8.0))) + 1.0) < 1e-15, "erf(-8)≠-1"
    g = np.array(derf1(jnp.linspace(-6.0, 6.0, 200)))
    assert np.all(np.diff(g) >= -1e-15), "erf not monotonic"
    print("  DERF1 identities (zero/odd/saturation/monotonic)  PASS")


def test_gamma_vs_scipy_positive():
    """GAMMA (Cody) == scipy.special.gamma over the positive range that Morrison actually uses
    (size-dist shape params, ~1..15), spanning all three positive branches: y<1 (reduction up),
    1≤y<12 (rational + integer reduction down), y≥12 (Stirling). Cody is a ~20-digit fit."""
    if not _HAVE_SCIPY:
        print("  GAMMA vs scipy: SKIPPED (no scipy)"); return
    x = np.concatenate([np.linspace(0.05, 15.0, 500),
                        np.array([1.0, 2.0, 11.999, 12.0, 12.001, 0.999, 1.001])])
    j = np.array(gamma(jnp.asarray(x)))
    s = _scipy_gamma(x)
    rel = np.abs(j - s) / np.abs(s)
    assert rel.max() < 1e-12, f"GAMMA vs scipy max rel {rel.max():.2e}"
    print(f"  GAMMA == scipy.special.gamma over [0.05,15] (max rel {rel.max():.1e}, all 3 branches)  PASS")


def test_gamma_known_and_reflection():
    """Exact known values + the negative-argument reflection branch (−π/sin(πr))."""
    sqrtpi = np.sqrt(np.pi)
    assert abs(float(gamma(jnp.array(1.0))) - 1.0) < 1e-14, "Γ(1)=1"
    assert abs(float(gamma(jnp.array(5.0))) - 24.0) < 1e-12, "Γ(5)=24"
    assert abs(float(gamma(jnp.array(0.5))) - sqrtpi) < 1e-13, "Γ(0.5)=√π"
    # Reflection: Γ(-0.5) = -2√π, Γ(-1.5) = 4√π/3
    if _HAVE_SCIPY:
        for v in (-0.5, -1.5, -2.5, -0.3, -3.7):
            assert abs(float(gamma(jnp.array(v))) - _scipy_gamma(v)) / abs(_scipy_gamma(v)) < 1e-11, \
                f"reflection Γ({v})"
    print("  GAMMA known values + negative-arg reflection  PASS")


def test_gamma_poles_and_overflow():
    """Poles at non-positive integers and overflow beyond XBIG return XINF=3.4e38 (Fortran behavior)."""
    XINF = 3.4e38
    for v in (0.0, -1.0, -2.0, -5.0):
        assert float(gamma(jnp.array(v))) == XINF, f"pole Γ({v}) not XINF"
    assert float(gamma(jnp.array(40.0))) == XINF, "overflow not XINF (y>XBIG)"
    print("  GAMMA poles + overflow → XINF  PASS")


if __name__ == "__main__":
    print("Morrison special-function validation:")
    test_polysvp_triple_point()
    test_polysvp_transcription_exact()
    test_polysvp_physical_plausibility()
    test_polysvp_floor_and_monotonic()
    test_derf1_vs_scipy()
    test_derf1_identities()
    test_gamma_vs_scipy_positive()
    test_gamma_known_and_reflection()
    test_gamma_poles_and_overflow()
    print("All Morrison special-function tests PASSED.")
