"""Validate the JAX saturation port (saturation.py) against saturation.F90 logic.

Focus: the `sat_vapor_press_liq` dispatcher (extracted iter 307) selects the correct
leaf SVP approximation by `saturation_formula`, and `sat_mixrat_liq` is consistent with
it. The leaf polynomials (Flatau/Bolton) are checked against well-known reference values
near 0 degC; the dispatcher is checked bit-exact against the leaves it routes to. The
analytic checks never SKIP; the f2py bit-shadow (test_sat_mixrat_f2py) SKIPs cleanly when
clubb_f2py is unbuilt.
"""
import os
import sys

import numpy as np
import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.CLUBB_core.saturation import (
    sat_vapor_press_liq, sat_vapor_press_liq_flatau, sat_vapor_press_liq_bolton,
    sat_mixrat_liq, sat_mixrat_ice, SATURATION_FLATAU, SATURATION_BOLTON,
)

_T_FREEZE = 273.15


def test_dispatcher_matches_leaves():
    """sat_vapor_press_liq routes bit-exactly to the leaf chosen by saturation_formula."""
    T = jnp.linspace(220.0, 310.0, 64)
    for formula, leaf in ((SATURATION_FLATAU, sat_vapor_press_liq_flatau),
                          (SATURATION_BOLTON, sat_vapor_press_liq_bolton)):
        d = float(jnp.max(jnp.abs(sat_vapor_press_liq(T, formula) - leaf(T))))
        assert d == 0.0, f"formula {formula}: dispatcher != leaf ({d})"
    # The documented not-target formulas (saturation.py: gfdl/lookup) must raise — pin the SPECIFIC
    # constants_clubb codes (gfdl=2, lookup=4), not just a generic out-of-range value, so a future change
    # that silently routes them somewhere is caught (mirrors saturation.F90's gfdl/lookup leaves being unported).
    from clubb_jax.src.CLUBB_core.model_flags import saturation_gfdl, saturation_lookup
    for bad in (saturation_gfdl, saturation_lookup, 999):
        try:
            sat_vapor_press_liq(T, bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"unsupported saturation_formula={bad} should raise ValueError")
    print("  dispatcher routes bit-exactly (flatau/bolton) + rejects gfdl/lookup/unknown formula  PASS")


def test_svp_reference_values():
    """SVP over liquid at 0 degC ~ 611 Pa; both formulas agree to ~1% near 0-20 degC."""
    T0 = jnp.array([_T_FREEZE])
    es_flatau = float(sat_vapor_press_liq_flatau(T0)[0])
    es_bolton = float(sat_vapor_press_liq_bolton(T0)[0])
    assert 605.0 < es_flatau < 615.0, f"Flatau SVP(0C)={es_flatau}"
    assert 605.0 < es_bolton < 615.0, f"Bolton SVP(0C)={es_bolton}"
    T = jnp.linspace(_T_FREEZE, _T_FREEZE + 20.0, 21)
    rel = float(jnp.max(jnp.abs(sat_vapor_press_liq_flatau(T) - sat_vapor_press_liq_bolton(T))
                        / sat_vapor_press_liq_flatau(T)))
    assert rel < 0.02, f"Flatau vs Bolton disagree by {rel:.3f} near 0-20C"
    print(f"  SVP(0C): flatau {es_flatau:.2f} Pa, bolton {es_bolton:.2f} Pa; agree to {rel*100:.2f}%  PASS")


def test_sat_mixrat_liq_consistency_and_grad():
    """sat_mixrat_liq uses the dispatcher esat and yields physical, finite, differentiable rsat."""
    T = jnp.linspace(240.0, 305.0, 40)
    p = jnp.full_like(T, 90000.0)
    rsat = sat_mixrat_liq(p, T, SATURATION_FLATAU)
    assert bool(jnp.all(jnp.isfinite(rsat))) and bool(jnp.all(rsat > 0.0))
    # rsat increases monotonically with T at fixed p
    assert bool(jnp.all(jnp.diff(rsat) > 0.0)), "rsat not monotonic in T"
    # ice rsat < liquid rsat below freezing
    Tc = jnp.full((10,), 260.0)
    pc = jnp.full((10,), 90000.0)
    assert float(sat_mixrat_ice(pc, Tc)[0]) < float(sat_mixrat_liq(pc, Tc, SATURATION_FLATAU)[0])
    g = jax.grad(lambda t: jnp.sum(sat_mixrat_liq(jnp.atleast_1d(90000.0),
                                                  jnp.atleast_1d(t), SATURATION_FLATAU)))(290.0)
    assert bool(jnp.isfinite(g)), "sat_mixrat_liq grad not finite"
    print("  sat_mixrat_liq monotonic + ice<liq + grad finite  PASS")


def test_sat_mixrat_f2py():
    """Bit-shadow the saturation mixing ratios vs the f2py Fortran oracle. The single elementwise JAX
    sat_mixrat_liq / sat_mixrat_ice cover both the Fortran rank-1 (_k) and rank-2 (_2D) generic-interface
    procedures; validate the 2D form against f2py_sat_mixrat_liq_2d (both saturation_formula leaves) and
    f2py_sat_mixrat_ice_2d. The leaf SVP polynomials match coefficient-for-coefficient, so this is bit-faithful
    (the ice path's slightly larger residual is the polynomial evaluation order). SKIPs if clubb_f2py is unbuilt.
    (iter 435)"""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py sat_mixrat oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(0)
    worst_liq = worst_ice = 0.0
    for formula in (SATURATION_FLATAU, SATURATION_BOLTON):
        for _ in range(15):
            ng, nz = 2, 12
            p = rng.uniform(2e4, 1.02e5, (ng, nz)); T = rng.uniform(220.0, 310.0, (ng, nz))
            j = np.asarray(sat_mixrat_liq(jnp.asarray(p), jnp.asarray(T), formula))
            f = np.asarray(clubb_f2py.f2py_sat_mixrat_liq_2d(p, T, formula))
            worst_liq = max(worst_liq, float(np.max(np.abs(j - f) / np.maximum(np.abs(f), 1e-30))))
    for _ in range(15):
        ng, nz = 2, 12
        p = rng.uniform(2e4, 1.02e5, (ng, nz)); T = rng.uniform(220.0, 273.0, (ng, nz))
        j = np.asarray(sat_mixrat_ice(jnp.asarray(p), jnp.asarray(T)))
        f = np.asarray(clubb_f2py.f2py_sat_mixrat_ice_2d(p, T, SATURATION_FLATAU))
        worst_ice = max(worst_ice, float(np.max(np.abs(j - f) / np.maximum(np.abs(f), 1e-30))))
    assert worst_liq < 1e-12, f"sat_mixrat_liq f2py rel mismatch {worst_liq:.2e}"
    assert worst_ice < 1e-11, f"sat_mixrat_ice f2py rel mismatch {worst_ice:.2e}"
    print(f"  f2py sat_mixrat_liq (flatau+bolton) + sat_mixrat_ice (30+15 cases): rel-match, "
          f"liq {worst_liq:.1e} / ice {worst_ice:.1e}  PASS")


def test_saturation_formula_enum_values():
    """The `saturation_<formula>` enum VALUES (BOLTON=1, FLATAU=3) select the SVP approximation, so a drifted value would
    silently use the wrong formula. Source-grounded: parses `saturation_<name> = <n>` straight from model_flags.F90 and
    checks the two the JAX ports (BOLTON, FLATAU) match. (GFDL=2 / LOOKUP=4 are unported.) SKIPs if the F90 is absent.
    (iter 471)"""
    import re
    f90 = os.path.join(_ROOT, "src", "CLUBB_core", "model_flags.F90")
    if not os.path.exists(f90):
        print("  saturation enum values vs Fortran: SKIP (model_flags.F90 absent)")
        return
    fort = {}
    for raw in open(f90):
        line = raw.split("!", 1)[0]
        m = re.match(r"\s*saturation_([A-Za-z0-9_]+)\s*=\s*([0-9]+)[\s,&]*$", line)
        if m:
            fort[m.group(1).lower()] = int(m.group(2))
    assert fort.get("bolton") and fort.get("flatau"), f"saturation enums not parsed from F90: {fort}"
    mism = []
    if SATURATION_BOLTON != fort["bolton"]:
        mism.append(f"BOLTON: JAX {SATURATION_BOLTON} vs Fortran {fort['bolton']}")
    if SATURATION_FLATAU != fort["flatau"]:
        mism.append(f"FLATAU: JAX {SATURATION_FLATAU} vs Fortran {fort['flatau']}")
    assert not mism, "saturation enum value(s) diverge from model_flags.F90:\n  " + "\n  ".join(mism)
    print(f"  saturation enum values match model_flags.F90 (BOLTON={SATURATION_BOLTON}, FLATAU={SATURATION_FLATAU})  PASS")


if __name__ == "__main__":
    print("saturation.py vs saturation.F90 logic:")
    test_dispatcher_matches_leaves()
    test_svp_reference_values()
    test_sat_mixrat_liq_consistency_and_grad()
    test_sat_mixrat_f2py()
    test_saturation_formula_enum_values()
    print("Done.")
