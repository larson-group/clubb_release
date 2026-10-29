"""Validate the JAX saturation port (saturation.py) against saturation.F90 logic."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os
import pytest

import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

_ROOT = str(_REPO_ROOT)

from clubb_jax.src.CLUBB_core.saturation import (
    sat_vapor_press_liq, sat_vapor_press_liq_flatau, sat_vapor_press_liq_bolton,
    sat_vapor_press_liq_gfdl, SATURATION_GFDL,
    sat_mixrat_liq, sat_mixrat_ice, SATURATION_FLATAU, SATURATION_BOLTON,
)

_T_FREEZE = 273.15


def test_dispatcher_matches_leaves():
    """sat_vapor_press_liq routes bit-exactly to the leaf chosen by saturation_formula."""
    T = jnp.linspace(220.0, 310.0, 64)
    for formula, leaf in ((SATURATION_FLATAU, sat_vapor_press_liq_flatau),
                          (SATURATION_BOLTON, sat_vapor_press_liq_bolton),
                          (SATURATION_GFDL, sat_vapor_press_liq_gfdl)):
        d = float(jnp.max(jnp.abs(sat_vapor_press_liq(T, formula) - leaf(T))))
        assert d == 0.0, f"formula {formula}: dispatcher != leaf ({d})"
    from clubb_jax.src.CLUBB_core.model_flags import saturation_lookup
    for bad in (saturation_lookup, 999):
        with pytest.raises(ValueError):
            sat_vapor_press_liq(T, bad)
    print("  dispatcher matches Flatau/Bolton/GFDL and rejects lookup/unknown formulas  PASS")


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


def test_saturation_formula_enum_values():
    """The `saturation_<formula>` enum VALUES (BOLTON=1, FLATAU=3) select the SVP approximation, so a drifted value would
    silently use the wrong formula. Source-grounded: parses `saturation_<name> = <n>` straight from model_flags.F90 and
    checks the two the JAX ports (BOLTON, FLATAU) match. (LOOKUP=4 remains unsupported.) The Fortran source is required.
    (iter 471)"""
    import re
    f90 = os.path.join(_ROOT, "src", "CLUBB_core", "model_flags.F90")
    if not os.path.exists(f90):
        raise AssertionError("  saturation enum values vs Fortran: required input missing (model_flags.F90 absent)")
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
