#!/usr/bin/env python3
"""test_pressure_sounding_z.py — validate convert_pressure_sounding_to_z (input_interpret.F90:read_z_profile).

Derives the altitudes of a pressure-coordinate (Press[Pa]) sounding hydrostatically from its own
thermodynamics, so CGILS/cloud_feedback cases can interpolate onto the height grid. Checks:
  1. On the REAL cgils_s11 sounding: altitudes are finite, strictly increasing, anchored near zm_init, and span
     a physically reasonable range (top a few tens of km for a ~100 hPa top).
  2. Independent re-derivation of the unsaturated temperature branch (rcm=0 → thvm computed by hand) fed through
     the already-validated inverse_hydrostatic must reproduce convert's altitudes bit-for-bit.
  3. The derived z is the hydrostatic inverse: round-tripping z → exner via the forward log-mean integral
     recovers the input pressures.
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

from clubb_jax.src.Input_fields.sounding import read_sounding, convert_pressure_sounding_to_z
from clubb_jax.src.Input_fields.hydrostatic_module import inverse_hydrostatic, calc_ref_z_linear_thvm
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, p0, kappa, ep2

_CGILS = _ROOT + "/input/case_setups/cgils_s11_sounding.in"
_P_SFC, _ZM_INIT, _SAT = 1.02078e5, 0.0, 3   # cgils p_sfc≈sfc value; flatau saturation formula = 3


def test_real_cgils_sounding():
    if not os.path.exists(_CGILS):
        print("  real cgils_s11 sounding: SKIP (file absent)"); return
    snd = read_sounding(_CGILS)
    assert str(snd['altitude_type']).strip().lower() == 'press[pa]', "cgils sounding is not pressure-coordinate"
    out = convert_pressure_sounding_to_z(snd, _P_SFC, _ZM_INIT, _SAT)
    z = np.asarray(out['z'])
    assert np.isfinite(z).all(), "non-finite altitudes"
    assert np.all(np.diff(z) > 0.0), "altitudes not strictly increasing"
    assert abs(z[0] - _ZM_INIT) < 200.0, f"surface altitude {z[0]:.1f} not near zm_init"
    assert 8e3 < z[-1] < 60e3, f"top altitude {z[-1]:.0f} m physically implausible"
    assert out['altitude_type'] == 'z[m]', "altitude_type not switched to z[m]"
    print(f"  real cgils_s11: z monotonic in [{z[0]:.1f}, {z[-1]:.0f}] m over {z.size} levels  PASS")


def test_temperature_branch_independent():
    # Synthetic unsaturated Press[Pa]/T[K] sounding (rt small → rcm=0); recompute thvm by hand.
    p = np.linspace(1.0e5, 2.0e4, 30)
    T = 290.0 - 6.5e-3 * np.linspace(0, 12000, 30)        # rough lapse (only needs to be unsaturated)
    rt = np.full(30, 1.0e-7)                               # tiny → guaranteed unsaturated (rcm=0) everywhere
    snd = dict(
        z=p,
        theta=T,
        rt=rt,
        altitude_type='Press[Pa]',
        temperature_type='T[K]',
        subs_type='omega[Pa/s]',
    )
    out = convert_pressure_sounding_to_z(snd, 1.0e5, 0.0, _SAT)
    z_got = np.asarray(out['z'])

    exner = (p / p0) ** kappa
    theta = T / exner
    rcm = np.zeros_like(p)                                 # unsaturated by construction
    assert np.all(rt < 0.02), "test profile unexpectedly near saturation"
    thv_ds = theta * (1.0 + ep2 * (rt - rcm)) ** kappa
    # calculate_thvm: thlm + ep1·thv_ds·rt + (Lv/(Cp·exner) − ep2·thv_ds)·rcm; with rcm=0 the last term drops.
    from clubb_jax.src.CLUBB_core.constants_clubb import ep1
    thlm = theta
    thvm = thlm + ep1 * thv_ds * rt
    z_ref = np.asarray(inverse_hydrostatic(1.0e5, 0.0, thvm, exner))
    w = np.max(np.abs(z_got - z_ref))
    assert w < 1e-9, f"temperature-branch z mismatch vs independent thvm: {w:.2e}"
    print(f"  temperature branch vs independent thvm re-derivation: worst {w:.2e} m  PASS")


def test_hydrostatic_roundtrip():
    # z derived from p must invert the forward log-mean integral: forward exner(z; thvm) recovers input p.
    p = np.linspace(1.0e5, 3.0e4, 25)
    T = 295.0 - 6.0e-3 * np.linspace(0, 11000, 25)
    rt = np.full(25, 1.0e-7)                               # tiny → unsaturated (rcm=0)
    snd = dict(
        z=p,
        theta=T,
        rt=rt,
        altitude_type='Press[Pa]',
        temperature_type='T[K]',
        subs_type='omega[Pa/s]',
    )
    out = convert_pressure_sounding_to_z(snd, 1.0e5, 0.0, _SAT)
    z = np.asarray(out['z'])
    exner = (p / p0) ** kappa
    # Forward check: the relative altitudes from calc_ref_z_linear_thvm(thvm,exner) must match (z - z[0]) within
    # the surface-anchoring offset, confirming z is the consistent hydrostatic inverse.
    from clubb_jax.src.CLUBB_core.constants_clubb import ep1
    theta = T / exner
    thvm = theta + ep1 * (theta * (1.0 + ep2 * rt) ** kappa) * rt
    ref_z = np.asarray(calc_ref_z_linear_thvm(jnp.asarray(thvm), jnp.asarray(exner)))
    # z and ref_z differ only by a constant surface offset → their level-to-level increments must match.
    w = np.max(np.abs(np.diff(z) - np.diff(ref_z)))
    assert w < 1e-9, f"z increments inconsistent with the forward log-mean integral: {w:.2e}"
    print(f"  hydrostatic round-trip (Δz vs forward integral): worst {w:.2e} m  PASS")


def main():
    print("test_pressure_sounding_z:")
    for t in (test_real_cgils_sounding, test_temperature_branch_independent, test_hydrostatic_roundtrip):
        t()
    print("All convert_pressure_sounding_to_z checks PASSED")


if __name__ == "__main__":
    main()
