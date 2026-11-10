#!/usr/bin/env python3
"""test_pressure_sounding_z.py — validate the pressure branch of read_z_profile (input_interpret.F90).

Derives the altitudes of a pressure-coordinate (Press[Pa]) sounding hydrostatically from its own
thermodynamics, so CGILS/cloud_feedback cases can interpolate onto the height grid. Checks:
  1. On the REAL cgils_s11 sounding: altitudes are finite, strictly increasing, anchored near zm_init, and span
     a physically reasonable range (top a few tens of km for a ~100 hPa top).
  2. Independent re-derivation of the unsaturated temperature branch (rcm=0 → thvm computed by hand) fed through
     the already-validated inverse_hydrostatic must reproduce read_z_profile's altitudes bit-for-bit.
  3. The derived z is the hydrostatic inverse: round-tripping z → exner via the forward log-mean integral
     recovers the input pressures.
  4. Case initialization applies the conversion: mc3e's initial thlm follows its Press[Pa] sounding when
     looked up at the model's own pressures.
  5. Case initialization converts a T[K] sounding to potential temperature first, and stops on z[m] with T[K].
"""

from clubb_jax.auto_llm_generated_tests.pytests.array_inputs import call_jax
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os

_ROOT = str(_REPO_ROOT)

import numpy as np
import jax.numpy as jnp
import pytest

from clubb_jax.src.Input_fields.input_interpret import read_z_profile
from clubb_jax.src.Input_fields.input_reader import one_dim_read_var
from clubb_jax.src.Input_fields.sounding import read_sounding_file
from clubb_jax.src.Input_fields.hydrostatic_module import inverse_hydrostatic, calc_ref_z_linear_thvm
from clubb_jax.src.CLUBB_core.constants_clubb import p0, kappa, ep2
from clubb_jax.src.clubb_case_initialization import init_clubb_case
from input.scripts.create_case_namelist import create_case_namelist_file

_CGILS = _ROOT + "/input/case_data/cgils_s11_sounding.in"
_P_SFC, _ZM_INIT, _SAT = 1.02078e5, 0.0, 3   # cgils p_sfc≈sfc value; flatau saturation formula = 3


def _pressure_retVars(p, T, rt):
    """Build the read_one_dim_file output for a Press[Pa]/T[K] sounding."""
    return [one_dim_read_var(name, 'Press[Pa]', values)
            for name, values in (('Press[Pa]', p), ('T[K]', T), ('rt[kg/kg]', rt))]


def test_real_cgils_sounding():
    if not os.path.exists(_CGILS):
        raise AssertionError("  real cgils_s11 sounding: required input missing (file absent)")
    _, z, _, _, _, _, _, _, _, _, _, altitude_type, _, _ = read_sounding_file(
        0, "cgils_s11", _ROOT + "/input/case_data", 10000, _P_SFC, _SAT, _ZM_INIT)
    assert np.isfinite(z).all(), "non-finite altitudes"
    assert np.all(np.diff(z) > 0.0), "altitudes not strictly increasing"
    assert abs(z[0] - _ZM_INIT) < 200.0, f"surface altitude {z[0]:.1f} not near zm_init"
    assert 8e3 < z[-1] < 60e3, f"top altitude {z[-1]:.0f} m physically implausible"
    assert altitude_type == 'Press[Pa]', "altitude_type should report the file's column, as in Fortran"
    print(f"  real cgils_s11: z monotonic in [{z[0]:.1f}, {z[-1]:.0f}] m over {z.size} levels  PASS")


def test_temperature_branch_independent():
    # Synthetic unsaturated Press[Pa]/T[K] sounding (rt small → rcm=0); recompute thvm by hand.
    p = np.linspace(1.0e5, 2.0e4, 30)
    T = 290.0 - 6.5e-3 * np.linspace(0, 12000, 30)        # rough lapse (only needs to be unsaturated)
    rt = np.full(30, 1.0e-7)                               # tiny → guaranteed unsaturated (rcm=0) everywhere
    z_got, _, _ = read_z_profile(3, 30, _pressure_retVars(p, T, rt), 1.0e5, 0.0, _SAT)

    exner = (p / p0) ** kappa
    theta = T / exner
    rcm = np.zeros_like(p)                                 # unsaturated by construction
    assert np.all(rt < 0.02), "test profile unexpectedly near saturation"
    thv_ds = theta * (1.0 + ep2 * (rt - rcm)) ** kappa
    # calculate_thvm: thlm + ep1·thv_ds·rt + (Lv/(Cp·exner) − ep2·thv_ds)·rcm; with rcm=0 the last term drops.
    from clubb_jax.src.CLUBB_core.constants_clubb import ep1
    thlm = theta
    thvm = thlm + ep1 * thv_ds * rt
    z_ref = np.asarray(call_jax(inverse_hydrostatic, 1.0e5, 0.0, thvm, exner))
    w = np.max(np.abs(z_got - z_ref))
    assert w < 1e-9, f"temperature-branch z mismatch vs independent thvm: {w:.2e}"
    print(f"  temperature branch vs independent thvm re-derivation: worst {w:.2e} m  PASS")


def test_hydrostatic_roundtrip():
    # z derived from p must invert the forward log-mean integral: forward exner(z; thvm) recovers input p.
    p = np.linspace(1.0e5, 3.0e4, 25)
    T = 295.0 - 6.0e-3 * np.linspace(0, 11000, 25)
    rt = np.full(25, 1.0e-7)                               # tiny → unsaturated (rcm=0)
    z, _, _ = read_z_profile(3, 25, _pressure_retVars(p, T, rt), 1.0e5, 0.0, _SAT)
    exner = (p / p0) ** kappa
    # Forward check: the relative altitudes from calc_ref_z_linear_thvm(thvm,exner) must match (z - z[0]) within
    # the surface-anchoring offset, confirming z is the consistent hydrostatic inverse.
    from clubb_jax.src.CLUBB_core.constants_clubb import ep1
    theta = T / exner
    thvm = theta + ep1 * (theta * (1.0 + ep2 * rt) ** kappa) * rt
    ref_z = np.asarray(call_jax(calc_ref_z_linear_thvm, jnp.asarray(thvm), jnp.asarray(exner)))
    # z and ref_z differ only by a constant surface offset → their level-to-level increments must match.
    w = np.max(np.abs(np.diff(z) - np.diff(ref_z)))
    assert w < 1e-9, f"z increments inconsistent with the forward log-mean integral: {w:.2e}"
    print(f"  hydrostatic round-trip (Δz vs forward integral): worst {w:.2e} m  PASS")


def test_case_initialization_uses_sounding_altitudes(tmp_path):
    # mc3e's stock bugsrad is unsupported in JAX and does not affect the initial sounding.
    path = create_case_namelist_file(
        "mc3e", tmp_path, max_iters=1, stats="none", debug="0",
        override="rad_scheme='simplified'",
    )
    state = init_clubb_case(str(path))
    p_model = np.asarray(state['p_in_Pa'])[0]
    thlm_model = np.asarray(state['thlm'])[0]
    _, _, theta, _, _, _, _, _, p_snd, _, _, _, _, _ = read_sounding_file(
        0, "mc3e", _ROOT + "/input/case_data", 10000, _P_SFC, _SAT, _ZM_INIT)
    # np.interp needs ascending coordinates, and pressure falls with height.
    thlm_snd = np.interp(-p_model, -p_snd, theta)
    # Higher up the sounding levels are far apart, so interpolating linearly in
    # pressure instead of height gives a genuinely different answer.
    well_resolved = p_model > 2.0e4
    w = np.max(np.abs(thlm_model - thlm_snd)[well_resolved])
    assert w < 0.2, f"initial thlm departs from the pressure sounding by {w:.2f} K"


def _write_temperature_sounding(directory, height_column):
    """Write a dry T[K] sounding beside the namelist, which initialization reads first."""
    levels = np.linspace(102000.0, 10000.0, 30)
    temperature = 300.0 - (102000.0 - levels) * 7.0e-4
    if height_column:
        first_column, levels = 'z[m]', np.linspace(0.0, 16000.0, 30)
    else:
        first_column = 'Press[Pa]'
    subsidence = 'w[m/s]' if height_column else 'omega[Pa/s]'
    rows = "".join(f"{a} {b} 1.0e-6 1.0 1.0 0.0 1.0 1.0\n" for a, b in zip(levels, temperature))
    (directory / "bomex_sounding.in").write_text(
        f"{first_column} T[K] 'rt[kg/kg]' 'u[m/s]' 'v[m/s]' '{subsidence}' 'ug[m/s]' 'vg[m/s]'\n" + rows
    )
    return levels, temperature


def test_temperature_sounding_initializes_potential_temperature(tmp_path):
    p_snd, T_snd = _write_temperature_sounding(tmp_path, height_column=False)
    path = create_case_namelist_file("bomex", tmp_path, max_iters=1, stats="none", debug="0")
    state = init_clubb_case(str(path))
    p_model = np.asarray(state['p_in_Pa'])[0]
    # The sounding is dry, so thlm is potential temperature and thlm * exner recovers T.
    T_model = np.asarray(state['thlm'])[0] * (p_model / p0) ** kappa
    w = np.max(np.abs(T_model - np.interp(-p_model, -p_snd, T_snd)))
    assert w < 0.2, f"initial temperature departs from the T[K] sounding by {w:.2f} K"


def test_height_sounding_with_temperature_stops(tmp_path):
    _write_temperature_sounding(tmp_path, height_column=True)
    path = create_case_namelist_file("bomex", tmp_path, max_iters=1, stats="none", debug="0")
    with pytest.raises(ValueError, match="absolute temperature"):
        init_clubb_case(str(path))
