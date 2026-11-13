#!/usr/bin/env python3
"""test_rad_extended_atmosphere.py — validate the case-specific radiation extended atmosphere (Iter92).

When l_use_default_std_atmosphere=.false. (CGILS/cloud_feedback/astex/twp_ice) the radiation extended atmosphere
above the model top is built from the case's OWN deep sounding + {case}_ozone_sounding.in, not the default
US-standard atmosphere (sounding.F90:convert_snd2extended_atm; Iter92 port). This guards that leaf code — it
affects all 12 CGILS/cloud_feedback cases and previously had only end-to-end validation. Checks:
  1. On the REAL cgils_s11 sounding + ozone sounding: structure matches convert_snd2extended_atm (alt ascending =
     the converted sounding altitudes; T_in_K = the T[K] column verbatim; sp_hmdty = rt/(1+rt); p_in_mb = p/100;
     o3l = the ozone column; all lengths equal; physical sanity).
  2. The thm[K] / thlm[K] branch: T_in_K = θ·exner with the level-1 exner from p_sfc (the Fortran quirk).
"""
from pathlib import Path
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os

_ROOT = str(_REPO_ROOT)

import numpy as np

from clubb_jax.src.Input_fields.input_reader import one_dim_read_var, read_one_dim_file
from clubb_jax.src.Input_fields.sounding import read_sounding_file, convert_snd2extended_atm
from clubb_jax.src.CLUBB_core.constants_clubb import p0, kappa

_CASE_DATA = _ROOT + "/input/case_data"
_P_SFC, _SAT = 1.02078e5, 3


def test_real_cgils_extended_atmosphere():
    snd_f, oz_f = _CASE_DATA + "/cgils_s11_sounding.in", _CASE_DATA + "/cgils_s11_ozone_sounding.in"
    if not (os.path.exists(snd_f) and os.path.exists(oz_f)):
        raise AssertionError("  real cgils extended atmosphere: required input missing (files absent)")
    n, z, theta, rt, _, _, _, _, p_in_Pa, _, _, _, _, retVars = read_sounding_file(
        0, Path(_CASE_DATA) / 'cgils_s11_sounding.in', 10000, _P_SFC, _SAT, 0.0)
    o3l = read_one_dim_file(0, 1, oz_f)[0].values
    assert len(o3l) == n, f"ozone levels {len(o3l)} != sounding levels {n}"

    ext = convert_snd2extended_atm(0, Path(_CASE_DATA) / 'cgils_s11_ozone_sounding.in', 8, _P_SFC, 0.0, retVars, _SAT)
    # Structure / convert_snd2extended_atm semantics
    for k in ('alt', 'T_in_K', 'sp_hmdty', 'p_in_mb', 'o3l'):
        assert k in ext and len(ext[k]) == n, f"ext['{k}'] missing or wrong length"
    assert np.all(np.diff(ext['alt']) > 0), "extended altitudes not ascending"
    assert np.allclose(ext['alt'], z), "extended altitudes != converted sounding altitudes"
    assert np.allclose(ext['T_in_K'], theta), "T[K] column not copied verbatim"
    assert np.allclose(ext['sp_hmdty'], rt / (1.0 + rt)), "sp_hmdty != rt/(1+rt)"
    assert np.allclose(ext['p_in_mb'], p_in_Pa / 100.0), "p_in_mb != p/100"
    assert np.allclose(ext['o3l'], o3l), "o3l != ozone sounding"
    # Physical sanity: reaches the stratosphere, pressure decreasing with altitude, ozone positive.
    assert ext['alt'][-1] > 25e3, "extension does not reach the stratosphere (>25 km)"
    assert np.all(np.diff(ext['p_in_mb']) < 0), "pressure not monotonically decreasing with altitude"
    assert np.all(ext['o3l'] > 0) and np.all(ext['T_in_K'] > 150), "non-physical ozone/temperature"
    print(f"  real cgils_s11 extended atmosphere: {n} levels, alt→{ext['alt'][-1]/1e3:.1f} km, structure OK  PASS")


def test_theta_branch_exner_conversion(tmp_path):
    # Synthetic thm[K] (potential-temperature) sounding: T_in_K must be θ·exner, exner[0] from p_sfc.
    p = np.linspace(1.0e5, 2.0e4, 12)
    theta = np.full(12, 300.0)                              # constant θ
    rt = np.full(12, 1.0e-3)
    p_sfc = 1.01e5
    retVars = [one_dim_read_var(name, 'Press[Pa]', values)
               for name, values in (('Press[Pa]', p), ('thm[K]', theta), ('rt[kg/kg]', rt))]
    (tmp_path / "synthetic_ozone_sounding.in").write_text(
        "'o3[kg/kg]'\n" + "5.0e-8\n" * 12
    )
    ext = convert_snd2extended_atm(0, Path(tmp_path) / 'synthetic_ozone_sounding.in', 3, p_sfc, 0.0, retVars, _SAT)
    exner = (p / p0) ** kappa
    exner[0] = (p_sfc / p0) ** kappa                        # Fortran uses p_sfc for level 1
    assert np.allclose(ext['T_in_K'], theta * exner), "thm[K] branch: T_in_K != θ·exner (or wrong level-1 exner)"
    # level-1 exner uses p_sfc, not p[0]: verify it differs from the naive (p[0]/p0)**kappa when p_sfc != p[0]
    assert abs(ext['T_in_K'][0] - theta[0] * (p[0] / p0) ** kappa) > 0, "level-1 exner did not use p_sfc"
    print("  thm[K] branch: T_in_K = θ·exner with p_sfc level-1 exner  PASS")
