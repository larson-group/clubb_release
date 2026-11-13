"""Check read_sounding (sounding.F90) on small sounding files written for each test.

Each test writes a sounding into a temporary directory and builds a minimal grid,
then checks one behavior of the port: the interpolation, the surface values, a
fatal check, blank filling, scalar soundings or the DYCOMS II RF02 profiles.
"""
from pathlib import Path

from types import SimpleNamespace

import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.constants_clubb import g_per_kg
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.Input_fields.sounding import read_sounding, read_sounding_file

_HEADER = "z[m] thm[K] 'rt[kg/kg]' 'u[m/s]' 'v[m/s]' 'w[m/s]' 'ug[m/s]' 'vg[m/s]'\n"
_NO_SCALARS = SclrIdx(-1, -1, -1, -1, -1, -1)


def _write_sounding(directory, runtype="testcase", z=(0.0, 1000.0, 2000.0, 3000.0, 4000.0)):
    """Write a sounding whose columns are linear in z, so results can be computed by hand."""
    rows = "".join(
        f"{level} {300.0 + 0.005 * level} {0.01 - 1.0e-6 * level} {0.002 * level} "
        f"{-0.001 * level} 0.0 {0.002 * level} {-0.001 * level}\n"
        for level in z
    )
    (directory / f"{runtype}_sounding.in").write_text(_HEADER + rows)


def _grid(zt, zm_sfc=0.0):
    """Build the grid fields read_sounding uses, for one column."""
    return SimpleNamespace(zt=np.array([zt]), zm=np.array([[zm_sfc]]))


def _read(directory, gr, runtype="testcase", sclr_dim=0, edsclr_dim=0, sclr_idx=_NO_SCALARS,
          l_cubic=False):
    """Call read_sounding with the arguments these tests do not vary."""
    return read_sounding(1, sclr_dim, edsclr_dim, sclr_idx,
                         gr, 0, runtype, Path(directory) / (runtype + '_sounding.in'), Path(directory) / (runtype + '_sclr_sounding.in'), Path(directory) / (runtype + '_edsclr_sounding.in'), Path(directory) / (runtype + '_ozone_sounding.in'), '', 0.0,
                         1,
                         False, l_cubic,
                         1.0e5)


def test_linear_interpolation_and_surface_values(tmp_path):
    """Linear columns come back exact at new heights and at the model surface."""
    _write_sounding(tmp_path)
    zt = np.array([50.0, 1500.0, 3999.0])
    thlm, temperature_type, rtm, um, vm, ugm, vgm, altitude_type, press, subs_type, wm, \
        rtm_sfc, thlm_sfc, sclrm, edsclrm = _read(tmp_path, _grid(zt, zm_sfc=10.0))

    assert np.allclose(thlm[0], 300.0 + 0.005 * zt)
    assert np.allclose(rtm[0], 0.01 - 1.0e-6 * zt)
    assert np.allclose(um[0], 0.002 * zt)
    assert np.allclose(vgm[0], -0.001 * zt)
    assert (temperature_type, altitude_type, subs_type) == ("thm[K]", "z[m]", "w[m/s]")
    assert np.all(press == -999.9)
    assert np.allclose(rtm_sfc, 0.01 - 1.0e-6 * 10.0)
    assert np.allclose(thlm_sfc, 300.0 + 0.005 * 10.0)
    assert sclrm.shape == (1, 3, 0) and edsclrm.shape == (1, 3, 0)


def test_level_at_first_sounding_level_keeps_its_value(tmp_path):
    """A model level exactly at z(1) takes the first sounding value, without interpolating."""
    _write_sounding(tmp_path, z=(100.0, 1000.0, 2000.0))
    thlm, *_ = _read(tmp_path, _grid(np.array([100.0, 500.0])))
    assert thlm[0, 0] == 300.0 + 0.005 * 100.0


def test_surface_values_flag_a_sounding_above_the_surface(tmp_path):
    """rtm_sfc and thlm_sfc are -999 when the sounding starts above the surface."""
    _write_sounding(tmp_path, z=(100.0, 1000.0, 2000.0))
    result = _read(tmp_path, _grid(np.array([500.0]), zm_sfc=0.0))
    assert result[11][0] == -999.0 and result[12][0] == -999.0


def test_lowest_grid_level_below_sounding_stops(tmp_path):
    """The Fortran error stop fires when gr%zt(1,1) is below z(1)."""
    _write_sounding(tmp_path, z=(100.0, 1000.0))
    with pytest.raises(ValueError, match="Lowest level of input sounding"):
        _read(tmp_path, _grid(np.array([50.0, 500.0])))


def test_sounding_below_grid_top_stops(tmp_path):
    """The Fortran error stop fires when the sounding does not reach the model top."""
    _write_sounding(tmp_path, z=(0.0, 1000.0))
    with pytest.raises(ValueError, match="Not enough sounding data"):
        _read(tmp_path, _grid(np.array([500.0, 1500.0])))


def test_missing_sounding_file_stops(tmp_path):
    """A case with no sounding file stops with a clear message."""
    with pytest.raises(ValueError, match="Cannot open .*no_such_case_sounding.in"):
        _read(tmp_path, _grid(np.array([500.0])), runtype="no_such_case")


def test_readiopdata_is_rejected(tmp_path):
    """The unported IOP path fails instead of silently reading the sounding."""
    _write_sounding(tmp_path)
    with pytest.raises(ValueError, match="l_readiopdata"):
        read_sounding(1, 0, 0, _NO_SCALARS, _grid(np.array([500.0])), 0, "testcase", Path(tmp_path) / ("testcase" + '_sounding.in'), Path(tmp_path) / ("testcase" + '_sclr_sounding.in'), Path(tmp_path) / ("testcase" + '_edsclr_sounding.in'), Path(tmp_path) / ("testcase" + '_ozone_sounding.in'), '',
                      0.0, 1, True, False, 1.0e5)


def test_blanks_are_filled_before_pressure_conversion(tmp_path):
    """A blank thm[K] on a pressure sounding gives the altitudes of the filled sounding.

    Fortran fills blanks before read_z_profile, which needs a value at every level.
    """
    header = "Press[Pa] thm[K] 'rt[kg/kg]' 'u[m/s]' 'v[m/s]' 'omega[Pa/s]' 'ug[m/s]' 'vg[m/s]'\n"
    template = "{p} {th} 1.0e-5 1.0 1.0 0.0 1.0 1.0\n"
    p = (100000.0, 90000.0, 80000.0, 70000.0)
    # theta is linear in pressure, so filling the blank recovers this exact value.
    theta = [300.0 + (100000.0 - level) * 1.0e-4 for level in p]
    (tmp_path / "full_sounding.in").write_text(
        header + "".join(template.format(p=a, th=b) for a, b in zip(p, theta))
    )
    (tmp_path / "blank_sounding.in").write_text(
        header + "".join(template.format(p=a, th=-999.9 if a == 80000.0 else b)
                         for a, b in zip(p, theta))
    )
    z_full = read_sounding_file(0, Path(tmp_path) / 'full_sounding.in', 10000, 1.0e5, 1, 0.0)[1]
    z_blank = read_sounding_file(0, Path(tmp_path) / 'blank_sounding.in', 10000, 1.0e5, 1, 0.0)[1]
    assert np.allclose(z_blank, z_full)
    assert np.all(np.diff(z_blank) > 0.0)


def test_scalar_sounding_reads_only_the_columns_in_use(tmp_path):
    """sclr_dim = 1 reads the first of two columns, and the index check still applies."""
    _write_sounding(tmp_path)
    rows = "".join(f"{level} 7.0\n" for level in (1.0, 2.0, 3.0, 4.0, 5.0))
    (tmp_path / "testcase_sclr_sounding.in").write_text("'rt[kg/kg]' thm[K]\n" + rows)
    gr = _grid(np.array([500.0]))

    sclrm = _read(tmp_path, gr, sclr_dim=1, sclr_idx=SclrIdx(1, -1, -1, -1, -1, -1))[13]
    assert sclrm.shape == (1, 1, 1)
    assert np.allclose(sclrm[0, 0, 0], 1.5)

    with pytest.raises(ValueError, match="iisclr_rt index does not match column"):
        _read(tmp_path, gr, sclr_dim=1, sclr_idx=SclrIdx(2, -1, -1, -1, -1, -1))


def test_cubic_interpolation_reproduces_linear_columns(tmp_path):
    """Steffen's method is exact on linear data, including near the top of the sounding."""
    _write_sounding(tmp_path)
    zt = np.array([50.0, 1500.0, 2500.0, 3999.0])
    thlm, _, rtm, *_ = _read(tmp_path, _grid(zt), l_cubic=True)
    assert np.allclose(thlm[0], 300.0 + 0.005 * zt)
    assert np.allclose(rtm[0], 0.01 - 1.0e-6 * zt)


def test_dycoms2_rf02_uses_the_analytic_profiles(tmp_path):
    """DYCOMS II RF02 ignores the sounding's thermodynamics and winds (Wyant et al. 2007)."""
    _write_sounding(tmp_path, runtype="dycoms2_rf02")
    zt = np.array([500.0, 1000.0])
    thlm, _, rtm, um, vm, ugm, vgm, _, _, _, wm, *_ = _read(tmp_path, _grid(zt), runtype="dycoms2_rf02")

    assert np.allclose(um[0], 3.0 + 4.3 * zt / 1000.0)
    assert np.allclose(vm[0], -9.0 + 5.6 * zt / 1000.0)
    assert np.array_equal(ugm, um) and np.array_equal(vgm, vm)
    assert np.allclose(thlm[0], [288.3, 295.0 + (1000.0 - 795.0)**(1.0 / 3.0)])
    assert np.allclose(rtm[0], [9.45 / g_per_kg,
                                (5.0 - 3.0 * (1.0 - np.exp((795.0 - 1000.0) / 500.0))) / g_per_kg])
    assert np.all(wm == 0.0)


def test_explicit_sounding_filename_is_independent_of_case(tmp_path):
    """Read a profile whose filename differs from the scientific runtype."""
    _write_sounding(tmp_path)
    selected = tmp_path / "observed profile.dat"
    (tmp_path / "testcase_sounding.in").rename(selected)
    result = read_sounding(1, 0, 0, _NO_SCALARS,
                           _grid(np.array([500.0])), 0, "bomex", selected, '', '', '', '',
                           0.0, 1, False, False, 1.0e5)
    assert np.allclose(result[0], [[302.5]])
