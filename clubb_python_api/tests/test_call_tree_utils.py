"""Test wrappers for utility routines in advance_clubb_core call tree."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=11, dz=500.0):
    """Create a standard test grid and return it."""
    clubb_api.init_err_info(ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(flags)

    zm = np.arange(nzmax, dtype=np.float64) * dz
    zt = 0.5 * (zm[:-1] + zm[1:])

    gr, _err_info = clubb_api.setup_grid(
        nzmax=nzmax, ngrdcol=ngrdcol,
        sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
        l_implemented=True, l_ascending_grid=True,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        zm_init=np.zeros(ngrdcol, dtype=np.float64),
        zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
        thermodynamic_heights=np.tile(zt, (ngrdcol, 1)),
        err_info=ErrInfo(ngrdcol=ngrdcol),
    )
    return gr


@pytest.fixture(scope="module")
def gr():
    """Shared test grid."""
    return _make_grid()


def test_julian2gregorian_date_known():
    """Julian day 2451545 corresponds to 2000-01-01."""
    day, month, year = clubb_api.julian2gregorian_date(2451545)
    assert int(day) == 1
    assert int(month) == 1
    assert int(year) == 2000


def test_gregorian2julian_day_and_leap_year_known():
    """Gregorian helper wrappers should match known leap-year behavior."""
    assert clubb_api.leap_year(2000) is True
    assert clubb_api.leap_year(1900) is False
    assert clubb_api.gregorian2julian_day(1, 1, 2000) == 1
    assert clubb_api.gregorian2julian_day(31, 12, 2000) == 366
    assert clubb_api.gregorian2julian_day(31, 12, 2001) == 365


def test_lin_interpolate_two_points_midpoint():
    """Linear interpolation should match expected analytic value."""
    val = clubb_api.lin_interpolate_two_points(
        height_int=2.5,
        height_high=10.0,
        height_low=0.0,
        var_high=40.0,
        var_low=20.0,
    )
    np.testing.assert_allclose(val, 25.0)


def test_calculate_thvm_basic(gr):
    """calculate_thvm should return finite zt-grid output with correct shape."""
    shape = (gr.ngrdcol, gr.nzt)
    thvm = clubb_api.calculate_thvm(
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        thlm=np.full(shape, 300.0),
        rtm=np.full(shape, 0.01),
        rcm=np.zeros(shape),
        exner=np.full(shape, 1.0),
        thv_ds_zt=np.full(shape, 300.0),
    )
    assert thvm.shape == shape
    assert np.all(np.isfinite(thvm))


def test_calc_ri_zm_matches_formula(gr):
    """calc_ri_zm should apply clipping and ratio formula elementwise."""
    shape = (gr.ngrdcol, gr.nzm)
    bv = np.linspace(0.0, 2.0, gr.nzm, dtype=np.float64).reshape(1, -1)
    shear = np.linspace(0.1, 1.0, gr.nzm, dtype=np.float64).reshape(1, -1)
    ri = clubb_api.calc_ri_zm(
        nzm=gr.nzm, ngrdcol=gr.ngrdcol,
        bv_freq_sqd=bv, shear=shear, lim_bv=0.25, lim_shear=0.2
    )
    expected = np.maximum(bv, 0.25) / np.maximum(shear, 0.2)
    assert ri.shape == shape
    np.testing.assert_allclose(ri, expected)


def test_vertical_avg_and_integral():
    """vertical_avg and vertical_integral should match direct computation."""
    rho_ds = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    field = np.array([10.0, 20.0, 30.0], dtype=np.float64)
    dz = np.array([2.0, 2.0, 2.0], dtype=np.float64)

    avg = clubb_api.vertical_avg(total_idx=rho_ds.shape[0], rho_ds=rho_ds, field=field, dz=dz)
    integ = clubb_api.vertical_integral(total_idx=rho_ds.shape[0], rho_ds=rho_ds, field=field, dz=dz)

    expected_integral = float(np.sum(rho_ds * field * dz))
    expected_avg = expected_integral / float(np.sum(rho_ds * dz))

    np.testing.assert_allclose(integ, expected_integral)
    np.testing.assert_allclose(avg, expected_avg)


def test_zm2zt2zm_preserves_constant_field(gr):
    """A constant profile should remain constant through zm->zt->zm smoothing."""
    azm = np.full((gr.ngrdcol, gr.nzm), 7.5, dtype=np.float64)
    azm_smoothed = clubb_api.zm2zt2zm(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, azm=azm, zm_min=0.0
    )
    np.testing.assert_allclose(azm_smoothed, azm)


def test_zt2zm2zt_preserves_constant_field(gr):
    """A constant profile should remain constant through zt->zm->zt smoothing."""
    azt = np.full((gr.ngrdcol, gr.nzt), 3.25, dtype=np.float64)
    azt_smoothed = clubb_api.zt2zm2zt(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, azt=azt, zt_min=0.0
    )
    np.testing.assert_allclose(azt_smoothed, azt)


def test_ddzm_matches_grid_spacing_formula(gr):
    """ddzm should compute first differences across zt with invrs_dzt scaling."""
    azm = np.tile(np.arange(gr.nzm, dtype=np.float64), (gr.ngrdcol, 1))
    out = clubb_api.ddzm(gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, azm=azm)
    expected = (azm[:, 1:] - azm[:, :-1]) * gr.invrs_dzt
    assert out.shape == (gr.ngrdcol, gr.nzt)
    np.testing.assert_allclose(out, expected)


def test_ddzt_matches_grid_spacing_formula(gr):
    """ddzt should compute first differences across zm with edge replication."""
    azt = np.tile(np.arange(gr.nzt, dtype=np.float64), (gr.ngrdcol, 1))
    out = clubb_api.ddzt(gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, azt=azt)
    expected = np.empty((gr.ngrdcol, gr.nzm), dtype=np.float64)
    expected[:, 1:-1] = (azt[:, 1:] - azt[:, :-1]) * gr.invrs_dzm[:, 1:-1]
    expected[:, 0] = expected[:, 1]
    expected[:, -1] = expected[:, -2]
    assert out.shape == (gr.ngrdcol, gr.nzm)
    np.testing.assert_allclose(out, expected)


def test_pvertinterp_hits_exact_level(gr):
    """Interpolating exactly at a native level should return that level value."""
    p_mid_1d = np.linspace(100000.0, 10000.0, gr.nzt, dtype=np.float64)
    p_mid = np.tile(p_mid_1d, (gr.ngrdcol, 1))
    input_var = np.tile(np.linspace(0.0, 9.0, gr.nzt, dtype=np.float64), (gr.ngrdcol, 1))

    target_k = 4
    p_out = float(p_mid[0, target_k])
    interp = clubb_api.pvertinterp(
        gr=gr, nzt=gr.nzt, ngrdcol=gr.ngrdcol, p_mid=p_mid, p_out=p_out, input_var=input_var
    )

    assert interp.shape == (gr.ngrdcol,)
    np.testing.assert_allclose(interp, input_var[:, target_k])


def test_smooth_max_supports_scalar_array_and_array_scalar(gr):
    """smooth_max should match the analytic smoothed max for supported 2D forms."""
    arr = np.array([[1.0, -1.0, 3.0]], dtype=np.float64)
    coef = 0.25
    ngrdcol, nz = arr.shape

    scalar_array = clubb_api.smooth_max(nz=nz, ngrdcol=ngrdcol, input_var1=0.5, input_var2=arr, smth_coef=coef)
    array_scalar = clubb_api.smooth_max(nz=nz, ngrdcol=ngrdcol, input_var1=arr, input_var2=0.5, smth_coef=coef)

    expected = 0.5 * ((arr + 0.5) + np.sqrt((arr - 0.5) ** 2 + coef ** 2))
    np.testing.assert_allclose(scalar_array, expected)
    np.testing.assert_allclose(array_scalar, expected)


def test_smooth_min_supports_scalar_array_and_array_scalar(gr):
    """smooth_min should match the analytic smoothed min for supported 2D forms."""
    arr = np.array([[1.0, -1.0, 3.0]], dtype=np.float64)
    coef = 0.25
    ngrdcol, nz = arr.shape

    scalar_array = clubb_api.smooth_min(nz=nz, ngrdcol=ngrdcol, input_var1=0.5, input_var2=arr, smth_coef=coef)
    array_scalar = clubb_api.smooth_min(nz=nz, ngrdcol=ngrdcol, input_var1=arr, input_var2=0.5, smth_coef=coef)

    expected = 0.5 * ((arr + 0.5) - np.sqrt((arr - 0.5) ** 2 + coef ** 2))
    np.testing.assert_allclose(scalar_array, expected)
    np.testing.assert_allclose(array_scalar, expected)


def test_read_grid_heights_reads_momentum_grid_file(tmp_path):
    """read_grid_heights should read momentum heights for grid_type=3."""
    ngrdcol = 1
    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    path = tmp_path / "zm.grd"
    path.write_text("0\n100\n250\n450\n")

    pulled_err, momentum_heights, thermodynamic_heights = clubb_api.read_grid_heights(
        nzmax=4,
        grid_type=3,
        zm_grid_fname=str(path),
        zt_grid_fname="",
        file_unit=98,
        err_info=err_info,
    )

    np.testing.assert_allclose(momentum_heights, np.array([0.0, 100.0, 250.0, 450.0], dtype=np.float64))
    np.testing.assert_allclose(thermodynamic_heights, np.zeros(3, dtype=np.float64))
    assert np.all(pulled_err.err_code == 0)


def test_setup_grid_heights_updates_stored_grid(tmp_path):
    """setup_grid_heights should recompute grid heights and interpolation weights."""
    gr = _make_grid(ngrdcol=1, nzmax=4, dz=100.0)
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)
    momentum = np.array([[0.0, 100.0, 250.0, 450.0]], dtype=np.float64)
    thermo = np.zeros((gr.ngrdcol, gr.nzt), dtype=np.float64)

    updated_gr, pulled_err = clubb_api.setup_grid_heights(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        l_implemented=False,
        l_ascending_grid=True,
        grid_type=3,
        deltaz=np.zeros(gr.ngrdcol, dtype=np.float64),
        zm_init=np.zeros(gr.ngrdcol, dtype=np.float64),
        momentum_heights=momentum,
        thermodynamic_heights=thermo,
        err_info=err_info,
    )

    np.testing.assert_allclose(updated_gr.zm, momentum)
    np.testing.assert_allclose(updated_gr.zt, np.array([[50.0, 175.0, 350.0]], dtype=np.float64))
    assert np.all(pulled_err.err_code == 0)


def test_calc_xpwp_matches_centered_difference_formula(gr):
    """calc_xpwp should apply Km * invrs_dzm * backward thermo difference on interior levels."""
    km_zm = np.full((gr.ngrdcol, gr.nzm), 2.0, dtype=np.float64)
    xm = np.tile(np.arange(gr.nzt, dtype=np.float64), (gr.ngrdcol, 1))

    out = clubb_api.calc_xpwp(gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, km_zm=km_zm, xm=xm)

    expected = np.zeros((gr.ngrdcol, gr.nzm), dtype=np.float64)
    expected[:, 1:-1] = km_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1] * (xm[:, 1:] - xm[:, :-1])

    assert out.shape == (gr.ngrdcol, gr.nzm)
    np.testing.assert_allclose(out, expected)
