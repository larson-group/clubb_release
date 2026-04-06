"""Test the stats system Python API.

Tests init/finalize, config, variable metadata, update/budget
operations, and data push/pull.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.clubb_api import (
    init_err_info,
    init_clubb_params,
    init_stats, finalize_stats, get_stats_config,
    get_param_names,
    get_stats_var_meta, get_stats_var_data, set_stats_var_data,
    stats_begin_timestep, stats_end_timestep,
    stats_update_scalar, stats_update_1d, stats_update_2d,
    var_on_stats_list,
)
from clubb_python.derived_types.err_info import ErrInfo

# Test grid dimensions
NCOL = 2
NZT = 5
NZM = 6  # nzm = nzt + 1

# Test timing
DT = 60.0        # 60s timestep
TSAMP = 60.0     # sample every timestep
TOUT = 120.0     # output every 2 timesteps

REGISTRY = str(Path(__file__).resolve().parent / "test_stats_registry.in")
PARAM_NAMES = get_param_names()


def _decode(val):
    """Decode bytes to str if needed, then strip."""
    if isinstance(val, bytes):
        return val.decode('ascii', errors='replace').strip()
    return str(val).strip()


@pytest.fixture
def stats_env():
    """Set up stats with the test registry, tear down afterward."""
    init_err_info(NCOL)
    err_info = ErrInfo(ngrdcol=NCOL)
    zt = np.linspace(50.0, 250.0, NZT)
    zm = np.linspace(0.0, 250.0, NZM)
    clubb_params = init_clubb_params(NCOL, iunit=10, filename="")
    init_stats(
        registry_path=REGISTRY,
        output_path="",         # no NetCDF output
        ncol=NCOL,
        stats_tsamp=TSAMP,
        stats_tout=TOUT,
        dt_main=DT,
        day_in=1, month_in=1, year_in=2000,
        time_initial=0.0,
        nzt=NZT, zt=zt, nzm=NZM, zm=zm, err_info=err_info,
        sclr_dim=0, edsclr_dim=0,
        clubb_params=clubb_params,
        param_names=PARAM_NAMES,
    )
    yield err_info
    finalize_stats(err_info=err_info)


class TestStatsConfig:
    """Test stats configuration."""

    def test_disabled_by_default(self):
        """Before init, stats should be disabled."""
        init_err_info(1)
        cfg = get_stats_config()
        assert cfg[0] == 0  # enabled

    def test_init_enables_stats(self, stats_env):
        """After init with registry, stats should be enabled."""
        cfg = get_stats_config()
        assert cfg[0] == 1  # enabled

    def test_config_fields(self, stats_env):
        """Config should have correct ncol and nvars."""
        cfg = get_stats_config()
        assert cfg[1] == NCOL   # ncol
        assert cfg[2] == 5      # nvars (5 entries in test registry)

    def test_sampling_params(self, stats_env):
        """Sampling parameters should match init args."""
        cfg = get_stats_config()
        # stats_nsamp = stats_tsamp / dt_main = 60/60 = 1
        assert cfg[3] == 1  # stats_nsamp
        # stats_nout = stats_tout / dt_main = 120/60 = 2
        assert cfg[4] == 2  # stats_nout

    def test_finalize_disables(self, stats_env):
        """After finalize, stats should be disabled."""
        finalize_stats(err_info=stats_env)
        cfg = get_stats_config()
        assert cfg[0] == 0  # enabled
        # Re-init for the fixture teardown to not double-finalize
        zt = np.linspace(50.0, 250.0, NZT)
        zm = np.linspace(0.0, 250.0, NZM)
        clubb_params = init_clubb_params(NCOL, iunit=10, filename="")
        init_stats(
            registry_path=REGISTRY, output_path="", ncol=NCOL,
            stats_tsamp=TSAMP, stats_tout=TOUT, dt_main=DT,
            day_in=1, month_in=1, year_in=2000, time_initial=0.0,
            nzt=NZT, zt=zt, nzm=NZM, zm=zm, sclr_dim=0, edsclr_dim=0, err_info=stats_env,
            clubb_params=clubb_params, param_names=PARAM_NAMES,
        )


class TestVarMeta:
    """Test stats variable metadata."""

    def test_first_var_is_thlm(self, stats_env):
        """First registered variable should be 'thlm' on zt grid."""
        name, grid, units, long_name, grid_id, nz = get_stats_var_meta(0)
        assert _decode(name) == "thlm"
        assert _decode(grid) == "zt"
        assert _decode(units) == "K"
        assert int(nz) == NZT

    def test_zm_var(self, stats_env):
        """wp2 should be on zm grid with nz=NZM."""
        name, grid, units, long_name, grid_id, nz = get_stats_var_meta(1)
        assert _decode(name) == "wp2"
        assert _decode(grid) == "zm"
        assert int(nz) == NZM

    def test_sfc_var(self, stats_env):
        """lwp should be on sfc grid with nz=1."""
        name, grid, units, long_name, grid_id, nz = get_stats_var_meta(2)
        assert _decode(name) == "lwp"
        assert _decode(grid) == "sfc"
        assert int(nz) == 1

    def test_var_long_name(self, stats_env):
        """Long names should be preserved."""
        name, grid, units, long_name, grid_id, nz = get_stats_var_meta(0)
        assert "Liquid water potential temperature" in _decode(long_name)


class TestVarOnStatsList:
    """Test variable name lookup."""

    def test_registered_var(self, stats_env):
        """Registered variables should be found."""
        assert var_on_stats_list("thlm") is True
        assert var_on_stats_list("wp2") is True
        assert var_on_stats_list("lwp") is True

    def test_unregistered_var(self, stats_env):
        """Unregistered variables should not be found."""
        assert var_on_stats_list("nonexistent_var") is False


class TestDataRoundtrip:
    """Test push/pull of variable buffer data."""

    def test_initial_buffer_zeros(self, stats_env):
        """Buffers should be zero after init."""
        buf, nsamp = get_stats_var_data(0, NCOL, NZT)
        assert buf.shape == (NCOL, NZT)
        assert nsamp.shape == (NCOL, NZT)
        np.testing.assert_array_equal(buf, 0.0)
        np.testing.assert_array_equal(nsamp, 0)

    def test_push_pull_roundtrip(self, stats_env):
        """Set buffer data and read it back."""
        buf = np.arange(NCOL * NZT, dtype=np.float64).reshape(NCOL, NZT)
        nsamp = np.ones((NCOL, NZT), dtype=np.int32) * 3
        set_stats_var_data(0, NCOL, NZT, buf, nsamp)
        buf2, nsamp2 = get_stats_var_data(0, NCOL, NZT)
        np.testing.assert_allclose(buf2, buf)
        np.testing.assert_array_equal(nsamp2, nsamp)


class TestTimestepAndUpdate:
    """Test begin/end timestep and update operations."""

    def test_begin_timestep_sets_l_sample(self, stats_env):
        """After begin_timestep, l_sample should be True."""
        stats_begin_timestep(0)
        cfg = get_stats_config()
        assert cfg[7] == 1  # l_sample

    def test_update_scalar(self, stats_env):
        """Update a scalar value and verify buffer accumulation."""
        stats_begin_timestep(0)
        stats_update_scalar("thlm", 300.0, icol=0, level=0)
        stats_end_timestep(DT, err_info=stats_env)

        buf, nsamp = get_stats_var_data(0, NCOL, NZT)
        assert nsamp[0, 0] == 1
        np.testing.assert_allclose(buf[0, 0], 300.0)

    def test_update_1d(self, stats_env):
        """Update a 1D profile and verify buffer."""
        stats_begin_timestep(0)
        profile = np.arange(1.0, NZT + 1.0)
        stats_update_1d("thlm", NZT, profile, icol=0)
        stats_end_timestep(DT, err_info=stats_env)

        buf, nsamp = get_stats_var_data(0, NCOL, NZT)
        np.testing.assert_allclose(buf[0, :], profile)

    def test_update_2d(self, stats_env):
        """Update a 2D field and verify buffer."""
        stats_begin_timestep(0)
        field = np.ones((NCOL, NZM)) * 42.0
        stats_update_2d("wp2", NCOL, NZM, field)
        stats_end_timestep(DT, err_info=stats_env)

        buf, nsamp = get_stats_var_data(1, NCOL, NZM)
        np.testing.assert_allclose(buf, 42.0)


class TestNetcdfMetadata:
    """Test metadata written to NetCDF output."""

    def test_init_writes_param_names_and_values(self, tmp_path):
        """Stats init should write CLUBB parameter names and values to output."""
        init_err_info(NCOL)
        err_info = ErrInfo(ngrdcol=NCOL)
        zt = np.linspace(50.0, 250.0, NZT)
        zm = np.linspace(0.0, 250.0, NZM)
        output_path = tmp_path / "stats_with_params.nc"
        clubb_params = init_clubb_params(NCOL, iunit=10, filename="")
        init_stats(
            registry_path=REGISTRY,
            output_path=str(output_path),
            ncol=NCOL,
            stats_tsamp=TSAMP,
            stats_tout=TOUT,
            dt_main=DT,
            day_in=1,
            month_in=1,
            year_in=2000,
            time_initial=0.0,
            nzt=NZT,
            zt=zt,
            nzm=NZM,
            zm=zm,
            sclr_dim=0,
            edsclr_dim=0,
            clubb_params=clubb_params,
            param_names=PARAM_NAMES,
            err_info=err_info,
        )
        finalize_stats(err_info=err_info)

        with Dataset(output_path) as ds:
            assert "param_name" in ds.variables
            assert "clubb_params" in ds.variables

            raw_names = np.ma.filled(ds.variables["param_name"][:], b" ")
            decoded_names = [b"".join(row).decode("ascii").strip() for row in raw_names]
            assert decoded_names == PARAM_NAMES

            stored_params = ds.variables["clubb_params"][:]
            np.testing.assert_allclose(stored_params, clubb_params.T)
