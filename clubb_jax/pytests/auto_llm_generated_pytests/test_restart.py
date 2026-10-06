"""Synthetic restart-reader and initialization boundary checks.

Exercise raw zeros, units, column/dimension order, date offsets, saved times,
interpolation, broadcasting and invalid input without advancing a model case.
The real restart workflow lives in tests/run_restart_test.py and compares
one selected final-timestep field, using the same contract as the native test.
"""

import jax.numpy as jnp
import numpy as np
import pytest
from netCDF4 import Dataset

from clubb_jax.src.Input_fields import input_fields
from clubb_jax.src.clubb_case_initalization import init_clubb_case
from utilities.create_case_namelist import create_case_namelist_file, set_stats_string


def _write_profile(path, data, units="K", heights=(0.0, 10.0), dims=("time", "zt", "col")):
    with Dataset(path, "w") as ds:
        ds.createDimension("time", 2)
        ds.createDimension("zt", len(heights))
        ds.createDimension("col", data.shape[-1])
        time = ds.createVariable("time", "f8", ("time",))
        time.units = "seconds since 2000-01-01 00:00:00.0"
        time[:] = (60.0, 120.0)
        ds.createVariable("zt", "f8", ("zt",))[:] = heights
        var = ds.createVariable("thlm", "f8", dims, fill_value=0.0)
        var.units = units
        # Test noncanonical dimension order as well as normal CLUBB layout.
        axes = [("time", "zt", "col").index(dim) for dim in dims]
        var[:] = data.transpose(axes)


@pytest.mark.parametrize("dims", [("time", "zt", "col"), ("col", "time", "zt")])
@pytest.mark.parametrize("units,divisor", [("K", 1.0), ("g/kg", 1000.0), ("K/day", 86400.0)])
def test_raw_values_units_and_distinct_columns(tmp_path, dims, units, divisor):
    path = tmp_path / "case_stats.nc"
    data = np.array([[[0.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 0.0]]])
    _write_profile(path, data, units=units, dims=dims)
    result, error = input_fields.get_clubb_variable_interpolated(
        True, path, "thlm", 2, 2,
        np.array([0.0, 10.0]), jnp.zeros((2, 2)),
    )
    assert not error
    np.testing.assert_array_equal(result, data[1].T / divisor)


def test_exact_record_selection_date_offset_and_rejection(tmp_path, monkeypatch):
    path = tmp_path / "case_stats.nc"
    _write_profile(path, np.ones((2, 2, 1)))
    for key, value in (("clubb_year", 2000), ("clubb_month", 1), ("clubb_day", 1)):
        monkeypatch.setattr(input_fields, key, value)
    assert input_fields.compute_timestep(path, True, 60.0) == 1
    assert input_fields.compute_timestep(path, True, 120.0) == 2
    with pytest.raises(ValueError, match="not a saved output time"):
        input_fields.compute_timestep(path, True, 90.0)
    monkeypatch.setattr(input_fields, "clubb_day", 2)
    assert input_fields.compute_timestep(path, True, 120.0 - 86400.0) == 2


def test_split_and_unified_filenames(tmp_path):
    prefix = tmp_path / "case"
    assert input_fields.set_filenames(prefix) == (tmp_path / "case_stats.nc",) * 3
    (tmp_path / "case_zt.nc").touch()
    assert input_fields.set_filenames(prefix) == tuple(
        tmp_path / f"case_{grid}.nc" for grid in ("zt", "zm", "sfc")
    )


def test_single_column_broadcast_interpolation_and_errors(tmp_path):
    path = tmp_path / "case_stats.nc"
    _write_profile(path, np.array([[[1.0], [3.0]], [[2.0], [4.0]]]))
    result, error = input_fields.get_clubb_variable_interpolated(
        True, path, "thlm", 3, 2,
        np.array([-5.0, 5.0, 10.0]), jnp.zeros((2, 3)),
    )
    assert not error
    np.testing.assert_array_equal(result, [[2.0, 3.0, 4.0]] * 2)
    for name, heights, ncol in (("missing", [0.0, 10.0], 1), ("thlm", [0.0, 20.0], 1)):
        original = jnp.zeros((ncol, 2))
        result, error = input_fields.get_clubb_variable_interpolated(
            True, path, name, 2, 1, np.array(heights), original,
        )
        assert error
        np.testing.assert_array_equal(result, original)
    result, error = input_fields.get_clubb_variable_interpolated(
        False, path, "missing", 2, 1, np.array([0.0, 10.0]), original,
    )
    assert not error


@pytest.mark.parametrize("time_restart,message", [
    (301.0, "not a multiple of dt_main"),
    (-60.0, "time_restart must lie"),
])
def test_initialization_rejects_invalid_restart_clock(tmp_path, time_restart, message):
    path = create_case_namelist_file(
        "bomex", tmp_path, max_iters=10, stats="none", debug="0",
        override=f"l_restart=.true.,time_restart={time_restart}",
    )
    with pytest.raises(ValueError, match=message):
        init_clubb_case(str(path))


@pytest.mark.parametrize('suffix', ['stats', 'zt', 'zm', 'sfc'])
@pytest.mark.parametrize('alias', ['direct', 'symlink', 'hardlink'])
def test_initialization_preserves_reference_when_output_path_collides(tmp_path, suffix, alias):
    reference = tmp_path / f"bomex_{suffix}.nc"
    _write_profile(reference, np.ones((2, 2, 1)))
    if suffix != 'stats':
        (tmp_path / 'bomex_zt.nc').touch(exist_ok=True)
    original = reference.read_bytes()
    output = reference
    if alias != 'direct':
        output = tmp_path / 'alias.nc'
        if alias == 'symlink':
            output.symlink_to(reference)
        else:
            output.hardlink_to(reference)
    path = create_case_namelist_file(
        "bomex", tmp_path, max_iters=10, debug="0",
        override=f"l_restart=.true.,time_restart=60.,restart_path_case='{tmp_path / 'bomex'}'",
    )
    path.write_text(set_stats_string(path.read_text(), 'stats_output_filename', output.name))
    with pytest.raises(ValueError, match="must use different paths"):
        init_clubb_case(str(path))
    assert reference.read_bytes() == original


def test_batched_restart_rejected_before_reading_reference(tmp_path):
    path = create_case_namelist_file(
        'bomex', tmp_path, multicol='4', batch_size=2, stats='none',
        debug='-1', max_iters=1, override='l_restart=.true.',
    )
    with pytest.raises(ValueError, match='restart does not yet support runtime batching'):
        init_clubb_case(str(path))
