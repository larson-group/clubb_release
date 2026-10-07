"""Check that a failed SILHS loss column preserves its healthy neighbor.

This extends the reusable loss scenarios in src/clubb_loss_driver_test.F90
with the error-boundary behavior of src/clubb_driver.F90
(advance_clubb_to_end) and src/clubb_loss_driver.F90
(clubb_get_loss_for_params and calculate_field_loss).

A four-step, 60-second RICO case uses KK microphysics, native JAX SILHS draws,
two distinct C8 columns and two loss windows against a synthetic benchmark.
A controlled fatal microphysics flag is injected into one column. Its metrics
must receive finite penalties while the healthy column exactly reproduces
its uninterrupted loss values. The per-column isolation assertion is a JAX
extension to the general native loss tests; no Fortran executable is needed.
"""

import jax
import jax.numpy as jnp
import numpy as np
from netCDF4 import Dataset

from clubb_jax.src import advance_clubb_to_end as driver
from clubb_jax.src import clubb_loss_driver as loss
from clubb_jax.src.CLUBB_core import error_code
from utilities.create_case_namelist import (
    create_case_namelist_file, prune_clubb_stats_namelist, set_stats_string, set_stats_value,
)


def test_real_silhs_loss_penalizes_only_failed_microphysics_column(tmp_path, monkeypatch):
    original_debug_level = error_code._debug_level
    benchmark = tmp_path / 'truth.nc'
    times = np.arange(0., 241., 60.)
    heights = np.arange(0., 6001., 20.)
    with Dataset(benchmark, 'w') as ds:
        for name, values in [('time', times), ('z', heights), ('y', [1.]), ('x', [1.])]:
            ds.createDimension(name, len(values))
            ds.createVariable(name, 'f8', (name,))[:] = values
        ds['time'].units = 'seconds since 2004-12-16 00:00:00'
        ds['z'].units = 'm'
        ds.createVariable('thlm', 'f8', ('time', 'z', 'y', 'x'))[:] = (
            300. + heights[None, :] * .004 + times[:, None] * .00001
        )[:, :, None, None]
        ds['thlm'].units = 'K'
    path = create_case_namelist_file(
        'rico_silhs', tmp_path, multicol='C8/0.2:0.8/2', debug='-1',
        max_iters=4, dt_main=60, dt_rad=60,
    )
    text = prune_clubb_stats_namelist(path.read_text(), ['thlm'])
    text = set_stats_string(text, 'stats_output_filename', '')
    for name, value in dict(stats_tstart=0, stats_tend=240, stats_tout=120, stats_tsamp=60).items():
        text = set_stats_value(text, name, str(value))
    text += f'''\n&tuner_loss_nl
 les_stats_file="{benchmark}"
 clubb_var_names(1)="thlm"
 benchmark_var_name(1)="thlm"
 altitude_comparison_range=20,2940
 time_average_range=0,240
 num_time_windows=2
/\n'''
    path.write_text(text)
    try:
        _, params = loss.init_clubb_loss(str(path), True)
        expected = loss.clubb_get_loss_for_params(params)
        assert all(np.isfinite(values).all() for values in expected)
        for name, error_index in [('pdf_hydromet_microphys_prep', 1),
                                  ('advance_microphys', 9)]:
            original = getattr(driver, name)

            def fail_column(*args, **kwargs):
                result = list(original(*args, **kwargs))
                result[error_index] = result[error_index].set_fatal(jnp.array([True, False]))
                return tuple(result)

            with monkeypatch.context() as patch:
                patch.setattr(driver, name, fail_column)
                result = loss.clubb_get_loss_for_params(params)
            for values, healthy, penalty in zip(result, expected, loss.set_invalid_field_metric_outputs()):
                np.testing.assert_array_equal(values[..., 0], penalty)
                np.testing.assert_array_equal(values[..., 1], healthy[..., 1])
    finally:
        loss.finalize_clubb_loss()
        error_code._debug_level = original_debug_level
        jax.clear_caches()
