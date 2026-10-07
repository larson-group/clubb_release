"""Check source-derived loss metrics, benchmark reads and input boundaries.

src/clubb_loss_driver.F90 defines the metric and request contracts. Analytic
profiles and synthetic NetCDF records exercise metric edge cases, units,
interpolation, dimensions and configuration errors without advancing a case.
Real reset/rerun checks adapted from src/clubb_loss_driver_test.F90 and the
native executable oracle live in clubb_jax/tests/test_loss_driver.py.
"""
import importlib.abc
from pathlib import Path
import subprocess
import sys

import jax
import jax.numpy as jnp
import numpy as np
from netCDF4 import Dataset
import pytest

from clubb_jax.src import clubb_loss_driver as loss
from clubb_jax.src.CLUBB_core.advance_clubb_core_module import advance_clubb_core
from clubb_jax.src.CLUBB_core.parameters_tunable import PNAME_IDX, init_clubb_params
from clubb_jax.src.CLUBB_core.stats_netcdf import StatsWriter
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from utilities.create_case_namelist import (
    create_case_namelist_file, prune_clubb_stats_namelist, set_stats_string, set_stats_value,
)

ROOT = Path(loss.__file__).resolve().parents[2]


@pytest.mark.parametrize('model,truth,expected', [
    ([1., 2., 3.], [1., 2., 3.], [1., 1., 0., 0.]),
    ([3., 2., 1.], [1., 2., 3.], [-1., 1., 2., 0.]),
    ([2., 2., 2.], [1., 2., 3.], [0., 0., 1., 0.]),
    ([2., 4., 6.], [1., 2., 3.], [1., 2., 1., np.sqrt(6.)]),
    ([1., 2., 3.], [2., 2., 2.], [1., 1., np.sqrt(2./3.), 0.]),
    ([4.], [2.], [1., 1., 0., 2.]),
])
def test_taylor_metrics_source_conventions(model, truth, expected):
    result = loss.calculate_taylor_metrics(jnp.array(model), jnp.array(truth))
    np.testing.assert_allclose(result, expected, rtol=1.e-13, atol=1.e-13)


@pytest.mark.parametrize('model,truth,message', [
    ([1.0], [1.0, 2.0], 'matching sizes'),
    ([], [], 'at least one level'),
])
def test_taylor_metrics_invalid_shapes(model, truth, message):
    with pytest.raises(ValueError, match=message):
        loss.calculate_taylor_metrics(jnp.array(model), jnp.array(truth))


@pytest.mark.parametrize('dimensions', [
    ('time', 'z', 'col'), ('col', 'time', 'z'),
])
@pytest.mark.parametrize('units,divisor', [
    ('K', 1.0), ('g/kg', 1000.0), ('K/day', 86400.0),
])
def test_benchmark_raw_values_interpolation_and_first_column(
    tmp_path, dimensions, units, divisor,
):
    registry = tmp_path / 'stats.in'
    registry.write_text(
        '&clubb_stats_nl\n'
        ' entry(1)="profile | zt | K | Compared profile"\n/\n'
    )
    writer = StatsWriter(
        registry_path=str(registry), output_path='', nzt=3, nzm=4, ngrdcol=2,
        zt=np.array([0.0, 5.0, 20.0]), zm=np.array([0.0, 5.0, 10.0, 20.0]),
        stats_tsamp=60.0, stats_tout=120.0, dt_main=60.0,
        day=1, month=1, year=2000, time_initial=0.0,
    )
    benchmark = tmp_path / 'benchmark.nc'
    data = np.array([
        [[99.0, 999.0], [99.0, 999.0], [99.0, 999.0]],
        [[0.0, 999.0], [4.0, 999.0], [8.0, 999.0]],
        [[4.0, 999.0], [8.0, 999.0], [12.0, 999.0]],
    ])
    with Dataset(benchmark, 'w') as dataset:
        for name, values in (
            ('time', [0.0, 60.0, 120.0]), ('z', [0.0, 10.0, 20.0]),
            ('col', [1.0, 2.0]),
        ):
            dataset.createDimension(name, len(values))
            dataset.createVariable(name, 'f8', (name,))[:] = values
        dataset['time'].units = 'seconds since 2000-01-01 00:00:00.0'
        dataset['z'].units = 'm'
        variable = dataset.createVariable('truth', 'f8', dimensions, fill_value=0.0)
        variable.units = units
        variable.scale_factor = 2.0
        variable.add_offset = 100.0
        variable.set_auto_scale(False)
        axes = [('time', 'z', 'col').index(name) for name in dimensions]
        variable[:] = data.transpose(axes)

    request = loss.loss_request_type(
        state={'ngrdcol': 2}, les_stats_file=str(benchmark),
        altitude_comparison_range=(0.0, 20.0), time_window_ranges=((0, 120),),
        fields=[loss.loss_field_type('profile', 'truth')],
    )
    try:
        loss.prepare_loss_request_for_scoring(request, writer)
        # The t=0 record and second column are excluded. Raw zero is valid;
        # scale_factor/add_offset are ignored, as with nf90_get_var.
        np.testing.assert_array_equal(
            request.fields[0].truth_profile, np.array([[2.0, 4.0, 10.0]]) / divisor,
        )
    finally:
        writer.finalize()


def test_loss_file_batches_keep_configured_column_slices(tmp_path, monkeypatch):
    parameter_file = tmp_path / 'params.in'
    parameter_file.write_text('&clubb_params_nl\n/\n')
    params = jnp.asarray(init_clubb_params(4, str(parameter_file)))
    state = {
        'ngrdcol': 2, 'stats_output_path': 'stats.nc', 'clubb_params': params[:2],
        '_initial_state': {'clubb_params': params[:2]},
        'err_info': ErrInfo.initialized(2), '_jax_stats': None,
    }
    request = loss.loss_request_type(
        l_initialized=True, total_param_sets=4, state=state,
        time_window_ranges=((0, 120),), dt_main_seconds=60.0,
    )
    monkeypatch.setattr(loss, 'active_request', request)
    batches = []

    def reset(state, params, batch_num=None):
        batches.append(batch_num)
        state['clubb_params'] = params

    monkeypatch.setattr(loss, 'set_case_initial_conditions', reset)
    monkeypatch.setattr(loss, 'advance_clubb_to_end', lambda *args, **kwargs: None)
    monkeypatch.setattr(
        loss, 'calculate_field_loss',
        lambda *args: tuple(jnp.zeros((1, 2)) for _ in range(5)),
    )
    result = loss.clubb_get_loss_for_params(params)
    assert batches == [1, 2]
    assert all(values.shape == (1, 1, 4) for values in result)
    with pytest.raises(ValueError, match='in-memory statistics'):
        loss.clubb_get_loss_for_params(params[:3])


@pytest.mark.parametrize('indexed_parser', [False, True])
def test_loss_namelist_preserves_first_blank_array_entry(monkeypatch, indexed_parser):
    class Namelist(dict):
        start_index = {'clubb_var_names': [2], 'benchmark_var_name': [2]}

    values = {
        'les_stats_file': 'truth.nc', 'altitude_comparison_range': [0.0, 100.0],
        'time_average_range': [0, 120],
    }
    if indexed_parser:
        values.update({
            'clubb_var_names(2)': 'thlm', 'benchmark_var_name(2)': 'thlm',
        })
    else:
        values = Namelist(values)
        values.update(clubb_var_names=['thlm'], benchmark_var_name=['thlm'])
    monkeypatch.setattr(
        loss, '_read_namelist_groups', lambda path: {'tuner_loss_nl': values}
    )
    with pytest.raises(ValueError, match='entries after the first blank'):
        loss.init_loss_request('unused.in', 1, loss.loss_request_type())


@pytest.mark.parametrize('option', ['-h', '-help', '--help'])
def test_loss_entry_point_accepts_help_without_initializing_model(monkeypatch, option):
    from clubb_jax.src import clubb_standalone_loss

    def unexpected_run(*args):
        pytest.fail('A help request must not initialize CLUBB')

    monkeypatch.setattr(clubb_standalone_loss, 'clubb_get_loss', unexpected_run)
    monkeypatch.setattr(sys, 'argv', ['loss_driver', option])
    assert clubb_standalone_loss.main() == 0


def test_loss_rejects_restart_before_initialization(tmp_path, monkeypatch):
    # Standalone restart support must not make the unsupported loss/restart
    # combination open output files or read a reference before its own gate.
    runfile = tmp_path / 'restart_loss.in'
    runfile.write_text('&model_setting\n l_restart=.true.\n/\n')
    monkeypatch.setattr(loss, 'active_request', loss.loss_request_type())

    def unexpected_initialization(*args):
        pytest.fail('Restart loss must be rejected before case initialization')

    monkeypatch.setattr(loss, 'init_clubb_case', unexpected_initialization)
    with pytest.raises(ValueError, match='does not support l_restart'):
        loss.init_clubb_loss(str(runfile))
