"""Numerical and I/O contracts for bounded, selection-first plot extraction."""

import os

import numpy as np
import pytest
from netCDF4 import Dataset

from dash_app.plot_tab.plot_types import shared
from dash_app.shared.array_cache import ArrayCache, value_bytes


@pytest.fixture(autouse=True)
def empty_caches():
    shared.clear_all_caches()
    yield
    shared.clear_all_caches()


@pytest.fixture(params=[('time', 'zt', 'col'), ('col', 'zt', 'time'), ('zt', 'time', 'col')])
def stats(request, tmp_path):
    path = str(tmp_path / 'bomex_stats.nc')
    cube = np.random.default_rng(123).normal(size=(8, 4, 3))
    cube[0, 0, 0] = 0.0  # CLUBB sometimes uses zero as _FillValue; it is real data.
    subcolumns = np.stack([cube, cube * 2], axis=-1)
    with Dataset(path, 'w') as ds:
        for dim, size in [('time', 8), ('zt', 4), ('col', 3), ('subcolumn', 2)]:
            ds.createDimension(dim, size)
        ds.createVariable('time', 'f8', ('time',))[:] = np.arange(1, 9) * 60
        ds.variables['time'].units = 'seconds'
        ds.createVariable('zt', 'f8', ('zt',))[:] = [0, 100, 200, 300]
        ds.variables['zt'].units = 'm'
        dims = request.param
        field = ds.createVariable('field', 'f8', dims, fill_value=0.0)
        field[:] = np.transpose(cube, [('time', 'zt', 'col').index(d) for d in dims])
        field.units = 'K'
        field.long_name = 'A test field'
        # Keep subcolumns before columns to exercise flattened line ordering.
        ds.createVariable('lh_nl_w', 'f8', ('col', 'subcolumn', 'zt', 'time'))[:] = subcolumns.transpose(2, 3, 1, 0)
        ds.createVariable('lwp', 'f8', ('col', 'time'))[:] = cube[:, 0, :].T
    return path, cube, subcolumns, request.param


@pytest.mark.parametrize('window', [(2, 4), (4, 2), (-20, 100), (7, 7)])
def test_profile_matches_reference_time_mean(stats, window):
    path, cube, _, _ = stats
    first, last = sorted(window)
    first, last = max(0, min(first, 7)), max(0, min(last, 7))
    for column in range(3):
        result = shared.extract_time_avg_profile_for_path(path, 'field', window, col_index=column)
        np.testing.assert_allclose(result['profiles'], cube[first:last + 1, :, column].mean(axis=0)[None, :], rtol=1e-14, atol=1e-14)
        assert result['labels'] == (f'col {column + 1}',)
        assert result['units'] == 'K'
        assert result['z_units'] == 'm'
        assert not result['profiles'].flags.writeable


def test_filtered_columns_and_plot_families_preserve_values_and_labels(stats):
    path, cube, subcolumns, _ = stats
    selection = dict(column_mode='all', column_filter_indices=[2, 0, 2, 99])
    profile = shared.extract_time_avg_profile_for_path(path, 'field', [2, 4], **selection)
    np.testing.assert_allclose(profile['profiles'], cube[2:5, :, [0, 2]].mean(axis=0).T)
    assert profile['labels'] == ('col 1', 'col 3')
    timeheight = shared.extract_time_height_for_path(path, 'field', **selection)
    np.testing.assert_allclose(timeheight[2], cube[:, :, [0, 2]].mean(axis=2))
    timeseries = shared.extract_timeseries_for_path(path, 'lwp', **selection)
    np.testing.assert_array_equal(timeseries[1], cube[:, 0, [0, 2]])
    assert timeseries[2] == ('col 1', 'col 3')
    z, values, labels, genuine = shared.extract_subcolumn_profiles_for_path(path, 'w', [2, 4], **selection)
    np.testing.assert_allclose(values, subcolumns[2:5, :, [0, 2], :].mean(axis=0).reshape(4, 4))
    assert genuine
    assert labels == ('col 1 / subcol 1', 'col 1 / subcol 2', 'col 3 / subcol 1', 'col 3 / subcol 2')
    assert shared.extract_time_avg_profile_for_path(path, 'field', [2, 4], column_mode='all', column_filter_indices=[]) is None


def test_single_subcolumn_and_non_subcolumn_fallback(stats):
    path, cube, subcolumns, _ = stats
    z, values, labels, genuine = shared.extract_subcolumn_profiles_for_path(path, 'w', [2, 4], col_index=1)
    np.testing.assert_allclose(values, subcolumns[2:5, :, 1, :].mean(axis=0))
    assert genuine and labels == ('col 2 / subcol 1', 'col 2 / subcol 2')
    _, values, labels, genuine = shared.extract_subcolumn_profiles_for_path(path, 'field', [2, 4], col_index=1)
    np.testing.assert_allclose(values, cube[2:5, :, 1].mean(axis=0)[:, None])
    assert not genuine and labels == ('col 2',)


def test_read_slices_before_loading_and_reuses_reduced_result(stats, monkeypatch):
    path, cube, _, dims = stats
    factory = shared.DatasetInfo
    reads = []

    class Variable:
        def __init__(self, variable):
            self.variable = variable

        def __getitem__(self, key):
            reads.append(key)
            return self.variable[key]

    class DatasetProxy:
        def __init__(self, dataset):
            self.dataset = dataset
            self.variables = dict(dataset.variables)
            self.variables['field'] = Variable(dataset.variables['field'])
            self.dimensions = dataset.dimensions

        def close(self):
            self.dataset.close()

    def instrumented(path):
        info = factory(path)
        info.ds = DatasetProxy(info.ds)
        return info

    monkeypatch.setattr(shared, 'DatasetInfo', instrumented)
    first = shared.extract_time_avg_profile_for_path(path, 'field', [2, 4], col_index=1)
    again = shared.extract_time_avg_profile_for_path(path, 'field', [2, 4], col_index=1)
    assert first is again
    assert len(reads) == 1
    selectors = dict(zip(dims, reads[0]))
    assert selectors['time'] == slice(2, 5)
    assert selectors['col'] == slice(1, 2)
    assert selectors['zt'] == slice(None)
    # No time × height × column cube survives in the result cache.
    assert all(value.ndim < 3 for item in shared._EXTRACT_CACHE.values()
               for value in (item.values() if isinstance(item, dict) else [item])
               if isinstance(value, np.ndarray))


def test_file_change_invalidates_cached_results(stats):
    path, _, _, _ = stats
    before = shared.extract_time_avg_profile_for_path(path, 'field', [2, 4])
    old_signature = shared._file_signature(path)
    with Dataset(path, 'a') as ds:
        ds.variables['field'][:] = 42.0
    os.utime(path, ns=(old_signature[1] + 1000000, old_signature[1] + 1000000))
    after = shared.extract_time_avg_profile_for_path(path, 'field', [2, 4])
    assert after is not before
    np.testing.assert_array_equal(after['profiles'], np.full((1, 4), 42.0))
    assert all(key[1] != old_signature for key in shared._EXTRACT_CACHE)


def test_cache_enforces_byte_and_entry_limits_and_lru():
    array = np.arange(100, dtype=float)
    size = value_bytes(array)
    cache = ArrayCache(max_bytes=2 * size, max_entries=3)
    cache['a'] = array
    cache['b'] = array
    cache.move_to_end('a')
    cache['c'] = array
    assert list(cache) == ['a', 'c']
    assert cache.total_bytes == 2 * size
    cache['huge'] = np.ones(1000)
    assert list(cache) == ['a', 'c']
    cache['a'] = array[:1]
    assert cache.total_bytes == size + value_bytes(array[:1])
    del cache['c']
    cache.clear()
    assert not cache and cache.total_bytes == 0
    cache = ArrayCache(max_bytes=100000, max_entries=2)
    for key in range(4):
        cache[key] = array
    assert list(cache) == [2, 3]


def test_profile_callback_resolves_compact_metadata_and_patches_playback(stats, monkeypatch):
    from types import SimpleNamespace
    from dash import Dash, Patch
    from dash_app.plot_tab.case_cache import compact_case_data
    from dash_app.plot_tab.plot_types import profile_plot

    path, cube, _, _ = stats
    case = {
        'name': 'bomex', 'files': [path], 'output_dirs': [os.path.dirname(path)],
        'time_controls_physical': True, 'time_len': 8,
        'time_seconds': (np.arange(1, 9) * 60).tolist(),
        'time_slider_duration_min_minutes': 1, 'default_time_duration_minutes': 2,
        'default_time_start_seconds': 0, 'default_height_range': [0, 300],
        'profile_vars': [{'label': 'field', 'value': 'field'}],
    }
    descriptor = compact_case_data(case)
    app = Dash(__name__)
    plot = profile_plot.ProfilePlotType()
    plot.register_callbacks(app)
    from dash_app.plot_tab.tasks import _HANDLERS
    update = _HANDLERS['profile']
    extracts = []
    extractor = shared.extract_time_avg_profile_for_path
    def counted(*args, **kwargs):
        extracts.append(args)
        return extractor(*args, **kwargs)
    monkeypatch.setattr(shared, 'extract_time_avg_profile_for_path', counted)
    monkeypatch.setattr(profile_plot, 'callback_context', SimpleNamespace(triggered_id='plots-case-data'))
    figure, _, _ = update('field', descriptor, [], 2, 120, None, [0, 300], 1,
                          'single', {}, 'dark', 'normal', None, {'index': 7})
    np.testing.assert_allclose(figure.data[0].x, cube[2:4, :, 1].mean(axis=0))
    assert len(extracts) == 1  # Error display and figure reuse the same traces.
    monkeypatch.setattr(profile_plot, 'callback_context', SimpleNamespace(triggered_id='plots-global-time-point'))
    patch, _, _ = update('field', descriptor, [], 2, 240, None, [0, 300], 1,
                         'single', {}, 'dark', 'normal', {'xaxis.range': [-2, 2]}, {'index': 7})
    assert isinstance(patch, Patch)
    operations = patch.to_plotly_json()['operations']
    update_x = next(op for op in operations if op['location'] == ['data', 0, 'x'])
    np.testing.assert_allclose(update_x['params']['value'], cube[4:6, :, 1].mean(axis=0))
    assert not any(op['location'] == ['layout', 'xaxis', 'range'] for op in operations)
    assert len(extracts) == 2


def test_budget_and_subcolumn_patches_support_physical_averaging_windows(stats, monkeypatch):
    from dash import Patch
    from dash_app.plot_tab.plot_types import budget_plot, subcolumn_plot

    path, cube, subcolumns, _ = stats
    context = {
        'case_data': {'files': [path], 'time_controls_physical': True, 'time_len': 8,
                      'time_seconds': (np.arange(1, 9) * 60).tolist(),
                      'time_slider_duration_min_minutes': 1},
        'time_range': 2, 'time_point': 120, 'selected_column': 1,
        'height_range': [0, 300],
    }
    monkeypatch.setattr(budget_plot, 'BUDGET_GROUPS', {'field': {'terms': ['field'], 'label': 'Field'}})
    budget = budget_plot.BudgetPlotType().build_patch({'var': 'field'}, context)
    subcolumn, _ = subcolumn_plot.SubcolumnPlotType().build_patch({'var': 'w'}, context)
    assert isinstance(budget, Patch) and isinstance(subcolumn, Patch)
    np.testing.assert_allclose(budget.to_plotly_json()['operations'][0]['params']['value'], cube[2:4, :, 1].mean(axis=0))
    np.testing.assert_allclose(subcolumn.to_plotly_json()['operations'][0]['params']['value'], subcolumns[2:4, :, 1, 0].mean(axis=0))
