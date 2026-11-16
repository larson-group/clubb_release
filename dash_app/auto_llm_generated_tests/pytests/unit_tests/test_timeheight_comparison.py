"""Check stacked directory comparisons using small, independent NetCDF fields."""

import numpy as np
import pytest
from dash import Dash
from netCDF4 import Dataset

from dash_app.plot_tab import callbacks_grid
from dash_app.plot_tab.plot_types import shared, timeheight_plot


@pytest.fixture(autouse=True)
def clear_caches():
    """Keep synthetic file metadata independent across tests."""
    shared.clear_all_caches()
    yield
    shared.clear_all_caches()


def _write_stats(path, steps, heights, offset):
    """Write known time/height/column values and a file-specific variable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    values = offset + np.arange(steps * len(heights) * 2).reshape(steps, len(heights), 2)
    with Dataset(path, 'w') as ds:
        for name, count in [('time', steps), ('zt', len(heights)), ('col', 2)]:
            ds.createDimension(name, count)
        time = ds.createVariable('time', 'f8', ('time',))
        time[:] = np.arange(1, steps + 1) * 60
        time.units = 'seconds'
        z = ds.createVariable('zt', 'f8', ('zt',))
        z[:] = heights
        z.units = 'm'
        field = ds.createVariable('field', 'f8', ('time', 'zt', 'col'))
        field[:] = values
        field.units = 'K'
        ds.createVariable(f'only_{offset}', 'f8', ('time', 'zt'))[:] = values[:, :, 0]
        # A name present everywhere is insufficient if one file lacks height.
        dims = ('time', 'zt') if offset == 0 else ('time',)
        ds.createVariable('incompatible', 'f8', dims)[:] = 1
    return values


def _case(paths):
    """Use the production comparison catalog with one-minute physical bins."""
    case = shared.build_case_data('synthetic', [str(path) for path in paths])
    case.update(time_slider_start_min_seconds=0, time_slider_duration_min_minutes=1,
                default_time_duration_minutes=1, time_slider_final_end_seconds=360)
    return case


@pytest.mark.parametrize('include_les', [False, True])
@pytest.mark.parametrize('mode,filters', [('single', None), ('all', None), ('all', [0])])
def test_comparison_stacks_each_directory_and_retains_les(tmp_path, monkeypatch, mode, filters, include_les):
    """Different grids and time extents retain their values and share color limits."""
    first, second = (tmp_path / name / 'synthetic_stats.nc' for name in ('first', 'second'))
    first_values = _write_stats(first, 4, [0, 100], 0)
    second_values = _write_stats(second, 6, [0, 80, 200], 100)
    case = _case([first, second])
    options = {option['value'] for option in case['timeheight_vars']}
    assert case['compare_mode'] and 'field' in options
    assert not {'only_0', 'only_100', 'incompatible'} & options
    assert timeheight_plot.PLOT.supports_case_data(case)
    collection = shared.DatasetCollection([str(first), str(second)])
    try:
        uncached_options = timeheight_plot.PLOT.get_variable_options(collection, {'compare_mode': True})
        assert {option['value'] for option in uncached_options} == options
    finally:
        collection.close()
    les_values = np.full((6, 2), -5.)
    monkeypatch.setattr(timeheight_plot, 'extract_benchmark_timeheight_panels',
                        lambda *_args: [{'title': 'SAM', 'data': (
                            np.arange(1, 7) * 60, [0, 100], les_values, 'K', '', 'm')}] if include_les else [])
    figure = timeheight_plot.PLOT.build_figure({'var': 'field'}, {
        'case_data': case, 'column_mode': mode, 'column_filter_indices': filters,
        'selected_column': 1, 'enabled_benchmarks': ['sam'], 'time_range': 1,
    })
    assert len(figure.data) == 2 + int(include_les)
    assert [a.text for a in figure.layout.annotations] == [
        label if mode == 'all' else f'{label} · CLUBB col 2'
        for label in case['source_labels']
    ] + (['SAM'] if include_les else [])
    for trace, values, heights in zip(figure.data[:2], (first_values, second_values), ([0, 100], [0, 80, 200])):
        expected = values[:, :, 1] if mode == 'single' else values[:, :, filters or [0, 1]].mean(axis=2)
        actual = np.asarray(trace.z)
        assert actual.shape == (len(heights), 6)
        np.testing.assert_allclose(actual[:, :len(values)], expected.T)
        assert np.isnan(actual[:, len(values):]).all()
        np.testing.assert_array_equal(trace.y, heights)
        np.testing.assert_array_equal(trace.x, np.arange(1, 7) * 60)
        assert trace.coloraxis == 'coloraxis'
    if include_les:
        np.testing.assert_array_equal(figure.data[-1].z, les_values.T)
    assert figure.layout.coloraxis.cmin == min(np.nanmin(trace.z) for trace in figure.data)
    assert figure.layout.coloraxis.cmax == max(np.nanmax(trace.z) for trace in figure.data)
    assert figure.layout.xaxis.matches == f'x{len(figure.data)}'
    assert figure.layout.yaxis.matches == f'y{len(figure.data)}'


def test_single_output_still_uses_the_single_heatmap(tmp_path, monkeypatch):
    """One directory without LES retains the original plot and unit colorbar."""
    path = tmp_path / 'one' / 'synthetic_stats.nc'
    values = _write_stats(path, 6, [0, 100], 0)
    monkeypatch.setattr(timeheight_plot, 'extract_benchmark_timeheight_panels', lambda *_args: [])
    figure = timeheight_plot.PLOT.build_figure({'var': 'field'}, {
        'case_data': _case([path]), 'time_range': 1, 'selected_column': 0,
    })
    assert len(figure.data) == 1
    np.testing.assert_array_equal(figure.data[0].z, values[:, :, 0].T)
    assert figure.data[0].colorbar.title.text == 'K'


def test_comparison_card_grows_with_directories_and_returns_to_normal():
    """The native card shell must leave room for every stacked panel."""
    app = Dash(__name__)
    callbacks_grid.register_grid_callbacks(app)
    callback = next(entry['callback'].__wrapped__ for key, entry in app.callback_map.items()
                    if 'plots-size-toggle' in key and '.disabled' in key)
    state = {'0': {'plot_type': 'timeheight'}}
    card, shell, _text, _class, locked = callback('normal', [], {
        'compare_mode': True, 'files': ['a', 'b', 'c', 'd'],
    }, {'index': 0}, state)
    assert locked and 'size-large' in card
    assert int(shell['height'].removesuffix('px')) >= 4 * 260
    card, shell, _text, _class, locked = callback('normal', [], {
        'compare_mode': False, 'files': ['a'],
    }, {'index': 0}, state)
    assert not locked and 'size-normal' in card
    assert shell == shared.plot_graph_shell_style('normal')
