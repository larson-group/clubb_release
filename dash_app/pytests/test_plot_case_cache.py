"""Compact case descriptors, process handoff, and refresh-safe user intent."""

import json
import multiprocessing
from types import SimpleNamespace

import diskcache
import pytest
from dash import Dash, dcc, html, no_update

from dash_app.persistence import enable_workspace_persistence
from dash_app.plot_tab import callbacks_case, case_cache
from dash_app.plot_tab.tasks import _HANDLERS
from dash_app.plot_tab.plot_types import shared


@pytest.fixture
def metadata_cache(tmp_path, monkeypatch):
    cache = diskcache.Cache(str(tmp_path / 'metadata'))
    monkeypatch.setattr(case_cache, '_cache', lambda: cache)
    case_cache._load_metadata.cache_clear()
    yield cache
    case_cache._load_metadata.cache_clear()
    cache.close()


def metadata(name='bomex'):
    return {
        'name': name, 'files': ['/tmp/bomex_stats.nc'], 'output_dirs': ['/tmp'],
        'time_seconds': list(range(10000)),
        'time_bounds_seconds': [[i, i + 1] for i in range(10000)],
        'profile_vars': [{'label': f'Field {i}', 'value': f'field_{i}'} for i in range(1000)],
        'time_len': 10000, 'columns_len': 3, 'preserve_plot_view': True,
        'available_cases': ['bomex', 'rico'],
    }


def test_compact_snapshot_round_trip_and_persistence_policy(metadata_cache):
    full = metadata()
    descriptor = case_cache.compact_case_data(full)
    assert len(json.dumps(descriptor)) < 1000
    assert not (case_cache._DERIVED_FIELDS & descriptor.keys())
    assert case_cache.resolve_case_data(descriptor) == full
    assert case_cache.resolve_case_data(full) is full
    root = html.Div([dcc.Store(id='plots-case-data', data=descriptor),
                     dcc.Store(id='plots-case-selection', data={'name': 'bomex'})])
    enable_workspace_persistence(root)
    assert getattr(root.children[0], 'storage_type', 'memory') == 'memory'
    assert root.children[1].storage_type == 'local'


def _publish_from_worker(directory, data, queue):
    cache = diskcache.Cache(directory)
    case_cache._cache = lambda: cache
    queue.put(case_cache.compact_case_data(data))
    cache.close()


def test_background_worker_snapshot_resolves_in_main_process(metadata_cache):
    context = multiprocessing.get_context('spawn')
    queue = context.Queue()
    worker = context.Process(target=_publish_from_worker, args=(metadata_cache.directory, metadata(), queue))
    worker.start()
    try:
        descriptor = queue.get(timeout=15)
        worker.join(timeout=15)
        assert worker.exitcode == 0
        assert case_cache.resolve_case_data(descriptor) == metadata()
    finally:
        if worker.is_alive():
            worker.terminate()
            worker.join(timeout=5)
        queue.close()


def test_evicted_snapshot_rebuilds_and_other_case_remains_independent(metadata_cache, monkeypatch):
    from dash_app.services import profiles

    first = metadata()
    second = {**metadata('rico'), 'time_seconds': [11, 12]}
    one = case_cache.compact_case_data(first)
    two = case_cache.compact_case_data(second)
    metadata_cache.delete(one['metadata_key'])
    case_cache._load_metadata.cache_clear()
    rebuilt = []
    monkeypatch.setattr(profiles, 'build_case_metadata', lambda *args: rebuilt.append(args) or first)
    assert case_cache.resolve_case_data(one) == first
    assert case_cache.resolve_case_data(two) == second
    assert len(rebuilt) == 1


def test_hydration_keeps_saved_view_and_scans_once(metadata_cache, monkeypatch):
    app = Dash(__name__)
    callbacks_case.register_case_callbacks(app)
    select = _HANDLERS['case']
    calls = []
    monkeypatch.setattr(callbacks_case, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    monkeypatch.setattr(callbacks_case, 'scan_output_cases', lambda dirs: calls.append(dirs) or {'bomex': ['/tmp/bomex_stats.nc'], 'rico': ['/tmp/rico_stats.nc']})
    full = {
        **metadata(), 'time_controls_physical': True,
        'time_slider_duration_min_minutes': 1, 'time_slider_duration_max_minutes': 120,
        'time_slider_duration_step_minutes': 1, 'time_slider_start_min_seconds': 0,
        'time_slider_start_max_seconds': 7200, 'time_slider_final_end_seconds': 10800,
        'height_slider_min': 0, 'height_slider_max': 3000, 'default_height_range': [0, 3000],
        'benchmarks': {'available_sources': []},
    }
    monkeypatch.setattr(callbacks_case, 'build_case_data', lambda name, *args: {**full, 'name': name})
    saved = {'name': 'rico', 'output_dirs': ['/tmp'], 'time_start_seconds': 1200,
             'average_minutes': 10, 'height_range': [100, 1500], 'selected_column': 2,
             'column_mode': 'all', 'time_override': None}
    result = select([], ['/tmp'], None, [9], {'9': {'plot_type': 'profile', 'var': 'wp2'}},
                    10, None, [], 0, 'single', 1, 0, [0, 1], saved, None)
    assert len(calls) == 1
    assert result[0]['name'] == 'rico'
    assert result[0]['preserve_plot_view']
    assert result[0]['available_cases'] == shared.ordered_case_names(['bomex', 'rico'])
    assert result[2] == [9] and result[4] == 10
    assert result[5:7] == (2, 'all')
    assert result[9] == 10 and result[14] == 1200
    assert result[19] == [100, 1500]
    assert 'profile_vars' not in result[0]

    remember = app.callback_map['plots-case-selection.data']['callback'].__wrapped__
    assert remember(None, 1, 0, [0, 1], 0, 'single', None, saved) is no_update
    assert remember(result[0], 10, 1200, [100, 1500], 2, 'all', None, saved) is no_update


def test_unchanged_catalog_returns_no_update(monkeypatch):
    app = Dash(__name__)
    callbacks_case.register_case_callbacks(app)
    refresh = _HANDLERS['catalog']
    catalog = [{'path': '/tmp/a'}]
    scans = []
    monkeypatch.setattr(callbacks_case, 'callback_context', SimpleNamespace(triggered_id='plots-output-menu-expanded'))
    monkeypatch.setattr(callbacks_case, 'discover_output_directories', lambda **kwargs: scans.append(kwargs) or catalog)
    assert refresh(False, 0, ['/tmp/a'], catalog) is no_update
    assert not scans
    assert refresh(True, 0, ['/tmp/a'], catalog) is no_update
    assert len(scans) == 1


def test_legacy_browser_workspace_migration_and_import():
    """Execute the real browser asset with a localStorage shim, without a server."""
    import shutil
    import subprocess
    from pathlib import Path

    node = shutil.which('node')
    if node is None:
        pytest.skip('Node is needed to exercise the browser workspace migration')
    script = r'''
const fs = require('fs');
const vm = require('vm');
const assert = require('assert').strict;
const source = fs.readFileSync(process.argv[1], 'utf8');
function workspace(entries) {
  const map = new Map(Object.entries(entries));
  const storage = {
    getItem: key => map.has(key) ? map.get(key) : null,
    setItem: (key, value) => map.set(key, value), removeItem: key => map.delete(key),
    key: index => Array.from(map.keys())[index], get length() { return map.size; }
  };
  const context = {document: {title: 'clubb_test'}, TextDecoder,
    window: {localStorage: storage, dash_clientside: {no_update: 'NO_UPDATE'},
      atob: value => Buffer.from(value, 'base64').toString('binary'),
      setTimeout: () => {}, location: {reload: () => {}}}};
  vm.runInNewContext(source, context);
  return {storage, actions: context.window.dash_clientside.dashboardWorkspace};
}
const token = 'clubb_test:dashboard-workspace:v1';
const key = id => '_dash_persistence.' + id + '.value.' + JSON.stringify(token);
const legacy = {'plots-case-data': JSON.stringify({name: 'arm', output_dirs: ['/tmp/run'], time_seconds: [1,2,3]}),
  [key('plots-global-time-range')]: JSON.stringify([60, 1]),
  [key('plots-global-time-point')]: JSON.stringify([7200, 0]),
  [key('plots-global-height-range')]: JSON.stringify([[100, 2000], [0, 1]]),
  [key('plots-column-mode')]: JSON.stringify(['all', 'single']),
  'plots-selected-column': '2'};
const migrated = workspace(legacy);
const selected = JSON.parse(migrated.storage.getItem('plots-case-selection'));
assert.deepEqual(selected, {name: 'arm', output_dirs: ['/tmp/run'], average_minutes: 60,
  time_start_seconds: 7200, height_range: [100,2000], column_mode: 'all', selected_column: 2});
assert.equal(migrated.storage.getItem('plots-case-data'), null);
const existing = workspace({...legacy, 'plots-case-selection': JSON.stringify({name:'rico'})});
assert.equal(JSON.parse(existing.storage.getItem('plots-case-selection')).name, 'rico');
const malformed = workspace({...legacy, [key('plots-global-time-range')]: 'invalid'});
assert.equal(JSON.parse(malformed.storage.getItem('plots-case-selection')).name, 'arm');
const meta = {schema_version: 1, schema_key: 'clubb_test:dashboard-workspace:schema', token,
  store_ids: ['plots-case-selection', 'plots-output-dirs', 'plots-selected-column'], extra_keys: []};
const payload = Buffer.from(JSON.stringify({schema_version:1, local_storage:legacy})).toString('base64');
const restored = workspace({});
assert.match(restored.actions.importWorkspace('data:application/json;base64,'+payload, 'old.json', meta), /^Restored/);
assert.deepEqual(JSON.parse(restored.storage.getItem('plots-case-selection')), selected);
assert.equal(restored.storage.getItem('plots-case-data'), null);
restored.actions.resetWorkspace(1, meta);
assert.equal(restored.storage.getItem('plots-case-selection'), null);
'''
    asset = Path(__file__).resolve().parents[1] / 'assets' / '37_workspace_persistence.js'
    result = subprocess.run([node, '-e', script, str(asset)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_initial_plot_layout_defers_discovery_and_netcdf(monkeypatch):
    from dash_app.plot_tab import state, tab
    from dash_app.services import profiles

    discoveries = []
    monkeypatch.setattr(profiles, 'discover_output_directories', lambda **kwargs: discoveries.append(kwargs) or [])
    def unexpected(*args, **kwargs):
        raise AssertionError('initial layout must not open NetCDF metadata')
    monkeypatch.setattr(state, 'build_case_data', unexpected)
    app = Dash(__name__)
    plot_tab = tab.build_tab(app)
    app.layout = html.Div([dcc.Tabs([plot_tab], id='dashboard-tabs', value='plots'),
                          dcc.Store(id='theme-store', data='dark'),
                          dcc.Store(id='dashboard-request'), dcc.Store(id='dashboard-broker-jobs')])
    enable_workspace_persistence(app.layout)
    client = app.server.test_client()
    response = client.get('/_dash-layout')
    assert response.status_code == 200
    assert not discoveries
    assert b'plots-output-refresh-interval' not in response.data
    assert client.get('/_dash-dependencies').status_code == 200


def test_identical_metadata_from_different_outputs_has_separate_identity(metadata_cache):
    first = metadata()
    second = {**first, 'files': ['/tmp/other/bomex_stats.nc']}
    assert case_cache.compact_case_data(first)['metadata_key'] != case_cache.compact_case_data(second)['metadata_key']


def test_external_output_is_inspected_without_recursing(tmp_path, monkeypatch):
    app = Dash(__name__)
    callbacks_case.register_case_callbacks(app)
    refresh = _HANDLERS['catalog']
    calls = []
    monkeypatch.setattr(callbacks_case, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    monkeypatch.setattr(callbacks_case, 'discover_output_directories',
                        lambda **kwargs: calls.append(kwargs) or [])
    result = refresh(True, 0, [str(tmp_path)], [])
    assert calls == [{'root': str(tmp_path), 'selected_dirs': [], 'recursive': False}]
    assert result[0]['path'] == str(tmp_path)
    assert not result[0]['available']
    assert refresh(True, 0, [str(tmp_path)], result) is no_update
    assert len(calls) == 1  # Known selections never trigger another scan.


def test_restored_external_output_gets_counts_from_worker_scan(monkeypatch):
    app = Dash(__name__)
    callbacks_case.register_case_callbacks(app)
    refresh = _HANDLERS['catalog']
    monkeypatch.setattr(callbacks_case, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    record = {'path': '/tmp/run', 'case_names': ['bomex'], 'case_count': 1, 'available': True}
    monkeypatch.setattr(callbacks_case, 'discover_output_directories', lambda **kw: [record])
    result = refresh(False, 0, ['/tmp/run'], [])
    assert result[0]['case_names'] == ['bomex']
    assert result[0]['available']
