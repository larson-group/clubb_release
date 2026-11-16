"""Check lazy LES acquisition through Dash's existing benchmark extraction owner.

The clubb_dash/-dash suite owns these checks. Small synthetic NetCDF files
exercise real conversion and plotting; network transfer is already covered by
input/scripts/confirm_les_data's generated tests. No model or external data is
required here.
"""

import os
import base64
from concurrent.futures import Future
import json
from pathlib import Path
import shutil
import subprocess

from netCDF4 import Dataset
import numpy as np
import pytest

from dash_app.plot_tab import benchmark_overlay as overlay, case_definitions
from dash_app.plot_tab import tasks
from dash_app.plot_tab.layout import benchmark_button
from dash_app.plot_tab.plot_types import shared, timeheight_plot
from input.scripts import confirm_les_data as acquisition

VALUES = np.arange(6, dtype=float).reshape(2, 3) + 290
RELATIVE = 'sam_benchmark_runs/synthetic/les.nc'


def _write_stats(path, sam=True):
    """Write a known two-record thermodynamic field on a three-level grid."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, 'w', format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('time', 2)
        ds.createDimension('z', 3)
        t = ds.createVariable('time', 'f8', ('time',))
        t[:] = [60, 120]
        t.units = 'seconds'
        z = ds.createVariable('z', 'f8', ('z',))
        z[:] = [10, 100, 200]
        z.units = 'm'
        field = ds.createVariable('THETAL' if sam else 'thlm', 'f8', ('time', 'z'))
        field[:] = VALUES
        field.units = 'K'
        if sam:
            ds.createVariable('CLD', 'f8', ('time', 'z'))[:] = 0


def _case(root):
    """Build metadata from canonical paths that need not exist locally."""
    return {
        'name': 'synthetic', 'time_len': 2, 'time_seconds': [60, 120],
        'time_bounds_seconds': [[0, 60], [60, 120]], 'time_controls_physical': True,
        'time_interval_semantics': '(start, end]', 'time_slider_start_min_seconds': 0,
        'time_slider_duration_min_minutes': 1,
        'benchmarks': overlay.build_case_benchmark_info({
            'sam_benchmark_file': {'sam_benchmark': str(root / RELATIVE)},
            'coamps_benchmark_file': {'sm': str(root / 'les_runs/other.nc')},
        }),
    }


@pytest.fixture
def roots(tmp_path, monkeypatch):
    """Isolate acquisition roots and metadata caches from the user's dataset."""
    data = tmp_path / 'input' / acquisition.COLLECTION
    source = tmp_path / 'shared'
    output = tmp_path / 'output'
    monkeypatch.setattr(acquisition, 'INPUT_ROOT', data)
    monkeypatch.setattr(acquisition, 'SHARED_ROOTS', (source,))
    monkeypatch.setattr(acquisition, 'urlopen', lambda *_a, **_k: pytest.fail('unexpected network'))
    monkeypatch.setattr(overlay.tempfile, 'gettempdir', lambda: str(tmp_path))
    overlay.clear_benchmark_caches()
    shared.clear_all_caches()
    yield data, source, output
    overlay.clear_benchmark_caches()
    shared.clear_all_caches()


def test_discovery_keeps_missing_sources_without_acquiring_data(tmp_path, monkeypatch, roots):
    """LES controls must remain selectable when no symlink or downloaded data exists."""
    monkeypatch.setattr(case_definitions, 'REPO_ROOT', str(tmp_path))
    legacy = tmp_path / 'pyplotgen' / acquisition.COLLECTION
    legacy.mkdir(parents=True)
    expected = str(roots[0]) + os.sep
    assert case_definitions._benchmark_env()['BENCHMARK_OUTPUT_ROOT'] == expected
    monkeypatch.setattr(overlay, 'confirm_les_data', lambda *_a: pytest.fail('discovery acquired data'))
    case = _case(roots[0])
    assert case['benchmarks']['available_sources'] == ['sam', 'coamps']
    assert overlay.resolve_enabled_sources(case, ['sam'], strict=True) == ['sam']
    assert overlay.toggle_enabled_source(case, [], 'coamps') == ['coamps']
    assert overlay.extract_benchmark_timeheight_panels(case, 'thlm', []) == []
    assert not roots[0].exists() and not roots[2].exists()


def test_shared_data_are_rendered_without_an_input_symlink(roots, tmp_path):
    """Profile, PDF/custom-profile extraction and time-height panels reuse shared data."""
    data, source, output = roots
    _write_stats(source / RELATIVE)
    case = _case(data)
    profile = overlay.extract_benchmark_profile(case, 'sam', 'thlm', 'range', 2, 0)
    assert profile['profile'].shape == (3,)
    np.testing.assert_array_equal(profile['profile'], VALUES.mean(axis=0))
    profiles = overlay.extract_normalized_benchmark_profiles(case, 'sam', ['thlm'], 'range', 2, 0)
    np.testing.assert_array_equal(profiles['thlm']['profile'], VALUES.mean(axis=0))
    assert not data.exists() and not output.exists()

    model = tmp_path / 'model_stats.nc'
    _write_stats(model, sam=False)
    case['files'] = [str(model)]
    figure = timeheight_plot.PLOT.build_figure({'var': 'thlm'}, {
        'case_data': case, 'time_range': 1, 'selected_column': 0, 'enabled_benchmarks': ['sam'],
    })
    assert len(figure.data) == 2
    for heatmap in figure.data:
        assert np.asarray(heatmap.z).shape == (3, 2)
        np.testing.assert_array_equal(heatmap.z, VALUES.T)


def test_missing_selected_source_is_acquired_once_into_input(roots, monkeypatch):
    """Enabled SAM downloads into input after replacing a broken link; COAMPS is untouched."""
    data, source, output = roots
    data.parent.mkdir(parents=True)
    data.symlink_to(source, target_is_directory=True)
    calls = []

    def download(relative, target, temporary, descriptor):
        """Stand in for HTTP transport while retaining the real find/lock/download workflow."""
        os.close(descriptor)
        calls.append(relative)
        _write_stats(temporary)
        os.replace(temporary, target)
        return target.resolve()

    monkeypatch.setattr(acquisition, '_download', download)
    case = _case(data)
    assert calls == []
    # Figure extraction cannot start network work or wait for acquisition.
    assert overlay.extract_benchmark_profile(case, 'sam', 'thlm', 'range', 2, 0) is None
    assert calls == []
    for _ in range(2):
        assert overlay.acquire_benchmark_source(case, 'sam')['state'] == 'ready'
        profile = overlay.extract_benchmark_profile(case, 'sam', 'thlm', 'range', 2, 0)
        np.testing.assert_array_equal(profile['profile'], VALUES.mean(axis=0))
    assert calls == [RELATIVE]
    assert data.is_dir() and not data.is_symlink()
    assert (data / RELATIVE).is_file()
    assert not (data / 'les_runs/other.nc').exists()
    assert not source.exists() and not output.exists()


def test_acquisition_failure_is_visible_and_can_be_retried(roots, monkeypatch):
    """An acquisition error reaches the button state and does not poison a later retry."""
    case = _case(roots[0])

    def fail(_path):
        """Model an unavailable published benchmark."""
        raise RuntimeError('LES download failed: requested file is not published')

    monkeypatch.setattr(overlay, 'confirm_les_data', fail)
    result = overlay.acquire_benchmark_source(case, 'sam')
    assert result['state'] == 'error' and 'LES download failed' in result['error']
    assert not overlay._SOURCE_CACHE
    actual = roots[1] / RELATIVE
    _write_stats(actual)
    monkeypatch.setattr(overlay, 'confirm_les_data', acquisition.confirm_les_data)
    assert overlay.acquire_benchmark_source(case, 'sam')['state'] == 'ready'
    profile = overlay.extract_benchmark_profile(case, 'sam', 'thlm', 'range', 2, 0)
    np.testing.assert_array_equal(profile['profile'], VALUES.mean(axis=0))


def test_readiness_cannot_enable_another_case_or_an_unselected_overlay(roots):
    """Render inputs accept only the current case's prepared, selected sources."""
    case = _case(roots[0])
    ready = {'name': case['name'], 'benchmarks': case['benchmarks'], 'ready': ['sam'], 'sources': []}
    assert overlay.ready_benchmark_sources(case, ready) == []
    assert overlay.ready_benchmark_sources(case, ready, include_unselected=True) == ['sam']
    ready['sources'] = ['sam']
    assert overlay.ready_benchmark_sources(case, ready) == ['sam']
    assert overlay.ready_benchmark_sources({**case, 'name': 'other'}, ready) == []
    assert overlay.ready_benchmark_sources({**case, 'benchmarks': {}}, ready) == []


def test_native_callbacks_add_les_after_preparation_without_blocking_model_plots(roots, tmp_path, monkeypatch):
    """Render real native callbacks before/after preparation, including an LES expression."""
    from dash_app.plot_tab import async_callbacks

    handlers = {}
    monkeypatch.setattr(tasks, '_HANDLERS', handlers)
    monkeypatch.setattr(async_callbacks, '_HANDLERS', handlers)
    model = tmp_path / 'model_stats.nc'
    _write_stats(model, sam=False)
    _write_stats(roots[1] / RELATIVE)
    case = {**_case(roots[0]), 'files': [str(model)], 'default_height_range': [0, 200]}

    def execute(task, args):
        """Use the same fixed callback dispatch and JSON boundary as the worker."""
        return json.loads(tasks.execute({'task': task, 'args': args}))

    args = ['thlm', case, None, 2, 0, None, [0, 200], 0, 'single', {}, 'dark', 'normal', None, {'index': 1}]
    assert len(execute('profile', args)[0]['data']) == 1
    custom_args = ['thlm - les.thlm', case, 2, 0, None, [0, 200], 0, 'single', {}, None, 'dark', 'normal', None]
    pending = execute('custom', custom_args)[0]
    assert any('Preparing LES' in item['text'] for item in pending['layout']['annotations'])
    status = execute('benchmark_sam', [case, ['sam'], {}])
    assert status['state'] == 'ready'
    packet = {'name': case['name'], 'benchmarks': case['benchmarks'], 'sources': ['sam'], 'ready': ['sam']}
    args[2] = packet
    assert len(execute('profile', args)[0]['data']) == 2
    custom_args[9] = packet
    figure = execute('custom', custom_args)[0]
    assert len(figure['data']) == 1
    encoded = figure['data'][0]['x']
    values = np.frombuffer(base64.b64decode(encoded['bdata']), dtype=encoded['dtype'])
    assert values.shape == (3,)
    np.testing.assert_array_equal(values, [0, 0, 0])


def test_les_worker_does_not_hold_up_figures_or_case_selection():
    """A blocked acquisition leaves the existing figure and metadata queues free."""
    pools = []

    class Executor:
        """Retain unresolved futures to model a deliberately slow transfer."""
        def __init__(self):
            """Start with no submitted work."""
            self.calls = []

        def submit(self, function, request):
            """Keep one work item pending until the test explicitly completes it."""
            future = Future()
            self.calls.append((request, future))
            return future

        def shutdown(self, **_kwargs):
            """No processes are used by this isolated queue check."""

    def factory():
        """Give each queue lane its own observable executor."""
        executor = Executor()
        pools.append(executor)
        return executor

    queue = tasks.PlotTasks(factory)
    try:
        for task in ('benchmark_sam', 'profile', 'case'):
            queue.submit(dict(task=task, scope='browser:' + task, revision=1, args=[]))
        assert len(pools) == 3
        for pool, expected in zip(pools[1:], ('profile', 'case')):
            assert pool.calls[0][0]['task'] == expected
            pool.calls[0][1].set_result(json.dumps({'usable': True}))
            assert queue.poll('browser:' + expected, 1)['state'] == 'ready'
        assert queue.poll('browser:benchmark_sam', 1)['state'] == 'running'
        pools[0].calls[0][1].set_result(json.dumps({'state': 'ready'}))
        assert queue.poll('browser:benchmark_sam', 1)['result'] == {'state': 'ready'}
    finally:
        queue.close()


def test_les_button_exposes_busy_state_and_retryable_error():
    """Loading keeps the source toggle usable; errors are visible and accessible."""
    button = benchmark_button('sam', 'SAM LES', True, selected=True, acquisition={'state': 'loading'})
    props = button.to_plotly_json()['props']
    assert props['aria-busy'] == 'true' and props['disabled'] is False
    assert any(getattr(child, 'className', None) == 'plots-les-spinner' for child in button.children)
    failure = benchmark_button('sam', 'SAM LES', True, selected=True,
                               acquisition={'state': 'error', 'error': 'offline'})
    assert failure.children[1].role == 'alert'
    assert 'offline' in failure.children[1].children


def test_http_completion_enables_overlays_independent_of_json_field_order(roots, monkeypatch):
    """Use the real HTTP response to verify completed preparation reaches browser readiness."""
    from dash import Dash, html
    from dash_app.plot_tab.async_callbacks import register_task_routes
    from dash_app.shared import broker_client

    node = shutil.which('node')
    if not node:
        pytest.skip('Node is needed for browser callback tests')
    _write_stats(roots[1] / RELATIVE)
    case = _case(roots[0])
    # Deliberately reverse file-label order, as a multi-file source can do.
    path = str(roots[0] / RELATIVE)
    case['benchmarks']['sam'] = {'z_label': path, 'a_label': path}
    status = overlay.acquire_benchmark_source(case, 'sam')
    assert status['state'] == 'ready'
    monkeypatch.setattr(broker_client, 'perform_action', lambda *_a, **_k: {
        'results': [{'state': 'ready', 'result': status}],
    })
    app = Dash(__name__)
    app.layout = html.Div()
    register_task_routes(app)
    response = app.server.test_client().post('/plots/tasks', json={'operation': 'poll', 'requests': []})
    assert response.status_code == 200
    asset = Path(__file__).resolve().parents[3] / 'assets' / '33_plot_tasks.js'
    script = r'''
const assert = require('assert'), fs = require('fs');
global.window = global;
window.dash_clientside = {no_update: null};
eval(fs.readFileSync(process.argv[1], 'utf8'));
const input = JSON.parse(fs.readFileSync(0, 'utf8'));
const status = JSON.parse(input.httpBody).results[0].result;
let readiness = clubbPlotTasks.benchmarkReadiness(input.caseData, ['sam'], status, null, null);
assert.deepEqual(readiness.sources, ['sam'], 'HTTP completion must enable the prepared overlay');
// Field-order independence must not weaken the stale-case and path checks.
input.caseData = {...input.caseData, name: 'another_case'};
assert.deepEqual(clubbPlotTasks.benchmarkReadiness(input.caseData, ['sam'], status, null, null).sources, []);
input.caseData.name = status.request.name;
input.caseData.benchmarks.sam.a_label = '/another/file.nc';
assert.deepEqual(clubbPlotTasks.benchmarkReadiness(input.caseData, ['sam'], status, null, null).sources, []);
'''
    result = subprocess.run([node, '-e', script, str(asset)],
                            input=json.dumps({'caseData': case, 'httpBody': response.get_data(as_text=True)}),
                            text=True, capture_output=True, timeout=5)
    assert result.returncode == 0, result.stdout + result.stderr


def test_browser_keeps_model_figures_usable_until_background_les_is_ready():
    """Exercise real task delivery for pending, ready, disabled, failed and stale sources."""
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is needed for browser callback tests')
    asset = Path(__file__).resolve().parents[3] / 'assets' / '33_plot_tasks.js'
    script = r'''
const assert = require('assert');
const fs = require('fs');
global.window = global;
const NO = {no_update: true}, updates = [], submissions = [], ready = new Map(), states = {};
window.dash_clientside = {no_update: NO, set_props(id, props) {
    updates.push([id, props]);
    if (typeof id === 'string' && id.startsWith('plots-benchmark-')) states[id.slice(16)] = props.data;
}};
global.fetch = async (_url, options) => {
    const payload = JSON.parse(options.body);
    if (payload.operation === 'submit') {
        submissions.push(payload);
        return {ok: true, json: async () => ({state: 'queued'})};
    }
    return {ok: true, json: async () => ({results: payload.requests.map(
        item => ready.get(item.scope + ':' + item.revision) || {state: 'running'})})};
};
eval(fs.readFileSync(process.argv[1], 'utf8'));
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
const api = clubbPlotTasks;
let caseData = {name: 'first', benchmarks: {available_sources: ['sam', 'coamps'], sam: {s: '/sam.nc'}, coamps: {s: '/coamps.nc'}}};
const context = {triggered: []};
function start(source, selected, plots = {}) {
    const value = api.run('benchmark_' + source, [caseData, selected, plots],
        [{id: 'plots-benchmark-' + source, property: 'data'}], null, context);
    if (value !== NO) states[source] = value;
}
function finish(task, result, state = 'ready') {
    const request = submissions.filter(item => item.task === task).at(-1);
    assert(request);
    ready.set(request.scope + ':' + request.revision, state === 'ready' ? {state, result} : {state, error: result});
}
const available = selected => api.benchmarkReadiness(caseData, selected, states.sam, states.coamps, null);
(async () => {
    start('sam', ['sam']);
    assert.equal(states.sam.state, 'loading', 'busy feedback must be immediate');
    assert.deepEqual(available(['sam']).sources, []);
    // Ordinary figures can finish while the LES transfer is still pending.
    const graph = {type: 'profile-graph', index: 1};
    api.run('profile', [], [{id: graph, property: 'figure'}], graph, context);
    await wait(65);
    finish('profile', {data: ['model only']});
    await wait(180);
    assert(updates.some(([id, props]) => id === graph || id?.index === 1 && props.figure?.data[0] === 'model only'));
    assert.equal(states.sam.state, 'loading');
    finish('benchmark_sam', {request: states.sam.request, state: 'ready'});
    await wait(180);
    assert.deepEqual(available(['sam']).sources, ['sam'], 'completion enables the overlay automatically');
    assert.equal(api.benchmarkReadiness(caseData, ['sam'], states.sam, states.coamps, available(['sam'])), NO);
    const count = submissions.length;
    start('sam', ['sam', 'coamps']);
    start('coamps', ['sam', 'coamps']);
    await wait(20);
    assert.equal(submissions.length, count + 1, 'another source must not restart a ready download');
    assert.deepEqual(available(['sam', 'coamps']).sources, ['sam']);
    finish('benchmark_coamps', {request: states.coamps.request, state: 'error', error: 'offline'});
    await wait(180);
    assert.equal(states.coamps.error, 'offline');
    start('coamps', ['sam']);
    start('coamps', ['sam', 'coamps']);
    await wait(20);
    assert.equal(states.coamps.state, 'loading', 'off/on retries a failed transfer');
    const old = submissions.at(-1);
    const oldRequest = states.coamps.request;
    start('coamps', ['sam']);
    ready.set(old.scope + ':' + old.revision, {state: 'ready', result: {request: oldRequest, state: 'ready'}});
    await wait(180);
    assert.equal(states.coamps.state, 'disabled', 'late completion cannot re-enable a deselected source');
    assert.deepEqual(available(['sam']).sources, ['sam']);
    // A late result from another case cannot supply readiness to this case.
    start('coamps', ['sam', 'coamps']);
    await wait(20);
    const stale = submissions.at(-1), staleRequest = states.coamps.request;
    caseData = {...caseData, name: 'second'};
    start('sam', ['sam']);
    start('coamps', []);
    ready.set(stale.scope + ':' + stale.revision, {state: 'ready', result: {request: staleRequest, state: 'ready'}});
    await wait(180);
    assert.deepEqual(available(['sam']).sources, []);
    assert.equal(states.coamps.state, 'disabled');
    finish('benchmark_sam', {request: states.sam.request, state: 'ready'});
    await wait(180);
    assert.deepEqual(available(['sam']).sources, ['sam']);
    // LES expressions prepare their source without adding an unsolicited overlay.
    start('sam', []);
    start('sam', [], {'1': {plot_type: 'custom', expression: 'thlm + les.thlm'}});
    await wait(20);
    assert.equal(states.sam.state, 'loading');
    finish('benchmark_sam', {request: states.sam.request, state: 'ready'});
    await wait(180);
    assert.deepEqual(available([]).ready, ['sam']);
    assert.deepEqual(available([]).sources, []);
    start('sam', []); // stop polling before the process exits
})().catch(error => {console.error(error); process.exitCode = 1;});
'''
    result = subprocess.run([node, '-e', script, str(asset)], text=True, capture_output=True, timeout=8)
    assert result.returncode == 0, result.stdout + result.stderr


def test_legacy_variable_extraction_reads_the_resolved_source(roots, monkeypatch):
    """Fields handled by the existing variable catalog also use the shared lookup."""
    _write_stats(roots[1] / RELATIVE)
    case = _case(roots[0])
    monkeypatch.setattr(overlay, 'load_case_benchmark_catalog', lambda _case: {
        'native_only': {'sam': {'entries': [{'kind': 'var', 'name': 'THETAL'}], 'conversion_factor': 1.0}},
    })
    profile = overlay.extract_benchmark_profile(case, 'sam', 'native_only', 'range', 2, 0)
    assert profile['profile'].shape == (3,)
    np.testing.assert_array_equal(profile['profile'], VALUES.mean(axis=0))
