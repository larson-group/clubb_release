"""Selection latency, deferred diagnostics and refresh races without model runs."""
import inspect
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest
from dash import Dash, no_update

from dash_app.compile_tab import runner_selector
from dash_app.plot_tab import tasks
from dash_app.run_tab.callbacks_selection import register_selection_callbacks
from dash_app.tune_tab import callbacks_display as display
from tuner.sample_history import SampleHistoryWriter
from utilities.loss_metrics import LOSS_METRIC_NAMES

ROOT = Path(__file__).resolve().parents[4]


def node(script):
    """Run a bounded browser-asset scenario with deterministic event scheduling."""
    executable = shutil.which('node')
    if not executable:
        pytest.skip('Node is required for browser asset checks')
    result = subprocess.run([executable, '-e', script], cwd=ROOT, text=True,
                            capture_output=True, timeout=15)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.fixture
def diagnostic_handlers():
    """Register the real deferred callbacks and obtain their worker handlers."""
    app = Dash(__name__)
    display.register_display_callbacks(app)
    return {name: tasks._HANDLERS[name] for name in tasks.TUNE_TASKS}


def call(function, **changes):
    """Supply named callback state while leaving unused controls empty."""
    args = {name: None for name in inspect.signature(function).parameters}
    args.update(enabled=[], tab='tune')
    args.update(changes)
    return function(**args)


def test_disabled_and_hidden_diagnostics_do_no_history_io(diagnostic_handlers, monkeypatch):
    """All disabled and hidden paths must return before reading or scanning history."""
    def forbidden(*args, **kwargs):
        raise AssertionError('Off/hidden diagnostics must not touch history')
    monkeypatch.setattr(display, '_sample_history_signature', forbidden)
    monkeypatch.setattr(display, '_read_job_request', forbidden)
    for handler in diagnostic_handlers.values():
        assert all(value is no_update for value in call(handler))
        assert all(value is no_update for value in call(handler, enabled=['landscape', 'taylor'], tab='run'))


def test_only_selected_landscape_builder_runs(diagnostic_handlers, monkeypatch):
    """Choosing correlations must not invoke scatter, sensitivity or interaction builders."""
    monkeypatch.setattr(display, '_sample_history_signature', lambda path: 'chunk1')
    monkeypatch.setattr(display, '_load_sample_history', lambda path: {})
    monkeypatch.setattr(display, '_read_job_request', lambda path: {})
    monkeypatch.setattr(display, '_landscape_metric_series', lambda *args: np.array([1., 2.]))
    def forbidden(*args, **kwargs):
        raise AssertionError('An unselected plot was computed')
    for name in ['build_landscape_figure', 'build_field_sensitivity_figure', 'build_field_interaction_figure']:
        monkeypatch.setattr(display, name, forbidden)
    monkeypatch.setattr(display, 'build_parameter_correlation_figure', lambda *args, **kwargs: {'data': []})
    values = call(diagnostic_handlers['tune_landscape'], enabled=['correlation'], theme='light')
    assert values[0] is no_update and values[2] is no_update and values[3] is no_update
    assert values[1] is not no_update
    assert 'correlation' in values[4]


def test_history_cache_reuses_chunks_and_invalidates_on_publish(tmp_path, monkeypatch):
    """Repeated plots reuse decompressed data; a new published chunk becomes visible."""
    writer = SampleHistoryWriter(job_dir=tmp_path, param_names=['C8'], case_names=['bomex'],
        field_names=['wp2'], metric_names=list(LOSS_METRIC_NAMES), case_configs=[],
        case_window_counts={'bomex': 1})
    def sample(index):
        """Construct one synthetic sample with complete labeled metric metadata."""
        return {'sample_id': index, 'batch_id': index, 'all_params': {'C8': .2 + index * .1},
            'case_loss': {'bomex': 1.}, 'field_metrics': {'bomex': {'wp2': {'subwindows': [
                {metric: 1. for metric in LOSS_METRIC_NAMES}]}}}}
    writer.append([sample(0)])
    first = display._load_sample_history(tmp_path)
    with monkeypatch.context() as patch:
        patch.setattr(display.np, 'load', lambda *a, **k: pytest.fail('Unchanged chunks reloaded'))
        assert display._load_sample_history(tmp_path) is first
    writer.append([sample(1)])
    writer.flush_if_due(force=True)
    second = display._load_sample_history(tmp_path)
    assert second is not first
    assert second['sample_id'].tolist() == [0, 1]


def test_runtime_inspection_only_listens_to_explicit_session_requests(monkeypatch):
    """Opening/selecting cannot reach inspection; the initial request runs in the worker."""
    app = Dash(__name__)
    runner_selector.register_runner_selector(app)
    entry = app.callback_map['compile-jax-runtime-info.data']
    assert 'callback' not in entry
    assert entry['inputs'] == [{'id': 'runner-runtime-request', 'property': 'data'}]
    calls = []
    monkeypatch.setattr(runner_selector, 'inspect_jax_runtime_profiles',
                        lambda **kwargs: calls.append(kwargs) or {'cpu': {'status': 'ready'}})
    handler = tasks._HANDLERS['runner_info']
    assert handler(None, {}) is no_update
    assert handler({'kind': 'inventory', 'revision': 1}, {})['cpu']['status'] == 'ready'
    assert calls == [{}]
    assert tasks.PlotTasks.lane({'task': 'runner_info'}) != tasks.PlotTasks.lane({'task': 'tune_landscape'})


def test_case_selection_runs_in_browser_without_http():
    """A case toggles immediately and zero-click mount events retain the selection."""
    app = Dash(__name__)
    register_selection_callbacks(app, {'small': ['arm', 'bomex']})
    callbacks = [entry for entry in app._callback_list if 'run-selected-cases.data' in entry['output']]
    assert len(callbacks) == 3 and all('clientside_function' in entry for entry in callbacks)
    source = next(script for script in app._inline_scripts if 'const next = Array.from' in script)
    node('''const assert = require('assert'); global.window = global;
      global.dash_clientside = {no_update: 'unchanged', callback_context: {triggered_id: {name: 'arm'}, triggered: [{value: 1}]}};
    ''' + source + '''
      const fn = Object.values(dash_clientside._dashprivate_clientside_funcs).at(-1);
      assert.deepEqual(fn([1], []), ['arm']);
      assert.deepEqual(fn([2], ['arm', 'bomex']), ['bomex']);
      dash_clientside.callback_context.triggered[0].value = 0;
      assert.equal(fn([0], ['arm']), 'unchanged');
    ''')


def test_optional_plot_updates_preserve_pending_work_and_cancel_stale_delivery():
    """Unchanged data ticks preserve pending work; disabling plots invalidates it."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      global.dash_clientside = {no_update: 'unchanged', set_props() {}};
      global.setTimeout = () => 1; global.clearTimeout = () => {};
      eval(fs.readFileSync('dash_app/assets/33_plot_tasks.js', 'utf8'));
      const api = clubbPlotTasks, specs = [{id: 'tune-field-interaction-plot', property: 'figure'}];
      const context = (trigger, enabled, tab='tune') => ({triggered_id: trigger, triggered: [], inputs_list: [
        {id: 'tune-diagnostic-plots', value: enabled}, {id: 'dashboard-tabs', value: tab}]});
      api.run('tune_landscape', [], specs, null, context('tune-diagnostic-plots', []));
      assert.equal(api.scopes.size, 0);
      api.run('tune_landscape', [], specs, null, context('tune-diagnostic-plots', ['interaction']));
      const first = [...api.scopes.values()][0], revision = first.revision;
      api.run('tune_landscape', [], specs, null, context('tune-status', ['interaction']));
      assert.equal(first.revision, revision); assert.equal(first.done, false);
      api.run('tune_landscape', [], specs, null, context('tune-diagnostic-plots', []));
      assert(first.revision > revision); assert(first.done);
      api.run('tune_landscape', [], specs, null, context('tune-diagnostic-plots', ['interaction']));
      assert([...api.scopes.values()][0].revision > first.revision);
    ''')


def test_run_acknowledgement_during_idle_fetch_preserves_immediate_refresh():
    """A run acknowledged during an older idle fetch must avoid the five-second idle delay."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      let nextTimer = 0, requests = 0, complete;
      const timers = new Map();
      const button = {disabled: false, style: {}};
      global.document = {readyState: 'complete', addEventListener() {}, querySelectorAll() {return [];},
        getElementById(id) {return id === 'run-tab-layout' ? {offsetParent: {}} :
          ['run-button', 'run-cancel'].includes(id) ? button : null;}};
      global.setTimeout = (fn, delay) => {const id = ++nextTimer; timers.set(id, {fn, delay}); return id;};
      global.clearTimeout = id => timers.delete(id); global.setInterval = () => 0;
      global.fetch = () => {requests++; return new Promise(resolve => {complete = resolve;});};
      function fire() {const [id, timer] = [...timers][0]; timers.delete(id); timer.fn();}
      eval(fs.readFileSync('dash_app/assets/37_run_console_state.js', 'utf8'));
      (async () => {
        fire(); assert.equal(requests, 1);
        dash_clientside.dashboardRunConsoleState.synchronizeInputs(['arm'], {action: 'run', at: 1}, '/output');
        fire(); assert.equal(requests, 1); // No concurrent telemetry fetch.
        complete({ok: true, json: async () => ({revision: 1, runs: [], progress: {}, logs: {}})});
        for (let i=0; i<8; i++) await Promise.resolve();
        assert.deepEqual([...timers.values()].map(timer => timer.delay), [0]);
        fire(); assert.equal(requests, 2);
        complete({ok: true, json: async () => ({revision: 2, runs: [], progress: {}, logs: {}})});
        for (let i=0; i<8; i++) await Promise.resolve();
        assert.deepEqual([...timers.values()].map(timer => timer.delay), [5000]);
      })().catch(error => {console.error(error); process.exitCode=1;});
    ''')


def test_runner_opens_and_closes_locally_without_metadata_response():
    """The anchored menu opens and dismisses synchronously, before any server response."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      const updates = {}; let requests = 0;
      global.fetch = () => {requests++;}; innerWidth = 1000; innerHeight = 800;
      global.document = {addEventListener() {}, getElementById() {return {getBoundingClientRect() {return {left: 900, top: 50, bottom: 90};}};}};
      global.dash_clientside = {no_update: 'unchanged', set_props(id, props) {updates[id] = props;},
        callback_context: {triggered_id: 'run-selected-build-badge', triggered: [{prop_id:'run-selected-build-badge.n_clicks',value: 1}]}};
      eval(fs.readFileSync('dash_app/assets/35_runner_selector.js', 'utf8'));
      const anchor = clubbRunnerSelector.position();
      assert.equal(anchor.trigger_id, 'run-selected-build-badge');
      assert.equal(updates['compile-build-selector-popover'].className, 'compile-build-selector-popover');
      assert.equal(updates['runner-surface'].style.left, '472px');
      dash_clientside.callback_context.triggered_id = 'compile-build-selector-dismiss';
      dash_clientside.callback_context.triggered[0].prop_id = 'compile-build-selector-dismiss.n_clicks';
      assert.equal(clubbRunnerSelector.position(), null);
      assert(updates['compile-build-selector-popover'].className.endsWith('-closed'));
      assert.equal(requests, 0);
    ''')


def test_uptime_advances_without_status_ticks_and_stops_at_completion():
    """A browser clock keeps running through quiet workers and freezes for terminal jobs."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      let now = 0, tick; const node = {textContent: ''};
      global.performance = {now: () => now}; global.document = {getElementById: () => node};
      global.setInterval = fn => {tick = fn;}; global.dash_clientside = {no_update: null};
      eval(fs.readFileSync('dash_app/assets/36_tune_display.js', 'utf8'));
      clubbTuneDisplay.clock({job_dir: 'one', state: 'running', elapsed_seconds: 2});
      now = 3200; tick(); assert.equal(node.textContent, '5s');
      clubbTuneDisplay.clock({job_dir: 'one', state: 'running', elapsed_seconds: 2});
      now = 4200; tick(); assert.equal(node.textContent, '6s');
      clubbTuneDisplay.clock({job_dir: 'one', state: 'finished', elapsed_seconds: 7});
      now = 9000; tick(); assert.equal(node.textContent, '7s');
      clubbTuneDisplay.clock({job_dir: 'two', state: 'running', elapsed_seconds: 0});
      tick(); assert.equal(node.textContent, '0s');
    ''')


def test_late_optional_plot_response_is_discarded_after_deselection():
    """Deselect during an in-flight poll; its eventual figure must never be applied."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      const updates = []; let resolvePoll;
      global.dash_clientside = {no_update: null, set_props(id, props) {updates.push([id, props]);}};
      global.fetch = async (url, options) => {
        const data = JSON.parse(options.body);
        if (data.operation === 'submit') return {ok: true, json: async () => ({state: 'queued'})};
        return new Promise(resolve => {resolvePoll = () => resolve({ok: true, json: async () =>
          ({results: [{state: 'ready', result: {data: ['OLD']}}]})});});
      };
      eval(fs.readFileSync('dash_app/assets/33_plot_tasks.js', 'utf8'));
      const context = enabled => ({triggered_id: 'tune-diagnostic-plots', triggered: [], inputs_list: [
        {id: 'tune-diagnostic-plots', value: enabled}, {id: 'dashboard-tabs', value: 'tune'}]});
      const specs = [{id: 'tune-landscape-plot', property: 'figure'}];
      (async () => {
        clubbPlotTasks.run('tune_landscape', [{}], specs, null, context(['landscape']));
        for (let i=0; !resolvePoll && i<30; i++) await new Promise(resolve => setTimeout(resolve, 10));
        assert(resolvePoll, 'The worker result poll must be in flight');
        clubbPlotTasks.run('tune_landscape', [{}], specs, null, context([]));
        resolvePoll(); await new Promise(resolve => setTimeout(resolve, 10));
        assert(!updates.some(([id, props]) => id === 'tune-landscape-plot' && props.figure));
      })().catch(error => {console.error(error); process.exitCode=1;});
    ''')


@pytest.mark.parametrize('mode', ['mean', 'median', 'min', 'max', 'p90', 'unknown'])
def test_vectorized_observation_reductions_preserve_missing_and_penalty_values(mode):
    """All reductions match finite per-row references, including empty and penalty-only rows."""
    values = np.array([[1., 2., 5., np.nan, np.inf], [np.nan, np.inf, -np.inf, 1e100, -1e100],
                       [-3., np.nan, 4., 0., 5.]])
    reducer = {'median': np.median, 'min': np.min, 'max': np.max,
               'p90': lambda x: np.percentile(x, 90)}.get(mode, np.mean)
    expected = [reducer([1., 2., 5.]), np.nan, reducer([-3., 4., 0., 5.])]
    np.testing.assert_allclose(display._aggregate_sample_matrix(values, mode), expected, equal_nan=True)
    assert display._aggregate_sample_matrix(np.empty((2, 0)), mode).shape == (2,)


def test_native_worker_registers_static_diagnostic_tasks():
    """A spawned native worker can initialize and serialize the deferred Tune handler."""
    import time
    queue = tasks.PlotTasks()
    app = Dash(__name__)
    display.register_display_callbacks(app)
    function = tasks._HANDLERS['tune_landscape']
    args = {name: None for name in inspect.signature(function).parameters}
    args.update(enabled=[], tab='tune')
    request = {'scope': 'native:tune_landscape', 'revision': 1, 'task': 'tune_landscape',
               'args': list(args.values()), 'triggered': [], 'selection': 0}
    try:
        queue.submit(request)
        deadline = time.monotonic() + 15
        while time.monotonic() < deadline:
            result = queue.poll(request['scope'], 1)
            if result['state'] in {'ready', 'error'}:
                break
            time.sleep(.02)
        assert result['state'] == 'ready', result
        assert len(result['result']) == 5
        assert all('_dash_no_update' in item for item in result['result'])
    finally:
        queue.close()


def test_execution_tab_loaders_do_not_mount_other_workflow_bodies():
    """Each workflow loads alone; the lightweight runner state is present before visiting tabs."""
    from dash._utils import to_json
    from dash_app.app import build_dashboard_tabs
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = build_dashboard_tabs(app)
    for tab in ['compile', 'run', 'profile', 'tune']:
        entry = next(item for key, item in app.callback_map.items()
                     if f'dashboard-tabs-lazy-{tab}-loaded.data' in key)
        assert [out.component_id for out in entry['output']] == [
            f'dashboard-tab-{tab}', f'dashboard-tabs-lazy-{tab}-loaded']
    assert 'runner-selected-build' in to_json(app.layout)
    assert 'compile-discovery' in to_json(app.layout)


def test_config_save_callbacks_only_depend_on_their_mounted_tab(monkeypatch):
    """Tune saves publish a shared catalog; opening/canceling requires no Run controls."""
    from types import SimpleNamespace
    from dash_app.tune_tab import callbacks_settings as settings
    app = Dash(__name__)
    settings.register_settings_callbacks(app)
    entry = next(item for item in app.callback_map.values()
                 if item.get('callback') and item['callback'].__name__ == 'save_tuned_config')
    assert all(not item['id'].startswith('run-') for item in entry['state'])
    assert all(not output.component_id.startswith('run-') for output in entry['output'])
    callback = entry['callback'].__wrapped__
    for trigger, clicks in [('tune-config-save-open', 1), ('tune-config-save-cancel', 0), ('other', 0)]:
        monkeypatch.setattr(settings, 'callback_context', SimpleNamespace(triggered_id=trigger))
        assert len(callback(clicks, 0, 0, None, '', None, [], [], {})) == len(entry['output'])


def test_runner_ignores_unrelated_jobs_and_compile_log_heartbeats():
    """Log traffic cannot repeatedly rebuild an open menu; lifecycle changes still propagate."""
    node('''const fs = require('fs'), assert = require('assert'); global.window = global;
      global.document = {addEventListener() {}}; global.dash_clientside = {no_update: 'unchanged'};
      eval(fs.readFileSync('dash_app/assets/35_runner_selector.js', 'utf8'));
      const fn = clubbRunnerSelector.compileState;
      const jobs = {compile: {state: 'running', job: {pid: 123}, log_tail: 'one'}};
      const first = fn(jobs, {});
      assert.equal(fn({...jobs, tune: {state: 'running'}}, first), 'unchanged');
      assert.equal(fn({compile: {...jobs.compile, log_tail: 'two', updated_at: 2}}, first), 'unchanged');
      assert.equal(fn({compile: {...jobs.compile, state: 'finished'}}, first).compile.state, 'finished');
    ''')
