"""Native Plot queue isolation, stale results, and browser-first interactions."""
from concurrent.futures import Future
import json
from pathlib import Path
import shutil
import subprocess
import time

import pytest
from dash import Dash, html

from dash_app.plot_tab import tasks
from dash_app.plot_tab.async_callbacks import register_task_routes


class Executor:
    def __init__(self):
        self.calls = []
    def submit(self, function, request):
        future = Future()
        self.calls.append((request, future))
        return future
    def shutdown(self, **kwargs):
        pass


def request(scope='browser:profile:1', revision=1, task='profile', args=None):
    return dict(scope=scope, revision=revision, task=task, args=args or [])


def test_only_newest_queued_request_runs_and_old_result_cannot_replace_it():
    executor = Executor()
    queue = tasks.PlotTasks(lambda: executor)
    queue.submit(request(revision=1))
    queue.submit(request(revision=2))
    queue.submit(request(revision=3))
    assert len(executor.calls) == 1
    executor.calls[0][1].set_result('["old"]')
    assert [call[0]['revision'] for call in executor.calls] == [1, 3]
    assert queue.poll('browser:profile:1', 3)['state'] == 'running'
    executor.calls[1][1].set_result('["new"]')
    assert queue.poll('browser:profile:1', 3)['result'] == ['new']
    assert queue.poll('browser:profile:1', 1)['state'] == 'superseded'
    queue.submit(request(revision=2))
    assert len(executor.calls) == 2


def test_metadata_and_figures_have_independent_lanes_and_browser_scopes():
    pools = []
    def factory():
        executor = Executor()
        pools.append(executor)
        return executor
    queue = tasks.PlotTasks(factory)
    queue.submit(request())
    queue.submit(request(scope='other:profile:1'))
    queue.submit(request(scope='browser:catalog', task='catalog'))
    assert len(pools) == 2
    assert len(pools[0].calls) == 1 and len(pools[1].calls) == 1
    pools[1].calls[0][1].set_result('[]')
    assert queue.poll('browser:catalog', 1)['state'] == 'ready'
    assert queue.poll('browser:profile:1', 1)['state'] == 'running'
    pools[0].calls[0][1].set_result('[]')
    assert pools[0].calls[1][0]['scope'] == 'other:profile:1'


def test_new_selection_discards_other_cards_queued_for_the_old_case():
    executor = Executor()
    queue = tasks.PlotTasks(lambda: executor)
    queue.submit({**request(scope='page:profile:1'), 'selection': 1})
    queue.submit({**request(scope='page:profile:2'), 'selection': 1})
    queue.submit({**request(scope='page:case', task='case'), 'selection': 2})
    assert 'page:profile:2' not in queue.pending
    assert queue.submit({**request(scope='page:profile:3'), 'selection': 1})['state'] == 'superseded'
    executor.calls[0][1].set_result('["old"]')
    assert queue.poll('page:profile:1', 1)['state'] == 'missing'


def test_expired_inflight_record_cannot_replace_a_resubmitted_request():
    executor = Executor()
    queue = tasks.PlotTasks(lambda: executor)
    queue.submit(request())
    queue.records['browser:profile:1']['touched'] -= 121
    queue.submit(request())
    executor.calls[0][1].set_result('["expired"]')
    assert queue.poll('browser:profile:1', 1)['state'] == 'running'
    executor.calls[1][1].set_result('["fresh"]')
    assert queue.poll('browser:profile:1', 1)['result'] == ['fresh']


def test_failure_and_result_budget_do_not_stall_queue():
    executor = Executor()
    queue = tasks.PlotTasks(lambda: executor, max_scopes=2, max_bytes=20)
    queue.submit(request())
    executor.calls[0][1].set_exception(ValueError('Unreadable NetCDF'))
    assert queue.poll('browser:profile:1', 1)['error'] == 'Unreadable NetCDF'
    queue.submit(request(revision=2))
    executor.calls[1][1].set_result(json.dumps('x' * 30))
    assert queue.poll('browser:profile:1', 2)['state'] == 'error'
    for scope in ['b', 'c']:
        queue.submit(request(scope=scope))
    assert len(queue.records) == 2
    assert queue.poll('browser:profile:1', 2)['state'] == 'missing'
    queue.close()
    with pytest.raises(RuntimeError, match='stopping'):
        queue.submit(request())


def test_one_poll_can_deliver_multiple_cards_without_losing_revision_checks():
    executor = Executor()
    queue = tasks.PlotTasks(lambda: executor)
    queue.submit(request())
    executor.calls[0][1].set_result('["figure"]')
    response = queue.poll_many([
        {'scope': 'browser:profile:1', 'revision': 1},
        {'scope': 'browser:profile:1', 'revision': 2},
        {'scope': 'another:profile:1', 'revision': 1},
    ])
    assert [item['state'] for item in response['results']] == ['ready', 'superseded', 'missing']
    assert response['results'][0]['result'] == ['figure']
    with pytest.raises(ValueError, match='batch'):
        queue.poll_many([{}] * 129)


def test_initialized_ui_can_poll_without_a_full_broker_status_probe(monkeypatch):
    from io import BytesIO
    from dash_app.shared import broker_client
    checks = []
    monkeypatch.setattr(broker_client, 'ensure_broker', lambda: checks.append(True))
    monkeypatch.setattr(broker_client, 'read_connection', lambda: {'url': 'http://localhost', 'token': 'test'})
    monkeypatch.setattr(broker_client, 'urlopen', lambda *a, **kw: BytesIO(b'{"results": []}'))
    assert broker_client.perform_action('plot_task_poll', {}, ensure_running=False) == {'results': []}
    assert not checks
    broker_client.perform_action('plot_task_poll', {})
    assert checks == [True]


def test_real_worker_inspects_missing_directory_and_serializes_no_update(tmp_path):
    # Spawn a real worker; no scientific build, live broker or dashboard needed.
    queue = tasks.PlotTasks()
    try:
        missing = str(tmp_path / 'missing')
        payload = request(scope='test:catalog', task='catalog', args=[True, 0, [missing], []])
        payload['triggered'] = [{'prop_id': 'plots-output-dirs.data', 'value': [missing]}]
        queue.submit(payload)
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            response = queue.poll('test:catalog', 1)
            if response['state'] in {'ready', 'error'}:
                break
            time.sleep(.02)
        assert response['state'] == 'ready', response
        assert response['result'][0]['available'] is False
        assert 'No readable' in response['result'][0]['error']
    finally:
        queue.close()


def test_native_route_only_proxies_typed_plot_actions(monkeypatch):
    from dash_app.shared import broker_client
    calls = []
    monkeypatch.setattr(broker_client, 'perform_action', lambda name, data, **kw: calls.append((name, data)) or {'state': 'queued'})
    app = Dash(__name__)
    app.layout = html.Div()
    register_task_routes(app)
    client = app.server.test_client()
    assert client.post('/plots/tasks', json={'operation': 'submit', **request()}).status_code == 200
    assert calls[0][0] == 'plot_task_submit'
    assert client.post('/plots/tasks', json={'operation': 'run_shell'}).status_code == 400
    assert client.post('/plots/tasks', json={'operation': 'submit'}, headers={'Origin': 'https://example.org'}).status_code == 403


def test_browser_updates_selection_immediately_and_ignores_stale_results():
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is needed for browser callback tests')
    assets = Path(__file__).resolve().parents[1] / 'assets'
    script = r'''
const assert = require('assert');
const fs = require('fs');
global.window = global;
const updates = [], submissions = [], pollBatches = [], ready = new Map();
const NO = {no_update: true};
window.dash_clientside = {no_update: NO, callback_context: {}, set_props: (id, props) => updates.push([id, props])};
global.fetch = async (_url, options) => {
    const data = JSON.parse(options.body);
    if (data.operation === 'submit') { submissions.push(data); return {ok: true, json: async () => ({state: 'queued'})}; }
    pollBatches.push(data.requests.length);
    return {ok: true, json: async () => ({results:data.requests.map(item => ready.get(item.revision) || {state: 'running'})})};
};
for (const path of process.argv.slice(1)) eval(fs.readFileSync(path, 'utf8'));
const wait = ms => new Promise(resolve => setTimeout(resolve, ms));
(async () => {
    const paths = {repo: '/repo', output: '/repo/output', home: '/home/test'};
    const chooser = dash_clientside.plotsChooser;
    dash_clientside.callback_context = {triggered_id: {type: 'plots-add-output-dir', path: '/b'}, triggered: [{value: 123}]};
    const chosen = chooser.choose([], 0, [], ['/a'], '', [{path: '/b'}], paths);
    assert.deepEqual(chosen[0], ['/a', '/b']);
    assert.equal(submissions.length, 0, 'selection must never wait for an HTTP request');
    const rendered = chooser.render([], ['/a'], false, null, paths);
    assert.equal(rendered[1][0].props.children[1].props.children, 'Loading…');
    dash_clientside.callback_context = {triggered_id: 'plots-add-extra-dir', triggered: [{value: 123}]};
    assert.deepEqual(chooser.choose([],1,[],[], '~/run/../data', [], paths)[0], ['/home/test/data']);
    const outputs = [{id: 'plots-case-data', property: 'data'}];
    clubbPlotTasks.run('case', [], outputs, null, {triggered_id: {type:'plots-case-button',name:'bomex'},triggered:[{value:0}]});
    assert.equal(clubbPlotTasks.scopes.size, 0, 'mounting a button must not create a load request');
    // Dash may batch a button mount with the real directory change.
    const context = {triggered_id: {type:'plots-case-button',name:'bomex'}, triggered: [
        {prop_id:'{"name":"bomex","type":"plots-case-button"}.n_clicks',value:0},
        {prop_id: 'plots-output-dirs.data', value: ['/a']}]};
    assert.equal(clubbPlotTasks.run('case', [], outputs, null, context), NO);
    await wait(10);
    assert.equal(submissions.at(-1).triggered[0].prop_id, 'plots-output-dirs.data');
    clubbPlotTasks.run('case', [], outputs, null, context);
    await wait(10);
    ready.set(1, {state: 'ready', result: {name: 'OLD'}});
    ready.set(2, {state: 'ready', result: {name: 'NEW'}});
    await wait(180);
    const cases = updates.filter(([id]) => id === 'plots-case-data');
    assert.deepEqual(cases.map(([,props]) => props.data.name), ['NEW']);
    // A late figure cannot overwrite the display after a directory change.
    const graph = {type: 'profile-graph', index: 9};
    ready.clear();
    clubbPlotTasks.run('profile', [], [{id: graph, property: 'figure'}], graph, {triggered: []});
    await wait(60);
    const count = submissions.length;
    clubbPlotTasks.run('profile', [], [{id: graph, property: 'figure'}], graph, {triggered: [{prop_id:'plots-selected-column.data',value:0}]});
    await wait(60);
    assert.equal(submissions.length, count, 'unchanged column feedback must not restart the figure');
    clubbPlotTasks.invalidateSelection();
    ready.set(1, {state: 'ready', result: {data: ['OLD FIGURE']}});
    await wait(160);
    assert.equal(updates.filter(([id,p]) => typeof id === 'object' && id.type === 'profile-graph' && p.figure).length, 0);
    ready.clear();
    clubbPlotTasks.run('params', [], [{id:'plots-param-data',property:'data'}], null, {triggered:[]});
    ready.set(1, {state:'error',error:'Unreadable parameters'});
    await wait(160);
    assert(updates.some(([id,p]) => id === 'plots-param-status' && p.children.includes('Unreadable parameters')));
    assert.equal(updates.filter(([id]) => id === 'plots-param-panel').at(-1)[1].className, 'plots-param-panel--loading');
    clubbPlotTasks.run('params', [], [{id:'plots-param-data',property:'data'}], null, {triggered:[]});
    ready.set(2, {state:'ready',result:{params:{}}});
    await wait(160);
    assert.equal(updates.filter(([id]) => id === 'plots-param-panel').at(-1)[1].className, '');
    ready.clear();
    clubbPlotTasks.run('case', [], outputs, null, context);
    ready.set(3, {state:'ready',result:{name:'THIRD'}});
    await wait(160);
    assert(updates.some(([id,p]) => id?.type === 'profile-graph' && p.relayoutData === null));
    const graphs = [10,11,12].map(index=>({type:'profile-graph',index}));
    ready.set(1, {state:'ready',result:{data:[]}});
    graphs.forEach(id=>clubbPlotTasks.run('profile', [], [{id,property:'figure'}], id, {triggered:[]}));
    await wait(220);
    assert(pollBatches.some(size=>size>=3), 'cards must share a polling request');
    const completeCount = submissions.length;
    clubbPlotTasks.run('profile', [], [{id:graphs[0],property:'figure'}], graphs[0], {triggered:[]});
    await wait(60);
    assert.equal(submissions.length, completeCount, 'completed figures ignore unchanged inputs');
    ready.set(2, {state:'ready',result:{data:[]}});
    clubbPlotTasks.run('profile', ['new variable'], [{id:graphs[0],property:'figure'}], graphs[0], {triggered:[]});
    await wait(220);
    assert.equal(submissions.length, completeCount+1, 'changed inputs must still render');
    assert.deepEqual(clubbPlotTasks.patchFigure({data:[{x:[1]}]}, {operations:[{operation:'Assign',location:['data',0,'x'],params:{value:[2]}}]}), {data:[{x:[2]}]});
})().catch(error => {console.error(error); process.exitCode = 1;});
'''
    result = subprocess.run([node, '-e', script, str(assets / '33_plot_tasks.js'), str(assets / '34_plot_chooser.js')], text=True, capture_output=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
