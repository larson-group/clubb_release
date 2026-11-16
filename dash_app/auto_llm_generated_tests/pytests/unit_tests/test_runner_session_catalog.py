"""Runner inventories persist until Compile, runtime setup or explicit refresh changes them."""
import hashlib
from pathlib import Path
import shutil
import subprocess

import pytest
from dash import Dash, no_update
from dash._utils import to_json

from dash_app.compile_tab import build_selector, callbacks, layout, runner_selector
from dash_app.plot_tab import tasks
from clubb_jax import runtime_info

ROOT = Path(__file__).resolve().parents[4]


def browser(script):
    """Exercise the actual component asset without HTTP or wall-clock sleeps."""
    executable = shutil.which('node')
    if not executable:
        pytest.skip('Node is needed for browser asset checks')
    setup = """const fs=require('fs'), assert=require('assert'); global.window=global;
    global.document={addEventListener(){}};
    global.dash_clientside={no_update:'unchanged'};
    eval(fs.readFileSync('dash_app/assets/35_runner_selector.js','utf8'));
    const api=clubbRunnerSelector;
    """
    result = subprocess.run([executable, '-e', setup + script], cwd=ROOT, text=True,
                            capture_output=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr


def test_runtime_requests_ignore_time_selection_and_initial_stamp():
    """Only startup, manual refresh and a changed setup revision request runtime work."""
    browser("""
    const first=api.runtimeRequest(0,null,null);
    assert.equal(first.kind,'inventory'); assert.equal(first.revision,1);
    Date.now=()=>9999999999999;
    for(let i=0;i<20;i++) assert.equal(api.runtimeRequest(0,null,first),'unchanged');
    assert.equal(api.runtimeRequest(0,{revision:0},first),'unchanged');
    const ready=api.runtimeRequest(0,{revision:1},first);
    assert.equal(ready.kind,'environment'); assert.equal(ready.revision,2);
    assert.equal(api.runtimeRequest(0,{revision:1},ready),'unchanged');
    const refresh=api.runtimeRequest(1,{revision:1},ready);
    assert.equal(refresh.kind,'inventory'); assert.equal(refresh.revision,3);
    """)


def test_compile_completion_reuses_mounted_compile_and_refreshes_when_unvisited():
    """Compile supplies its own updates when mounted; other completions refresh once."""
    browser("""
    const first=api.catalogRequest(0,{},null,null);
    assert.equal(first.revision,1);
    const done={compile:{state:'finished',job:{pid:8,start_time:10}}};
    assert.equal(api.catalogRequest(0,done,'compile-tab-layout',first),'unchanged');
    const completed=api.catalogRequest(0,done,null,first);
    assert.equal(completed.revision,2);
    assert.equal(api.catalogRequest(0,done,null,completed),'unchanged');
    assert.equal(api.catalogRequest(0,{compile:{state:'running'}},null,completed),'unchanged');
    assert.equal(api.catalogRequest(1,done,'compile-tab-layout',completed).revision,3);
    """)


def test_compile_and_chooser_share_one_catalog_and_compile_mount_does_not_scan(monkeypatch):
    """One pair of root stores serves both consumers; mounting Compile adds no discovery."""
    from dash_app.compile_tab import tab
    def forbidden(*args, **kwargs):
        raise AssertionError('Compile mount must reuse the shared catalog')
    monkeypatch.setattr(runner_selector, 'discover_compile_state', forbidden)
    root = to_json(build_selector.build_selector_overlay())
    body = to_json(layout.build_layout())
    assert root.count('"id":"compile-discovery"') == 1
    assert root.count('"id":"compile-build-statuses"') == 1
    assert 'compile-discovery' not in body
    assert 'runner-catalog"' not in root
    assert 'compile-tab-layout' in to_json(tab.build_tab(Dash(__name__)))
    app = Dash(__name__)
    runner_selector.register_runner_selector(app)
    assert tasks._HANDLERS['runner_catalog'](None) == (no_update, no_update)
    for entry in app.callback_map.values():
        if entry['inputs'] == [{'id': 'runner-runtime-request', 'property': 'data'}]:
            break
    else:
        pytest.fail('Runtime inspection must be driven only by its request store')


def cpu_report(tmp_path, monkeypatch):
    """Create a real CPU inventory report for an absent synthetic managed environment."""
    requirements = tmp_path / 'requirements.txt'
    requirements.write_text('jax==0.11.0\n')
    venv = tmp_path / 'venv'
    monkeypatch.setattr(runtime_info, '_installed_runtime', lambda path: {})
    return runtime_info.inspect_runtime('cpu', 'cpu', requirements, venv, '0.11.0'), requirements, venv


def test_setup_refreshes_readiness_without_hardware_queries_or_input_mutation(tmp_path, monkeypatch):
    """A newly prepared environment changes readiness while preserving the original inventory."""
    original, requirements, venv = cpu_report(tmp_path, monkeypatch)
    def forbidden(*args, **kwargs):
        raise AssertionError('A readiness refresh must never rediscover hardware')
    monkeypatch.setattr(runtime_info.backends.cpu, 'hardware_info', forbidden)
    monkeypatch.setattr(runtime_info.backends.cpu, 'inspect_devices', forbidden)
    venv.mkdir()
    (venv / '.clubb-jax-requirements.sha256').write_text(hashlib.sha256(requirements.read_bytes()).hexdigest())
    monkeypatch.setattr(runtime_info, '_installed_runtime', lambda path: {
        'python': '3.12.0', 'packages': runtime_info.backends.cpu.expected_packages('0.11.0')})
    updated = runtime_info.refresh_environment(original, requirements)
    assert updated['status'] == 'ready'
    assert updated['runtime']['environment_ready'] is True
    assert updated['hardware'] == original['hardware']
    assert original['status'] == 'setup_required'
    monkeypatch.setattr('clubb_jax.run_jax.runtime_paths', lambda accelerator: (requirements, venv, 'cpu'))
    monkeypatch.setattr(runner_selector, 'inspect_jax_runtime_profiles', forbidden)
    runner_selector.register_runner_selector(Dash(__name__))
    assert tasks._HANDLERS['runner_info']({'kind': 'environment'}, {'cpu': original})['cpu']['status'] == 'ready'


def test_environment_stamp_only_invalidates_after_a_change(tmp_path, monkeypatch):
    """The existing heartbeat notices setup completion once, using only managed file stats."""
    report, requirements, venv = cpu_report(tmp_path, monkeypatch)
    first = runner_selector.environment_change({'cpu': report}, None)
    assert first['revision'] == 0
    assert runner_selector.environment_change({'cpu': report}, first) is no_update
    venv.mkdir()
    (venv / '.clubb-jax-requirements.sha256').write_text('prepared')
    next_value = runner_selector.environment_change({'cpu': report}, first)
    assert next_value['revision'] == 1
    assert runner_selector.environment_change({'cpu': report}, next_value) is no_update


def test_gpu_selection_uses_cached_compatibility_and_common_environment(monkeypatch):
    """Switching GPUs uses static compatibility checks and enables preallocation immediately."""
    info = {'status': 'unavailable', 'selectable': False, 'reason': 'old default unavailable',
        'runtime': {'accelerator': 'cuda13', 'environment_ready': True},
        'hardware': {'gpus': [
            {'uuid': 'GPU-old', 'index': 0, 'name': 'Old', 'driver_version': '580.126', 'compute_capability': '5.0'},
            {'uuid': 'GPU-new', 'index': 1, 'name': 'New', 'driver_version': '580.126', 'compute_capability': '8.0'}]}}
    def forbidden(*args, **kwargs):
        raise AssertionError('Selecting a GPU must not probe hardware')
    monkeypatch.setattr(runner_selector, 'inspect_jax_runtime_profiles', forbidden)
    old = build_selector.jax_device_info(info, 'GPU-old')
    new = build_selector.jax_device_info(info, 'GPU-new')
    assert old['status'] == 'unavailable'
    assert new['status'] == 'ready' and new['hardware']['selected_gpu']['uuid'] == 'GPU-new'
    assert callbacks.jax_preallocation_available({'gpu': info}, 'gpu', 'GPU-new')
    assert info['hardware'].get('selected_gpu') is None
    assert build_selector.jax_device_info(info, 'GPU-missing')['selectable'] is False


def test_shared_catalog_writers_tolerate_unmounted_compile_controls():
    """Root stores exist before Compile; its callbacks must ignore absent local controls."""
    app = Dash(__name__)
    callbacks.register_compile_callbacks(app)
    functions = {entry['callback'].__wrapped__.__name__: entry for entry in app.callback_map.values()
                 if entry.get('callback')}
    status = functions['update_build_statuses']
    assert all(item.get('allow_optional') for item in status['inputs'] if item['id'] != 'compile-discovery')
    assert status['callback'].__wrapped__({}, None, None) is no_update
    selected = functions['select_existing_build']
    assert all(value is no_update for value in selected['callback'].__wrapped__([], None, {}, {}, None, None, None))
    assert functions['refresh_discovery']['callback'].__wrapped__(None) is no_update
    assert functions['refresh_discovery_after_job']['callback'].__wrapped__(None) is no_update


def test_stopped_compile_refreshes_partial_catalog_once():
    """A stopped rebuild can leave new directories, so unvisited Compile still invalidates once."""
    browser("""
    const initial=api.catalogRequest(0,{},null,null);
    const stopped={compile:{state:'stopped',job:{pid:123,start_time:10}}};
    const refreshed=api.catalogRequest(0,stopped,null,initial);
    assert.equal(refreshed.revision,initial.revision+1);
    assert.equal(api.catalogRequest(0,stopped,null,refreshed),'unchanged');
    """)


def test_compile_catalog_update_rejects_an_older_in_flight_scan():
    """Publishing Compile's newer catalog prevents an old native scan from overwriting it."""
    browser("""
    let finishOld; const updates=[];
    dash_clientside.set_props=(id,props)=>updates.push([id,props]);
    global.fetch=async (url,options)=>{
      const data=JSON.parse(options.body);
      if(data.operation==='submit')return {ok:true,json:async()=>({state:'queued'})};
      return new Promise(resolve=>{finishOld=()=>resolve({ok:true,json:async()=>({results:[
        {state:'ready',result:[{builds:['OLD']},{}]}]})});});
    };
    eval(fs.readFileSync('dash_app/assets/33_plot_tasks.js','utf8'));
    (async()=>{
      clubbPlotTasks.run('runner_catalog',[{revision:1}],
        [{id:'compile-discovery',property:'data'},{id:'compile-build-statuses',property:'data'}],
        null,{triggered:[]});
      for(let i=0;!finishOld&&i<30;i++)await new Promise(resolve=>setTimeout(resolve,10));
      assert(finishOld);
      dash_clientside.set_props('compile-discovery',{data:{builds:['NEW']}});
      clubbPlotTasks.cancel('runner_catalog');
      finishOld(); await new Promise(resolve=>setTimeout(resolve,10));
      assert.deepEqual(updates.filter(([id])=>id==='compile-discovery').at(-1)[1].data.builds,['NEW']);
    })().catch(error=>{console.error(error);process.exitCode=1;});
    """)
