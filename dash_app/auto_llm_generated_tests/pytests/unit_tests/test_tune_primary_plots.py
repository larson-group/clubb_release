"""Permanent Tune diagnostics follow retained results and published sample history."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest
from dash import Dash
from dash.development.base_component import Component

from dash_app.tune_tab import callbacks_display, runtime, tab
from dash_app.tune_tab.layout import build_layout
from dash_app.tune_tab.state import PRIMARY_DIAGNOSTIC_PLOTS

ROOT = Path(__file__).resolve().parents[4]


def nodes(component):
    """Walk the real layout to check user-visible controls and initial visibility."""
    if not isinstance(component, Component):
        return
    yield component
    children = getattr(component, 'children', None)
    for child in children if isinstance(children, (list, tuple)) else [children]:
        yield from nodes(child)


def browser(script):
    """Run the actual task transport/publisher without browser or worker processes."""
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is required for browser update checks')
    result = subprocess.run([node, '-e', script], cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr


def test_primary_plots_are_visible_and_have_no_hide_controls():
    """The three required graphs are visible initially; only extra plots have checkboxes."""
    layout = build_layout(tab.build_initial_tune_state())
    by_id = {item.id: item for item in nodes(layout) if isinstance(getattr(item, 'id', None), str)}
    assert by_id['tune-diagnostic-plots'].data == list(PRIMARY_DIAGNOSTIC_PLOTS)
    assert {item['value'] for item in by_id['tune-postprocessing'].options} == {
        'correlation', 'sensitivity', 'interaction'}
    for name in ['tune-taylor-diagram', 'tune-parameter-box-plot', 'tune-landscape-plot',
                 'tune-taylor-container', 'tune-landscape-container']:
        assert getattr(by_id[name], 'style', {}).get('display') != 'none'
    for name in ['tune-parameter-correlation-plot', 'tune-field-sensitivity-plot', 'tune-field-interaction-plot']:
        assert by_id[name].style['display'] == 'none'


def test_optional_selection_cannot_remove_primary_plots():
    """The real client publisher always includes the required plots, even for no extras."""
    app = Dash(__name__)
    callbacks_display.register_display_callbacks(app)
    script = next(item for item in app._inline_scripts if 'Array.from(new Set' in item)
    browser("global.window=global;global.dash_clientside={};const assert=require('assert');" + script + '''
      const fn=Object.values(dash_clientside._dashprivate_clientside_funcs)[0];
      assert.deepEqual(fn([]),['taylor','parameters','landscape']);
      assert.deepEqual(fn(null),['taylor','parameters','landscape']);
      assert.deepEqual(fn(['sensitivity']),['taylor','parameters','landscape','sensitivity']);
      assert.deepEqual(fn(['junk','taylor']),['taylor','parameters','landscape']);
    ''')


def test_status_detects_chunk_publication_with_unchanged_sample_count(tmp_path, monkeypatch):
    """A timed history flush invalidates plots even when the displayed sample count is unchanged."""
    path = tmp_path / 'status.json'
    path.write_text(json.dumps({'job_dir': str(tmp_path), 'state': 'running', 'samples_evaluated': 42}))
    history = tmp_path / 'sample_history'
    history.mkdir()
    # Invalid chunk contents prove status only stats files; it never loads arrays.
    (history / 'sample_history_000001.npz').write_bytes(b'not opened')
    monkeypatch.setattr(callbacks_display.np, 'load', lambda *a, **k: pytest.fail('Status loaded history arrays'))
    first = runtime.read_tuning_status(path)
    assert len(first['sample_history_signature']) == 64
    path.write_text(json.dumps({**first, 'elapsed_seconds': 99}))
    assert runtime.read_tuning_status(path)['sample_history_signature'] == first['sample_history_signature']
    (history / 'sample_history_000002.npz').write_bytes(b'new publication')
    second = runtime.read_tuning_status(path)
    assert second['samples_evaluated'] == first['samples_evaluated']
    assert second['sample_history_signature'] != first['sample_history_signature']


def test_automatic_refresh_is_data_driven_and_does_not_starve_pending_work():
    """New chunks/results refresh; heartbeat-only ticks skip and pending work can complete."""
    browser('''const fs=require('fs'),assert=require('assert');global.window=global;
      global.dash_clientside={no_update:'unchanged',set_props(){}};
      global.setTimeout=()=>1;global.clearTimeout=()=>{};
      eval(fs.readFileSync('dash_app/assets/33_plot_tasks.js','utf8'));
      const api=clubbPlotTasks;
      const context=trigger=>({triggered_id:trigger,triggered:[],inputs_list:[
        {id:'tune-diagnostic-plots',value:['taylor','parameters','landscape']},
        {id:'dashboard-tabs',value:'tune'}]});
      const spec=[{id:'tune-landscape-plot',property:'figure'}];
      const status=(stamp,elapsed=1)=>({job_dir:'/job',sample_history_signature:stamp,elapsed_seconds:elapsed});
      function run(stamp,trigger='tune-status') {api.run('tune_landscape',[status(stamp)],spec,null,context(trigger));}
      run('one','tune-diagnostic-plots');
      const first=[...api.scopes.values()][0];first.done=first.success=true;
      run('one');assert.strictEqual([...api.scopes.values()][0],first);
      run('two');const next=[...api.scopes.values()][0];assert(next.revision>first.revision);
      run('three');assert.strictEqual([...api.scopes.values()][0],next,'Do not discard pending data work');
      next.done=next.success=true;
      run('three');assert([...api.scopes.values()][0].revision>next.revision,'The next tick must catch up');
      const args=best=>[best,{},'aggregate','normalized','average','dark',
        ['taylor','parameters','landscape'],0,'tune',status('three'),''];
      const diag=[{id:'tune-taylor-diagram',property:'figure'}];
      api.run('tune_diagnostics',args([{total_loss:2}]),diag,null,context('tune-diagnostic-plots'));
      const before=[...api.scopes.values()].find(item=>item.task==='tune_diagnostics');before.done=before.success=true;
      api.run('tune_diagnostics',args([{total_loss:2}]),diag,null,context('tune-status'));
      assert.strictEqual([...api.scopes.values()].find(item=>item.task==='tune_diagnostics'),before);
      api.run('tune_diagnostics',args([{total_loss:1}]),diag,null,context('tune-best-results'));
      assert([...api.scopes.values()].find(item=>item.task==='tune_diagnostics').revision>before.revision);
    ''')


def test_last_data_update_is_rendered_after_pending_work_without_another_tick():
    """A final publication during a slow plot refresh must be delivered even after polling stops."""
    browser('''const fs=require('fs'),assert=require('assert');global.window=global;
      const updates=[],submissions=[];let finishFirst;
      global.dash_clientside={no_update:null,set_props(id,props){updates.push([id,props]);}};
      global.fetch=async(url,options)=>{
        const data=JSON.parse(options.body);
        if(data.operation==='submit'){submissions.push(data);return {ok:true,json:async()=>({state:'queued'})};}
        if(data.requests[0].revision===1)return new Promise(resolve=>{
          finishFirst=()=>resolve({ok:true,json:async()=>({results:[{state:'ready',result:{data:['first']}}]})});
        });
        return {ok:true,json:async()=>({results:[{state:'ready',result:{data:['last']}}]})};
      };
      eval(fs.readFileSync('dash_app/assets/33_plot_tasks.js','utf8'));
      const context=trigger=>({triggered_id:trigger,triggered:[],inputs_list:[
        {id:'tune-diagnostic-plots',value:['taylor','parameters','landscape']},
        {id:'dashboard-tabs',value:'tune'}]});
      const spec=[{id:'tune-landscape-plot',property:'figure'}];
      (async()=>{
        clubbPlotTasks.run('tune_landscape',[{job_dir:'/job',sample_history_signature:'first'}],spec,null,context('tune-diagnostic-plots'));
        for(let i=0;!finishFirst&&i<100;i++)await new Promise(r=>setTimeout(r,10));
        assert(finishFirst);
        clubbPlotTasks.run('tune_landscape',[{job_dir:'/job',sample_history_signature:'last'}],spec,null,context('tune-status'));
        assert.equal(submissions.length,1,'Keep the in-flight figure instead of continuously superseding it');
        finishFirst();
        for(let i=0;!updates.some(([id,p])=>id==='tune-landscape-plot'&&p.figure?.data[0]==='last')&&i<150;i++)
          await new Promise(r=>setTimeout(r,10));
        assert.equal(submissions.length,2);
        assert(updates.some(([id,p])=>id==='tune-landscape-plot'&&p.figure?.data[0]==='last'));
      })().catch(error=>{console.error(error);process.exitCode=1;});
    ''')
