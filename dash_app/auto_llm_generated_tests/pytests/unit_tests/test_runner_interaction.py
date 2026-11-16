"""Runner clicks paint locally and publish the exact intended target before Launch."""
import json
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import pytest
from dash import Dash

from dash_app.compile_tab import callbacks, runner_selector

ROOT = Path(__file__).resolve().parents[4]


def test_batched_mount_and_click_resolve_the_positive_event_id(monkeypatch):
    """Reproduce the live request: a Fortran zero mount precedes the real JAX click."""
    first = {'type': 'compile-run-implementation-choice', 'index': 'fortran'}
    clicked = {'type': 'compile-run-implementation-choice', 'index': 'jax'}
    monkeypatch.setattr(callbacks, 'callback_context', SimpleNamespace(
        triggered_id=first, triggered=[
            {'prop_id': json.dumps(first) + '.n_clicks', 'value': 0},
            {'prop_id': json.dumps(clicked) + '.n_clicks', 'value': 1},
        ]))
    assert callbacks.clicked_trigger_id() == clicked


def test_metadata_rendering_has_no_selection_or_anchor_inputs():
    """Preference changes cannot queue a server render or a server selection mutation."""
    app = Dash(__name__)
    runner_selector.register_runner_selector(app)
    assert all(not any(key.startswith(name) for name in ['compile-run-implementation.data', 'compile-run-jax-profile.data', 'compile-run-jax-xla-prealloc.data']) for key in app.callback_map)
    entry = next(item for item in app.callback_map.values()
                 if item.get('callback') and item['callback'].__name__ == 'update_build_selector')
    assert {item['id'] for item in entry['inputs']} == {
        'compile-discovery', 'compile-build-statuses', 'runner-failures',
        'runner-compile-state', 'compile-jax-runtime-info'}
    panels, views = entry['callback'].__wrapped__({}, {}, {}, {}, {}, None, 'fortran', 'cpu', '', False, 'run')
    assert panels[1].id == 'runner-fortran-options'
    assert panels[2].id == 'runner-jax-options'
    assert set(views) == {'fortran', 'jax'}


def test_native_click_paints_before_store_publication_and_launch_flushes_it():
    """The actual delegated event skips Dash's click queue and prevents stale Launch state."""
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is needed for the native click lifecycle')
    script = r'''
    const fs=require('fs'),assert=require('assert');global.window=global;
    const frames=[],updates=[],listeners={};
    function button(index,type='compile-run-implementation-choice') {
      const attrs={},classes=new Set();
      return {id:JSON.stringify({index,type}),attrs,dataset:{},style:{},
        classList:{toggle(name,on){if(on)classes.add(name);else classes.delete(name);}},
        setAttribute(key,value){attrs[key]=value;},querySelector(){return null;}};
    }
    const fortran=button('fortran'),jax=button('jax');
    const panels={'runner-fortran-options':{style:{}},'runner-jax-options':{style:{}},
      'runner-context-note':{textContent:''}};
    global.document={
      addEventListener(name,fn){listeners[name]=fn;},
      getElementById(id){return panels[id]||null;},querySelector(){return null;},
      querySelectorAll(selector){return selector.includes('implementation-choice')?[fortran,jax]:[];}
    };
    global.requestAnimationFrame=fn=>{frames.push(fn);return frames.length;};
    global.dash_clientside={no_update:'unchanged',set_props(id,props){updates.push([id,props]);}};
    global.fetch=()=>{throw new Error('Selection must not use HTTP');};
    eval(fs.readFileSync('dash_app/assets/35_runner_selector.js','utf8'));
    const api=clubbRunnerSelector,views={fortran:{},jax:{cpu:{},gpu:{can_preallocate:true},'GPU-B':{can_preallocate:true}}};
    api.synchronize(views,null,'fortran','cpu','',false);
    function click(target) {
      let prevented=false,stopped=false;
      listeners.click({target:{closest:()=>target},preventDefault(){prevented=true;},
        stopImmediatePropagation(){stopped=true;}});
      return {prevented,stopped};
    }
    assert.deepEqual(click(jax),{prevented:true,stopped:true});
    assert.equal(jax.attrs['aria-pressed'],'true');
    assert.equal(panels['runner-jax-options'].style.display,'');
    assert.equal(panels['runner-fortran-options'].style.display,'none');
    assert.deepEqual(updates,[],'The first paint must precede global callback work');
    api.synchronize(views,null,'fortran','cpu','',false);
    assert.equal(jax.attrs['aria-pressed'],'true','An older callback cannot undo the click');
    click({id:'run-button'});
    assert.deepEqual(updates.at(-1),['compile-run-implementation',{data:'jax'}],
      'Launch must see JAX even when clicked before the next frame');
    api.synchronize(views,null,'jax','cpu','',false);
    click(fortran);
    while(frames.length)frames.shift()();
    assert.deepEqual(updates.at(-1),['compile-run-implementation',{data:'fortran'}]);
    assert.equal(fortran.attrs['aria-pressed'],'true');
    api.chooseCompute('GPU-B');
    while(frames.length)frames.shift()();
    assert.deepEqual(updates.slice(-2),[
      ['compile-run-jax-gpu',{data:'GPU-B'}],['compile-run-jax-profile',{data:'gpu'}]]);
    api.synchronize(views,null,'fortran','gpu','GPU-B',false);
    api.chooseCompute('cpu'); while(frames.length)frames.shift()();
    assert.deepEqual(updates.at(-1),['compile-run-jax-profile',{data:'cpu'}]);
    api.synchronize(views,null,'fortran','cpu','GPU-B',false);
    api.chooseCompute('GPU-B'); while(frames.length)frames.shift()();
    assert.deepEqual(updates.at(-1),['compile-run-jax-profile',{data:'gpu'}]);
    api.synchronize(views,null,'fortran','gpu','GPU-B',false);
    api.chooseAllocation(true); while(frames.length)frames.shift()();
    assert.deepEqual(updates.at(-1),['compile-run-jax-xla-prealloc',{data:true}]);
    api.synchronize(views,null,'fortran','gpu','GPU-B',true);
    api.chooseAllocation(false); while(frames.length)frames.shift()();
    assert.deepEqual(updates.at(-1),['compile-run-jax-xla-prealloc',{data:false}]);
    const count=updates.length;
    api.chooseCompute('missing'); while(frames.length)frames.shift()();
    assert.equal(updates.length,count);
    api.synchronize(views,null,'fortran','cpu','GPU-B',false);
    api.chooseAllocation(true); while(frames.length)frames.shift()();
    assert.equal(updates.length,count);
    '''
    result = subprocess.run([node, '-e', script], cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr


def test_menu_opener_resolves_a_real_click_after_other_badge_mounts():
    """A zero-valued Run mount cannot swallow the Profile opener's actual click."""
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is required for browser callback checks')
    script = r'''
    const fs=require('fs'),assert=require('assert');global.window=global;
    global.innerWidth=1000;global.innerHeight=800;
    global.document={addEventListener(){},getElementById(){return {
      getBoundingClientRect(){return {left:10,top:10,bottom:40};}};}};
    global.dash_clientside={no_update:'unchanged',set_props(){},callback_context:{
      triggered_id:'run-selected-build-badge',triggered:[
        {prop_id:'run-selected-build-badge.n_clicks',value:0},
        {prop_id:'profile-selected-build-badge.n_clicks',value:1}]}};
    eval(fs.readFileSync('dash_app/assets/35_runner_selector.js','utf8'));
    assert.equal(clubbRunnerSelector.position().trigger_id,'profile-selected-build-badge');
    '''
    result = subprocess.run([node, '-e', script], cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
