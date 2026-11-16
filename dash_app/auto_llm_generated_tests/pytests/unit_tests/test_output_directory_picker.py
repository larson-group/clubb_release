"""Check the two-sided chooser and its public selected-folder list.

Owned by clubb_dash/-dash. Browser transitions run under Node without a Plot
transport; filesystem actions use temporary folders and no model runs.
"""

from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

from dash import Dash, no_update
import pytest

from dash_app.plot_tab import output_chooser
from dash_app.plot_tab.callbacks_case import register_case_callbacks
from dash_app.plot_tab.tasks import _HANDLERS


def _walk(component):
    """Visit the chooser tree independently of component positions."""
    yield component
    children = getattr(component, 'children', None)
    for child in children if isinstance(children, (list, tuple)) else [children]:
        if child is not None:
            yield from _walk(child)


def test_both_sides_and_their_state_belong_to_one_chooser():
    """Build one component containing its divider, controls and public list store."""
    chooser = output_chooser.build_output_chooser(['/tmp/one', '/tmp/two'])
    components = list(_walk(chooser))
    by_id = {item.id: item for item in components if isinstance(getattr(item, 'id', None), str)}
    assert chooser.id == 'plots-output-chooser'
    assert by_id['plots-output-dirs'].data == ['/tmp/one', '/tmp/two']
    assert {'plots-output-catalog', 'plots-output-refresh', 'plots-active-output-list',
            'plots-available-output-list', 'plots-output-delete-confirm', 'plots-output-delete-expiry'} <= by_id.keys()
    divider = next(item for item in components if getattr(item, 'className', '') == 'plots-output-chooser-divider')
    assert divider.role == 'separator' and divider.to_plotly_json()['props']['aria-orientation'] == 'vertical'
    assert by_id['plots-output-expanded'].data is False
    toolbar = next(item for item in components if getattr(item, 'className', '') == 'plots-output-chooser-heading')
    toolbar_ids = {getattr(item, 'id', None) for item in _walk(toolbar) if isinstance(getattr(item, 'id', None), str)}
    assert {'plots-output-refresh', 'plots-show-extra-dir', 'plots-output-expand'} <= toolbar_ids


def test_plots_listens_to_selection_without_owning_chooser_callbacks():
    """Case loading reads selected folders; discovery/rendering stays in the chooser."""
    app = Dash(__name__)
    register_case_callbacks(app)
    assert 'plots-output-catalog.data' not in app.callback_map
    case_callback = next(value for key, value in app.callback_map.items() if 'plots-case-data.data' in key)
    assert 'plots-output-dirs' in {item['id'] for item in case_callback['inputs']}
    assert 'plots-output-catalog' not in {item['id'] for item in case_callback['inputs'] + case_callback['state']}
    output_chooser.register_output_chooser(app)
    assert 'plots-output-catalog.data' in app.callback_map


def test_empty_catalog_hydration_discovers_all_choices_without_pruning_selection(monkeypatch):
    """A restored selection cannot replace initial discovery with a partial catalog."""
    app = Dash(__name__)
    output_chooser.register_output_chooser(app)
    selected = ['/tmp/missing']
    records = [{'path': '/tmp/new', 'available': True}, {'path': selected[0], 'available': False}]
    scans = []
    monkeypatch.setattr(output_chooser, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    monkeypatch.setattr(output_chooser, 'discover_output_directories', lambda **kwargs: scans.append(kwargs) or records)
    result = _HANDLERS['catalog']('plots-output-chooser', 0, selected, [])
    assert scans == [{'selected_dirs': selected}]
    assert result == records and selected == ['/tmp/missing']


def test_known_selection_changes_do_not_rescan_and_extra_paths_are_inspected_once(monkeypatch):
    """Moving known folders is local; an unknown path only needs one direct inspection."""
    app = Dash(__name__)
    output_chooser.register_output_chooser(app)
    catalog = [{'path': '/tmp/one', 'label': '/tmp/one', 'available': True}]
    scans = []
    monkeypatch.setattr(output_chooser, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    monkeypatch.setattr(output_chooser, 'discover_output_directories', lambda **kwargs: scans.append(kwargs) or [])
    refresh = _HANDLERS['catalog']
    assert refresh('plots-output-chooser', 0, ['/tmp/one'], catalog) is no_update
    assert not scans
    result = refresh('plots-output-chooser', 0, ['/tmp/one', '/tmp/extra'], catalog)
    assert scans == [{'root': '/tmp/extra', 'selected_dirs': [], 'recursive': False}]
    record = next(item for item in result if item['path'] == '/tmp/extra')
    assert record['catalog_origin'] == 'external' and record['available'] is False
    assert refresh('plots-output-chooser', 0, ['/tmp/one', '/tmp/extra'], result) is no_update
    assert len(scans) == 1


def test_confirmed_delete_keeps_other_active_outputs(tmp_path, monkeypatch):
    """Deletion stays guarded and preserves every other selected folder."""
    root = tmp_path / 'output'
    folder, other = root / 'one', root / 'two'
    folder.mkdir(parents=True)
    other.mkdir()
    monkeypatch.setattr(output_chooser, 'OUTPUT_ROOT', root)
    monkeypatch.setattr(output_chooser.time, 'time', lambda: 100)
    app = Dash(__name__)
    output_chooser.register_output_chooser(app)
    remove = next(item['callback'].__wrapped__ for item in app.callback_map.values()
                  if item.get('callback') and item['callback'].__wrapped__.__name__ == 'delete_output_directory')
    context = SimpleNamespace(triggered_id={'type': 'plots-delete-output-dir', 'path': str(folder)},
                              triggered=[{'value': 0}])
    monkeypatch.setattr(output_chooser, 'callback_context', context)
    selected = [str(folder), str(other)]
    catalog = [{'path': value} for value in selected]
    assert remove([], 0, selected, catalog, None)[0] is no_update
    assert folder.is_dir()
    context.triggered = [{'value': 1}]
    armed = remove([], 0, selected, catalog, None)
    assert folder.is_dir() and armed[3]['path'] == str(folder)
    result = remove([], 0, selected, catalog, armed[3])
    assert not folder.exists() and other.is_dir()
    assert result[0] == [str(other)] and [item['path'] for item in result[1]] == [str(other)]


def test_browser_moves_folders_between_sides_without_a_plot_dependency():
    """Run real chooser rendering and transitions for known, missing and external paths."""
    node = shutil.which('node')
    if not node:
        pytest.skip('Node is needed for browser chooser transitions')
    asset = Path(__file__).resolve().parents[3] / 'assets' / '34_plot_chooser.js'
    script = r'''
const assert = require('assert'), fs = require('fs');
global.window = global;
const NO = {no_update: true};
window.dash_clientside = {no_update: NO, callback_context: {}};
global.fetch = () => {throw new Error('selection must not wait for HTTP')};
eval(fs.readFileSync(process.argv[1], 'utf8'));
const chooser = dash_clientside.plotsChooser;
const paths = {repo: '/repo', output: '/repo/output', home: '/home/test'};
const records = ['newest', 'second', 'third', 'oldest'].map(name => ({
    path: '/repo/output/' + name, label: 'output/' + name, case_names: ['arm','bomex'], case_count: 2, available: true}));
const rowPath = row => row.props.children[0].props.id.path;
const choose = (trigger, selected, extra='', value=1) => {
    dash_clientside.callback_context = {triggered_id: trigger, triggered: [{value}]};
    return chooser.choose([], 0, [], selected, extra, records, paths);
};
let selected = [], rendered = chooser.render(records, selected, null, paths);
assert.equal(rendered[0].length, 4, 'all available folders must be accessible');
assert.equal(rendered[0][0].props.children[0].props.children, 'newest · 2');
assert.equal(rendered[0][0].props.children[0].props.title, records[0].path + '\nCases: arm, bomex');
assert.deepEqual(choose({type:'plots-add-output-dir',path:records[0].path}, [], '', 0), [NO,NO,NO]);
selected = choose({type:'plots-add-output-dir',path:records[0].path}, selected)[0];
rendered = chooser.render(records, selected, null, paths);
assert.deepEqual(rendered[0].map(rowPath), records.slice(1).map(row => row.path));
assert.equal(rendered[1][0].props.className, 'plots-output-active-item');
assert.equal(rendered[1][0].props.children[0].props.children, 'newest · 2');
selected = choose({type:'plots-add-output-dir',path:records[1].path}, selected)[0];
assert.deepEqual(selected, records.slice(0,2).map(row => row.path));
assert.deepEqual(choose({type:'plots-add-output-dir',path:records[1].path}, selected), [NO,NO,NO]);
selected = choose({type:'plots-remove-output-dir',path:records[0].path}, selected)[0];
assert.deepEqual(selected, [records[1].path]);
rendered = chooser.render(records, selected, null, paths);
assert.equal(rowPath(rendered[0][0]), records[0].path, 'removing a bubble returns the folder to Available');
const confirming = chooser.render(records, selected, {path:records[0].path}, paths, true);
assert.equal(confirming[0][0].props.children[1].props.children, 'Confirm');
const missing = {...records[1], available:false, error:'temporarily unavailable'};
rendered = chooser.render([missing], selected, null, paths);
assert(rendered[1][0].props.className.includes('--unavailable'));
assert.equal(rendered[1][0].props.children[1].props.children, 'Unavailable');
rendered = chooser.render([], ['/tmp/unknown'], null, paths);
assert.equal(rendered[1][0].props.children[1].props.children, 'Loading…');
const extra = choose('plots-add-extra-dir', selected, '~/run/../data');
assert.deepEqual(extra[0], [...selected, '/home/test/data']);
assert.equal(extra[1], '');
assert.equal(choose('plots-add-extra-dir', extra[0], '~/data')[0], NO);
const external = {path:'/home/test/data', label:'/home/test/data', available:true};
assert.equal(chooser.render([external], [], null, paths, true)[0][0].props.children[1].props.disabled, true);

const moreRecords = Array.from({length:8}, (_, i) => ({...records[0], path:'/repo/output/run' + i}));
const compact = chooser.render(moreRecords, [], null, paths, false);
const expanded = chooser.render(moreRecords, [], null, paths, true);
assert.equal(compact[0].length, 5, 'compact view shows five available entries');
assert.equal(expanded[0].length, 8, 'expanded view includes the remaining entries');
assert.equal(compact[3], 'false');
assert.equal(expanded[3], 'true');
assert(expanded[2].includes('--expanded'));
assert(!compact[2].includes('--expanded'));
assert.equal(compact[0][0].props.children.length, 1, 'delete controls stay in the expanded view');
dash_clientside.callback_context = {triggered_id:'plots-output-expand',triggered:[{value:0}]};
assert.deepEqual(chooser.toggle(0,0,false,{display:'none'}), [NO,NO]);
dash_clientside.callback_context.triggered = [{value:1}];
assert.deepEqual(chooser.toggle(1,0,false,{display:'none'}), [true,NO]);
assert.deepEqual(chooser.toggle(2,0,true,{display:'flex'}), [false,{display:'none'}]);
dash_clientside.callback_context.triggered_id = 'plots-show-extra-dir';
assert.deepEqual(chooser.toggle(0,1,false,{display:'none'}), [true,{display:'flex'}]);

assert.equal(global.clubbPlotTasks, undefined, 'chooser must work without a plot runtime');
'''
    result = subprocess.run([node, '-e', script, str(asset)], capture_output=True, text=True, timeout=5)
    assert result.returncode == 0, result.stdout + result.stderr


def test_restored_external_output_gets_counts_from_worker_scan(monkeypatch):
    """Restored external folders retain their case counts after the initial worker scan."""
    app = Dash(__name__)
    output_chooser.register_output_chooser(app)
    refresh = _HANDLERS['catalog']
    monkeypatch.setattr(output_chooser, 'callback_context', SimpleNamespace(triggered_id='plots-output-dirs'))
    record = {'path': '/tmp/run', 'case_names': ['bomex'], 'case_count': 1, 'available': True}
    monkeypatch.setattr(output_chooser, 'discover_output_directories', lambda **kw: [record])
    result = refresh('plots-output-chooser', 0, ['/tmp/run'], [])
    assert result[0]['case_names'] == ['bomex']
    assert result[0]['available']
