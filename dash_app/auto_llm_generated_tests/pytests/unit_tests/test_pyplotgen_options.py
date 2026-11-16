"""Check help-derived PyPlotGen controls and argument forwarding without generating plots.

The clubb_dash/-dash suite owns these checks. One bounded real -h invocation
checks discovery against the current CLI; synthetic help checks new options,
value handling, the help notecard and the existing broker launch boundary.
"""

from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace

from dash import Dash, dcc, html, no_update
import pytest

from dash_app.plot_tab import callbacks_pyplotgen as callbacks, pyplotgen_runtime as runtime, tasks
from dash_app.plot_tab.layout import _directory_case_selector
from dash_app.shared import actions


HELP = '''usage: pyplotgen.py [options]

options:
  -h, -help             Show help.
  -clubb CLUBB [CLUBB ...]
                        Input folders.
  -output_dir DIR       Export destination.
  -les                  Compare LES.
  -pdf                  Generate a PDF.
  -pdf_filesize_limit PDF_FILESIZE_LIMIT
                        Maximum PDF size.
  -future_choice {one,two}
                        A new option Dash has never seen.
                        With a wrapped description.
  -movies [{mp4,avi}]    Optional movie format.
  -cases CASES [CASES ...]
                        Cases to plot.
'''


@pytest.fixture
def help_data(monkeypatch):
    """Supply an independent CLI help sample without launching PyPlotGen."""
    data = {'text': HELP, 'options': runtime._parse_help_options(HELP)}
    monkeypatch.setattr(runtime, 'pyplotgen_help', lambda: data)
    return data


def _walk(component):
    """Visit components in a header or notecard without depending on child positions."""
    yield component
    children = getattr(component, 'children', None)
    for child in children if isinstance(children, (list, tuple)) else [children]:
        if child is not None:
            yield from _walk(child)


def _callbacks():
    """Expose the real callback bodies while retaining their registered dependencies."""
    app = Dash(__name__)
    callbacks.register_pyplotgen_callbacks(app)
    return {item['callback'].__wrapped__.__name__: item['callback'].__wrapped__
            for item in app.callback_map.values() if 'callback' in item}


def test_actual_help_is_discovered_without_generating_a_gallery(tmp_path, monkeypatch):
    """Read the current CLI's actual help and derive the PDF options from it."""
    monkeypatch.setattr(runtime, 'PYPLOTGEN_OUTPUT_ROOT', tmp_path)
    runtime._read_help.cache_clear()
    try:
        data = runtime.pyplotgen_help()
        options = {item['flag']: item for item in data['options']}
        assert 'usage:' in data['text'] and '-pdf' in data['text']
        assert options['-pdf']['metavar'] == ''
        assert options['-pdf_filesize_limit']['metavar'] == 'PDF_FILESIZE_LIMIT'
        assert '-clubb' not in options and '-output_dir' not in options and '-h' not in options
        assert not list(tmp_path.rglob('index.html'))
    finally:
        runtime._read_help.cache_clear()


def test_unseen_options_and_required_optional_and_multiple_values(help_data):
    """Flags and quoted values come from help, including an unknown-to-Dash option."""
    assert runtime.pyplotgen_arguments(['-pdf', '-future_choice', '-cases', '-movies'], {
        '-future_choice': 'two', '-cases': 'bomex arm', '-movies': '',
    }) == ['-pdf', '-future_choice', 'two', '-cases', 'bomex', 'arm', '-movies']
    assert runtime.pyplotgen_arguments(['-cases'], {'-cases': '"case with spaces"'}) == ['-cases', 'case with spaces']
    assert runtime.pyplotgen_arguments([]) == []
    assert runtime.pyplotgen_arguments(None) == ['-les']


@pytest.mark.parametrize('selected,values', [
    (['-not_in_help'], {}), (['-clubb'], {}), (['-output_dir'], {}),
    (['-pdf_filesize_limit'], {}), (['-future_choice'], {'-future_choice': 'one two'}),
    (['-cases'], {'-cases': '-output_dir /elsewhere'}),
])
def test_invalid_or_dash_managed_options_are_rejected_before_launch(help_data, selected, values):
    """Keep managed folders authoritative and reject unusable selections/value fields."""
    with pytest.raises(ValueError):
        runtime.pyplotgen_arguments(selected, values)


def test_help_command_is_cached_until_the_entry_point_changes(tmp_path, monkeypatch):
    """Avoid repeating CLI imports for an unchanged script; edits refresh the catalog."""
    script = tmp_path / 'pyplotgen.py'
    script.write_text('original')
    monkeypatch.setattr(runtime, 'PYPLOTGEN_SCRIPT', script)
    monkeypatch.setattr(runtime, 'PYPLOTGEN_OUTPUT_ROOT', tmp_path)
    calls = []

    def run(command, **kwargs):
        """Record help-only subprocess arguments without executing a script."""
        calls.append(command)
        return SimpleNamespace(stdout=HELP)

    monkeypatch.setattr(runtime.subprocess, 'run', run)
    runtime._read_help.cache_clear()
    try:
        assert runtime.pyplotgen_help() == runtime.pyplotgen_help()
        assert calls == [[runtime.sys.executable, str(script), '-h']]
        script.write_text('updated entry point')
        runtime.pyplotgen_help()
        assert len(calls) == 2
    finally:
        runtime._read_help.cache_clear()


def test_header_separates_export_and_supports_searchable_multiple_selection():
    """The export stays in the top header after its own vertical divider."""
    header = _directory_case_selector({'case_data': {}, 'enabled_benchmarks': []})
    sections = header.children[0].children
    export = sections[-1]
    assert export.className == 'plots-pyplotgen-control'
    assert sections[-2].className == 'plots-output-benchmark-divider'
    selector = next(item for item in _walk(export) if isinstance(item, dcc.Dropdown))
    assert selector.multi and selector.searchable and selector.closeOnSelect is False
    assert any(getattr(item, 'id', None) == 'plots-pyplotgen-help-open' for item in _walk(export))


def test_controls_help_and_values_follow_cli_text(help_data, monkeypatch):
    """Render a new CLI option/value field and preserve the verbatim help in its notecard."""
    functions = _callbacks()
    options, error = functions['update_pyplotgen_options'](help_data)
    assert not error
    assert next(item for item in options if item['value'] == '-future_choice')['title'] == (
        'A new option Dash has never seen. With a wrapped description.')
    fields = functions['update_pyplotgen_value_fields'](['-pdf', '-future_choice'], help_data,
                                                      ['two'], [{'flag': '-future_choice'}])
    value = next(item for item in _walk(html.Div(fields)) if isinstance(item, dcc.Input))
    assert value.value == 'two' and value.placeholder == '{one,two}'
    monkeypatch.setattr(callbacks, 'callback_context', SimpleNamespace(triggered_id='plots-pyplotgen-help-open'))
    dialog = functions['show_pyplotgen_help'](1, [], help_data, None)
    assert next(item for item in _walk(dialog) if isinstance(item, html.Pre)).children == HELP
    monkeypatch.setattr(callbacks, 'callback_context', SimpleNamespace(
        triggered_id={'type': 'plots-pyplotgen-help-close', 'index': 0}))
    assert functions['show_pyplotgen_help'](1, [0], help_data, dialog) is no_update
    assert functions['show_pyplotgen_help'](1, [1], help_data, dialog) == ''
    monkeypatch.setattr(callbacks, 'callback_context', SimpleNamespace(triggered_id='plots-pyplotgen-help'))
    assert functions['show_pyplotgen_help'](1, [], help_data, '') is no_update


def test_selected_options_reach_the_broker_owned_literal_command(tmp_path, monkeypatch, help_data):
    """Follow UI selection through the broker owner into argv without starting an export."""
    source = tmp_path / 'model output'
    source.mkdir()
    monkeypatch.setattr(runtime, 'PYPLOTGEN_OUTPUT_ROOT', tmp_path / 'galleries')
    calls = []

    def popen(command, **kwargs):
        """Capture argv and return a process stub without touching model data."""
        calls.append(command)
        return SimpleNamespace(pid=1234)

    monkeypatch.setattr(runtime.subprocess, 'Popen', popen)
    monkeypatch.setattr(actions, 'broker_jobs', lambda: {})
    monkeypatch.setattr(actions, 'set_broker_job', lambda *_a, **_k: None)
    monkeypatch.setattr(actions, 'publish_event', lambda *_a, **_k: None)
    monkeypatch.setattr(actions, '_background', lambda *_a, **_k: None)
    process = []
    start = runtime.start_pyplotgen

    def start_and_capture(*args, **kwargs):
        """Retain the fake process so its real log handle is closed after the check."""
        job, handle = start(*args, **kwargs)
        process.append(handle)
        return job, handle

    monkeypatch.setattr(actions, 'start_pyplotgen', start_and_capture)
    monkeypatch.setattr(callbacks, 'perform_action', lambda _name, payload, **_k: actions.launch_pyplotgen_request(**payload))
    functions = _callbacks()
    try:
        result = functions['launch_pyplotgen'](1, [str(source)], ['-pdf', '-pdf_filesize_limit', '-future_choice'],
                                              ['8', 'two'], [{'flag': '-pdf_filesize_limit'}, {'flag': '-future_choice'}], {}, {})
        assert result['kind'] == 'started'
        assert calls[0][3:8] == ['-pdf', '-pdf_filesize_limit', '8', '-future_choice', 'two']
        assert calls[0][8:10] == ['-clubb', str(source)]
        assert calls[0][10] == '-output_dir'
        assert Path(calls[0][11]).parent == tmp_path / 'galleries'
    finally:
        for handle in process:
            runtime.release_pyplotgen(handle)


def test_case_selection_cannot_discard_help_pending_for_the_whole_page():
    """Help survives a case generation change while sharing the metadata queue."""
    calls = []

    class Executor:
        """Hold tasks pending while the simulated browser changes case."""
        def submit(self, _function, request):
            """Return an unresolved future for the test's chosen completion order."""
            future = Future()
            calls.append((request, future))
            return future

        def shutdown(self, **_kwargs):
            """The test uses no child process."""

    queue = tasks.PlotTasks(Executor)
    try:
        queue.submit(dict(scope='page:pyplotgen_help', task='pyplotgen_help', revision=1, selection=0, args=[]))
        queue.submit(dict(scope='page:case', task='case', revision=1, selection=1, args=[]))
        calls[0][1].set_result('{"text":"help"}')
        assert queue.poll('page:pyplotgen_help', 1)['result'] == {'text': 'help'}
        calls[1][1].set_result('{}')
    finally:
        queue.close()
