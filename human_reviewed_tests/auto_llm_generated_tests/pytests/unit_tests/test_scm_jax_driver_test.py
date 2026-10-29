"""Driver-test selection reuses native SCM setup and managed JAX dispatch."""
import sys

import pytest

from clubb_jax import run_jax
from run_scripts import run_scm
from utilities.clubb_case_runner import choose_run_command


@pytest.mark.parametrize('selection', ['-jax', '-jax=cpu', '-jax=gpu'])
def test_scm_driver_test_uses_managed_launcher(monkeypatch, tmp_path, selection):
    """Require driver-test selection to reach the managed JAX module through the shared runner."""
    monkeypatch.setattr(sys, 'argv', ['run_scm.py', 'bomex', selection, '-driver_test'])
    monkeypatch.setattr(run_scm, 'resolve_output_dir', lambda value: tmp_path)
    monkeypatch.setattr(run_scm, 'create_case_namelist_from_args', lambda *args: 'case.in')
    calls = []
    class Runner:
        def __init__(self, options, **kwargs):
            """Capture the runtime options passed into the SCM runner."""
            calls.append(choose_run_command(options, report=lambda *_: None)[0])
        def __enter__(self):
            """Return the synthetic runner for the SCM context."""
            return self
        def __exit__(self, *args):
            """Accept context cleanup without launching a backend."""
            pass
        def run_and_wait(self, *args, **kwargs):
            """Record the namelist submission and report synthetic success."""
            return {'state': 'completed'}
    monkeypatch.setattr(run_scm, 'clubb_case_runner', Runner)
    assert run_scm.main() == 0
    assert calls[0][-1] == '-module=clubb_jax.src.clubb_driver_test'
    expected = [] if selection == '-jax' else ['-options=' + selection.split('=', 1)[1]]
    assert calls[0][1:-1] == expected


@pytest.mark.parametrize('arguments', [
    ['-python', '-driver_test'], ['-exe', '/model', '-driver_test'],
    ['-jax', '-python'], ['-jax', '-exe', '/model'],
])
def test_scm_rejects_conflicting_runtimes(monkeypatch, arguments):
    monkeypatch.setattr(sys, 'argv', ['run_scm.py', 'bomex', *arguments])
    with pytest.raises(SystemExit) as result:
        run_scm.main()
    assert result.value.code == 2


@pytest.mark.parametrize('module', [
    'clubb_jax.src.clubb_standalone', 'clubb_jax.src.clubb_driver_test',
])
def test_launcher_dispatches_driver_module(monkeypatch, tmp_path, module):
    python = tmp_path / 'python'
    calls = []
    monkeypatch.setattr(run_jax, '_prepare_environment', lambda *args: python)
    monkeypatch.setattr(run_jax, '_print_runtime_summary', lambda *args: None)
    monkeypatch.setattr(run_jax.os, 'chdir', lambda *args: None)
    monkeypatch.setattr(run_jax.os, 'execvpe', lambda *args: calls.append(args))
    assert run_jax.main(['-options=cpu', f'-module={module}', 'case.in']) == 0
    assert calls[0][1] == [str(python), '-m', module, 'case.in']


@pytest.mark.parametrize('options', [
    ['-module=unknown'],
    ['-module=clubb_jax.src.clubb_driver_test', '-module=clubb_jax.src.clubb_standalone'],
])
def test_launcher_rejects_unowned_or_duplicate_modules(options):
    with pytest.raises(run_jax.LauncherError, match='module'):
        run_jax.parse_launcher_args(options)
