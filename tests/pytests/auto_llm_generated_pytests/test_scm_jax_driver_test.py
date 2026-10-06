"""Driver-test selection reuses native SCM setup and managed JAX dispatch."""
import sys

import pytest

from clubb_jax import run_jax
from run_scripts import run_scm


@pytest.mark.parametrize('selection', ['-jax', '-jax=cpu', '-jax=gpu'])
def test_scm_driver_test_uses_managed_launcher(monkeypatch, tmp_path, selection):
    monkeypatch.setattr(sys, 'argv', ['run_scm.py', 'bomex', selection, '-driver_test'])
    monkeypatch.setattr(run_scm, 'resolve_output_dir', lambda value: tmp_path)
    monkeypatch.setattr(run_scm, 'create_case_namelist', lambda *args: 'case.in')
    calls = []
    monkeypatch.setattr(run_scm, 'run_case', lambda command, *args, **kwargs: calls.append(command) or 0)
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
