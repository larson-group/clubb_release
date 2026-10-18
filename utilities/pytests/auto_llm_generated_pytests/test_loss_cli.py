"""Check loss-frontend selection through the shared case runner.

Usage:
    tests/run_pytests.sh -unit -include_generated -k test_loss_cli
    pytest runs command-selection, incompatible-option and one-shot CLI checks.

Understanding:
    Synthetic installs and processes check native, Python and JAX selection
    without running model cases. Assertions preserve opaque JAX modifiers,
    installed Python paths and native exit status 6, while rejecting persistent
    loss requests: repeated retained evaluations belong to the existing loss API.
"""

import os
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from run_scripts import run_scm_loss
from utilities import clubb_case_runner as runtime


@pytest.mark.parametrize("options", [None, "cpu", "gpu"])
def test_shared_runner_selects_jax_loss_without_a_fortran_install(monkeypatch, options):
    """Resolve JAX loss through its launcher without inspecting any Fortran install."""
    def reject_fortran_install(*args):
        """Fail if a JAX-only request tries to locate the native runtime."""
        raise AssertionError("JAX loss evaluation must not resolve a Fortran install")

    monkeypatch.setattr(runtime, "choose_install_dir", reject_fortran_install)
    args = SimpleNamespace(jax=True, jax_options=options, python=False, driver_test=False)
    command, cwd, environment = runtime.choose_run_command(args, loss=True)
    assert command == [str(Path(runtime.CLUBB_ROOT) / 'clubb_jax/run_jax.py')] + (
        [] if options is None else [f'-options={options}']
    ) + ['-module=clubb_jax.src.clubb_standalone_loss']
    assert cwd == runtime.RUN_SCRIPTS
    assert environment is None


@pytest.mark.parametrize("python_backend,driver_test", [(True, False), (False, True)])
def test_shared_runner_rejects_conflicting_loss_backends(python_backend, driver_test):
    """Reject ambiguous backend selections and unsupported JAX driver-test combinations."""
    args = SimpleNamespace(jax=True, python=python_backend, driver_test=driver_test)
    with pytest.raises(ValueError):
        runtime.choose_run_command(args, loss=True)


@pytest.mark.parametrize('driver_test', [False, True])
def test_shared_runner_selects_native_loss_frontend(tmp_path, driver_test):
    """Select the native one-shot loss executable or its driver-test variant."""
    name = 'clubb_loss_driver_test' if driver_test else 'clubb_standalone_loss'
    (tmp_path / name).touch()
    command, _, environment = runtime.choose_run_command(
        SimpleNamespace(install_dir=str(tmp_path), driver_test=driver_test), loss=True,
    )
    assert command == [str(tmp_path / name)]
    assert environment is None


@pytest.mark.parametrize('driver_test', [False, True])
def test_shared_runner_selects_python_loss_frontend(monkeypatch, tmp_path, driver_test):
    """Select the Python loss frontend and preserve its installed runtime path."""
    from utilities import setup_python_venv

    install = tmp_path / 'install'
    python_runtime = install / 'python'
    driver = (tmp_path / 'clubb_python_api/tests/clubb_loss_driver_test.py' if driver_test
              else python_runtime / 'clubb_python/clubb_standalone_loss.py')
    driver.parent.mkdir(parents=True)
    driver.touch()
    monkeypatch.setattr(runtime, 'CLUBB_ROOT', str(tmp_path))
    monkeypatch.setattr(runtime, 'python_runtime_dir_from_install', lambda _: str(python_runtime))
    monkeypatch.setattr(setup_python_venv, 'prepare_environment', lambda _: Path(sys.executable))
    monkeypatch.setenv('PYTHONPATH', 'existing-path')
    command, _, environment = runtime.choose_run_command(
        SimpleNamespace(python=True, install_dir=str(install), driver_test=driver_test),
        loss=True,
    )
    module = ('clubb_python_api.tests.clubb_loss_driver_test' if driver_test
              else 'clubb_python.clubb_standalone_loss')
    assert command == [sys.executable, '-m', module]
    assert environment['PYTHONPATH'].split(os.pathsep) == [str(python_runtime), str(tmp_path), 'existing-path']


def test_loss_frontends_cannot_claim_api_persistence(monkeypatch):
    """Reject persistent loss frontend requests before resolving a model runtime."""
    def unexpected_runtime(*args, **kwargs):
        """Fail if unsupported persistence reaches runtime setup."""
        raise AssertionError('Unsupported worker must be rejected before resolving its runtime')

    monkeypatch.setattr(runtime, 'choose_run_command', unexpected_runtime)
    with runtime.clubb_case_runner({'exe': '/unused'}, loss=True, persistent=True) as runner:
        with pytest.raises(ValueError, match='reusable evaluations use the loss API'):
            runner.init()
        assert runner.status()['pid'] is None


def test_loss_cli_submits_runtime_options_to_the_case_runner(monkeypatch, tmp_path):
    """Require the CLI to prepare its namelist and submit it as a one-shot loss run."""
    args = SimpleNamespace(disable_stats_storage=True, driver_test=False, jax=True, jax_options='cpu')
    namelist = tmp_path / 'bomex.in'
    monkeypatch.setattr(run_scm_loss, 'build_loss_case_namelist', lambda *a: (namelist, ['thlm']))
    captured = []

    class Runner:
        def __init__(self, runtime_options, **kwargs):
            """Check the one-shot loss options supplied by the CLI."""
            assert runtime_options is args
            assert kwargs == {'persistent': False, 'loss': True, 'echo': True}
        def __enter__(self):
            """Return the synthetic runner for the CLI context."""
            return self
        def __exit__(self, *args):
            """Accept context exit without creating a model process."""
            pass
        def run_and_wait(self, path, **kwargs):
            """Record the submitted namelist and report synthetic completion."""
            captured.append(path)
            return {'state': 'completed'}

    monkeypatch.setattr(run_scm_loss, 'clubb_case_runner', Runner)
    assert run_scm_loss.run_loss_case(args, 'bomex', ['thlm'], tmp_path) == 0
    assert captured == [namelist]


def test_loss_driver_test_defaults_to_one_shot_execution(tmp_path):
    """Require a driver-test executable to receive only its namelist and return native success."""
    executable = tmp_path / 'clubb_loss_driver_test'
    executable.write_text(
        f'#!{sys.executable}\n'
        'import sys\n'
        'from pathlib import Path\n'
        'Path(sys.argv[1]+".args").write_text(repr(sys.argv[1:]))\n'
        'sys.exit(6)\n'
    )
    executable.chmod(0o755)
    namelist = tmp_path / 'bomex.in'
    namelist.touch()
    with runtime.clubb_case_runner({'install_dir': str(tmp_path), 'driver_test': True},
                                    loss=True, timeout=3) as runner:
        result = runner.run_and_wait(namelist)
        assert result['state'] == 'completed'
        assert result['return_code'] == 6
    assert Path(str(namelist)+'.args').read_text() == repr([str(namelist)])
