"""Loss runner backend selection and option forwarding."""

import importlib
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
from run_scripts import run_scm


@pytest.fixture
def loss_cli(monkeypatch):
    # The existing script imports run_scm by its executable module name.
    # Supply the same script directory when testing its command chooser directly.
    script_dir = Path(run_scm.__file__).resolve().parent
    monkeypatch.syspath_prepend(str(script_dir))
    return importlib.import_module("run_scripts.run_scm_loss")


@pytest.mark.parametrize("options", [None, "cpu", "gpu"])
def test_loss_cli_selects_jax_without_a_fortran_install(loss_cli, monkeypatch, options):
    def reject_fortran_install():
        raise AssertionError("JAX loss evaluation must not resolve a Fortran install")

    monkeypatch.setattr(loss_cli, "choose_install_dir", reject_fortran_install)
    args = SimpleNamespace(
        jax=True, jax_options=options, python=False, driver_test=False,
    )
    command, cwd, environment = loss_cli.choose_run_command(args)
    assert command[0] == sys.executable
    assert Path(command[1]).name == "run_jax.py"
    assert command[-1] == "-module=clubb_jax.src.clubb_standalone_loss"
    assert cwd == loss_cli.RUN_SCRIPTS
    assert environment is None
    if options is None:
        assert not any(token.startswith("-options=") for token in command)
    else:
        assert f"-options={options}" in command


@pytest.mark.parametrize(
    "python_backend,driver_test",
    [(True, False), (False, True)],
)
def test_loss_cli_rejects_conflicting_backends(loss_cli, python_backend, driver_test):
    args = SimpleNamespace(
        jax=True, jax_options=None,
        python=python_backend, driver_test=driver_test,
    )
    with pytest.raises(SystemExit, match="-jax cannot be combined"):
        loss_cli.choose_run_command(args)
