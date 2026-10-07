"""Keep loss-consistency runs on the same selected model backend."""

import sys

import pytest

from tests import run_loss_output_consistency as harness


@pytest.mark.parametrize("jax_option", [None, "-jax", "-jax=cpu", "-jax=gpu"])
def test_both_runs_receive_selected_backend(monkeypatch, tmp_path, jax_option):
    argv = ["run_loss_output_consistency.py", "bomex"]
    if jax_option is not None:
        argv.append(jax_option)
    monkeypatch.setattr(sys, "argv", argv)
    args = harness.parse_args()
    commands = [
        harness.normal_run_command(args, "params.in", str(tmp_path), 10800, 21600),
        harness.loss_run_command(args, "params.in", str(tmp_path), ["thlm", "wp2"]),
    ]
    for command in commands:
        if jax_option is None:
            assert not any(token.startswith("-jax") for token in command)
        else:
            assert command.count(jax_option) == 1


def test_duplicate_backend_options_are_rejected(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_loss_output_consistency.py", "bomex", "-jax", "-jax=cpu",
    ])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def test_loss_comparison_requests_one_full_window(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, 'argv', [
        'run_loss_output_consistency.py', 'bomex', '-jax', '-fields', 'thlm',
    ])
    args = harness.parse_args()
    command = harness.loss_run_command(args, 'params.in', str(tmp_path), ['thlm'])
    assert command[command.index('-num_time_windows') + 1] == '1'
