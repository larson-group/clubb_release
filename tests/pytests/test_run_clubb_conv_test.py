"""Focused command-construction checks for the convergence runner."""

from types import SimpleNamespace

from tests import run_clubb_conv_test


def test_user_override_is_merged_with_convergence_overrides(monkeypatch):
    commands = []
    monkeypatch.setattr(
        run_clubb_conv_test.subprocess,
        "run",
        lambda command, **_kwargs: commands.append(command) or SimpleNamespace(returncode=0),
    )

    run_clubb_conv_test.run_scm(
        case_name="bomex",
        dt=1,
        time_initial=0.0,
        time_final=5400.0,
        l_restart=False,
        run_scm_args=["-config", "default"],
        user_override="l_diag_Lscale_from_tau=.true.",
    )

    command = commands[0]
    assert command[command.index("-out_dir") + 1] == str(run_clubb_conv_test.OUTPUT_DIR)
    assert command.count("-override") == 1
    override = command[command.index("-override") + 1]
    assert override.startswith("l_diag_Lscale_from_tau=.true.,")
    assert "time_final=5400.0" in override
    assert override.endswith("l_restart=.false.")


def test_override_is_consumed_by_convergence_cli(monkeypatch):
    monkeypatch.setattr(
        run_clubb_conv_test.sys,
        "argv",
        ["run_clubb_conv_test.py", "-p", "-override", "l_diag_Lscale_from_tau=.true."],
    )

    args, forwarded = run_clubb_conv_test.parse_args()

    assert args.plot_result is True
    assert args.override == "l_diag_Lscale_from_tau=.true."
    assert forwarded == []
