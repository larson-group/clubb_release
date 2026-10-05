"""Synthetic comparison commands shared by ordinary and provisional SILHS checks."""

from tests import run_jax_vs_fortran_cases as harness


def run_comparison(tmp_path, config):
    commands = []

    class Supervisor:
        def run_and_log(self, cmd, _cwd, log_path):
            commands.append(cmd)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            if "-jax" in cmd:
                log_path.write_text("Completed 3 timesteps\n")
            elif cmd[1].endswith("run_scm.py"):
                log_path.write_text("iteration: 3 / 3 -- time = 180\n")
            return 0

    task = harness.TaskCtx(
        config=config,
        flag_data=harness.FlagData("default", None),
        repo_root=tmp_path,
        run_scm_args=["-max_iters", "3", "-dt_main", "60", "-debug", "0"],
        bindiff_threshold=1e-7,
        bindiff_percent_threshold=config.percent_threshold,
        run_output_root=tmp_path / "results",
        supervisor=Supervisor(),
    )
    return harness._run_case_w_flags(task), commands, task
