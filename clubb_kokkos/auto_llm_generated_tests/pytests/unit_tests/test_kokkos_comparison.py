"""Check executable selection, shared inputs and failure reporting without a GPU.

Actual SCM runs belong to integration_tests/run_kokkos_vs_fortran_cases.py.
These bounded checks use its command/report contracts to catch accidental
backend selection, lost options and false passes after incomplete model runs.
"""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from clubb_kokkos.auto_llm_generated_tests.integration_tests import (
    run_kokkos_vs_fortran_cases as harness,
)


def run_pair(
    tmp_path,
    *,
    steps=(3, 3),
    codes=(0, 0),
    diff_rc=0,
    missing_report=False,
    comparison=None,
    thresholds=(0, 0),
):
    """Supply completed-step logs and structured reports for one case pair."""
    commands = []
    kokkos_exe = tmp_path / "build/clubb_standalone"
    native_exe = tmp_path / "build/fortran/clubb_standalone"
    config = harness.CaseConfig(
        "warm_rain", input_case="lba", overrides={"l_local_kk": True}
    )

    class Supervisor:
        def run_and_log(self, command, _cwd, log_path):
            """Record commands, distinguish both executables and write owner-format evidence."""
            commands.append(command)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            if command[1].endswith("run_scm.py"):
                side = 0 if command[command.index("-exe") + 1] == str(kokkos_exe) else 1
                log_path.write_text(
                    "".join(
                        f"iteration: {i} / {steps[side]} -- time = {i * 60} s\n"
                        for i in range(1, steps[side] + 1)
                    )
                )
                return codes[side]
            log_path.write_text("comparison output\n")
            if not missing_report:
                report = Path(command[command.index("-result_json") + 1])
                report.write_text(
                    json.dumps(
                        {
                            "cases": {
                                "lba": {
                                    "files": [
                                        {
                                            "comparison": comparison
                                            or (
                                                "compared"
                                                if diff_rc
                                                else "byte_identical"
                                            ),
                                            "first_failing_prefix": None,
                                        }
                                    ]
                                }
                            }
                        }
                    )
                )
            return diff_rc

    task = harness.TaskCtx(
        config=config,
        flag_data=harness.FlagData("default", None),
        repo_root=tmp_path,
        run_scm_args=[
            "-max_iters",
            "3",
            "-dt_main",
            "60",
            "-stats",
            "multi_col_stats",
            "-override",
            '{"all":{"l_local_kk":false}}',
        ],
        bindiff_threshold=thresholds[0],
        bindiff_percent_threshold=thresholds[1],
        run_output_root=tmp_path / "results",
        supervisor=Supervisor(),
        kokkos_exe=kokkos_exe,
        fortran_exe=native_exe,
    )
    return harness._run_case_w_flags(task), commands, task


def test_kokkos_comparison_forwards_same_options_to_explicit_executables(tmp_path):
    """Both models receive the same variant, overrides and timing, through -exe."""
    result, commands, task = run_pair(tmp_path)
    assert result.status == "match" and result.byte_identical is True
    kokkos, native, diff = commands
    for command, executable in ((kokkos, task.kokkos_exe), (native, task.fortran_exe)):
        assert command[command.index("-exe") + 1] == str(executable)
        assert "-jax" not in command
        assert command.count("-max_iters") == command.count("-dt_main") == 1
        assert command[-1] == "lba"
        assert command[command.index("-override") + 1] == "l_local_kk=.true."
        assert command[-len(task.run_scm_args) - 5 : -5] == task.run_scm_args
    assert "-strict" in diff
    assert diff[diff.index("-threshold") + 1] == "0"
    assert diff[diff.index("-percent_threshold") + 1] == "0"


@pytest.mark.parametrize(
    "steps,codes,missing_report,expected",
    [
        ((2, 3), (0, 0), False, "step_mismatch"),
        ((3, 3), (1, 0), False, "kokkos_failed"),
        ((3, 3), (0, 1), False, "fortran_failed"),
        ((3, 3), (1, 1), False, "both_failed"),
        ((3, 3), (0, 0), True, "diff"),
    ],
)
def test_kokkos_comparison_rejects_incomplete_or_failed_evidence(
    tmp_path, steps, codes, missing_report, expected
):
    """A successful diff cannot mask a failed/incomplete model or missing report."""
    result, commands, _ = run_pair(
        tmp_path, steps=steps, codes=codes, missing_report=missing_report
    )
    assert result.status == expected
    assert sum(command[1].endswith("run_scm.py") for command in commands) == 2


def test_kokkos_comparison_preserves_numeric_failure(tmp_path):
    """The shared comparator's numeric failure remains a failing case result."""
    result, _, _ = run_pair(tmp_path, diff_rc=1)
    assert result.status == "diff" and result.bindiff_rc == 1
    assert result.byte_identical is False


@pytest.mark.parametrize(
    "thresholds,expected",
    [((0, 0), "diff"), ((1e-10, 0), "match"), ((0, 1e-10), "match")],
)
def test_kokkos_comparison_zero_limits_require_bfb(tmp_path, thresholds, expected):
    """A clipped tiny-field match cannot pass the zero-limit BFB gate.

    Both model logs finish and the shared comparator returns success, but its
    structured report records nonidentical files. Explicit diagnostic tolerances
    retain the comparator's numeric-match behavior.
    """
    result, _, _ = run_pair(
        tmp_path,
        comparison="compared",
        thresholds=thresholds,
    )
    assert result.status == expected
    assert result.byte_identical is False
    assert result.kokkos_timesteps == result.fortran_timesteps == 3


def test_kokkos_comparison_sparse_output_step_mapping(tmp_path):
    """Map a saved output at 19,980 seconds to step three from a 19,800-second start."""
    report = tmp_path / "diff.json"
    report.write_text(
        json.dumps(
            {
                "cases": {
                    "arm": {
                        "files": [
                            {
                                "first_failing_prefix": {
                                    "record": 0,
                                    "saved_records": 2,
                                    "output_time": 19980,
                                },
                            }
                        ]
                    }
                }
            }
        )
    )
    assert harness._first_failing_timestep(report, "arm", 6, (19800, 60)) == 3


def test_kokkos_comparison_rejoins_native_stdout_stderr_progress(tmp_path):
    """An observed native stderr insertion must not hide the completed final step."""
    log = tmp_path / "native.log"
    log.write_text(
        "iteration: 359 / 360 -- time = 21540.0s / 21600.0s\n"
        "iteration: 3 Program exited normally\n"
        "60 / 360 -- time = 21600.0s / 21600.0s\n"
    )
    assert harness._parse_completed_timesteps(log) == 360
    assert harness._parse_model_timing(log) == (0, 60)


@pytest.mark.parametrize(
    "arguments",
    [
        ["-exe", "other"],
        ["-jax"],
        ["-multicol", "2"],
        ["-output_dir", "other"],
        ["-bindiff_threshold", "nan"],
        ["-bindiff_percent_threshold", "-1"],
        ["-cases", "../bomex"],
        ["-cases", "/tmp/bomex"],
    ],
)
def test_kokkos_comparison_rejects_conflicting_controls(arguments):
    """Runtime selection, output ownership and limits are explicit."""
    with pytest.raises(SystemExit) as error:
        harness.parse_args(arguments)
    assert error.value.code == 2


@pytest.mark.parametrize(
    "devices,workers,expected",
    [
        (["SERIAL"], None, 3),
        (["SERIAL"], 2, 2),
        (["CUDA", "SERIAL"], None, 1),
        (["CUDA", "SERIAL"], 2, None),
    ],
)
def test_kokkos_comparison_workers_follow_build_backend(
    tmp_path, monkeypatch, devices, workers, expected
):
    """Allow CPU concurrency while a CUDA-enabled package always stays serialized."""
    (tmp_path / "kokkos_backends.json").write_text(json.dumps({"devices": devices, "kokkos_bfb_requested": False, "kokkos_bfb": False}))
    monkeypatch.setattr(harness, "default_max_workers", lambda: 3)
    arguments = ["-build_dir", str(tmp_path)]
    if workers is not None:
        arguments += ["-workers", str(workers)]
    if expected is None:
        with pytest.raises(SystemExit) as error:
            harness.parse_args(arguments)
        assert error.value.code == 2
    else:
        args = harness.parse_args(arguments)
        assert args.jobs == expected
        assert args.kokkos_backends == devices


def test_kokkos_comparison_requires_configured_backend_metadata(tmp_path, capsys):
    """A missing build manifest must report which build needs reconfiguration."""
    with pytest.raises(SystemExit) as error:
        harness.parse_args(["-build_dir", str(tmp_path)])
    assert error.value.code == 2
    message = capsys.readouterr().err
    assert str(tmp_path / "kokkos_backends.json") in message
    assert "Reconfigure" in message


def test_kokkos_comparison_rejects_duplicate_case_output_paths(tmp_path, capsys):
    """Concurrent tasks must not write the same case logs and model outputs."""
    (tmp_path / "kokkos_backends.json").write_text(json.dumps({"devices": ["SERIAL"], "kokkos_bfb_requested": False, "kokkos_bfb": False}))
    with pytest.raises(SystemExit) as error:
        harness.parse_args(["-build_dir", str(tmp_path), "-cases", "bomex", "bomex"])
    assert error.value.code == 2
    assert "share output paths" in capsys.readouterr().err


@pytest.mark.parametrize("requested,effective", [(False, False), (True, True), (True, False)])
def test_kokkos_comparison_records_requested_and_effective_bfb(tmp_path, requested, effective):
    """Report an unavailable requested BFB backend as custom, not as BFB-enabled."""
    (tmp_path / "kokkos_backends.json").write_text(json.dumps({
        "devices": ["SERIAL"], "kokkos_bfb_requested": requested, "kokkos_bfb": effective,
    }))
    args = harness.parse_args(["-build_dir", str(tmp_path)])
    assert args.kokkos_bfb_requested is requested
    assert args.kokkos_bfb is effective


@pytest.mark.parametrize("metadata", [
    {"devices": ["SERIAL"]},
    {"devices": ["SERIAL"], "kokkos_bfb_requested": True, "kokkos_bfb": "OFF"},
    {"devices": ["SERIAL"], "kokkos_bfb_requested": False, "kokkos_bfb": True},
])
def test_kokkos_comparison_rejects_unknown_or_inconsistent_bfb(tmp_path, metadata):
    """Stale or inconsistent build metadata cannot mislabel a numerical comparison."""
    (tmp_path / "kokkos_backends.json").write_text(json.dumps(metadata))
    with pytest.raises(SystemExit) as error:
        harness.parse_args(["-build_dir", str(tmp_path)])
    assert error.value.code == 2
