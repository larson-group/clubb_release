import argparse
import csv
import json
import multiprocessing
import os
import signal
import sys
import threading
import time
from pathlib import Path

import pytest

import utilities.time_clubb as time_clubb
from utilities.time_clubb import (
    BenchmarkOptions,
    TimerParseError,
    _restore_sigint_handler,
    parse_arguments,
    parse_positive_int_list,
    parse_timer_text,
    run_benchmark,
    run_process_group,
)
from utilities.timing_profiles import load_profiles, write_profile_manifest


CPU_TIMING = """
========================= CLUBB TIMER SUMMARY =========================
Backend: cpu_time
                                           Timer      Calls   Inclusive(s)   Exclusive(s)     Avg Inc(s)   Excl Pct
------------------------------------------------ ---------- -------------- -------------- -------------- ----------
                            initialize_clubb          1       0.250000       0.100000       0.250000      20.00
                        advance_clubb_to_end          4       1.500000       1.200000       0.375000      80.00
----------------------------------------------------------------------
Total exclusive timed cpu_time seconds:       1.300000
=======================================================================
"""

GPTL_TIMING = """
GPTL timing results
Stats for thread 0:
                                  Called  Recurse     Wall      max      min
  initialize_clubb                     1        -    0.200000 0.200000 0.200000
    advance_clubb_to_end               4        -     1.25e+0  1.25e+0  1.25e+0

Stats for thread 1:
                                  Called  Recurse     Wall      max      min
advance_clubb_to_end                   4        0    1.100000 1.100000 1.100000
"""


FAKE_RUN_SCM = r'''#!/usr/bin/env python3
import json
import os
import sys
import time
from pathlib import Path


def option(name):
    for index, token in enumerate(sys.argv[2:], start=2):
        if token == name:
            return sys.argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


case_name = sys.argv[1]
output_dir = Path(option("-out_dir"))
columns = int(option("-multicol"))
output_dir.mkdir(parents=True, exist_ok=True)
print(f" - using executable: {sys.executable}")
print(" - arguments:", json.dumps(sys.argv[1:]))
(output_dir / f"{case_name}.in").write_text(
    """&model
grid_type = 1
nzmax = 110
deltaz_nl = 40.0
zm_init_nl = 0.0
zm_top_nl = 400.0
time_initial = 0.0
time_final = 180.0
dt_main = 60.0
/
""",
    encoding="utf-8",
)
(output_dir / f"{case_name}_setup.txt").write_text("setup", encoding="utf-8")
(output_dir / f"{case_name}_log").write_text("model log", encoding="utf-8")
started = time.monotonic()
(output_dir / "invocation.json").write_text(
    json.dumps({"argv": sys.argv[1:], "pid": os.getpid(), "started": started}), encoding="utf-8"
)

override = option("-override") or ""
if "SLEEP=" in override:
    time.sleep(float(override.split("SLEEP=", 1)[1].split(",", 1)[0]))
else:
    time.sleep(0.15)

finished = time.monotonic()
payload = json.loads((output_dir / "invocation.json").read_text(encoding="utf-8"))
payload["finished"] = finished
(output_dir / "invocation.json").write_text(json.dumps(payload), encoding="utf-8")
(Path.cwd() / f"trace_{output_dir.name}.json").write_text(json.dumps(payload), encoding="utf-8")

if f"FAIL_COLUMNS={columns}" in override:
    raise SystemExit(9)

(output_dir / f"{case_name}.timing").write_text(
    """========================= CLUBB TIMER SUMMARY =========================
Backend: cpu_time
                                           Timer      Calls   Inclusive(s)   Exclusive(s)     Avg Inc(s)   Excl Pct
------------------------------------------------ ---------- -------------- -------------- -------------- ----------
                        advance_clubb_to_end          3       0.300000       0.250000       0.100000      90.00
----------------------------------------------------------------------
Total exclusive timed cpu_time seconds:       0.250000
=======================================================================
""",
    encoding="utf-8",
)
'''


def make_fake_run_scm(tmp_path: Path) -> Path:
    path = tmp_path / "fake_run_scm.py"
    path.write_text(FAKE_RUN_SCM, encoding="utf-8")
    return path


def make_options(tmp_path: Path, **changes) -> BenchmarkOptions:
    values = {
        "case_name": "arm",
        "processes": (1,),
        "batch_sizes": (1,),
        "warmups": 0,
        "repetitions": 1,
        "output": tmp_path / "timing",
        "name": "test run",
        "timeout": None,
        "continue_on_error": False,
    }
    values.update(changes)
    return BenchmarkOptions(**values)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def run_cancellable_group_from_ignoring_parent(
    options: BenchmarkOptions,
    group_dir: Path,
    fake_runner: Path,
    invocation_cwd: Path,
) -> None:
    """Exercise the signal state inherited by a wrapper launched by Dash."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    _restore_sigint_handler()
    try:
        run_process_group(
            options=options,
            forwarded=["-override", "SLEEP=30", "-stats", "none", "-debug", "0"],
            process_count=1,
            batch_size=1,
            group_dir=group_dir,
            run_scm_path=fake_runner,
            invocation_cwd=invocation_cwd,
        )
    except KeyboardInterrupt:
        raise SystemExit(130) from None


def process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def test_parse_positive_int_list():
    assert parse_positive_int_list("1,4,16") == (1, 4, 16)
    with pytest.raises(argparse.ArgumentTypeError):
        parse_positive_int_list("1,0")
    with pytest.raises(argparse.ArgumentTypeError):
        parse_positive_int_list("2,2")
    with pytest.raises(argparse.ArgumentTypeError):
        parse_positive_int_list("1,")


def test_parse_arguments_forwards_run_scm_options_and_adds_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options, forwarded = parse_arguments(
        [
            "arm",
            "-processes",
            "2,4",
            "-columns",
            "8,16",
            "-output",
            "results",
            "-overwrite",
            "-config",
            "paper",
            "-override",
            "C2=2",
            "-exe",
            "/tmp/clubb",
        ]
    )

    assert options.processes == (2, 4)
    assert options.batch_sizes == (8, 16)
    assert options.output == tmp_path / "results"
    assert options.overwrite
    assert forwarded == [
        "-config",
        "paper",
        "-override",
        "C2=2",
        "-exe",
        "/tmp/clubb",
        "-stats",
        "none",
        "-debug",
        "0",
    ]


@pytest.mark.parametrize("option", ["-multicol", "-batch_size", "-out_dir"])
def test_parse_arguments_rejects_managed_run_scm_options(option):
    with pytest.raises(SystemExit):
        parse_arguments(["arm", option, "2"])


def test_parse_arguments_rejects_configurable_required_timer():
    with pytest.raises(SystemExit):
        parse_arguments(["arm", "-primary_timer", "some_timer"])


def test_parse_arguments_preserves_explicit_stats_and_debug():
    _, forwarded = parse_arguments(["arm", "-stats", "custom.in", "-debug", "2"])
    assert forwarded == ["-stats", "custom.in", "-debug", "2"]


def test_parse_arguments_rejects_timing_disabled_debug_level():
    with pytest.raises(SystemExit):
        parse_arguments(["arm", "-debug", "-1"])


def test_main_restores_sigint_handler_before_running(monkeypatch):
    observed: list[object] = []
    previous = signal.signal(signal.SIGINT, signal.SIG_IGN)
    monkeypatch.setattr(time_clubb, "parse_arguments", lambda _argv: (object(), []))
    monkeypatch.setattr(
        time_clubb,
        "run_benchmark",
        lambda _options, _forwarded: observed.append(signal.getsignal(signal.SIGINT)) or 0,
    )
    try:
        assert time_clubb.main([]) == 0
        assert observed == [signal.default_int_handler]
    finally:
        signal.signal(signal.SIGINT, previous)


def test_parse_builtin_cpu_timer_summary():
    records = parse_timer_text(CPU_TIMING)

    assert [(record.timer_name, record.calls) for record in records] == [
        ("initialize_clubb", 1),
        ("advance_clubb_to_end", 4),
    ]
    assert records[1].backend == "cpu_time"
    assert records[1].time_basis == "cpu_time"
    assert records[1].inclusive_seconds == pytest.approx(1.5)
    assert records[1].exclusive_seconds == pytest.approx(1.2)


def test_parse_gptl_timer_summary_with_threads():
    records = parse_timer_text(GPTL_TIMING)

    advance = [record for record in records if record.timer_name == "advance_clubb_to_end"]
    assert [record.thread for record in advance] == [0, 1]
    assert [record.inclusive_seconds for record in advance] == pytest.approx([1.25, 1.1])
    assert all(record.backend == "gptl" for record in records)
    assert [record.exclusive_seconds for record in advance] == pytest.approx([1.25, 1.1])


def test_parse_gptl_derives_exclusive_costs_and_deduplicates_multiple_parents():
    text = """
Stats for thread 0:
                                  Called  Recurse     Wall      max      min
  parent_one                         1     -       2.0      2.0      2.0
*   shared_child                     2     -       0.5      0.3      0.2
* shared_child                       2     -       0.5      0.3      0.2
  parent_two                         1     -       1.0      1.0      1.0
Overhead sum = 0 wallclock seconds
"""

    records = parse_timer_text(text)

    assert [record.timer_name for record in records].count("shared_child") == 1
    parent = next(record for record in records if record.timer_name == "parent_one")
    child = next(record for record in records if record.timer_name == "shared_child")
    assert parent.exclusive_seconds == pytest.approx(1.5)
    assert child.exclusive_seconds == pytest.approx(0.5)


def test_parse_timer_rejects_unknown_text():
    with pytest.raises(TimerParseError):
        parse_timer_text("not a timing file")


def test_fake_runner_is_concurrent_and_writes_isolated_csv_rows(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, processes=(2,), batch_sizes=(3,))

    result = run_benchmark(
        options,
        ["-config", "paper", "-stats", "none", "-debug", "0"],
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert result == 0
    run_dir = next(options.output.iterdir())
    manifest = json.loads((run_dir / "profile.json").read_text(encoding="utf-8"))
    assert manifest["format"] == "clubb-timing-profile"
    assert manifest["format_version"] == 2
    assert manifest["status"] == "complete"
    assert manifest["benchmark"]["vertical_levels"] == 11
    assert manifest["build"]["resolved_executable"] == str(Path(sys.executable).resolve())
    batches = read_csv(run_dir / "batches.csv")
    timings = read_csv(run_dir / "timings.csv")
    assert batches[0]["batch_size"] == "3"
    assert batches[0]["total_batch_size"] == "6"
    assert {row["process"] for row in timings if row["timer"] == "advance_clubb_to_end"} == {"0", "1"}
    summaries, processes = load_profiles(options.output, [manifest["run_id"]])
    assert summaries[0]["total_columns"] == 6
    assert summaries[0]["timer_max_seconds"] == pytest.approx(0.3)
    assert summaries[0]["throughput_columns_per_second"] == pytest.approx(20.0)
    assert len(processes) == 2
    log_dir = run_dir / "logs" / "p0002_b000003"
    assert {path.name for path in log_dir.iterdir()} == {"arm.in", "arm_setup.txt", "arm.log", "arm.timing"}
    assert (log_dir / "arm.log").read_text(encoding="utf-8") == "model log"
    assert not (options.output / "profile_catalog.csv").exists()
    assert not list(run_dir.glob("**/invocation.json"))
    invocations = [json.loads(path.read_text(encoding="utf-8")) for path in tmp_path.glob("trace_process_*.json")]
    assert len(invocations) == 2
    assert max(item["started"] for item in invocations) < min(item["finished"] for item in invocations)
    assert all(item["argv"][item["argv"].index("-multicol") + 1] == "3" for item in invocations)
    assert all(item["argv"][item["argv"].index("-batch_size") + 1] == "3" for item in invocations)


def test_benchmark_overwrites_named_profile_without_creating_a_version(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, name="replace me", overwrite=True)
    old_profile = options.output / "replace_me"
    old_profile.mkdir(parents=True)
    write_profile_manifest(
        old_profile,
        {"run_id": "replace_me", "label": "replace me", "case_name": "arm"},
    )
    (old_profile / "old-result.txt").write_text("old", encoding="utf-8")

    result = run_benchmark(
        options,
        ["-stats", "none", "-debug", "0"],
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert result == 0
    assert [path.name for path in options.output.iterdir()] == ["replace_me"]
    assert not (old_profile / "old-result.txt").exists()
    manifest = json.loads((old_profile / "profile.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete"


def test_live_results_have_workload_metadata_before_benchmark_finishes(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, repetitions=3)
    outcome: dict[str, int] = {}

    def run() -> None:
        outcome["returncode"] = run_benchmark(
            options,
            ["-stats", "none", "-debug", "0"],
            run_scm_path=fake_runner,
            invocation_cwd=tmp_path,
        )

    worker = threading.Thread(target=run)
    worker.start()
    live_rows = []
    observed_while_running = False
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline and worker.is_alive():
        profiles = list(options.output.iterdir()) if options.output.exists() else []
        if profiles:
            run_id = profiles[0].name
            live_rows, _process_rows = load_profiles(options.output, [run_id])
            if live_rows:
                observed_while_running = worker.is_alive()
                break
        time.sleep(0.01)

    worker.join(timeout=5)

    assert not worker.is_alive()
    assert outcome == {"returncode": 0}
    assert observed_while_running
    assert live_rows
    assert live_rows[0]["process_count"] == 1
    assert live_rows[0]["columns_per_process"] == 1
    assert live_rows[0]["total_columns"] == 1


def test_failure_preserves_completed_and_failed_points(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, batch_sizes=(1, 2))

    result = run_benchmark(
        options,
        ["-override", "FAIL_COLUMNS=2", "-stats", "none", "-debug", "0"],
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert result == 1
    run_dir = next(options.output.iterdir())
    batches = read_csv(run_dir / "batches.csv")
    timings = read_csv(run_dir / "timings.csv")
    assert [row["status"] for row in batches] == ["complete", "partial"]
    assert [row["status"] for row in timings if row["timer"] != "__process_group_wall__"] == ["success", "failed"]
    failed = next(row for row in timings if row["status"] == "failed" and not row["timer"])
    assert failed["return_code"] == "9"
    assert (run_dir / "logs" / "p0001_b000002" / "arm.log").is_file()


def test_continue_on_error_completes_remaining_sweep_points(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, batch_sizes=(1, 2, 3), continue_on_error=True)

    result = run_benchmark(
        options,
        ["-override", "FAIL_COLUMNS=2", "-stats", "none", "-debug", "0"],
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert result == 1
    run_dir = next(options.output.iterdir())
    batches = read_csv(run_dir / "batches.csv")
    assert [(row["batch_size"], row["status"]) for row in batches] == [
        ("1", "complete"),
        ("2", "partial"),
        ("3", "complete"),
    ]


def test_warmups_are_aggregated_and_dash_ignores_them(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, warmups=1)

    result = run_benchmark(
        options,
        ["-stats", "none", "-debug", "0"],
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert result == 0
    run_dir = next(options.output.iterdir())
    timing_rows = read_csv(run_dir / "timings.csv")
    assert {row["phase"] for row in timing_rows} == {"warmup", "measured"}
    assert sum(row["timer"] == "advance_clubb_to_end" for row in timing_rows) == 2
    summaries, processes = load_profiles(options.output, [run_dir.name])
    assert len(summaries) == 2
    assert len(processes) == 1


def test_timeout_terminates_group_and_reports_status(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path, timeout=0.05)

    results, elapsed = run_process_group(
        options=options,
        forwarded=["-override", "SLEEP=5", "-stats", "none", "-debug", "0"],
        process_count=1,
        batch_size=1,
        group_dir=tmp_path / "timeout",
        run_scm_path=fake_runner,
        invocation_cwd=tmp_path,
    )

    assert elapsed < 2
    assert results[0].status == "timeout"
    assert results[0].return_code is not None


def test_sigint_from_ignoring_parent_terminates_child_group(tmp_path):
    fake_runner = make_fake_run_scm(tmp_path)
    options = make_options(tmp_path)
    group_dir = tmp_path / "cancelled"
    invocation_path = group_dir / "process_000" / "invocation.json"
    context = multiprocessing.get_context("spawn")
    wrapper = context.Process(
        target=run_cancellable_group_from_ignoring_parent,
        args=(options, group_dir, fake_runner, tmp_path),
    )
    child_pid = 0
    wrapper.start()
    try:
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline and wrapper.is_alive():
            if invocation_path.is_file():
                child_pid = int(json.loads(invocation_path.read_text(encoding="utf-8"))["pid"])
                break
            time.sleep(0.02)

        assert child_pid
        assert wrapper.pid is not None
        os.kill(wrapper.pid, signal.SIGINT)
        wrapper.join(timeout=5)

        assert not wrapper.is_alive()
        assert wrapper.exitcode == 130
        assert not process_exists(child_pid)
    finally:
        if wrapper.is_alive():
            wrapper.kill()
            wrapper.join(timeout=2)
        if child_pid and process_exists(child_pid):
            try:
                os.killpg(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
