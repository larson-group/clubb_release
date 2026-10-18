"""Check Fortran-derived timer accounting, Dash parsing and trace cleanup.

The stack and summary follow src/CLUBB_core/code_timer_module.F90, with
host wall-clock timing as an intentional JAX adaptation.
JAX-specific checks cover optional trace export and host exception handling.
"""

import gzip
import json

import jax
import jax.numpy as jnp
import pytest

from clubb_jax.src.CLUBB_core import code_timer_module as timer
from utilities.time_clubb import parse_timer_file


@pytest.fixture(autouse=True)
def reset_timer():
    """Release process-owned timer state before and after every check."""
    timer.timer_disable()
    yield
    timer.timer_disable()


@pytest.fixture
def capture(monkeypatch, tmp_path):
    """Enable real trace exports in a fresh directory with Python tracing disabled."""
    monkeypatch.setenv("CLUBB_JAX_TRACE", "1")
    monkeypatch.setenv("CLUBB_JAX_TRACE_DIR", str(tmp_path))
    monkeypatch.setenv("CLUBB_JAX_TRACE_PYTHON_LEVEL", "0")
    yield tmp_path
    timer.timer_disable()


def test_basic_timing_works_without_trace_or_device_barrier(monkeypatch, tmp_path):
    """Check nested wall-time accounting and Dash parsing against controlled intervals
    without device barriers.
    """
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")
    monkeypatch.setenv("CLUBB_JAX_TRACE_DIR", str(tmp_path))
    monkeypatch.setattr(jax.profiler, "start_trace", lambda *a, **k: pytest.fail("unexpected trace"))
    monkeypatch.setattr(jax, "block_until_ready", lambda *a: pytest.fail("unexpected barrier"))
    monkeypatch.setattr(jax, "effects_barrier", lambda: pytest.fail("unexpected barrier"))
    monkeypatch.setattr(timer.time, "process_time", lambda: pytest.fail("unexpected CPU clock"))
    timer.timer_initialize()
    clock = iter([0.0, 1.0, 3.0, 5.0, 6.0, 9.0])
    monkeypatch.setattr(timer.time, "perf_counter", lambda: next(clock))
    timer.timer_start("advance_clubb_to_end")
    timer.timer_start("inner")
    timer.timer_stop("inner")
    timer.timer_start("inner")
    timer.timer_stop("inner")
    timer.timer_stop("advance_clubb_to_end")
    output = tmp_path / "case.timing"
    assert timer.timer_finalize(output, state=object()) is None
    # Independent elapsed intervals: outer=9 s, children=2+1 s, self=6 s.
    outer, inner = parse_timer_file(output)
    assert (outer.calls, outer.inclusive_seconds, outer.exclusive_seconds) == (1, 9.0, 6.0)
    assert (inner.calls, inner.inclusive_seconds, inner.exclusive_seconds) == (2, 3.0, 3.0)
    assert outer.backend == "wall_clock"
    assert outer.time_basis == "wallclock"
    text = output.read_text()
    assert "Total exclusive timed wall_clock seconds:       9.000000" in text
    inner_line = next(line for line in text.splitlines() if line.startswith("inner"))
    assert inner_line.split()[4:] == ["1.500000", "33.33"]
    assert list(tmp_path.iterdir()) == [output]


def test_disable_discards_timers_until_reinitialized(monkeypatch, tmp_path):
    """Discard partial results on disable, then record only the newly initialized case."""
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")
    timer.timer_initialize()
    timer.timer_start("discarded")
    timer.timer_disable()
    timer.timer_start("disabled")
    timer.timer_stop("disabled")
    output = tmp_path / "case.timing"
    timer.timer_finalize(output)
    assert not output.exists()
    timer.timer_initialize()
    timer.timer_start("new_case")
    timer.timer_stop("new_case")
    timer.timer_finalize(output)
    assert [row.timer_name for row in parse_timer_file(output)] == ["new_case"]


def test_reinitialization_requires_balanced_stack_and_resets_counts(monkeypatch, tmp_path):
    """Reject unfinished regions at initialization/finalization and reset a completed
    registry.
    """
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")
    timer.timer_initialize()
    timer.timer_start("first_case")
    with pytest.raises(RuntimeError, match="still active"):
        timer.timer_initialize()
    with pytest.raises(RuntimeError, match="still active"):
        timer.timer_finalize(tmp_path / "first.timing")
    timer.timer_stop("first_case")
    timer.timer_initialize()
    timer.timer_start("second_case")
    timer.timer_stop("second_case")
    output = tmp_path / "second.timing"
    timer.timer_finalize(output)
    assert [row.timer_name for row in parse_timer_file(output)] == ["second_case"]


def test_source_name_validation_and_display_truncation(monkeypatch, tmp_path):
    """Reject empty/unknown names and preserve source-format truncation of long display
    names.
    """
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")
    timer.timer_initialize()
    with pytest.raises(RuntimeError, match="must not be empty"):
        timer.timer_start("   ")
    with pytest.raises(RuntimeError, match="Unknown timer name"):
        timer.timer_stop("unknown")
    name = "a_long_region_" + "x" * 50
    timer.timer_start(name + "   ")
    timer.timer_stop(name)
    output = tmp_path / "case.timing"
    timer.timer_finalize(output)
    assert parse_timer_file(output)[0].timer_name == name[:45] + "..."


def test_basic_driver_failure_discards_state_and_can_restart(monkeypatch, tmp_path):
    """Preserve a model exception, discard partial basic timing and successfully time
    another case.
    """
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")

    @timer.timer_guard
    def fail():
        """Raise a model error to exercise the guarded driver failure path."""
        timer.timer_start("failed_driver")
        raise ValueError("model failure")

    timer.timer_initialize()
    with pytest.raises(ValueError, match="model failure"):
        fail()
    output = tmp_path / "case.timing"
    timer.timer_finalize(output)
    assert not output.exists()
    timer.timer_initialize()
    timer.timer_start("recovered")
    timer.timer_stop("recovered")
    timer.timer_finalize(output)
    assert [row.timer_name for row in parse_timer_file(output)] == ["recovered"]


@pytest.mark.parametrize("variable,value", [
    ("CLUBB_JAX_TRACE", "yes"),
    ("CLUBB_JAX_TRACE_HOST_LEVEL", "4"),
    ("CLUBB_JAX_TRACE_PYTHON_LEVEL", "-1"),
    ("CLUBB_JAX_TRACE_HOST_LEVEL", "bad"),
])
def test_invalid_capture_settings_fail_before_start(capture, monkeypatch, variable, value):
    """Reject each invalid capture setting before creating an export directory."""
    monkeypatch.setenv(variable, value)
    with pytest.raises(ValueError, match=variable):
        timer.timer_initialize()
    assert not list(capture.iterdir())


def test_real_trace_contains_regions_and_both_formats(capture):
    """Export a real compiled kernel with nested regions, a timestep marker and both
    viewer formats.
    """
    timer.timer_initialize()
    timer.timer_start("outer")
    timer.timer_start("timestep", step_num=3)
    timer.timer_start("inner")
    result = jax.jit(lambda x: jnp.sin(x) + x * x)(jnp.arange(32.))
    timer.timer_stop("inner")
    timer.timer_stop("timestep")
    timer.timer_stop("outer")
    summary = capture / "case.timing"
    capture_path = timer.timer_finalize(summary, state=result)
    assert capture_path is not None
    assert capture_path.parent == capture
    assert (capture_path / "trace_info.json").is_file()
    assert {row.timer_name for row in parse_timer_file(summary)} == {"outer", "timestep", "inner"}
    paths = list(capture.rglob("perfetto_trace.json.gz"))
    assert len(paths) == 1
    with gzip.open(paths[0], "rt") as file:
        events = json.load(file)["traceEvents"]
    names = {event.get("name") for event in events}
    assert {"outer", "inner"} <= names
    steps = [event for event in events if event.get("name") == "timestep"]
    assert steps and any(str(event.get("args", {}).get("step_num")) == "3" for event in steps)
    assert list(capture.rglob("*.xplane.pb"))
    info = json.loads(next(capture.rglob("trace_info.json")).read_text())
    assert info["options"]["python_tracer_level"] == 0
    assert not info["interrupted"]
    assert info["export_seconds"] >= 0


def test_wrong_nesting_does_not_silently_close_another_region(capture):
    """Reject an out-of-order stop while keeping the active trace regions available for
    proper closure.
    """
    timer.timer_initialize()
    timer.timer_start("outer")
    timer.timer_start("inner")
    with pytest.raises(RuntimeError, match="nesting mismatch"):
        timer.timer_stop("outer")
    timer.timer_stop("inner")
    timer.timer_stop("outer")
    timer.timer_finalize(capture / "case.timing")


def test_driver_failure_closes_capture_and_allows_reinitialization(capture):
    """Mark a failed driver capture interrupted and permit another capture in the same
    process.
    """
    @timer.timer_guard
    def fail():
        """Raise a model error to exercise the guarded driver failure path."""
        timer.timer_start("failed_driver")
        raise ValueError("model failure")

    timer.timer_initialize()
    with pytest.raises(ValueError, match="model failure"):
        fail()
    first = json.loads(next(capture.rglob("trace_info.json")).read_text())
    assert first["interrupted"]
    timer.timer_initialize()
    timer.timer_finalize(capture / "case.timing")
    assert len(list(capture.rglob("trace_info.json"))) == 2


def test_failed_readiness_still_stops_collector(capture, monkeypatch):
    """Release the trace collector after an asynchronous readiness failure so capture
    can restart.
    """
    timer.timer_initialize()
    block_until_ready = jax.block_until_ready

    def fail_readiness(state):
        """Simulate an outstanding device computation that fails during finalization."""
        raise ValueError("device failure")

    monkeypatch.setattr(jax, "block_until_ready", fail_readiness)
    with pytest.raises(ValueError, match="device failure"):
        timer.timer_finalize(capture / "case.timing", state=object())
    monkeypatch.setattr(jax, "block_until_ready", block_until_ready)
    timer.timer_initialize()
    timer.timer_finalize(capture / "case.timing")
    assert len(list(capture.rglob("trace_info.json"))) == 2


def test_export_failure_does_not_replace_model_error(monkeypatch):
    """Preserve the original model error and attach an export failure as a diagnostic note."""
    def fail_export():
        """Simulate a failed trace export while the model exception is being propagated."""
        raise OSError("export failure")

    monkeypatch.setattr(timer, "timer_disable", fail_export)

    @timer.timer_guard
    def fail():
        """Raise a model error to exercise the guarded driver failure path."""
        raise ValueError("model failure")

    with pytest.raises(ValueError, match="model failure") as exc:
        fail()
    assert any("export failure" in note for note in exc.value.__notes__)


def test_initialization_failure_does_not_leave_a_live_profiler(capture):
    """Close a trace after a missing input file and permit a subsequent capture."""
    from clubb_jax.src.clubb_case_initalization import init_clubb_case

    with pytest.raises(FileNotFoundError):
        init_clubb_case(str(capture / "missing.in"))
    info = json.loads(next(capture.rglob("trace_info.json")).read_text())
    assert info["interrupted"]
    timer.timer_initialize()
    timer.timer_finalize(capture / "case.timing")
    assert len(list(capture.rglob("trace_info.json"))) == 2
