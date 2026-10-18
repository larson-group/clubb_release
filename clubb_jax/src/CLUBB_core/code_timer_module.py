"""Record named driver timings and optional full JAX execution traces.

Mirrors the built-in registry, nested stack and summary table in
src/CLUBB_core/code_timer_module.F90. JAX measures host wall-clock seconds and
can record a full profile through the same start/stop calls. Driver code owns
instrumentation; compiled physics functions do not call this module.

Usage:
    - timer_initialize(): reset the registry and start any requested trace
      before initializing a case. Initialization follows the same enable rule
      as Fortran; the driver calls timer_disable() after reading negative debug.
    - timer_start(name), timer_stop(name): surround a region, stopping nested
      regions before their parent. Repeated names accumulate time and calls.
      JAX's optional step_num argument identifies timesteps in full profiles.
    - timer_finalize(output_file, state=None): write the .timing summary and
      stop/export any full trace. Return its capture directory, or None when
      no trace was made. The caller presents the location and viewer commands.
    - timer_disable(): discard basic timing and disable calls until the next
      initialization; export an active trace as interrupted.
    - timer_guard(routine): decorate host driver entry points to discard basic
      timing, stop an interrupted trace and re-raise the original exception.

Understanding:
    Timing state belongs to the process and is reset for each case. Start/stop
    names must match and be nonempty. Inclusive seconds cover a region's whole
    interval; exclusive seconds subtract its nested regions. The summary also
    reports call counts, average inclusive seconds and exclusive percentages.
    The cpu_timer helper names are retained to match the Fortran source;
    Python lists replace its registry/stack allocation helpers. GPTL is not
    a JAX backend.

    JAX adaptation: perf_counter measures elapsed host wall time instead of
    CPU consumption. Basic timing adds no device synchronization: a region can
    end before its submitted work completes. Use completed process-group wall
    time from time_clubb.py for comparisons across implementations.

    CLUBB_JAX_TRACE enables optional full profiling; the JAX launcher sets it
    for -jax=profile. Each initialization uses a separate output/jax_traces
    directory unless CLUBB_JAX_TRACE_DIR overrides it. Finalization waits for
    outstanding case arrays/callbacks once and exports XProf and Perfetto data,
    plus trace_info.json with recording/export durations. The timer owns no
    viewer logic. timer_guard also handles failures before case state exists.
"""

from dataclasses import dataclass
from functools import wraps
import json
import os
from pathlib import Path
import tempfile
import time

import jax


display_name_length = 48


@dataclass
class timer_entry_t:
    """One named timer's accumulated host wall-clock results."""
    name: str
    inclusive_time: float = 0.0
    exclusive_time: float = 0.0
    call_count: int = 0


@dataclass
class timer_stack_entry_t:
    """One active timer and the time spent in its nested regions."""
    timer_id: int
    start_time: float
    child_time: float = 0.0
    trace_region: object | None = None  # JAX adaptation: optional profiler annotation.


l_timer_initialized = False
l_timer_enabled = True
timer_entries = []
timer_stack = []

_trace_dir = None
_capture_start = 0.0
_trace_options = {}


def timer_initialize():
    """Initialize the basic timing backend and any requested deep trace."""
    global l_timer_initialized, l_timer_enabled

    if l_timer_initialized:
        require_empty_stack("timer_initialize")
    if _trace_dir is not None:
        _finish_trace()
    reset_cpu_timer_state()
    l_timer_initialized = False
    l_timer_enabled = True
    _start_trace()
    l_timer_initialized = True


def timer_disable():
    """Disable timing and discard its registry; save an interrupted deep trace."""
    global l_timer_initialized, l_timer_enabled

    try:
        if _trace_dir is not None:
            _finish_trace(interrupted=True)
    finally:
        reset_cpu_timer_state()
        l_timer_initialized = False
        l_timer_enabled = False


def timer_start(timer_name, *, step_num=None):
    """Start or enter a named timer region.

    step_num is optional JAX profiling metadata identifying a model timestep.
    Ordinary timer calls retain the source's one-name interface.
    """
    if not l_timer_enabled:
        return
    if not l_timer_initialized:
        timer_initialize()
    cpu_timer_start(timer_name)
    if _trace_dir is not None:
        name = timer_name.rstrip()
        if step_num is None:
            region = jax.profiler.TraceAnnotation(name)
        else:
            region = jax.profiler.StepTraceAnnotation(name, step_num=int(step_num))
        region.__enter__()
        timer_stack[-1].trace_region = region


def timer_stop(timer_name):
    """Stop the most recent matching named timer region."""
    if not l_timer_enabled:
        return
    if not l_timer_initialized:
        timer_fail("timer_stop called before timer_initialize or timer_start.")
    cpu_timer_stop(timer_name)


def timer_finalize(output_file, *, state=None):
    """Write the wall-clock summary and finalize any enabled deep trace.

    JAX adaptations: state supplies outstanding arrays for trace completion;
    return the exported capture directory, or None when no capture was made.
    Basic timing does not wait for arrays, and neither backend synchronizes
    individual timer regions. The caller owns presentation of the capture.
    """
    global l_timer_initialized

    if not l_timer_enabled or not l_timer_initialized:
        return
    if not str(output_file).strip():
        timer_fail("Timer output file must not be empty.")
    require_empty_stack("timer_finalize")
    output_path = timer_absolute_path(output_file)
    print(f"Writing timing results to: {output_path}", flush=True)
    try:
        with output_path.open("w", encoding="utf-8") as timer_unit:
            print_cpu_timer_summary(timer_unit)
        if _trace_dir is not None:
            return _finish_trace(state=state)
    except BaseException as exc:
        try:
            timer_disable()
        except BaseException as cleanup_error:
            exc.add_note(f"JAX timer cleanup also failed: {cleanup_error}")
        raise
    finally:
        reset_cpu_timer_state()
        l_timer_initialized = False


def cpu_timer_start(timer_name):
    """Start or enter a named timer region using host wall-clock time."""
    start_time = time.perf_counter()
    timer_id = get_timer_id(timer_name, l_create=True)
    # JAX/Python adaptation: lists grow dynamically, replacing the source's
    # ensure_registry_capacity/ensure_stack_capacity allocation routines.
    timer_stack.append(timer_stack_entry_t(timer_id, start_time))
    timer_entries[timer_id].call_count += 1


def cpu_timer_stop(timer_name):
    """Stop the current timer; stops must match the most recent active start."""
    stop_time = time.perf_counter()
    timer_id = get_timer_id(timer_name, l_create=False)
    if not timer_stack:
        timer_fail(f"timer_stop called with no active timers: {timer_name.rstrip()}")
    active = timer_stack[-1]
    if active.timer_id != timer_id:
        timer_fail(
            "timer_stop nesting mismatch. Expected "
            f"{timer_entries[active.timer_id].name} but received {timer_name.rstrip()}."
        )
    elapsed = stop_time - active.start_time
    entry = timer_entries[timer_id]
    entry.inclusive_time += elapsed
    entry.exclusive_time += elapsed - active.child_time
    timer_stack.pop()
    if timer_stack:
        timer_stack[-1].child_time += elapsed
    if active.trace_region is not None:
        active.trace_region.__exit__(None, None, None)


def get_timer_id(timer_name, l_create):
    """Return a name's registry index, optionally creating a new entry."""
    trimmed_name = timer_name.rstrip()
    if not trimmed_name:
        timer_fail("Timer name must not be empty")
    for timer_id, entry in enumerate(timer_entries):
        if entry.name == trimmed_name:
            return timer_id
    if not l_create:
        timer_fail(f"Unknown timer name: {trimmed_name}")
    timer_entries.append(timer_entry_t(trimmed_name))
    return len(timer_entries) - 1


def print_cpu_timer_summary(timer_unit):
    """Print the source-format registry table with JAX wall-clock labels."""
    timer_unit.write("========================= CLUBB TIMER SUMMARY =========================\n")
    timer_unit.write("Backend: wall_clock\n")
    if not timer_entries:
        timer_unit.write("No timers were recorded.\n")
        timer_unit.write("=======================================================================\n")
        return
    total_exclusive_time = sum(entry.exclusive_time for entry in timer_entries)
    headings = ("Timer", "Calls", "Inclusive(s)", "Exclusive(s)", "Avg Inc(s)", "Excl Pct")
    timer_unit.write(
        f"{headings[0]:>48} {headings[1]:>10} {headings[2]:>14} "
        f"{headings[3]:>14} {headings[4]:>14} {headings[5]:>10}\n"
    )
    timer_unit.write(
        f"{'-' * display_name_length} {'-' * 10} "
        f"{'-' * 14} {'-' * 14} {'-' * 14} {'-' * 10}\n"
    )
    for entry in timer_entries:
        average_time = entry.inclusive_time / entry.call_count
        exclusive_percent = (
            100.0 * entry.exclusive_time / total_exclusive_time
            if total_exclusive_time > 0.0 else 0.0
        )
        timer_unit.write(
            f"{timer_display_name(entry.name)} {entry.call_count:10d} "
            f"{entry.inclusive_time:14.6f} {entry.exclusive_time:14.6f} "
            f"{average_time:14.6f} {exclusive_percent:10.2f}\n"
        )
    timer_unit.write("----------------------------------------------------------------------\n")
    timer_unit.write(f"Total exclusive timed wall_clock seconds: {total_exclusive_time:14.6f}\n")
    timer_unit.write("=======================================================================\n")


def timer_display_name(timer_name):
    """Return a fixed-width timer name for the timer summary table."""
    name = timer_name.rstrip()
    if len(name) > display_name_length:
        name = name[:display_name_length - 3] + "..."
    return name.ljust(display_name_length)


def require_empty_stack(caller_name):
    """Verify that no timer regions are active."""
    if timer_stack:
        name = timer_entries[timer_stack[-1].timer_id].name
        timer_fail(f"{caller_name.rstrip()} called while timer {name} is still active.")


def reset_cpu_timer_state():
    """Clear the timer registry and active stack."""
    timer_entries.clear()
    timer_stack.clear()


def timer_absolute_path(file_path):
    """Return an absolute output path, resolving relative paths from the cwd."""
    return Path(str(file_path).rstrip()).resolve()


def timer_fail(message):
    """Report an unrecoverable timer API error through the host error guard."""
    # JAX adaptation: preserve driver exception handling instead of ERROR STOP.
    raise RuntimeError(f"CLUBB timer error: {message}")


def _start_trace():
    """Start an opt-in capture in a separate directory for each initialization."""
    global _trace_dir, _capture_start, _trace_options
    enabled = os.environ.get("CLUBB_JAX_TRACE", "0")
    if enabled not in ("0", "1"):
        raise ValueError("CLUBB_JAX_TRACE must be 0 or 1")
    if enabled == "0":
        return

    options = jax.profiler.ProfileOptions()
    levels = {
        "host_tracer_level": ("CLUBB_JAX_TRACE_HOST_LEVEL", 2, 3),
        "python_tracer_level": ("CLUBB_JAX_TRACE_PYTHON_LEVEL", 1, 1),
    }
    configured = {}
    for name, (variable, default, maximum) in levels.items():
        try:
            level = int(os.environ.get(variable, str(default)))
        except ValueError as exc:
            raise ValueError(f"{variable} must be an integer from 0 to {maximum}") from exc
        if not 0 <= level <= maximum:
            raise ValueError(f"{variable} must be an integer from 0 to {maximum}")
        setattr(options, name, level)
        configured[name] = level

    # JAX needs the destination when capture starts, before driver output paths
    # have been resolved. Keep this experimental default beneath repository output.
    default_dir = Path(__file__).resolve().parents[3] / "output" / "jax_traces"
    output_dir = Path(os.environ.get("CLUBB_JAX_TRACE_DIR", str(default_dir))).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    trace_dir = Path(tempfile.mkdtemp(prefix="trace_", dir=output_dir)).resolve()
    jax.profiler.start_trace(
        trace_dir,
        create_perfetto_link=False,
        create_perfetto_trace=True,
        profiler_options=options,
    )
    _trace_dir = trace_dir
    _capture_start = time.perf_counter()
    _trace_options = configured


def _finish_trace(state=None, interrupted=False):
    """Export the capture and return its directory; always release profiler state."""
    global _trace_dir, _capture_start, _trace_options

    trace_dir = _trace_dir
    try:
        readiness_error = None
        try:
            if state is not None:
                jax.block_until_ready(state)
            jax.effects_barrier()
        except BaseException as exc:
            # A failed asynchronous operation must still release the collector.
            readiness_error = exc
        while timer_stack:
            active = timer_stack.pop()
            if active.trace_region is not None:
                active.trace_region.__exit__(None, None, None)
        capture_seconds = time.perf_counter() - _capture_start
        export_start = time.perf_counter()
        try:
            jax.profiler.stop_trace()
        except BaseException as exc:
            if readiness_error is not None:
                readiness_error.add_note(f"JAX trace export also failed: {exc}")
                raise readiness_error
            raise
        export_seconds = time.perf_counter() - export_start
        info = {
            "jax_version": jax.__version__,
            "backend": jax.default_backend(),
            "options": _trace_options,
            "interrupted": interrupted or readiness_error is not None,
            "capture_seconds": capture_seconds,
            "export_seconds": export_seconds,
        }
        (trace_dir / "trace_info.json").write_text(json.dumps(info, indent=2) + "\n")
        if readiness_error is not None:
            raise readiness_error
        return trace_dir
    finally:
        _trace_dir = None
        _capture_start = 0.0
        _trace_options = {}


def timer_guard(routine):
    """Preserve driver exceptions while discarding timers and stopping traces."""
    @wraps(routine)
    def guarded(*args, **kwargs):
        """Run the driver, releasing timer state if it raises."""
        try:
            return routine(*args, **kwargs)
        except BaseException as exc:
            try:
                timer_disable()
            except BaseException as cleanup_error:
                exc.add_note(f"JAX trace cleanup also failed: {cleanup_error}")
            raise
    return guarded
