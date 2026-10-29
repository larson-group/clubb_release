"""Select a CLUBB runtime and run complete cases from generated namelists.

Scripts use this module to choose native Fortran, the Python SCM or JAX once,
then submit cases without duplicating command or environment setup. Ordinary
standalones can keep their process alive between cases. Every request still
initializes, advances and cleans up a fresh model; JAX can retain compiled
programs in that process without retaining model allocations.

Usage:
    Run from a repository Python script with existing aggregate namelists:

        from utilities.clubb_case_runner import clubb_case_runner

        with clubb_case_runner({"jax": True, "jax_options": "cpu"}) as runner:
            for path in ("output/first/bomex.in", "output/second/arm.in"):
                result = runner.run_and_wait(path)
                if result["state"] != "completed":
                    raise RuntimeError(result["error"])

    Runtime selection and parsing:
    - extract_jax_options(argv): separate attached -jax modifiers from argv.
    - choose_install_dir(args): resolve explicit or selected/latest installs.
    - python_runtime_dir_from_install(path): require a complete F2PY runtime.
    - validate_runtime_options(args, loss=False): reject conflicting selections.
    - choose_run_command(args, loss=False): return command, cwd and environment.
    - add_runner_arguments(parser): declare the shared runtime CLI options.
    - parse_case_arguments(argv): parse SCM runtime and model options together.

    clubb_case_runner lifecycle:
    - Construct with runtime options; init(runtime_options=None) optionally prepares the
      runtime before the first case. Entering a context does not initialize it.
    - run(namelist) submits a case; wait(timeout=None) waits for its result.
      run_and_wait(namelist) combines these operations.
    - status() returns a snapshot without waiting. close() releases the process
      and log; a context manager calls it even when the caller raises.

Understanding:
    Runtime options are an argparse Namespace or a mapping using parser field
    names: for example, {"jax": True, "jax_options": "cpu"}. Selection defaults
    to the native install/selected tree if present, otherwise install/latest.
    Namelist creation belongs to input/scripts/create_case_namelist.py; run() needs
    a file path, not a case name or parameter array.

    Persistent workers exchange one absolute namelist path and a completion
    acknowledgment per request. init() waits for their readiness; run() may do
    that startup on first use, then returns after submission. Explicit commands,
    -exe, debugger/driver tests and loss frontends default to one-shot execution:
    init() resolves the runtime and each run() launches a fresh process.

    One case may be active per runner; callers must serialize submissions.
    A reader thread drains non-debugger output while the caller works or polls status.
    Logs append, and protocol markers are consumed rather than logged. wait()
    returns failed-case status; setup/submission errors raise. A wait timeout
    raises TimeoutError and leaves the case active for another wait or close.
    Failures are never retried. Replace a failed or closed runner.

    close() requests EOF shutdown for an idle worker and terminates an active
    process group. JAX setup stays in its launcher; importing this module does
    not import JAX or F2PY. loss=True chooses a one-shot loss frontend; retained
    loss allocations remain owned by the tuner's existing loss APIs.
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
from types import SimpleNamespace

CLUBB_ROOT = str(Path(__file__).resolve().parents[1])
RUN_SCRIPTS = os.path.join(CLUBB_ROOT, "run_scripts")
INSTALL_DIR = os.path.join(CLUBB_ROOT, "install")
SELECTED_INSTALL = os.path.join(INSTALL_DIR, "selected")
LATEST_INSTALL = os.path.join(INSTALL_DIR, "latest")
READY = "CLUBB_WORKER_READY"
DONE = "CLUBB_WORKER_DONE"


def extract_jax_options(argv):
    """Return argv with bare -jax tokens, its last attached value and its count.

    A bare -jax leaves the following case name untouched. Values are not parsed
    here; the JAX launcher validates their contents. Callers use the count to
    reject repeated -jax selections.
    """
    # argparse's nargs="?" would consume CASE in "-jax CASE". Extract only
    # attached values (-jax=VALUE); the JAX launcher interprets their contents.
    normalized = []
    value = None
    occurrences = 0
    for token in argv:
        option, separator, attached = token.partition("=")
        if option == "-jax":
            occurrences += 1
            value = attached if separator else None
            normalized.append("-jax")
        else:
            normalized.append(token)
    return normalized, value, occurrences


def choose_install_dir(args=None):
    """Return (install_path, source) from args.install_dir or selected/latest.

    source is explicit, selected or latest. This selects a path without checking
    its contents; command selection reports incomplete or broken installs.
    """
    install_dir = getattr(args, "install_dir", None)
    if install_dir:
        return os.path.abspath(install_dir), "explicit"
    if os.path.lexists(SELECTED_INSTALL):
        return SELECTED_INSTALL, "selected"
    return LATEST_INSTALL, "latest"


def python_runtime_dir_from_install(install_dir):
    """Return install_dir/python after checking its F2PY runtime files.

    Raise RuntimeError when the extension, backend library or Python package is
    missing. This checks file presence without importing the compiled backend.
    """
    runtime_dir = os.path.join(install_dir, "python")
    missing_reason = None
    backend_candidates = [
        os.path.join(runtime_dir, "libclubb_f2py_backend.so"),
        os.path.join(runtime_dir, "libclubb_f2py_backend.dylib"),
        os.path.join(runtime_dir, "clubb_f2py_backend.dll"),
    ]

    if not os.path.isdir(runtime_dir):
        missing_reason = "directory is missing"
    elif not glob.glob(os.path.join(runtime_dir, "clubb_f2py*.so")):
        missing_reason = "clubb_f2py extension is missing"
    elif not any(os.path.isfile(candidate) for candidate in backend_candidates):
        missing_reason = "libclubb_f2py_backend shared library is missing"
    elif not os.path.isdir(os.path.join(runtime_dir, "clubb_python")):
        missing_reason = "clubb_python package is missing"

    if missing_reason:
        raise RuntimeError(
            f"{runtime_dir} is not a complete Python runtime ({missing_reason}). "
            "Rebuild that install with ./compile.py -python or choose a Python-enabled install."
        )

    return runtime_dir


def validate_runtime_options(args, *, loss=False):
    """Raise ValueError for incompatible runtime fields on args.

    loss selects the loss-driver rules, which allow Python driver tests but
    exclude JAX driver tests. This validation does not prepare an environment.
    """
    if sum(bool(getattr(args, name, False)) for name in ("exe", "python", "jax")) > 1:
        raise ValueError("Only one of -exe, -python, or -jax may be specified.")
    if getattr(args, "driver_test", False) and getattr(args, "exe", False):
        raise ValueError("-driver_test cannot combine with -exe.")
    if getattr(args, "driver_test", False) and getattr(args, "python", False) and not loss:
        raise ValueError("-driver_test cannot combine with -python for ordinary runs.")
    if loss and getattr(args, "jax", False) and getattr(args, "driver_test", False):
        raise ValueError("-jax cannot be combined with -driver_test for loss runs.")
    if getattr(args, "gdb", False) and any(getattr(args, name, False) for name in ("python", "jax")):
        raise ValueError("-gdb can only be used with compiled executables, not -python or -jax.")


def choose_run_command(args, *, loss=False, report=print):
    """Return (command, cwd, env) for args, selecting a loss frontend if loss.

    command is an argument list; cwd is run_scripts; env is None to inherit the
    parent environment or a mapping for the installed Python frontend. Python
    selection prepares its managed environment; JAX setup happens later in its
    launcher. report receives selection messages. Invalid paths raise
    RuntimeError and conflicting selections raise ValueError.
    """
    validate_runtime_options(args, loss=loss)
    run_cwd = RUN_SCRIPTS
    run_env = None
    install_dir = None
    install_source = None
    show_install_dir = False
    f2py_runtime_dir = None

    if getattr(args, "exe", None):
        executable = os.path.abspath(args.exe)
        if not os.path.isfile(executable):
            raise RuntimeError(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]
    elif getattr(args, "python", False):
        install_dir, install_source = choose_install_dir(args)
        f2py_runtime_dir = python_runtime_dir_from_install(install_dir)
        if loss and getattr(args, "driver_test", False):
            module_name = "clubb_python_api.auto_llm_generated_tests.integration_tests.clubb_loss_driver_test"
            python_driver = os.path.join(CLUBB_ROOT, "clubb_python_api", "auto_llm_generated_tests",
                                         "integration_tests", "clubb_loss_driver_test.py")
        elif loss:
            module_name = "clubb_python.clubb_standalone_loss"
            python_driver = os.path.join(f2py_runtime_dir, "clubb_python", "clubb_standalone_loss.py")
        else:
            module_name = "clubb_python_driver.clubb_standalone"
            python_driver = os.path.join(CLUBB_ROOT, "clubb_python_driver", "clubb_standalone.py")
        if not os.path.isfile(python_driver):
            raise RuntimeError(f"Python driver not found: {python_driver} (rebuild with ./compile.py -python)")

        from utilities.setup_python_venv import prepare_environment

        python = str(prepare_environment("python"))
        executable = f"{python} -m {module_name}"
        run_cmd = [python, "-m", module_name]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        pythonpath_entries = [f2py_runtime_dir, CLUBB_ROOT]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    elif getattr(args, "jax", False):
        jax_launcher = os.path.join(CLUBB_ROOT, "clubb_jax", "run_jax.py")
        if not os.path.isfile(jax_launcher):
            raise RuntimeError(f"JAX launcher not found: {jax_launcher}")
        executable = jax_launcher
        run_cmd = [jax_launcher]
        if getattr(args, "jax_options", None) is not None:
            # Keep profile/modifier parsing in the launcher so CLI and Dash
            # share its validation and environment setup rules.
            run_cmd.append(f"-options={args.jax_options}")
        if loss:
            run_cmd.append("-module=clubb_jax.src.clubb_standalone_loss")
        elif getattr(args, 'driver_test', False):
            run_cmd.append("-module=clubb_jax.src.clubb_driver_test")
    else:
        install_dir, install_source = choose_install_dir(args)
        show_install_dir = True
        if os.path.islink(install_dir) and not os.path.exists(install_dir):
            raise RuntimeError(f"{install_dir} points to a missing install directory.")
        if not os.path.isdir(install_dir):
            raise RuntimeError(f"{install_dir} is not an install directory (did you re-compile?)")
        if loss:
            executable_name = "clubb_loss_driver_test" if getattr(args, "driver_test", False) else "clubb_standalone_loss"
        else:
            executable_name = "clubb_driver_test" if getattr(args, "driver_test", False) else "clubb_standalone"
        executable = os.path.join(install_dir, executable_name)
        if not os.path.isfile(executable):
            raise RuntimeError(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]

    if getattr(args, "gdb", False):
        gdb_path = shutil.which("gdb")
        if gdb_path is None:
            raise RuntimeError("gdb not found on PATH")
        run_cmd = [gdb_path, "--args", *run_cmd]

    if install_dir and (show_install_dir or install_source == "explicit"):
        report(f" - using install dir ({install_source}): {os.path.realpath(install_dir)}")
    if f2py_runtime_dir:
        report(f" - using Python runtime dir ({install_source}): {os.path.realpath(f2py_runtime_dir)}")
    report(f" - using executable: {executable}")
    if getattr(args, "gdb", False):
        report(f" - launching with debugger: {gdb_path}")
    return run_cmd, run_cwd, run_env


def add_runner_arguments(parser):
    """Add runtime options to an ArgumentParser or argument group in place.

    These declare ordinary SCM options. Loss scripts declare their narrower
    surface themselves and reuse validate_runtime_options/choose_run_command.
    Attached -jax values require extract_jax_options before argparse runs.
    """
    parser.add_argument("-exe", metavar="[EXECUTABLE]",
        help="CLUBB executable to use. Overrides -install_dir and selected/latest install dirs.")

    parser.add_argument("-install_dir", metavar="[DIR]",
        help="Install directory containing CLUBB executables.\nDefault: install/selected if present, otherwise install/latest.")

    parser.add_argument("-driver_test", action="store_true",
        help="Run the driver lifecycle test instead of the standalone; combines with -jax")

    parser.add_argument("-python", action="store_true",
        help="Run the Python standalone driver (python -m clubb_python_driver.clubb_standalone)")

    parser.add_argument("-jax", action="store_true",
        help=("Run through the JAX wrapper. An optional attached -jax=VALUE is "
              "forwarded unchanged; see clubb_jax/run_jax.py -launcher_help."))

    parser.add_argument(
        "-gdb",
        action="store_true",
        help="Launch the selected compiled executable with gdb."
    )


def parse_case_arguments(argv):
    """Return a Namespace for argv's case name, runtime and namelist options.

    Combine the shared runtime and model parsers, preserve attached -jax values,
    and reject incompatible options through argparse. No case is generated and
    no runtime starts here; help or malformed input exits through argparse.
    """
    from input.scripts.create_case_namelist import add_namelist_arguments, normalize_case_selection

    parser = argparse.ArgumentParser(description="Run the standalone CLUBB model",
                                     add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    add_runner_arguments(parser.add_argument_group("Runtime options"))
    add_namelist_arguments(parser.add_argument_group("Model settings"))
    parser.add_argument("case_name", nargs="?", help="Case preset or JSON override file")
    normalized, value, occurrences = extract_jax_options(argv)
    args = normalize_case_selection(parser.parse_args(normalized), parser)
    args.jax_options = value
    if occurrences > 1:
        parser.error("-jax may be specified only once.")
    try:
        validate_runtime_options(args)
    except ValueError as exc:
        parser.error(str(exc))
    if args.batch_size is not None and args.multicol is None:
        parser.error("-batch_size requires -multicol.")
    if (args.stats_tstart is None) != (args.stats_tend is None):
        parser.error("-stats_tstart and -stats_tend must be provided together.")
    return args


class clubb_case_runner:
    """One configured runtime with at most one active case.

    Use runtime_options for normal SCM/loss entry points. command/cwd/env can
    instead select an explicit program; such commands default to one-shot runs.
    A custom worker must set persistent=True and implement READY/path/DONE.
    Construction and context entry are lazy; close() owns process/log cleanup.
    """

    def __init__(self, runtime_options=None, *, command=None, cwd=None, env=None,
                 persistent=None, loss=False,
                 log_file=None, echo=False, timeout=None):
        """Store runtime options and optional logging without launching a model.

        runtime_options is a Namespace or mapping; init() can also supply it.
        persistent=None reuses ordinary standalone processes, except explicit
        commands, -exe, -gdb, driver tests and loss frontends. False launches one
        process per case. True requires a frontend implementing the protocol.
        loss selects the loss frontend and cannot combine with persistence.

        command overrides runtime selection; cwd defaults to run_scripts. With
        runtime selection, env overlays the selected/inherited environment; for
        explicit commands it is the complete subprocess environment (None means
        inherit). log_file is opened for append now; echo also prints output.
        timeout is the default startup/wait limit in seconds, or None for none.
        """
        self.options = runtime_options or SimpleNamespace()
        if isinstance(self.options, dict):
            self.options = SimpleNamespace(**self.options)
        self.command = list(command) if command is not None else None
        self.cwd = cwd or RUN_SCRIPTS
        self.env = env
        self.persistent = False if command is not None and persistent is None else persistent
        self.loss = loss
        self.echo = echo
        self.timeout = timeout
        self._process = None
        self._reader = None
        self._condition = threading.Condition()
        self._state = "new"
        self._closing = False
        self._namelist = None
        self._return_code = None
        self._error = ""
        self._log = None
        self._log_path = None
        self._finished_at = None
        if log_file is not None:
            self._set_log(log_file)

    def _set_log(self, path):
        """Open the supplied log path, creating its directory and closing the previous log."""
        with self._condition:
            if self._log is not None:
                self._log.close()
            self._log_path = Path(path)
            self._log_path.parent.mkdir(parents=True, exist_ok=True)
            self._log = self._log_path.open("a", encoding="utf-8")

    def _write(self, text):
        """Append output to the active log and optionally echo it, under the reader lock."""
        with self._condition:
            if self._log is not None:
                self._log.write(text)
                self._log.flush()
            if self.echo:
                print(text, end="", flush=True)

    def init(self, runtime_options=None):
        """Resolve the runtime and return self, starting a worker when persistent.

        runtime_options may replace construction options only before the first
        init(). Persistent startup waits for READY using the constructor timeout;
        one-shot startup only resolves the command. Repeated init() is a no-op.
        Setup failure closes resources, marks the runner failed and re-raises.
        """
        if self._state == "closed":
            raise RuntimeError("CLUBB case runner is closed")
        if self._state == "failed":
            raise RuntimeError(self._error)
        if self._state != "new":
            if runtime_options is not None:
                raise RuntimeError("Runtime options cannot change after initialization")
            return self
        if runtime_options is not None:
            self.options = (SimpleNamespace(**runtime_options)
                            if isinstance(runtime_options, dict) else runtime_options)
        try:
            # Select process reuse independently of the model's allocations.
            if self.persistent is None:
                self.persistent = not self.loss and not any(
                    getattr(self.options, name, False) for name in ("exe", "gdb", "driver_test")
                )
            if self.persistent and any(getattr(self.options, name, False)
                                       for name in ("gdb", "driver_test")):
                raise ValueError("Persistent cases require a standalone, without -gdb or -driver_test")
            if self.persistent and self.loss:
                raise ValueError('Loss frontends run once; reusable evaluations use the loss API')
            # Resolve the command once; backend-specific setup stays with its owner.
            if self.command is None:
                command, cwd, env = choose_run_command(
                    self.options, loss=self.loss,
                    report=lambda text: self._write(text + "\n"),
                )
                self.command, self.cwd = command, cwd
                self.env = (env or os.environ.copy()) | (self.env or {})
            else:
                self._write(f" - using executable: {self.command[0]}\n")
            # A worker must acknowledge startup before any case can be submitted.
            if self.persistent:
                self._launch([*self.command, "-worker"], "starting")
                with self._condition:
                    ready = self._condition.wait_for(
                        lambda: self._state != "starting", timeout=self.timeout,
                    )
                if not ready:
                    raise TimeoutError("CLUBB runtime startup timed out")
                if self._state != "ready":
                    raise RuntimeError(self._error)
            else:
                self._state = "ready"
        except BaseException as exc:
            self.close()
            self._state, self._error = "failed", str(exc)
            raise
        return self

    def _launch(self, command, state):
        """Start command in state and attach the background output/status reader."""
        # Debugger runs retain the terminal. Other children have their own
        # process group so closing a launcher also terminates its descendants.
        terminal = bool(getattr(self.options, "gdb", False))
        self._state = state
        self._process = subprocess.Popen(
            command, cwd=self.cwd, env=self.env,
            stdin=subprocess.PIPE if self.persistent else None,
            stdout=None if terminal else subprocess.PIPE,
            stderr=None if terminal else subprocess.STDOUT,
            text=True, errors="replace", bufsize=1,
            start_new_session=not terminal,
        )
        self._reader = threading.Thread(target=self._read_output, daemon=True)
        self._reader.start()

    def _read_output(self):
        """Drain the owned process pipe, updating status and waking waiting callers.

        Worker markers change readiness/completion; other lines go to the log.
        Unexpected markers or worker exit fail the runner. One-shot exit codes
        0 and 6 mean success; 6 preserves the native standalone's convention.
        """
        process = self._process
        try:
            if process.stdout is not None:
                for line in process.stdout:
                    marker = line.rstrip("\r\n")
                    with self._condition:
                        if self.persistent and marker == READY:
                            if self._state != "starting":
                                raise RuntimeError("Unexpected standalone readiness marker")
                            self._state = "ready"
                        elif self.persistent and marker == DONE:
                            if self._state != "running":
                                raise RuntimeError("Unexpected standalone completion marker")
                            self._state = "completed"
                            self._return_code = 0
                            self._finished_at = time.perf_counter()
                        else:
                            self._write(line)
                        self._condition.notify_all()
            code = process.wait()
            with self._condition:
                if not self._closing:
                    self._return_code = code
                    self._finished_at = time.perf_counter()
                    if not self.persistent and code in (0, 6):
                        self._state = "completed"
                    else:
                        self._state = "failed"
                        self._error = f"CLUBB runner exited with status {code}"
                    self._condition.notify_all()
        except BaseException as exc:
            with self._condition:
                self._state, self._error = "failed", str(exc)
                self._condition.notify_all()
        finally:
            if process.stdout is not None:
                process.stdout.close()

    def run(self, namelist, *, log_file=None):
        """Submit namelist's complete case and return self without waiting for it.

        Initialize the runtime on first use, resolve the existing file path and
        optionally switch to an append log_file. A persistent worker receives
        the path over stdin; a one-shot process receives it as an argument.
        Every request initializes and cleans up its model. Missing/invalid paths,
        an active case or a failed/closed runner raise instead of submitting.
        """
        if namelist is None:
            raise ValueError('run() requires a namelist')
        self.init()
        # Check before joining a one-shot reader: an active case must be rejected,
        # not silently waited out before launching the next request.
        with self._condition:
            if self._state not in ('ready', 'completed'):
                raise RuntimeError('CLUBB case runner already has an active case')
        path = Path(namelist).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        if '\n' in str(path) or '\r' in str(path):
            raise ValueError('Worker paths cannot contain newlines')
        if not self.persistent and self._reader is not None:
            self._reader.join()
        # Publish the new case under the same lock used by the output reader.
        with self._condition:
            if self._state not in ('ready', 'completed'):
                raise RuntimeError('CLUBB case runner already has an active case')
            if log_file is not None:
                self._set_log(log_file)
            self._namelist = str(path)
            self._return_code, self._error, self._finished_at = None, '', None
            self._state = 'running'
            try:
                if self.persistent:
                    self._process.stdin.write(str(path) + '\n')
                    self._process.stdin.flush()
                else:
                    self._launch([*self.command, str(path)], 'running')
            except BaseException as exc:
                self._state, self._error = 'failed', str(exc)
                self._condition.notify_all()
                raise
        return self

    def wait(self, timeout=None):
        """Return the submitted case's status once it stops running.

        timeout is seconds, defaulting to the constructor limit when None. A
        model failure is returned as state=failed; an elapsed limit raises
        TimeoutError without canceling the case. Calling before submission raises.
        """
        with self._condition:
            if self._state in ("new", "starting", "ready"):
                raise RuntimeError("No case has been submitted")
            done = self._condition.wait_for(
                lambda: self._state != "running",
                timeout=self.timeout if timeout is None else timeout,
            )
        if not done:
            raise TimeoutError("CLUBB case completion timed out")
        return self.status()

    def run_and_wait(self, namelist, *, log_file=None, timeout=None):
        """Submit namelist with an optional log_file and return wait(timeout)'s status.

        Setup/submission errors and timeouts raise; model failures are returned
        in the status dictionary just as for separate run() and wait() calls.
        """
        self.run(namelist, log_file=log_file)
        return self.wait(timeout)

    def status(self):
        """Return a nonblocking dictionary describing the runtime and latest case.

        state is new, starting, ready, running, completed, failed or closed.
        pid identifies the child; namelist/log_file identify its current inputs
        and log. return_code is the one-shot exit status or 0 after worker DONE.
        error describes failure. finished_at is a perf_counter timestamp in
        seconds on completion/exit, not an elapsed duration. Unset fields are None
        (error is an empty string). A completed worker can still be alive.
        """
        with self._condition:
            return {
                "state": self._state,
                "pid": self._process.pid if self._process is not None else None,
                "namelist": self._namelist,
                "return_code": self._return_code,
                "error": self._error,
                "log_file": str(self._log_path) if self._log_path is not None else None,
                "finished_at": self._finished_at,
            }

    def close(self):
        """Close the runner permanently and release its process and log.

        Idle persistent workers receive EOF. Active processes receive SIGTERM;
        shutdown escalates to SIGKILL after two seconds. Process groups include
        launcher descendants; debugger signals target the debugger itself.
        Repeated close() calls are harmless. This does not return case results.
        """
        if self._state == "closed":
            return
        process = self._process
        self._closing = True
        if process is not None:
            if self.persistent and self._state in ("ready", "completed"):
                process.stdin.close()  # EOF requests normal standalone shutdown.
            if process.poll() is None:
                if self._state not in ("ready", "completed") or not self.persistent:
                    self._signal(signal.SIGTERM)
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    self._signal(signal.SIGKILL)
                    process.wait()
            if self._reader is not None:
                self._reader.join(timeout=2)
            if process.stdin is not None and not process.stdin.closed:
                process.stdin.close()
        with self._condition:
            if self._return_code is None and process is not None:
                self._return_code = process.returncode
            self._state = "closed"
            if self._log is not None:
                self._log.close()
                self._log = None
            self._condition.notify_all()

    def _signal(self, sig):
        """Send a signal to the owned process group, or directly to an interactive debugger."""
        try:
            if getattr(self.options, "gdb", False):
                self._process.send_signal(sig)
            else:
                os.killpg(self._process.pid, sig)
        except ProcessLookupError:
            pass

    def __enter__(self):
        """Return this runner for context-managed use without starting it."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close the owned process and log on normal or exceptional context exit."""
        self.close()
