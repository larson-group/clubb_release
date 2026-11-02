"""Persistent Dash executors with explicit Python runtime selection.

RuntimeExecutor(runtime='python' or 'jax') provides submit(), shutdown() and
terminate_workers() for a single retained process. JAX startup goes through
clubb_jax.run_jax so backend selection, dependencies and environment remain
owned by the launcher. The broker owns executors and the workers retain their
ordinary module state between calls.

Tasks and results travel through a private inherited multiprocessing pipe.
The child imports only the submitted function's module, so numerical workers
need not install Dash or Plotly. stdout/stderr are retained for failure reports.
main() is the internal process entry point; it accepts an inherited connection
file descriptor and serves calls until EOF or shutdown. No public socket or
per-request input/output files are created.
"""
from concurrent.futures import Executor, Future
from multiprocessing import Pipe
from multiprocessing.connection import Connection
from pathlib import Path
import argparse
import os
import queue
import signal
import subprocess
import sys
import tempfile
import threading
import traceback

REPO_ROOT = Path(__file__).resolve().parents[1]


class RuntimeExecutor(Executor):
    """Serialize calls in one persistent process using the selected runtime."""

    def __init__(self, runtime='python', *, jax_options=None):
        """Select the runtime; start its process lazily on the first call."""
        if runtime not in {'python', 'jax'}:
            raise ValueError(f'Unsupported Dash worker runtime: {runtime}')
        self.runtime, self.jax_options = runtime, jax_options
        self._calls = queue.Queue()
        self._lock = threading.Lock()
        self._closed = False
        self._broken = None
        self._process = None
        self._connection = None
        self._log = tempfile.TemporaryFile()
        self._thread = threading.Thread(target=self._work, name=f'dash-{runtime}-worker', daemon=True)
        self._thread.start()

    def submit(self, fn, /, *args, **kwargs):
        """Return a Future for one call; retain the worker after completion."""
        with self._lock:
            if self._closed or self._broken:
                raise RuntimeError(self._broken or 'Dash worker is stopping')
            future = Future()
            self._calls.put((future, fn, args, kwargs))
            return future

    def _start(self):
        """Launch a worker with a private pipe in its own process group."""
        if self.runtime == 'jax':
            from clubb_jax.run_jax import runtime_arguments
            command = [sys.executable, str(REPO_ROOT / 'clubb_jax/run_jax.py'),
                       *runtime_arguments(self.jax_options), '-module=dash_app.workers']
        else:
            command = [sys.executable, '-m', 'dash_app.workers']
        parent, child = Pipe()
        try:
            with self._lock:
                if self._closed:
                    raise RuntimeError('Dash worker is stopping')
                self._process = subprocess.Popen(
                    [*command, '-connection_fd', str(child.fileno())], cwd=REPO_ROOT,
                    pass_fds=(child.fileno(),), start_new_session=True,
                    stdout=self._log, stderr=subprocess.STDOUT)
                self._connection = parent
        except BaseException:
            parent.close()
            raise
        finally:
            child.close()
        if parent.recv() != 'ready':
            raise RuntimeError('Dash worker did not acknowledge startup')

    def _work(self):
        """Run queued calls and resolve every Future, including transport failures."""
        try:
            while True:
                call = self._calls.get()
                if call is None:
                    break
                future, fn, args, kwargs = call
                if not future.set_running_or_notify_cancel():
                    continue
                try:
                    if self._process is None:
                        self._start()
                    self._connection.send((fn, args, kwargs))
                    success, result = self._connection.recv()
                except Exception as exc:
                    self._log.seek(0, os.SEEK_END)
                    self._log.seek(max(0, self._log.tell() - 8000))
                    error = RuntimeError(f'Dash {self.runtime} worker exited: {exc}\n'
                                         + self._log.read().decode(errors='replace'))
                    with self._lock:
                        self._broken = str(error)
                    future.set_exception(error)
                    self._cancel_pending(error)
                    break
                if success:
                    future.set_result(result)
                else:
                    future.set_exception(RuntimeError(result))
        finally:
            if self._connection is not None:
                try:
                    self._connection.send(None)
                except (EOFError, OSError):
                    pass
                self._connection.close()
            if self._process is not None:
                self._process.wait()

    def _cancel_pending(self, error=None):
        """Resolve queued calls during shutdown or a broken process connection."""
        while True:
            try:
                call = self._calls.get_nowait()
            except queue.Empty:
                return
            if call is not None:
                future = call[0]
                if error is None:
                    future.cancel()
                elif not future.done():
                    future.set_exception(error)

    def shutdown(self, wait=True, *, cancel_futures=False):
        """Stop after accepted calls, optionally canceling those not yet started."""
        with self._lock:
            self._closed = True
            if cancel_futures:
                self._cancel_pending()
            self._calls.put(None)
        if wait and threading.current_thread() is not self._thread:
            self._thread.join()
            self._log.close()

    def terminate_workers(self):
        """Cancel pending calls and stop the owned runtime process group."""
        self.shutdown(wait=False, cancel_futures=True)
        with self._lock:
            process = self._process
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        self.shutdown(wait=True, cancel_futures=True)


def _stop(_signum, _frame):
    """Exit through Python cleanup when the broker terminates an idle worker."""
    raise SystemExit('Dash worker stopped')


def main():
    """Serve calls on the broker's inherited pipe without importing UI packages."""
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-connection_fd', type=int, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, _stop)
    connection = Connection(args.connection_fd)
    try:
        connection.send('ready')
        while True:
            try:
                call = connection.recv()
            except EOFError:
                break
            if call is None:
                break
            fn, positional, keywords = call
            try:
                result = fn(*positional, **keywords)
            except Exception:
                connection.send((False, traceback.format_exc()))
            else:
                connection.send((True, result))
    finally:
        connection.close()


if __name__ == '__main__':
    main()
