#!/usr/bin/env python3
"""Open a saved JAX profile with XProf or Perfetto.

Use this command after a run with -jax=profile. It selects the latest completed
capture beneath output/jax_traces unless a capture directory is supplied.
XProf is installed on demand through the JAX launcher's managed uv; Perfetto
needs only Python's standard library and a browser.

Usage:
    ./clubb_jax/view_profile.py
    ./clubb_jax/view_profile.py -viewer perfetto output/jax_traces/trace_NAME
    main(argv=None): select a capture and run the requested viewer.
    print_profile_commands(capture_path): print commands for an exported capture;
        drivers pass None when profiling was disabled, which prints nothing.

Understanding:
    The model never starts a viewer. This command opens the browser and prints
    its URL for remote/headless use. XProf stays in the foreground until Ctrl+C.
    Perfetto opens a local page with an Open in Perfetto button. The page hands
    the selected trace to Perfetto in the browser; the loopback server exits
    once the file has been fetched. No JAX import or model runtime is needed.
    An explicit directory can also open an interrupted capture for diagnosis.

Reference:
    https://perfetto.dev/docs/visualization/deep-linking-to-perfetto-ui
"""

from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
import webbrowser

if __package__:
    from . import run_jax
else:
    import run_jax


REPO_ROOT = Path(__file__).resolve().parent.parent
TRACE_ROOT = REPO_ROOT / "output" / "jax_traces"

# Perfetto's documented postMessage interface avoids cross-origin requests from
# its public UI into loopback. Open its window directly on the user's click so
# the browser permits it, then fetch locally and transfer the compressed bytes.
_PERFETTO_PAGE = b"""<!doctype html>
<meta charset="utf-8">
<title>CLUBB JAX profile</title>
<h1>CLUBB JAX profile</h1>
<button id="open">Open in Perfetto</button>
<p id="status">The trace stays in your browser; it is not uploaded.</p>
<script>
document.getElementById('open').onclick = async () => {
  const origin = 'https://ui.perfetto.dev';
  const viewer = window.open(origin);
  const status = document.getElementById('status');
  if (!viewer) {
    status.textContent = 'The browser blocked the new tab. Allow it and click again.';
    return;
  }
  try {
    status.textContent = 'Loading the saved trace...';
    const response = await fetch('/perfetto_trace.json.gz');
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const buffer = await response.arrayBuffer();
    const ready = (event) => {
      if (event.source !== viewer || event.origin !== origin || event.data !== 'PONG') return;
      clearInterval(ping);
      window.removeEventListener('message', ready);
      viewer.postMessage({perfetto: {buffer, title: 'CLUBB JAX profile',
        fileName: 'perfetto_trace.json.gz'}}, origin, [buffer]);
      status.textContent = 'Profile opened in Perfetto. You can close this page.';
    };
    window.addEventListener('message', ready);
    const ping = setInterval(() => viewer.postMessage('PING', origin), 100);
  } catch (error) {
    status.textContent = `Could not load the trace: ${error}`;
  }
};
</script>
"""


def print_profile_commands(capture_path):
    """Print both viewer commands for a saved capture; do nothing for None."""
    if capture_path is None:
        return

    capture_path = Path(capture_path).resolve()
    script = str(Path(__file__).resolve())
    capture_arg = str(capture_path)
    if Path.cwd() == REPO_ROOT:
        script = "./clubb_jax/view_profile.py"
        capture_arg = os.path.relpath(capture_path, REPO_ROOT)

    print(f"JAX profile saved in: {capture_path}")
    print("View this profile with either command:")
    for viewer in ("xprof", "perfetto"):
        print("  " + shlex.join([script, "-viewer", viewer, capture_arg]))
    sys.stdout.flush()


def _find_capture(capture_path):
    """Resolve an explicit directory or the most recently exported normal run."""
    if capture_path is not None:
        capture = Path(capture_path).expanduser().resolve()
        if not capture.is_dir():
            raise ValueError(f"Profile directory does not exist: {capture}")
        return capture

    completed = [
        path for path in TRACE_ROOT.glob("*/trace_info.json")
        if not json.loads(path.read_text())["interrupted"]
    ]
    if not completed:
        raise ValueError(f"No completed JAX profiles found beneath {TRACE_ROOT}")
    return max(completed, key=lambda path: path.stat().st_mtime_ns).parent


def _view_perfetto(capture, port):
    """Serve one exported Perfetto file until the browser has received it."""
    traces = list(capture.rglob("perfetto_trace.json.gz"))
    if len(traces) != 1:
        raise ValueError(f"Expected one perfetto_trace.json.gz beneath {capture}; found {len(traces)}")
    trace = traces[0]

    class TraceHandler(BaseHTTPRequestHandler):
        """Expose the browser handoff page and its selected trace on loopback."""

        def do_GET(self):
            """Send the selected gzip file; reject other paths and directory access."""
            if self.path == "/":
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(_PERFETTO_PAGE)))
                self.end_headers()
                self.wfile.write(_PERFETTO_PAGE)
                return

            if self.path != "/perfetto_trace.json.gz":
                self.send_error(404)
                return

            with trace.open("rb") as source:
                self.send_response(200)
                self.send_header("Content-Type", "application/gzip")
                self.send_header("Content-Length", str(trace.stat().st_size))
                self.end_headers()
                shutil.copyfileobj(source, self.wfile)
            self.server.trace_delivered = True

    # Keep browser waiting outside the model process. Closing this server after
    # the transfer leaves Perfetto's loaded trace available in the browser.
    with HTTPServer(("127.0.0.1", port), TraceHandler) as server:
        server.trace_delivered = False
        url = f"http://127.0.0.1:{server.server_port}/"
        print(f"Open URL in browser: {url}", flush=True)
        webbrowser.open(url)
        while not server.trace_delivered:
            server.handle_request()
    return 0


def _view_xprof(capture, port):
    """Prepare XProf independently of JAX and keep its server in the foreground."""
    if not any(capture.rglob("*.xplane.pb")):
        raise ValueError(f"No XPlane profile found beneath {capture}")

    # XProf also runs a local worker service. Keep its port distinct for each
    # viewer, and reject occupied ports before readiness could find another app.
    worker_port = port + 1 if port < 65535 else port - 1
    for requested_port in (port, worker_port):
        with socket.socket() as probe:
            probe.bind(("127.0.0.1", requested_port))

    # Reuse the launcher's uv bootstrap and repository-local package/Python cache.
    tools_dir = Path(os.environ.get("CLUBB_JAX_TOOLS_DIR", str(REPO_ROOT / ".clubb-jax-tools")))
    if not tools_dir.is_absolute():
        tools_dir = REPO_ROOT / tools_dir
    env = os.environ | {
        "UV_CACHE_DIR": str(tools_dir / "cache"),
        "UV_PYTHON_INSTALL_DIR": str(tools_dir / "python"),
    }
    uv = run_jax._ensure_uv(tools_dir, env)
    # Keep the verified viewer dependencies separate from the model environment.
    command = [
        str(uv), "tool", "run", "--python", "3.12",
        "--from", "xprof==2.23.2", "--with", "setuptools<70",
        "xprof", "--logdir", str(capture), "--port", str(port),
        "--grpc_port", str(worker_port),
    ]

    # Wait for readiness before opening the browser; first use may install tools.
    url = f"http://127.0.0.1:{port}/"
    with subprocess.Popen(command, env=env) as process:
        try:
            deadline = time.monotonic() + 300
            while True:
                if process.poll() is not None:
                    raise RuntimeError(f"XProf exited before its viewer was ready (exit {process.returncode})")
                try:
                    with urllib.request.urlopen(url, timeout=1):
                        break
                except urllib.error.URLError:
                    if time.monotonic() >= deadline:
                        raise RuntimeError("XProf did not become ready within 300 seconds") from None
                    time.sleep(0.2)

            print(f"Open URL in browser: {url} (Ctrl+C stops XProf)", flush=True)
            webbrowser.open(url)
            return process.wait()
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def main(argv=None):
    """Parse viewer options, resolve a capture and return the viewer's exit code."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("capture", nargs="?", help="saved capture directory; default: latest completed run")
    parser.add_argument("-viewer", choices=("xprof", "perfetto"), default="xprof")
    parser.add_argument("-port", type=int, help="local server port; default: 8791 for XProf, 9001 for Perfetto")
    args = parser.parse_args(argv)
    port = args.port if args.port is not None else (8791 if args.viewer == "xprof" else 9001)
    if not 1 <= port <= 65535:
        parser.error("-port must be between 1 and 65535")

    try:
        capture = _find_capture(args.capture)
        print(f"Opening {capture} with {args.viewer}", flush=True)
        if args.viewer == "perfetto":
            return _view_perfetto(capture, port)
        return _view_xprof(capture, port)
    except KeyboardInterrupt:
        return 130
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
