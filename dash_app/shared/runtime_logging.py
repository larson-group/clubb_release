"""Private, cross-session logging for the local Dash launcher."""

from __future__ import annotations

import argparse
import fcntl
import os
import sys
from pathlib import Path

from .runtime import private_path

APP_LOG_NAME = "dash.log"
APP_LOG_MAX_BYTES = 5 * 1024 * 1024
APP_LOG_BACKUPS = 3


def dashboard_log_path(repo_root: Path | str) -> Path:
    return private_path(repo_root, APP_LOG_NAME)


def _rotate_locked(path: Path) -> None:
    if not path.exists() or path.stat().st_size < APP_LOG_MAX_BYTES:
        return
    for index in range(APP_LOG_BACKUPS - 1, 0, -1):
        source = path.with_name(f"{path.name}.{index}")
        target = path.with_name(f"{path.name}.{index + 1}")
        if source.exists():
            source.replace(target)
    path.replace(path.with_name(f"{path.name}.1"))


def _write_chunk(path: Path, chunk: bytes) -> None:
    if not chunk:
        return
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    path.parent.chmod(0o700)
    lock_path = path.with_name(f".{path.name}.lock")
    with lock_path.open("a+b") as lock_file:
        os.chmod(lock_path, 0o600)
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        _rotate_locked(path)
        with path.open("ab") as log_file:
            os.chmod(path, 0o600)
            log_file.write(chunk)
            log_file.flush()
            os.fsync(log_file.fileno())
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def relay(path: Path, stream_name: str) -> None:
    terminal = sys.stdout.buffer if stream_name == "stdout" else sys.stderr.buffer
    while True:
        # BufferedReader.read(size) may wait for the entire 64 KiB request,
        # hiding low-volume tracebacks until process shutdown. os.read returns
        # as soon as any pipe data is available, keeping terminal and file live.
        chunk = os.read(sys.stdin.fileno(), 64 * 1024)
        if not chunk:
            return
        terminal.write(chunk)
        terminal.flush()
        _write_chunk(path, chunk)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("path", "prepare", "relay"))
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--path", type=Path)
    parser.add_argument("--stream", choices=("stdout", "stderr"))
    args = parser.parse_args()
    path = args.path or dashboard_log_path(args.repo)
    if args.command == "path":
        print(path)
    elif args.command == "prepare":
        path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        path.parent.chmod(0o700)
        with path.open("ab"):
            os.chmod(path, 0o600)
        print(path)
    else:
        if not args.stream:
            parser.error("--stream is required for relay")
        relay(path, args.stream)


if __name__ == "__main__":
    main()
