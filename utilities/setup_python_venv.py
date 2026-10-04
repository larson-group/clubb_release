#!/usr/bin/env python3
"""Prepare CLUBB's shared Python environment and restart Python entry points in it."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / ".clubb-python-tools"
UV_VERSION = "0.11.32"  # Match clubb_jax/run_jax.py.


def venv_dir() -> Path:
    """Resolve custom virtualenv paths relative to the repository root."""
    path = Path(os.environ.get("CLUBB_PYTHON_VENV", ROOT / ".venv-python"))
    return path if path.is_absolute() else ROOT / path


def python_version(executable: str | Path) -> tuple[int, int] | None:
    try:
        result = subprocess.run(
            [str(executable), "-c", "import sys; print(*sys.version_info[:2])"],
            check=True, capture_output=True, text=True,
        )
        major, minor = map(int, result.stdout.split())
        return major, minor
    except (OSError, ValueError, subprocess.CalledProcessError):
        return None


def ensure_uv(env: dict[str, str]) -> Path:
    local_uv = TOOLS / "bin" / "uv"
    if local_uv.is_file() and os.access(local_uv, os.X_OK):
        return local_uv
    installed = shutil.which("uv")
    if installed:
        return Path(installed)
    local_uv.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://astral.sh/uv/{UV_VERSION}/install.sh"
    print(f"Installing uv {UV_VERSION} in {local_uv.parent}", file=sys.stderr, flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": "CLUBB-python-setup/1"})
    with urllib.request.urlopen(request, timeout=30) as response:
        installer = response.read()
    install_env = dict(env, UV_UNMANAGED_INSTALL=str(local_uv.parent))
    subprocess.run(["/bin/sh"], input=installer, env=install_env, stdout=sys.stderr, check=True)
    if not local_uv.is_file():
        raise RuntimeError(f"uv installation did not create {local_uv}")
    return local_uv


def selected_python() -> str:
    explicit = os.environ.get("PYTHON")
    if explicit:
        resolved = shutil.which(explicit) or explicit
        if (version := python_version(resolved)) is None or version < (3, 12):
            raise RuntimeError(f"PYTHON must be Python 3.12 or newer: {explicit}")
        return str(resolved)
    for name in ("python3.12", "python3.13", "python3.14", "python3", "python"):
        resolved = shutil.which(name)
        if resolved and (version := python_version(resolved)) and version >= (3, 12):
            return resolved
    return "3.12"  # uv downloads a managed interpreter when needed.


def prepare_profile(profile: str, uv: Path, env: dict[str, str]) -> None:
    venv = venv_dir()
    # Every profile includes the Python/F2PY tools; Dash and plots add only
    # the packages needed by those entry points.
    requirements = [ROOT / "clubb_python_api" / "requirements.txt"]
    packages = ["numpy", "pytest", "netCDF4", "meson", "ninja", "tabulate"]
    if profile == "dash":
        requirements.append(ROOT / "dash_app" / "requirements.txt")
        packages.extend((
            "dash", "plotly", "scipy", "matplotlib", "Pillow", "mcp",
            "pydantic", "diskcache", "multiprocess", "psutil",
        ))
    elif profile == "plot":
        requirements.append(ROOT / "postprocessing" / "pyplotgen" / "requirements.txt")
        packages.extend(("matplotlib", "cycler", "seaborn", "Pillow", "fpdf", "opencv-python-headless"))
    for requirement in requirements:
        if not requirement.is_file():
            raise RuntimeError(f"Missing requirements file: {requirement}")
    venv_python = venv / "bin" / "python"
    version = python_version(venv_python)
    if version is None or version < (3, 12):
        if venv.exists() and os.environ.get("CLUBB_PYTHON_VENV"):
            raise RuntimeError(f"Custom virtualenv needs Python 3.12 or newer: {venv}")
        command = [str(uv), "venv"]
        if venv.exists():
            command.append("--clear")
        command += ["--python", selected_python(), str(venv)]
        print(f"Creating CLUBB virtualenv: {venv}", file=sys.stderr, flush=True)
        subprocess.run(command, env=env, stdout=sys.stderr, check=True)
    digest = hashlib.sha256(b"\0".join(path.read_bytes() for path in requirements)).hexdigest()
    stamp = venv / f".clubb-{profile}-requirements.sha256"
    probe = "from importlib.metadata import version\n" + "\n".join(
        f"version({name!r})" for name in packages
    )
    # A matching stamp avoids resolution; the probe catches a damaged venv.
    ready = stamp.is_file() and stamp.read_text(encoding="utf-8").strip() == digest
    if ready:
        ready = subprocess.run([str(venv_python), "-c", probe], capture_output=True).returncode == 0
    if ready:
        return
    # Setup output goes to stderr so stdio tools (including MCP) keep stdout clean.
    print(f"Installing {profile} requirements into {venv}", file=sys.stderr, flush=True)
    command = [str(uv), "pip", "install", "--python", str(venv_python)]
    for requirement in requirements:
        command.extend(("-r", str(requirement)))
    subprocess.run(command, env=env, stdout=sys.stderr, check=True)
    stamp.write_text(digest + "\n", encoding="utf-8")
    if profile != "python":
        # Installing an extra profile also satisfies the base requirements.
        base_digest = hashlib.sha256(requirements[0].read_bytes()).hexdigest()
        (venv / ".clubb-python-requirements.sha256").write_text(base_digest + "\n", encoding="utf-8")
    print(f"Ready: {venv_python}", file=sys.stderr, flush=True)


def prepare_environment(profile: str = "python") -> Path:
    """Create or refresh the shared venv when this profile needs packages."""
    if profile not in ("python", "dash", "plot"):
        raise ValueError(f"Unknown Python environment profile: {profile}")
    TOOLS.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, UV_CACHE_DIR=str(TOOLS / "cache"), UV_PYTHON_INSTALL_DIR=str(TOOLS / "python"))
    # Concurrent Jenkins stages and launchers must not modify the venv together.
    with (TOOLS / "environment.lock").open("a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        uv = ensure_uv(env)
        prepare_profile(profile, uv, env)
    return venv_dir() / "bin" / "python"


def ensure_python_venv(profile: str = "python") -> None:
    """Prepare the venv, then restart this command under its Python if needed."""
    try:
        python = prepare_environment(profile)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"Python setup failed: {exc}") from exc
    if Path(sys.prefix).resolve() == venv_dir().resolve():
        return
    # Keep `-m module` and interpreter flags when Python provides orig_argv.
    original = getattr(sys, "orig_argv", None)
    if original is not None:
        arguments = original[1:]
    else:
        module = getattr(sys.modules["__main__"], "__spec__", None)
        arguments = ["-m", module.name, *sys.argv[1:]] if module else sys.argv
    os.execv(str(python), [str(python), *arguments])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument("profile", nargs="?", choices=("python", "dash", "plot"), default="python")
    args = parser.parse_args()
    prepare_environment(args.profile)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"Python setup failed: {exc}") from exc
