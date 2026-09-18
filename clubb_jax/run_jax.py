#!/usr/bin/env python3
"""Prepare a repository-local JAX runtime and launch CLUBB-JAX.

This file intentionally uses only the Python standard library.  It may run
under the host Python, install a supported Python with ``uv`` when necessary,
create the selected runtime environment, and finally replace itself with the
managed interpreter so signals and exit status pass through unchanged.
"""

from __future__ import annotations

import fcntl
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Sequence


UV_VERSION = "0.11.32"
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
RUNTIME_INFO = SCRIPT_DIR / "runtime_info.py"


class LauncherError(RuntimeError):
    """A user-facing launcher configuration or setup error."""


def fail(message: str) -> None:
    raise LauncherError(message)


def usage() -> str:
    return """Usage: clubb_jax/run_jax.py [--profile=cpu|gpu] [namelist-path] [driver-args...]
       clubb_jax/run_jax.py [--profile=cpu|gpu] --init_env
       clubb_jax/run_jax.py [--profile=cpu|gpu] --info[=json]

Creates or updates the JAX environment, then runs the CLUBB-JAX standalone.
Missing copies of uv and a supported Python are downloaded automatically.

Options:
  --options=VALUE        Parse the value forwarded by run_scm.py -jax=VALUE:
                         cpu or gpu, optionally followed by ,xla_prealloc
  --profile=cpu|gpu      Select CPU or the host-native GPU backend
  --accelerator=VALUE    Select an explicit backend: cpu, cuda13, or metal
  --xla-prealloc         Enable CUDA memory preallocation (CUDA only)
  --init_env             Prepare the environment without running a case
  --info[=json]          Inspect hardware and runtime readiness without setup
  --launcher-help        Show this help

Environment:
  CLUBB_JAX_ACCELERATOR  Backend used when no option is given: cpu, cuda13, metal
  CLUBB_JAX_VENV         Override the profile virtualenv path
  CLUBB_JAX_TOOLS_DIR    Managed uv/Python path (default: .clubb-jax-tools)
  PYTHON                 Python used when creating a new virtualenv
  XLA_PYTHON_CLIENT_PREALLOCATE  CUDA preallocation (default: false)
"""


def _split_options(value: str) -> tuple[str, bool]:
    if not value:
        fail("--options requires a profile; supported profiles are cpu and gpu")
    if "\n" in value or "\r" in value:
        fail("JAX options must be on one line")
    pieces = value.split(",")
    if any(not piece for piece in pieces):
        fail("Empty JAX option; expected cpu or gpu, optionally followed by ,xla_prealloc")
    profile = pieces[0].lower()
    xla_prealloc = False
    for modifier in pieces[1:]:
        if modifier.lower() != "xla_prealloc":
            fail(f"Unknown JAX option: {modifier}; supported modifier is xla_prealloc")
        if xla_prealloc:
            fail("xla_prealloc may be specified only once")
        xla_prealloc = True
    return profile, xla_prealloc


def parse_launcher_args(argv: Sequence[str]) -> tuple[dict[str, object], list[str]]:
    """Consume launcher options until the first driver argument."""
    values: dict[str, object] = {
        "profile": None,
        "accelerator": None,
        "xla_prealloc": False,
        "init_env": False,
        "info_format": None,
        "help": False,
    }
    options_seen = False
    prealloc_seen = False
    index = 0
    while index < len(argv):
        token = argv[index]
        if token.startswith("--options="):
            if options_seen:
                fail("--options may be specified only once")
            if values["profile"] is not None or values["accelerator"] is not None:
                fail("--options cannot be combined with --profile or --accelerator")
            options_seen = True
            profile, attached_prealloc = _split_options(token.split("=", 1)[1])
            values["profile"] = profile
            if attached_prealloc:
                if prealloc_seen:
                    fail("xla_prealloc may be specified only once")
                values["xla_prealloc"] = True
                prealloc_seen = True
        elif token.startswith("--profile="):
            if options_seen or values["profile"] is not None:
                fail("--profile may be specified only once and cannot follow --options")
            if values["accelerator"] is not None:
                fail("--profile and --accelerator cannot be combined")
            profile = token.split("=", 1)[1]
            if not profile:
                fail("--profile requires cpu or gpu")
            values["profile"] = profile.lower()
        elif token.startswith("--accelerator="):
            if options_seen or values["profile"] is not None:
                fail("--accelerator cannot be combined with --options or --profile")
            if values["accelerator"] is not None:
                fail("--accelerator may be specified only once")
            accelerator = token.split("=", 1)[1]
            if not accelerator:
                fail("--accelerator requires cpu, cuda13, or metal")
            values["accelerator"] = accelerator.lower()
        elif token == "--xla-prealloc":
            if prealloc_seen:
                fail("xla_prealloc may be specified only once")
            values["xla_prealloc"] = True
            prealloc_seen = True
        elif token == "--init_env":
            if values["init_env"]:
                fail("--init_env may be specified only once")
            values["init_env"] = True
        elif token in ("--info", "--info=json"):
            if values["info_format"] is not None:
                fail("--info may be specified only once")
            values["info_format"] = "json" if token.endswith("=json") else "human"
        elif token.startswith("--info="):
            fail("--info supports only the optional '=json' format")
        elif token == "--launcher-help":
            values["help"] = True
        else:
            break
        index += 1

    if values["init_env"] and values["info_format"] is not None:
        fail("--init_env and --info cannot be used together")
    return values, list(argv[index:])


def _native_gpu_accelerator() -> str:
    return "metal" if platform.system() == "Darwin" else "cuda13"


def resolve_accelerator(values: dict[str, object]) -> tuple[str, str]:
    accelerator = values["accelerator"]
    profile = values["profile"]
    if accelerator is None and profile is not None:
        if profile == "cpu":
            accelerator = "cpu"
        elif profile == "gpu":
            accelerator = _native_gpu_accelerator()
        elif profile in {"cuda13", "metal"}:
            accelerator = profile
        else:
            fail(f"JAX profile must be 'cpu' or 'gpu'; got: {profile}")
    if accelerator is None:
        accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").strip().lower()
        if accelerator == "gpu":
            accelerator = _native_gpu_accelerator()
    if accelerator not in {"cpu", "cuda13", "metal"}:
        fail(f"CLUBB_JAX_ACCELERATOR must be cpu, cuda13, or metal; got: {accelerator}")
    resolved_profile = "cpu" if accelerator == "cpu" else "gpu"
    return str(accelerator), resolved_profile


def runtime_paths(accelerator: str) -> tuple[Path, Path, str]:
    if accelerator == "cpu":
        return SCRIPT_DIR / "requirements.txt", REPO_ROOT / ".venv-jax", "cpu"
    if accelerator == "cuda13":
        return SCRIPT_DIR / "requirements-cuda13.txt", REPO_ROOT / ".venv-jax-cuda13", "cuda"
    return SCRIPT_DIR / "requirements-metal.txt", REPO_ROOT / ".venv-jax-metal", "METAL"


def _python_version(executable: Path | str) -> tuple[int, int]:
    result = subprocess.run(
        [
            str(executable),
            "-c",
            "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    major, minor = result.stdout.strip().split(".", 1)
    return int(major), int(minor)


def _python_is_supported(executable: Path | str) -> bool:
    try:
        return _python_version(executable) >= (3, 11)
    except (OSError, ValueError, subprocess.CalledProcessError):
        return False


def _python_is_compatible(executable: Path | str, accelerator: str) -> bool:
    if not _python_is_supported(executable):
        return False
    version = _python_version(executable)
    return accelerator != "metal" or version < (3, 13)


def _jax_version_for(executable: Path | str, accelerator: str) -> str:
    if accelerator == "metal":
        # jax-metal 0.1.1 is Apple's latest plugin and is qualified against
        # the JAX 0.4.34 release line. Keep this runtime isolated from the
        # newer CPU/CUDA environments.
        return "0.4.34"
    return "0.11.0" if _python_version(executable) >= (3, 12) else "0.10.0"


def _find_python(accelerator: str = "cpu") -> Path | None:
    explicit = os.environ.get("PYTHON")
    if explicit:
        found = shutil.which(explicit) or (explicit if Path(explicit).is_file() else None)
        if not found:
            fail(f"PYTHON does not exist: {explicit}")
        if not _python_is_compatible(found, accelerator):
            requirement = "Python 3.11 or 3.12" if accelerator == "metal" else "Python 3.11+"
            fail(f"JAX {accelerator} requires {requirement}; selected: {explicit}")
        return Path(found)
    names = (
        ("python3.12", "python3.11")
        if accelerator == "metal"
        else ("python3.14", "python3.13", "python3.12", "python3.11", "python3", "python")
    )
    for name in names:
        found = shutil.which(name)
        if found and _python_is_compatible(found, accelerator):
            return Path(found)
    return None


def _ensure_uv(tools_dir: Path, env: dict[str, str]) -> Path:
    local_uv = tools_dir / "bin" / "uv"
    if local_uv.is_file() and os.access(local_uv, os.X_OK):
        return local_uv
    path_uv = shutil.which("uv")
    if path_uv:
        return Path(path_uv)

    installer_url = f"https://astral.sh/uv/{UV_VERSION}/install.sh"
    (tools_dir / "bin").mkdir(parents=True, exist_ok=True)
    print(f"==> Installing uv {UV_VERSION} in {tools_dir / 'bin'}", flush=True)
    try:
        request = urllib.request.Request(
            installer_url,
            headers={"User-Agent": "CLUBB-JAX-launcher/1"},
        )
        # nosec B310: this is a pinned release from uv's official installer host.
        with urllib.request.urlopen(request, timeout=30) as response:
            installer = response.read()
    except OSError as exc:
        fail(f"could not download uv {UV_VERSION}: {exc}")
    install_env = dict(env)
    install_env["UV_UNMANAGED_INSTALL"] = str(tools_dir / "bin")
    result = subprocess.run(["/bin/sh"], input=installer, env=install_env)
    if result.returncode != 0 or not local_uv.is_file():
        fail(f"uv installation did not create {local_uv}")
    print(
        f"==> uv is ready: {subprocess.check_output([str(local_uv), '--version'], text=True).strip()}",
        flush=True,
    )
    return local_uv


def _inspection_python(venv: Path, accelerator: str) -> tuple[Path, str, str]:
    """Use the host Python for stdlib-only discovery, regardless of JAX readiness.

    The target interpreter is consulted only to plan package versions. Missing
    environments are reported as setup_required and are never created here.
    """
    venv_python = venv / "bin" / "python"
    candidate = None
    if venv_python.is_file() and _python_is_compatible(venv_python, accelerator):
        candidate = venv_python
    else:
        try:
            candidate = _find_python(accelerator)
        except LauncherError:
            # An unsuitable PYTHON override must not hide the hardware.
            # Environment creation validates that override before installing.
            pass
    if candidate is not None:
        planned_python = ".".join(map(str, _python_version(candidate)))
        required_jax = _jax_version_for(candidate, accelerator)
    else:
        # Match _prepare_environment's managed Python fallback.
        planned_python = "3.12"
        required_jax = "0.4.34" if accelerator == "metal" else "0.11.0"
    return Path(sys.executable), planned_python, required_jax


def _run_inspection(
    profile: str,
    accelerator: str,
    requirements: Path,
    venv: Path,
    info_format: str,
    env: dict[str, str],
    require_selectable: bool = False,
) -> int:
    info_python, planned_python, required_jax = _inspection_python(venv, accelerator)
    command = [
        str(info_python),
        str(RUNTIME_INFO),
        "--profile",
        profile,
        "--accelerator",
        accelerator,
        "--requirements",
        str(requirements),
        "--venv",
        str(venv),
        "--required-jax",
        required_jax,
        "--python-version",
        planned_python,
        "--format",
        info_format,
    ]
    if require_selectable:
        command.append("--require-selectable")
    return subprocess.run(command, env=env).returncode


def _normalize_cuda_visibility(env: dict[str, str], venv: Path) -> None:
    if "CUDA_VISIBLE_DEVICES" not in env:
        return
    info_python, _, _ = _inspection_python(venv, "cuda13")
    result = subprocess.run(
        [
            str(info_python),
            str(RUNTIME_INFO),
            "--resolve-visible-devices",
            env["CUDA_VISIBLE_DEVICES"],
        ],
        env=env,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or "unknown selection error"
        fail(f"Could not resolve CUDA_VISIBLE_DEVICES against nvidia-smi: {detail}")
    env["CUDA_VISIBLE_DEVICES"] = result.stdout.strip()


def _packages_are_ready(venv_python: Path, required_jax: str, accelerator: str) -> bool:
    script = """
from importlib.metadata import version
import sys
assert version('jax') == sys.argv[1]
assert version('jaxlib') == sys.argv[1]
for package in ('netCDF4', 'pytest', 'tabulate'):
    version(package)
if sys.argv[2] == 'cuda13':
    assert version('jax-cuda13-plugin') == sys.argv[1]
    assert version('jax-cuda13-pjrt') == sys.argv[1]
elif sys.argv[2] == 'metal':
    assert version('jax-metal') == '0.1.1'
"""
    result = subprocess.run(
        [str(venv_python), "-c", script, required_jax, accelerator],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def _prepare_environment(
    accelerator: str,
    requirements: Path,
    venv: Path,
    tools_dir: Path,
    env: dict[str, str],
) -> Path:
    if not requirements.is_file():
        fail(f"Missing requirements file: {requirements}")
    tools_dir.mkdir(parents=True, exist_ok=True)
    env["UV_CACHE_DIR"] = str(tools_dir / "cache")
    env["UV_PYTHON_INSTALL_DIR"] = str(tools_dir / "python")

    lock_path = tools_dir / "environment.lock"
    with lock_path.open("a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        venv_python = venv / "bin" / "python"
        if venv_python.is_file() and _python_is_compatible(venv_python, accelerator):
            print(
                f"==> Reusing JAX virtualenv: {venv} "
                f"(Python {'.'.join(map(str, _python_version(venv_python)))})",
                flush=True,
            )
        else:
            if venv.exists() and os.environ.get("CLUBB_JAX_VENV"):
                fail(f"Custom virtualenv is incomplete or uses Python older than 3.11: {venv}")
            python_cmd = _find_python(accelerator)
            if python_cmd is not None:
                python_spec = str(python_cmd)
                print(
                    f"==> Creating JAX virtualenv with existing Python "
                    f"{'.'.join(map(str, _python_version(python_cmd)))}",
                    flush=True,
                )
            else:
                python_spec = "3.12"
                print(
                    f"==> No compatible Python found for {accelerator}; "
                    "uv will download Python 3.12",
                    flush=True,
                )
            uv = _ensure_uv(tools_dir, env)
            command = [str(uv), "venv"]
            if venv.exists():
                command.append("--clear")
            command.extend(["--python", python_spec, str(venv)])
            subprocess.run(command, env=env, check=True)

        venv_python = venv / "bin" / "python"
        required_jax = _jax_version_for(venv_python, accelerator)
        requirements_hash = hashlib.sha256(requirements.read_bytes()).hexdigest()
        stamp = venv / ".clubb-jax-requirements.sha256"
        try:
            installed_hash = stamp.read_text(encoding="utf-8").strip()
        except OSError:
            installed_hash = ""
        if installed_hash != requirements_hash or not _packages_are_ready(
            venv_python, required_jax, accelerator
        ):
            uv = _ensure_uv(tools_dir, env)
            print(
                f"==> Installing JAX {required_jax} ({accelerator}) and "
                f"requirements from {requirements}",
                flush=True,
            )
            subprocess.run(
                [str(uv), "pip", "install", "--python", str(venv_python), "-r", str(requirements)],
                env=env,
                check=True,
            )
            if not _packages_are_ready(venv_python, required_jax, accelerator):
                fail(f"JAX {accelerator} packages failed validation after installation")
            stamp.write_text(requirements_hash + "\n", encoding="utf-8")
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
    return venv_python


def _verify_backend(venv_python: Path, accelerator: str, env: dict[str, str]) -> None:
    expected = {"cpu": "cpu", "cuda13": "gpu", "metal": "metal"}[accelerator]
    script = """
import jax, jaxlib, sys
backend = jax.default_backend().lower()
expected = sys.argv[1]
if expected == 'metal':
    ok = backend == 'metal' or any(str(d.platform).lower() == 'metal' for d in jax.devices())
else:
    ok = backend == expected
assert ok, f'requested {expected} but JAX initialized {backend}: {jax.devices()}'
print(
    f'JAX environment ready: jax={jax.__version__} '
    f'jaxlib={jaxlib.__version__} backend={backend} devices={jax.devices()}'
)
"""
    subprocess.run([str(venv_python), "-c", script, expected], env=env, check=True)


def ensure_environment() -> None:
    """Prepare the runtime and restart the calling script in it once."""
    marker = "_CLUBB_JAX_ENVIRONMENT_PYTHON"
    if os.environ.get(marker) == sys.executable:
        return
    result = subprocess.run([str(SCRIPT_DIR / "run_jax.py"), "--init_env"])
    if result.returncode:
        raise SystemExit(result.returncode)

    values, _ = parse_launcher_args([])
    accelerator, _ = resolve_accelerator(values)
    _, default_venv, _ = runtime_paths(accelerator)
    venv = Path(os.environ.get("CLUBB_JAX_VENV", str(default_venv)))
    python = str((REPO_ROOT / venv / "bin" / "python").absolute())
    os.environ["CLUBB_JAX_ACCELERATOR"] = accelerator
    os.environ[marker] = python
    if sys.executable != python:
        os.execve(python, [python, *sys.argv], os.environ.copy())


def main(argv: Sequence[str] | None = None) -> int:
    values, driver_args = parse_launcher_args(sys.argv[1:] if argv is None else argv)
    if values["help"]:
        print(usage())
        return 0

    accelerator, profile = resolve_accelerator(values)
    requirements, default_venv, jax_platform = runtime_paths(accelerator)
    venv = Path(os.environ.get("CLUBB_JAX_VENV", str(default_venv)))
    tools_dir = Path(os.environ.get("CLUBB_JAX_TOOLS_DIR", str(REPO_ROOT / ".clubb-jax-tools")))
    if not venv.is_absolute():
        venv = REPO_ROOT / venv
    if not tools_dir.is_absolute():
        tools_dir = REPO_ROOT / tools_dir

    env = dict(os.environ)
    env["CLUBB_JAX_PROFILE"] = profile
    env["CLUBB_JAX_ACCELERATOR"] = accelerator
    env["JAX_PLATFORMS"] = jax_platform
    xla_prealloc = bool(values["xla_prealloc"])
    if xla_prealloc and accelerator != "cuda13":
        fail("--xla-prealloc is a CUDA-only option")
    if accelerator == "cuda13":
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = (
            "true" if xla_prealloc else env.get("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
        )
        if not env.get("CUDA_DEVICE_ORDER"):
            env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
            env["CLUBB_JAX_CUDA_DEVICE_ORDER_SOURCE"] = "launcher default"
        else:
            env["CLUBB_JAX_CUDA_DEVICE_ORDER_SOURCE"] = "user setting"
    elif accelerator == "metal":
        env.setdefault("ENABLE_PJRT_COMPATIBILITY", "1")
        # Apple's Metal plugin supports float32 but not float64. Users can
        # still override this explicitly to diagnose future plugin releases.
        env.setdefault("CLUBB_JAX_PRECISION", "single")
        env.setdefault(
            "PYTHONWARNINGS",
            "ignore:Explicitly requested dtype:UserWarning",
        )

    info_format = values["info_format"]
    if info_format is not None:
        return _run_inspection(profile, accelerator, requirements, venv, str(info_format), env)

    if accelerator != "cpu":
        if _run_inspection(
            profile,
            accelerator,
            requirements,
            venv,
            "human",
            env,
            require_selectable=True,
        ) != 0:
            backend_name = "CUDA" if accelerator == "cuda13" else "Metal"
            fail(
                f"{accelerator} profile is unavailable; "
                f"no {backend_name} environment was created"
            )
    if accelerator == "cuda13":
        _normalize_cuda_visibility(env, venv)

    print(f"==> CLUBB JAX profile: {profile} ({accelerator})", flush=True)
    venv_python = _prepare_environment(accelerator, requirements, venv, tools_dir, env)
    if values["init_env"]:
        _verify_backend(venv_python, accelerator, env)
        return 0

    os.chdir(REPO_ROOT)
    os.execvpe(
        str(venv_python),
        [str(venv_python), "-m", "clubb_jax.src.clubb_standalone", *driver_args],
        env,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except LauncherError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except subprocess.CalledProcessError as exc:
        print(
            f"ERROR: command failed with exit code {exc.returncode}: "
            f"{' '.join(map(str, exc.cmd))}",
            file=sys.stderr,
        )
        raise SystemExit(exc.returncode or 1) from exc
