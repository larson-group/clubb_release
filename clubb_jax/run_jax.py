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
import json
import os
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
if __package__:
    from . import backends
    from .backends.common import LauncherError, fail
else:
    import backends
    from backends.common import LauncherError, fail


def usage() -> str:
    return """Usage: clubb_jax/run_jax.py [-profile=cpu|gpu] [namelist-path] [driver-args...]
       clubb_jax/run_jax.py [-profile=cpu|gpu] -init_env
       clubb_jax/run_jax.py [-profile=cpu|gpu] -info[=json]

Creates or updates the JAX environment, then runs the CLUBB-JAX standalone.
Missing copies of uv and a supported Python are downloaded automatically.

Options:
  -options=VALUE        Parse the value forwarded by run_scm.py -jax=VALUE:
                         cpu or gpu, with comma-separated device=UUID,
                         prealloc_gpu_mem=true|false, or xla_prealloc modifiers
  -profile=cpu|gpu      Select CPU or the host-native GPU backend
  -accelerator=VALUE    Select an explicit backend: cpu, cuda13, rocm, or metal
  -xla_prealloc         Enable CUDA memory preallocation (CUDA only)
  -init_env             Prepare the environment without running a case
  -info[=json]          Inspect hardware and runtime readiness without setup
  -module=NAME         Entry module: CLUBB standalone, CLUBB driver test,
                        tuner.tune_clubb,
                        or clubb_jax.src.clubb_standalone_loss
  -launcher_help        Show this help

Environment:
  CLUBB_JAX_ACCELERATOR  Backend used when no option is given: cpu, cuda13, rocm, metal
  CLUBB_JAX_VENV         Override the profile virtualenv path
  CLUBB_JAX_TOOLS_DIR    Managed uv/Python path (default: .clubb-jax-tools)
  PYTHON                 Python used when creating a new virtualenv
  XLA_PYTHON_CLIENT_PREALLOCATE  CUDA memory preallocation (default: false)
"""


def _parse_options(value: str) -> dict:
    """Interpret the opaque selection string forwarded by SCM callers."""
    if not value:
        fail("-options requires a profile; supported profiles are cpu and gpu")
    if "\n" in value or "\r" in value:
        fail("JAX options must be on one line")
    pieces = value.split(",")
    if any(not piece for piece in pieces):
        fail("Empty JAX option; expected a profile followed by selection modifiers")
    selection = {"profile": pieces[0].lower(), "device": "", "prealloc_gpu_mem": None}
    seen = set()
    for modifier in pieces[1:]:
        key, separator, setting = modifier.partition("=")
        key = key.lower()
        if key == "xla_prealloc" and not separator:
            key, setting = "prealloc_gpu_mem", "true"
        elif key not in {"device", "prealloc_gpu_mem"} or not separator:
            fail(
                f"Unknown JAX option: {modifier}; supported modifiers are "
                "device=UUID, prealloc_gpu_mem=true|false, and xla_prealloc"
            )
        if key in seen:
            fail(f"JAX {key} may be specified only once")
        seen.add(key)
        if key == "device":
            selection[key] = backends.cuda.normalize_device(setting)
        else:
            if setting not in {"true", "false"}:
                fail("JAX prealloc_gpu_mem requires true or false")
            selection[key] = setting == "true"
    return selection


def parse_launcher_args(argv: Sequence[str]) -> tuple[dict[str, object], list[str]]:
    """Consume launcher options until the first driver argument."""
    values: dict[str, object] = {
        "profile": None,
        "accelerator": None,
        "xla_prealloc": False,
        "device": "",
        "prealloc_gpu_mem": None,
        "init_env": False,
        "info_format": None,
        "help": False,
        "module": "clubb_jax.src.clubb_standalone",
    }
    module_seen = False
    options_seen = False
    prealloc_seen = False
    index = 0
    while index < len(argv):
        token = argv[index]
        if token.startswith("-module="):
            if module_seen:
                fail("-module may be specified only once")
            module_seen = True
            values["module"] = token.split("=", 1)[1]
            if values["module"] not in {
                "clubb_jax.src.clubb_standalone", "clubb_jax.src.clubb_driver_test",
                "tuner.tune_clubb", "clubb_jax.src.clubb_standalone_loss",
            }:
                fail("Unsupported JAX entry module")
        elif token.startswith("-options="):
            if options_seen:
                fail("-options may be specified only once")
            if values["profile"] is not None or values["accelerator"] is not None:
                fail("-options cannot be combined with -profile or -accelerator")
            options_seen = True
            selection = _parse_options(token.split("=", 1)[1])
            values.update(selection)
            if selection["prealloc_gpu_mem"] is True:
                if prealloc_seen:
                    fail("xla_prealloc may be specified only once")
                values["xla_prealloc"] = True
                prealloc_seen = True
        elif token.startswith("-profile="):
            if options_seen or values["profile"] is not None:
                fail("-profile may be specified only once and cannot follow -options")
            if values["accelerator"] is not None:
                fail("-profile and -accelerator cannot be combined")
            profile = token.split("=", 1)[1]
            if not profile:
                fail("-profile requires cpu or gpu")
            values["profile"] = profile.lower()
        elif token.startswith("-accelerator="):
            if options_seen or values["profile"] is not None:
                fail("-accelerator cannot be combined with -options or -profile")
            if values["accelerator"] is not None:
                fail("-accelerator may be specified only once")
            accelerator = token.split("=", 1)[1]
            if not accelerator:
                fail("-accelerator requires cpu, cuda13, rocm, or metal")
            values["accelerator"] = accelerator.lower()
        elif token == "-xla_prealloc":
            if prealloc_seen:
                fail("xla_prealloc may be specified only once")
            values["xla_prealloc"] = True
            prealloc_seen = True
        elif token == "-init_env":
            if values["init_env"]:
                fail("-init_env may be specified only once")
            values["init_env"] = True
        elif token in ("-info", "-info=json"):
            if values["info_format"] is not None:
                fail("-info may be specified only once")
            values["info_format"] = "json" if token.endswith("=json") else "human"
        elif token.startswith("-info="):
            fail("-info supports only the optional '=json' format")
        elif token == "-launcher_help":
            values["help"] = True
        else:
            break
        index += 1

    if values["init_env"] and values["info_format"] is not None:
        fail("-init_env and -info cannot be used together")
    if values["xla_prealloc"] and values["prealloc_gpu_mem"] is False:
        fail("xla_prealloc conflicts with prealloc_gpu_mem=false")
    return values, list(argv[index:])


def runtime_selection(
    options: str | None = None, *,
    device: str = "", prealloc_gpu_mem: bool | None = None,
) -> dict:
    """Validate a requested runtime without installing packages or initializing JAX.

    A missing profile preserves the launcher's environment-based default.
    Device identifiers are opaque to callers; their format belongs to CUDA.
    """
    try:
        selection = (
            _parse_options(options) if options is not None
            else {"profile": None, "device": "", "prealloc_gpu_mem": None}
        )
        profile = selection["profile"]
        if profile is not None and profile not in {"cpu", "gpu"}:
            fail("JAX profile must be cpu or gpu")
        if prealloc_gpu_mem is not None and not isinstance(prealloc_gpu_mem, bool):
            fail("JAX prealloc_gpu_mem must be a boolean")
        if selection["prealloc_gpu_mem"] is not None and prealloc_gpu_mem is not None:
            if selection["prealloc_gpu_mem"] != prealloc_gpu_mem:
                fail("JAX option conflicts with the saved prealloc_gpu_mem setting")
        if prealloc_gpu_mem is not None:
            selection["prealloc_gpu_mem"] = prealloc_gpu_mem
        device = backends.cuda.normalize_device(device)
        if device and selection["device"] and device != selection["device"]:
            fail("JAX option conflicts with the saved device selection")
        if device:
            selection["device"] = device
        if selection["device"] and profile != "gpu":
            fail("JAX device requires the GPU profile")
        if selection["prealloc_gpu_mem"] is True and profile != "gpu":
            fail("JAX prealloc_gpu_mem requires the GPU profile")
    except LauncherError as exc:
        raise ValueError(str(exc)) from exc
    return selection


def runtime_arguments(
    options: str | None = None, *,
    device: str = "", prealloc_gpu_mem: bool | None = None, scm: bool = False,
) -> list[str]:
    """Encode a runtime selection for this launcher or an SCM forwarding script."""
    selection = runtime_selection(options, device=device, prealloc_gpu_mem=prealloc_gpu_mem)
    profile = selection["profile"]
    modifiers = []
    if selection["device"]:
        modifiers.append(f"device={selection['device']}")
    if selection["prealloc_gpu_mem"] is not None:
        modifiers.append(f"prealloc_gpu_mem={str(selection['prealloc_gpu_mem']).lower()}")
    if profile is None and modifiers:
        # Preserve an inherited explicit backend when encoding a memory setting.
        values, _ = parse_launcher_args([])
        profile, _ = resolve_accelerator(values)
    if profile is None:
        return ["-jax"] if scm else []
    value = ",".join([profile, *modifiers])
    return [f"-jax={value}" if scm else f"-options={value}"]


def resolve_accelerator(values: dict[str, object]) -> tuple[str, str]:
    accelerator = values["accelerator"]
    profile = values["profile"]
    if accelerator is None and profile is not None:
        if profile == "cpu":
            accelerator = "cpu"
        elif profile == "gpu":
            accelerator = backends.native_gpu_accelerator()
        elif profile in backends.BACKENDS:
            accelerator = profile
        else:
            fail(f"JAX profile must be 'cpu' or 'gpu'; got: {profile}")
    if accelerator is None:
        accelerator = os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").strip().lower()
        if accelerator == "gpu":
            accelerator = backends.native_gpu_accelerator()
    if accelerator not in backends.BACKENDS:
        fail(f"CLUBB_JAX_ACCELERATOR must be cpu, cuda13, rocm, or metal; got: {accelerator}")
    resolved_profile = backends.get_backend(str(accelerator)).PROFILE
    return str(accelerator), resolved_profile


def runtime_paths(accelerator: str) -> tuple[Path, Path, str]:
    backend = backends.get_backend(accelerator)
    return SCRIPT_DIR / backend.REQUIREMENTS, REPO_ROOT / backend.VENV, backend.PLATFORM


def _runtime_environment(
    values: dict, accelerator: str, profile: str, tools_dir: Path,
) -> dict[str, str]:
    """Build the child environment; backend details stay behind this boundary."""
    backend = backends.get_backend(accelerator)
    env = dict(os.environ)
    env["CLUBB_JAX_PROFILE"] = profile
    env["CLUBB_JAX_ACCELERATOR"] = accelerator
    env["JAX_PLATFORMS"] = backend.PLATFORM
    device = str(values["device"])
    prealloc_gpu_mem = values["prealloc_gpu_mem"]
    if device and accelerator != "cuda13":
        fail("Explicit device selection requires the CUDA GPU backend")
    if values["xla_prealloc"] and not getattr(backend, "SUPPORTS_PREALLOCATION", False):
        fail("-xla_prealloc is a CUDA-only option")
    if prealloc_gpu_mem is True and not getattr(backend, "SUPPORTS_PREALLOCATION", False):
        fail("prealloc_gpu_mem=true is a CUDA-only option")
    configure_selection = getattr(backend, "configure_selection", None)
    if configure_selection is not None:
        configure_selection(env, device, prealloc_gpu_mem)
    configure_environment = getattr(backend, "configure_environment", None)
    if configure_environment is not None:
        configure_environment(env, tools_dir, bool(values["xla_prealloc"]))
    return env


def runtime_configuration(
    options: str | None = None, *,
    device: str = "", prealloc_gpu_mem: bool | None = None,
) -> dict:
    """Describe the launch configuration for provenance, without setup or JAX.

    This describes requested/inherited settings. Use -info=json for hardware
    discovery and backend readiness; it does not confirm JAX initialization.
    """
    values, _ = parse_launcher_args(runtime_arguments(
        options, device=device, prealloc_gpu_mem=prealloc_gpu_mem,
    ))
    accelerator, profile = resolve_accelerator(values)
    requirements, default_venv, _ = runtime_paths(accelerator)
    tools_dir = Path(os.environ.get("CLUBB_JAX_TOOLS_DIR", str(REPO_ROOT / ".clubb-jax-tools")))
    if not tools_dir.is_absolute():
        tools_dir = REPO_ROOT / tools_dir
    env = _runtime_environment(values, accelerator, profile, tools_dir)
    precision_value = env.get("CLUBB_JAX_PRECISION", "double").strip().lower()
    precision = (
        "single" if precision_value in {"single", "float32", "f32", "32", "real4", "sp"}
        else "double"
    )
    environment_keys = (
        "CUDA_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES",
        "CLUBB_JAX_ACCELERATOR", "CLUBB_JAX_PRECISION", "CLUBB_JAX_VENV",
        "CLUBB_JAX_TOOLS_DIR", "XLA_PYTHON_CLIENT_PREALLOCATE", "JAX_PLATFORMS",
    )
    return {
        "profile": profile, "accelerator": accelerator, "precision": precision,
        "requirements": str(requirements),
        "venv": env.get("CLUBB_JAX_VENV", str(default_venv)),
        "environment": {key: env.get(key, "") for key in environment_keys},
    }


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
    backend = backends.get_backend(accelerator)
    maximum = getattr(backend, "MAX_PYTHON", None)
    return version >= backend.MIN_PYTHON and (maximum is None or version <= maximum)


def _jax_version_for(executable: Path | str, accelerator: str) -> str:
    backend = backends.get_backend(accelerator)
    if hasattr(backend, "JAX_VERSION"):
        return backend.JAX_VERSION
    return backend.jax_version(_python_version(executable))


def _find_python(accelerator: str = "cpu") -> Path | None:
    explicit = os.environ.get("PYTHON")
    if explicit:
        found = shutil.which(explicit) or (explicit if Path(explicit).is_file() else None)
        if not found:
            fail(f"PYTHON does not exist: {explicit}")
        if not _python_is_compatible(found, accelerator):
            requirement = backends.get_backend(accelerator).PYTHON_REQUIREMENT
            fail(f"JAX {accelerator} requires {requirement}; selected: {explicit}")
        return Path(found)
    names = backends.get_backend(accelerator).PYTHON_NAMES
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
        backend = backends.get_backend(accelerator)
        required_jax = (
            backend.JAX_VERSION if hasattr(backend, "JAX_VERSION")
            else backend.jax_version((3, 12))
        )
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
        "-profile",
        profile,
        "-accelerator",
        accelerator,
        "-requirements",
        str(requirements),
        "-venv",
        str(venv),
        "-required_jax",
        required_jax,
        "-python_version",
        planned_python,
        "-format",
        info_format,
    ]
    if require_selectable:
        command.append("-require_selectable")
    return subprocess.run(command, env=env).returncode


def _packages_are_ready(venv_python: Path, required_jax: str, accelerator: str) -> bool:
    expected = backends.get_backend(accelerator).expected_packages(required_jax)
    script = """
from importlib.metadata import version
import json
import sys
for package, required in json.loads(sys.argv[1]).items():
    assert version(package) == required
for package in ('netCDF4', 'pytest', 'tabulate'):
    version(package)
"""
    result = subprocess.run(
        [str(venv_python), "-c", script, json.dumps(expected)],
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
        if not (venv_python.is_file() and _python_is_compatible(venv_python, accelerator)):
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
            command = [str(uv), "pip", "install", "--python", str(venv_python), "-r", str(requirements)]
            backend = backends.get_backend(accelerator)
            install_arguments = getattr(backend, "install_arguments", None)
            if install_arguments is not None:
                command.extend(install_arguments(_python_version(venv_python), tools_dir))
            subprocess.run(
                command,
                env=env,
                check=True,
            )
            if not _packages_are_ready(venv_python, required_jax, accelerator):
                fail(f"JAX {accelerator} packages failed validation after installation")
            stamp.write_text(requirements_hash + "\n", encoding="utf-8")
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
    return venv_python


def _print_runtime_summary(venv_python: Path, accelerator: str) -> None:
    python_version = ".".join(map(str, _python_version(venv_python)))
    print("CLUBB JAX runtime:")
    print(f"  Environment: {venv_python.parent.parent} (Python {python_version})")
    print(f"  Accelerator: {accelerator}", flush=True)


def _verify_backend(venv_python: Path, accelerator: str, env: dict[str, str]) -> None:
    backend = backends.get_backend(accelerator)
    script = """
import jax, jaxlib, sys
backend = jax.default_backend().lower()
expected = sys.argv[1]
""" + backend.VERIFY_SCRIPT + "\n" + """
assert ok, f'requested {expected} but JAX initialized {backend}: {jax.devices()}'
devices = jax.devices()
labels = []
for device in devices:
    location = f'{device.platform.lower()}:{device.id}'
    kind = str(getattr(device, 'device_kind', '')).strip()
    labels.append(f'{location} ({kind})' if kind and kind.lower() != device.platform.lower()
                  else location)
print(f'  JAX: {jax.__version__} (jaxlib {jaxlib.__version__})')
print(f'  {"Device" if len(labels) == 1 else "Devices"}: {", ".join(labels)}')
"""
    subprocess.run([str(venv_python), "-c", script, backend.EXPECTED_BACKEND], env=env, check=True)


def ensure_environment() -> None:
    """Prepare the runtime and restart the calling script in it once."""
    marker = "_CLUBB_JAX_ENVIRONMENT_PYTHON"
    if os.environ.get(marker) == sys.executable:
        return
    result = subprocess.run([str(SCRIPT_DIR / "run_jax.py"), "-init_env"])
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
    requirements, default_venv, _ = runtime_paths(accelerator)
    venv = Path(os.environ.get("CLUBB_JAX_VENV", str(default_venv)))
    tools_dir = Path(os.environ.get("CLUBB_JAX_TOOLS_DIR", str(REPO_ROOT / ".clubb-jax-tools")))
    if not venv.is_absolute():
        venv = REPO_ROOT / venv
    if not tools_dir.is_absolute():
        tools_dir = REPO_ROOT / tools_dir

    backend = backends.get_backend(accelerator)
    env = _runtime_environment(values, accelerator, profile, tools_dir)

    info_format = values["info_format"]
    if info_format is not None:
        return _run_inspection(profile, accelerator, requirements, venv, str(info_format), env)

    if backend.PROFILE == "gpu":
        if _run_inspection(
            profile,
            accelerator,
            requirements,
            venv,
            "human",
            env,
            require_selectable=True,
        ) != 0:
            backend_name = backend.LABEL
            fail(
                f"{accelerator} profile is unavailable; "
                f"no {backend_name} environment was created"
            )
    prepare_run = getattr(backend, "prepare_run", None)
    if prepare_run is not None:
        prepare_run(env)

    venv_python = _prepare_environment(accelerator, requirements, venv, tools_dir, env)
    env["_CLUBB_JAX_ENVIRONMENT_PYTHON"] = str(venv_python)
    _print_runtime_summary(venv_python, accelerator)
    if values["init_env"]:
        _verify_backend(venv_python, accelerator, env)
        return 0

    os.chdir(REPO_ROOT)
    os.execvpe(
        str(venv_python),
        [str(venv_python), "-m", str(values["module"]), *driver_args],
        env,
    )
    return 0


def _record_tuner_launch_error(argv: Sequence[str], error_message: str) -> None:
    """Persist setup errors that occur before the tuner module can run."""
    if "-module=tuner.tune_clubb" not in argv:
        return
    job_dir = None
    for index, token in enumerate(argv):
        if token == "-job_dir" and index + 1 < len(argv):
            job_dir = Path(argv[index + 1])
        elif token.startswith("-job_dir="):
            job_dir = Path(token.split("=", 1)[1])
    if job_dir is not None:
        sys.path.insert(0, str(REPO_ROOT))
        from tuner.status import write_job_error

        write_job_error(job_dir / "status.json", job_dir / "results.json", error_message)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (LauncherError, ValueError) as exc:
        _record_tuner_launch_error(sys.argv[1:], str(exc))
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    except subprocess.CalledProcessError as exc:
        _record_tuner_launch_error(sys.argv[1:], str(exc))
        print(
            f"ERROR: command failed with exit code {exc.returncode}: "
            f"{' '.join(map(str, exc.cmd))}",
            file=sys.stderr,
        )
        raise SystemExit(exc.returncode or 1) from exc
    except OSError as exc:
        _record_tuner_launch_error(sys.argv[1:], str(exc))
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
