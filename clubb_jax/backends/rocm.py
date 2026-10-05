"""ROCm 7.2 runtime policy, wheel provisioning, and GPU crash workarounds."""

from __future__ import annotations

import hashlib
import os
import platform
import re
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path

from .common import GpuSelection, fail, first_error_line
from .common import expected_packages as base_packages

NAME = "rocm"
LABEL = "ROCm"
DISPLAY_NAME = "ROCM"
PROFILE = "gpu"
REQUIREMENTS = "requirements-rocm.txt"
VENV = ".venv-jax-rocm"
PLATFORM = "rocm,cpu"
EXPECTED_BACKEND = "rocm"
MIN_PYTHON = (3, 12)
MAX_PYTHON = (3, 14)
PYTHON_NAMES = ("python3.14", "python3.13", "python3.12", "python3.11", "python3", "python")
PYTHON_REQUIREMENT = "Python 3.12 through 3.14"
JAX_VERSION = "0.11.0"
ROCM_WHEELHOUSE_URL = "https://github.com/ROCm/jax/releases/download/rocm-jax-v0.11.0/wheelhouse_legacy_rocm7.2.0.zip"
ROCM_WHEELHOUSE_SHA256 = "6ce7c9267682ed78fc889b6e2db8d1be7af42c42cd3977838ad971e747cd66b1"
VERIFY_SCRIPT = """ok = backend in ('rocm', 'gpu') and bool(jax.devices('rocm'))
import importlib
for component in ('_linalg', '_solver', '_sparse'):
    try:
        importlib.import_module('jax_rocm7_plugin.' + component)
    except ImportError as exc:
        raise RuntimeError(
            'ROCm JAX requires loadable hipBLAS, hipSOLVER/rocSOLVER, and hipSPARSE libraries; '
            'install the complete ROCm HIP SDK or set ROCM_PATH to a complete runtime'
        ) from exc"""


def expected_packages(required_jax: str) -> dict[str, str]:
    return {
        **base_packages(required_jax),
        "jax-rocm7-plugin": required_jax,
        "jax-rocm7-pjrt": required_jax,
    }


def has_gpu(gfx_target_version: int | None = None) -> bool:
    for properties in Path("/sys/class/kfd/kfd/topology/nodes").glob("*/properties"):
        try:
            values = dict(line.split(maxsplit=1) for line in properties.read_text().splitlines())
            if (int(values.get("simd_count", "0")) > 0 and
                    (gfx_target_version is None or int(values.get("gfx_target_version", "0")) == gfx_target_version)):
                return True
        except (OSError, ValueError, IndexError):
            continue
    return False


def configure_environment(env: dict[str, str], tools_dir: Path, xla_prealloc: bool = False) -> None:
    local_roots = (
        tools_dir / "rocm724-root" / "opt" / "rocm-7.2.4",
        tools_dir / "rocm72-root" / "opt" / "rocm-7.2.0",
    )
    local_root = next((candidate for candidate in local_roots if candidate.is_dir()), local_roots[-1])
    root = Path(env.get("ROCM_PATH", "/opt/rocm")).absolute()
    if (local_root.is_dir() and (not env.get("ROCM_PATH")
            or (root == Path("/opt/rocm") and not (root / "lib" / "libamdhip64.so").exists()))):
        root = local_root.absolute()
    env["ROCM_PATH"] = str(root)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(entry for entry in (
        str(root / "lib"), str(root / "lib64"), env.get("LD_LIBRARY_PATH", "")
    ) if entry)
    env.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    flags = env.get("XLA_FLAGS", "")
    if has_gpu(110501):
        env.setdefault("CLUBB_JAX_PORTABLE_CHOLESKY", "1")
        for flag, value in (("xla_gpu_autotune_level", "0"),
                            ("xla_gpu_enable_command_buffer", "")):
            if flag not in flags:
                flags = (flags + f" --{flag}={value}").strip()
        env["XLA_FLAGS"] = flags


def prepare_wheelhouse(tools_dir: Path) -> Path:
    wheelhouse = tools_dir / "wheelhouse-legacy-rocm7.2.0"
    ready = wheelhouse / ".ready"
    if (ready.is_file() and ready.read_text().strip() == ROCM_WHEELHOUSE_SHA256
            and len(list(wheelhouse.glob("*.whl"))) == 4):
        return wheelhouse
    archive_path = tools_dir / "wheelhouse_legacy_rocm7.2.0.zip"
    if not archive_path.is_file():
        print("==> Downloading pinned AMD JAX wheels for ROCm 7.2", flush=True)
        request = urllib.request.Request(ROCM_WHEELHOUSE_URL, headers={"User-Agent": "CLUBB-JAX-launcher/1"})
        temporary = archive_path.with_suffix(".tmp")
        try:
            with urllib.request.urlopen(request, timeout=60) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
            temporary.replace(archive_path)
        except OSError as exc:
            temporary.unlink(missing_ok=True)
            fail(f"Could not download ROCm JAX wheels: {exc}")
    with archive_path.open("rb") as archive_file:
        digest = hashlib.sha256()
        while chunk := archive_file.read(1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != ROCM_WHEELHOUSE_SHA256:
        fail(f"ROCm wheel archive checksum mismatch: {archive_path}")
    wheelhouse.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as archive:
        wheels = [entry for entry in archive.infolist()
                  if Path(entry.filename).name.startswith(("jax_rocm7_plugin-", "jax_rocm7_pjrt-"))
                  and entry.filename.endswith(".whl")]
        if len(wheels) != 4:
            fail("Pinned ROCm archive does not contain the expected JAX wheels")
        for entry in wheels:
            with archive.open(entry) as source, (wheelhouse / Path(entry.filename).name).open("wb") as output:
                shutil.copyfileobj(source, output)
    ready.write_text(ROCM_WHEELHOUSE_SHA256 + "\n")
    return wheelhouse


def install_arguments(python_version: tuple[int, int], tools_dir: Path) -> list[str]:
    wheelhouse = prepare_wheelhouse(tools_dir)
    python_tag = "cp" + "".join(map(str, python_version))
    return [str(wheel) for wheel in wheelhouse.glob("*.whl")
            if "-py3-" in wheel.name or f"-{python_tag}-" in wheel.name]


def query_gpus() -> tuple[list[dict[str, object]], str]:
    if platform.system() != "Linux":
        return [], "The ROCm backend requires Linux"
    root = Path(os.environ.get("ROCM_PATH", "/opt/rocm"))
    try:
        version = (root / ".info" / "version").read_text().strip()
    except OSError:
        return [], f"ROCm 7.2.x was not found at {root}; set ROCM_PATH if installed elsewhere"
    if not version.startswith("7.2."):
        return [], f"The pinned JAX ROCm wheels require ROCm 7.2.x; found {version}"
    if not os.access("/dev/kfd", os.R_OK | os.W_OK):
        return [], "ROCm requires read/write access to /dev/kfd (including outside a sandbox)"
    executable = shutil.which("rocminfo") or str(root / "bin" / "rocminfo")
    try:
        result = subprocess.run([executable], capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [], f"rocminfo could not inspect the GPU: {exc}"
    if result.returncode:
        return [], "ROCm GPU inspection failed: " + first_error_line(result.stderr or result.stdout)
    gpus = []
    for agent in re.split(r"(?m)^\s*Agent \d+\s*$", result.stdout):
        target = re.search(r"(?m)^\s*Name:\s+(gfx\w+)\s*$", agent)
        if target is None:
            continue
        marketing = re.search(r"(?m)^\s*Marketing Name:\s+(.+)$", agent)
        gpus.append({
            "index": str(len(gpus)), "uuid": "",
            "name": marketing.group(1).strip() if marketing else target.group(1),
            "driver_version": version, "memory_mib": None,
            "compute_capability": "", "pci_bus_id": "",
            "backend": "rocm", "gfx_target": target.group(1),
        })
    return gpus, "" if gpus else "rocminfo did not report an AMD GPU"


def inspect_devices():
    gpus, error = query_gpus()
    selectable = bool(gpus) and not error
    selection = GpuSelection(tuple(gpus), None, order_known=False)
    if not selectable:
        selection = GpuSelection((), None, order_known=False)
    return gpus, selection, selectable, error


def format_devices(info: dict[str, object], section_divider: str) -> list[str]:
    runtime = info["runtime"]
    hardware = info["hardware"]
    lines = []
    lines.extend([section_divider, " COMPUTE DEVICE", section_divider])
    for gpu in hardware["gpus"]:
        lines.append(f" AMD GPU: {gpu['name']} ({gpu['gfx_target']})")
        lines.append(f" ROCm: {gpu['driver_version']}")
    for variable in ("ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES"):
        lines.append(f" {variable}: {os.environ.get(variable, 'all GPUs')}")
    lines.append(" Device selection and ordering: verified when the ROCm backend initializes")
    if runtime.get("xla_flags"):
        lines.append(f" XLA_FLAGS: {runtime['xla_flags']}")
    return lines
