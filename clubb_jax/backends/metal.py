"""Apple Metal runtime policy and hardware compatibility."""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from .common import GpuSelection, first_error_line
from .common import expected_packages as base_packages

NAME = "metal"
LABEL = "Metal"
DISPLAY_NAME = "METAL"
PROFILE = "gpu"
REQUIREMENTS = "requirements-metal.txt"
VENV = ".venv-jax-metal"
PLATFORM = "METAL,cpu"
EXPECTED_BACKEND = "metal"
MIN_PYTHON = (3, 11)
MAX_PYTHON = (3, 12)
PYTHON_NAMES = ("python3.12", "python3.11")
PYTHON_REQUIREMENT = "Python 3.11 or 3.12"
JAX_VERSION = "0.4.34"
METAL_MIN_MACOS = (14, 4)
VERIFY_SCRIPT = """ok = backend == 'metal' or any(
    str(device.platform).lower() == 'metal' for device in jax.devices()
)"""


def is_native_host() -> bool:
    return platform.system() == "Darwin"


def expected_packages(required_jax: str) -> dict[str, str]:
    return {**base_packages(required_jax), "jax-metal": "0.1.1"}


def configure_environment(env: dict[str, str], tools_dir: Path, xla_prealloc: bool = False) -> None:
    env.setdefault("ENABLE_PJRT_COMPATIBILITY", "1")
    env.setdefault("CLUBB_JAX_PRECISION", "single")
    env.setdefault("PYTHONWARNINGS", "ignore:Explicitly requested dtype:UserWarning")


def query_gpus() -> tuple[list[dict[str, object]], str]:
    """Return Apple GPUs advertised as Metal-capable by macOS."""
    if platform.system() != "Darwin":
        return [], "The Metal backend requires macOS"
    if platform.machine().lower() not in {"arm64", "aarch64"}:
        return [], "The JAX Metal plugin requires an Apple Silicon Mac"

    executable = shutil.which("system_profiler") or "/usr/sbin/system_profiler"
    try:
        result = subprocess.run(
            [executable, "SPDisplaysDataType", "-json"],
            capture_output=True,
            text=True,
            timeout=8,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [], f"system_profiler could not inspect the GPU: {exc}"
    if result.returncode != 0:
        detail = first_error_line(result.stderr or result.stdout)
        return [], "Metal GPU inspection failed" + (f": {detail}" if detail else "")
    try:
        displays = json.loads(result.stdout).get("SPDisplaysDataType", [])
    except (AttributeError, json.JSONDecodeError):
        return [], "system_profiler returned invalid display information"

    gpus: list[dict[str, object]] = []
    for item in displays:
        if not isinstance(item, dict):
            continue
        metal_supported = item.get("spdisplays_metal") == "spdisplays_supported"
        metal_family = str(item.get("spdisplays_mtlgpufamilysupport", ""))
        if not metal_supported and not metal_family.startswith("spdisplays_metal"):
            continue
        cores = item.get("sppci_cores")
        try:
            core_count = int(str(cores))
        except (TypeError, ValueError):
            core_count = None
        index = str(len(gpus))
        gpus.append(
            {
                "index": index,
                # CUDA UUIDs are deliberately not invented for a unified-memory
                # Apple GPU; the dashboard treats an empty UUID as host-native.
                "uuid": "",
                "name": item.get("sppci_model") or item.get("_name") or "Apple GPU",
                "driver_version": platform.mac_ver()[0],
                "memory_mib": None,
                "compute_capability": "",
                "pci_bus_id": "",
                "metal_cores": core_count,
                "backend": "metal",
            }
        )
    return gpus, "" if gpus else "system_profiler did not report a Metal-capable GPU"


def compatibility(
    gpus: Sequence[dict[str, object]], query_error: str
) -> tuple[bool, str]:
    if query_error:
        return False, query_error
    if not gpus:
        return False, "No Metal-capable Apple GPU was detected."
    version = tuple(
        int(part) for part in platform.mac_ver()[0].split(".")[:2] if part.isdigit()
    )
    if version and version < METAL_MIN_MACOS:
        required = ".".join(map(str, METAL_MIN_MACOS))
        return False, f"jax-metal requires macOS {required} or newer"
    return True, ""


def inspect_devices():
    gpus, error = query_gpus()
    selectable, reason = compatibility(gpus, error)
    selection = (
        GpuSelection(tuple(gpus), None)
        if selectable else GpuSelection((), None, order_known=False)
    )
    return gpus, selection, selectable, reason


def format_devices(info: dict[str, object], section_divider: str) -> list[str]:
    hardware = info["hardware"]
    lines = []
    lines.extend(
        [
            section_divider,
            " COMPUTE DEVICE",
            section_divider,
        ]
    )
    for gpu in hardware["gpus"]:
        core_count = gpu.get("metal_cores")
        core_label = f" ({core_count} GPU cores)" if core_count else ""
        lines.append(f" Apple GPU: {gpu.get('name') or 'Unknown GPU'}{core_label}")
    lines.extend(
        [
            " Unified memory: shared with the CPU",
            " Precision: float32 (jax-metal does not support float64)",
        ]
    )
    return lines
