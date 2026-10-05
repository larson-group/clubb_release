"""CPU runtime policy and host inspection."""

from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path

from .common import (
    GpuSelection,
    MIN_PYTHON,
    PYTHON_NAMES,
    PYTHON_REQUIREMENT,
    expected_packages,
    jax_version,
)

NAME = "cpu"
LABEL = "CPU"
DISPLAY_NAME = "CPU"
PROFILE = "cpu"
REQUIREMENTS = "requirements.txt"
VENV = ".venv-jax"
PLATFORM = "cpu"
EXPECTED_BACKEND = "cpu"
VERIFY_SCRIPT = "ok = backend == 'cpu'"


def hardware_info() -> dict[str, object]:
    model = ""
    try:
        for line in Path("/proc/cpuinfo").read_text(
            encoding="utf-8", errors="replace"
        ).splitlines():
            if line.lower().startswith(("model name", "hardware")) and ":" in line:
                model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    if not model and platform.system() == "Darwin":
        try:
            model = subprocess.check_output(
                ["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"],
                text=True,
                stderr=subprocess.DEVNULL,
                timeout=2,
            ).strip()
        except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass
    model = model or platform.processor().strip()
    return {
        "model": model or platform.machine() or "Unknown CPU",
        "logical_cpus": os.cpu_count(),
    }


def inspect_devices():
    return [], GpuSelection((), None, order_known=False), True, ""


def format_devices(info: dict[str, object], section_divider: str) -> list[str]:
    hardware = info["hardware"]
    lines = []
    cpu = hardware["cpu"]
    count = cpu.get("logical_cpus")
    suffix = f" ({count} logical CPUs)" if count else ""
    lines.extend(
        [
            section_divider,
            " COMPUTE DEVICE",
            section_divider,
            f" CPU: {cpu['model']}{suffix}",
        ]
    )
    return lines
