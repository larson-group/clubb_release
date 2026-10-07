"""CUDA 13 runtime policy, GPU inventory, and device selection."""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from .common import (
    GpuSelection,
    MIN_PYTHON,
    PYTHON_NAMES,
    PYTHON_REQUIREMENT,
    fail,
    first_error_line,
    jax_version,
)
from .common import expected_packages as base_packages

NAME = "cuda13"
LABEL = "CUDA"
DISPLAY_NAME = "CUDA 13"
PROFILE = "gpu"
REQUIREMENTS = "requirements-cuda13.txt"
VENV = ".venv-jax-cuda13"
PLATFORM = "cuda,cpu"
EXPECTED_BACKEND = "gpu"
SUPPORTS_PREALLOCATION = True
CUDA_ORDINALS = True
VISIBLE_DEVICES_VARIABLE = "CUDA_VISIBLE_DEVICES"
CUDA_MAJOR = 13
CUDA_MIN_DRIVER = 580
CUDA_MIN_COMPUTE_CAPABILITY = 7.5
VERIFY_SCRIPT = "ok = backend == 'gpu'"


def normalize_device(value: str | None) -> str:
    """An explicit job selection names one whole GPU; empty inherits visibility."""
    device = str(value or "").strip()
    if device and not re.fullmatch(
        r"GPU-[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-"
        r"[0-9a-fA-F]{4}-[0-9a-fA-F]{12}", device,
    ):
        raise ValueError("JAX device must be a full GPU UUID or empty for Default")
    return device


def configure_selection(
    env: dict[str, str], device: str, prealloc_gpu_mem: bool | None,
) -> None:
    """Apply this job's overrides before CUDA inspection or initialization."""
    if device:
        env[VISIBLE_DEVICES_VARIABLE] = normalize_device(device)
    if prealloc_gpu_mem is not None:
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = str(prealloc_gpu_mem).lower()


def has_gpu() -> bool:
    return bool(shutil.which("nvidia-smi"))


def expected_packages(required_jax: str) -> dict[str, str]:
    return {
        **base_packages(required_jax),
        "jax-cuda13-plugin": required_jax,
        "jax-cuda13-pjrt": required_jax,
    }


def configure_environment(env: dict[str, str], tools_dir: Path, xla_prealloc: bool = False) -> None:
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = (
        "true" if xla_prealloc else env.get("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    )
    if not env.get("CUDA_DEVICE_ORDER"):
        env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        env["CLUBB_JAX_CUDA_DEVICE_ORDER_SOURCE"] = "launcher default"
    else:
        env["CLUBB_JAX_CUDA_DEVICE_ORDER_SOURCE"] = "user setting"


def prepare_run(env: dict[str, str]) -> None:
    if "CUDA_VISIBLE_DEVICES" not in env:
        return
    try:
        env["CUDA_VISIBLE_DEVICES"] = resolve_visible_devices(env["CUDA_VISIBLE_DEVICES"])
    except ValueError as exc:
        fail(f"Could not resolve CUDA_VISIBLE_DEVICES against nvidia-smi: {exc}")


def query_gpus() -> tuple[list[dict[str, object]], str]:
    executable = shutil.which("nvidia-smi")
    if not executable:
        return [], "nvidia-smi was not found"

    fields = (
        "index",
        "uuid",
        "name",
        "driver_version",
        "memory.total",
        "compute_cap",
        "pci.bus_id",
    )
    command = [
        executable,
        f"--query-gpu={','.join(fields)}",
        "--format=csv,noheader,nounits",
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=3)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [], f"nvidia-smi could not inspect the GPU: {exc}"

    # Older nvidia-smi versions may not expose compute_cap. Retrying without it
    # still gives enough information to check the CUDA 13 driver requirement.
    if result.returncode != 0:
        fields = tuple(field for field in fields if field != "compute_cap")
        command = [
            executable,
            f"--query-gpu={','.join(fields)}",
            "--format=csv,noheader,nounits",
        ]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=3)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return [], f"nvidia-smi could not inspect the GPU: {exc}"

    if result.returncode != 0:
        detail = first_error_line(result.stderr or result.stdout)
        return [], "NVIDIA driver is unavailable" + (f": {detail}" if detail else "")

    has_compute_capability = "compute_cap" in fields
    gpus = []
    for row in csv.reader(result.stdout.splitlines(), skipinitialspace=True):
        if len(row) < len(fields):
            continue
        compute_capability = row[5].strip() if has_compute_capability else ""
        pci_bus_id = row[6 if has_compute_capability else 5].strip()
        try:
            memory_mib = int(float(row[4].strip()))
        except ValueError:
            memory_mib = None
        gpus.append(
            {
                "index": row[0].strip(),
                "uuid": row[1].strip(),
                "name": row[2].strip(),
                "driver_version": row[3].strip(),
                "memory_mib": memory_mib,
                "compute_capability": compute_capability,
                "pci_bus_id": pci_bus_id,
            }
        )
    return gpus, "" if gpus else "nvidia-smi did not report a GPU"


def resolve_visible_devices(visible_devices: str) -> str:
    """Use the same selection policy as inspection before CUDA initializes."""
    gpus, error = query_gpus()
    if error:
        raise ValueError(error)
    return select_gpus(gpus, visible_devices).resolved_devices


def _numeric_prefix(value: object) -> float | None:
    match = re.match(r"\d+(?:\.\d+)?", str(value).strip())
    return float(match.group()) if match else None


def select_gpus(
    gpus: list[dict[str, object]],
    requested: str | None,
    device_order: str = "FASTEST_FIRST",
) -> GpuSelection:
    """Resolve physical indices/unique UUID prefixes without importing JAX.

    Explicit lists retain their order and are converted to full UUIDs. With no
    list, all GPUs remain visible; only PCI ordering can be inferred here.
    """
    if requested is None:
        if not gpus:
            raise ValueError("No NVIDIA GPUs were detected.")
        if device_order == "PCI_BUS_ID":
            try:
                ordered = sorted(
                    gpus,
                    key=lambda gpu: tuple(
                        int(part, 16)
                        for part in re.split(r"[:.]", str(gpu["pci_bus_id"]))
                    ),
                )
                return GpuSelection(tuple(ordered), None)
            except (KeyError, ValueError):
                pass
        return GpuSelection(
            tuple(gpus), None, order_known=False,
            note="All GPUs are visible; the default device will be chosen by CUDA. "
            "Specify a GPU index or UUID to select one explicitly.",
        )

    if requested.strip() in {"", "-1"}:
        raise ValueError("CUDA_VISIBLE_DEVICES exposes no GPUs.")
    selected = []
    for selector in (item.strip() for item in requested.split(",")):
        if selector.isascii() and selector.isdecimal():
            matches = [gpu for gpu in gpus if str(gpu.get("index")) == str(int(selector))]
        elif selector.startswith("GPU-"):
            matches = [gpu for gpu in gpus if str(gpu.get("uuid", "")).startswith(selector)]
        else:
            raise ValueError(
                f"Unsupported GPU selector {selector!r}; use an nvidia-smi index "
                "or a unique GPU UUID (MIG selection is not supported by this preflight)."
            )
        if len(matches) != 1:
            detail = "ambiguous" if matches else "not found"
            raise ValueError(f"GPU selector {selector!r} is {detail} in nvidia-smi inventory.")
        gpu = matches[0]
        if not str(gpu.get("uuid", "")).startswith("GPU-"):
            raise ValueError(f"GPU {selector!r} has no usable UUID.")
        if gpu in selected:
            raise ValueError(f"GPU selector {selector!r} selects the same GPU twice.")
        selected.append(gpu)
    return GpuSelection(tuple(selected), ",".join(str(gpu["uuid"]) for gpu in selected))


def compatibility(
    gpus: Sequence[dict[str, object]], query_error: str
) -> tuple[bool, str]:
    if query_error:
        return False, query_error

    if not gpus:
        return False, "No GPUs are selected."
    # Every exposed GPU must be compatible, including non-default devices that
    # the CUDA backend may initialize. Hidden GPUs do not influence this check.
    for gpu in gpus:
        label = f"GPU {gpu.get('index', '?')} ({gpu.get('name', 'unknown')})"
        driver = gpu.get('driver_version') or 'unknown'
        if (_numeric_prefix(driver) or 0) < CUDA_MIN_DRIVER:
            return False, (
                f"{label}: CUDA {CUDA_MAJOR} requires NVIDIA driver "
                f"{CUDA_MIN_DRIVER} or newer; detected {driver}"
            )
        compute = gpu.get('compute_capability')
        # Older nvidia-smi may omit this field; backend initialization remains
        # the final check when it cannot be inspected in advance.
        if compute and (_numeric_prefix(compute) or 0) < CUDA_MIN_COMPUTE_CAPABILITY:
            return False, (
                f"{label}: CUDA {CUDA_MAJOR} requires compute capability "
                f"{CUDA_MIN_COMPUTE_CAPABILITY:.1f} or newer; detected {compute}"
            )
    return True, ""


def inspect_devices():
    gpus, error = query_gpus()
    requested = os.environ.get("CUDA_VISIBLE_DEVICES")
    selection = GpuSelection((), requested, order_known=False)
    try:
        if error:
            raise ValueError(error)
        selection = select_gpus(
            gpus, requested, os.environ.get("CUDA_DEVICE_ORDER", "FASTEST_FIRST")
        )
        selectable, reason = compatibility(selection.devices, "")
    except ValueError as exc:
        selectable, reason = False, str(exc)
    return gpus, selection, selectable, reason


def runtime_metadata() -> dict[str, object]:
    return {
        "xla_preallocate": os.environ.get("XLA_PYTHON_CLIENT_PREALLOCATE", "false"),
        "cuda_major": CUDA_MAJOR,
        "minimum_driver": CUDA_MIN_DRIVER,
        "minimum_compute_capability": CUDA_MIN_COMPUTE_CAPABILITY,
    }


def _memory_label(memory_mib: object) -> str:
    if not isinstance(memory_mib, int):
        return ""
    return f"{memory_mib / 1024:.1f} GiB"


def format_devices(info: dict[str, object], section_divider: str) -> list[str]:
    hardware = info["hardware"]
    lines = []
    visible_devices = hardware["cuda_visible_devices"]
    requested_visible_devices = hardware.get("requested_cuda_visible_devices")
    if requested_visible_devices and requested_visible_devices != visible_devices:
        visibility_lines = [
            f" Requested CUDA_VISIBLE_DEVICES: {requested_visible_devices}",
            " Resolved CUDA_VISIBLE_DEVICES: "
            + (str(visible_devices) if visible_devices is not None else "all GPUs")
            + " (GPU UUIDs)",
        ]
    else:
        visibility_lines = [
            " CUDA_VISIBLE_DEVICES: "
            + (
                str(visible_devices)
                if visible_devices is not None
                else "all GPUs"
            )
        ]
    lines.extend(
        [
            section_divider,
            " DEVICE SELECTION",
            section_divider,
            *visibility_lines,
            " CUDA_DEVICE_ORDER: "
            f"{hardware['cuda_device_order']} "
            f"({hardware['cuda_device_order_source']})",
            "",
            " PHYSICAL GPU INVENTORY",
            " GPU   JAX*  ACCESS        DEVICE",
            " ----  ----  ------------  --------------------------------------",
        ]
    )
    for gpu in hardware["gpus"]:
        physical_index = str(gpu.get("index") or "?")
        jax_device_index = gpu.get("jax_device_index")
        jax_label = str(jax_device_index) if isinstance(jax_device_index, int) else "--"
        if isinstance(jax_device_index, int):
            access = "SELECTED" if jax_device_index == 0 else "visible"
        else:
            access = "visible" if gpu.get("visible") else "hidden"
        lines.append(
            f" {physical_index:>3}  {jax_label:>4}  "
            f"{access:<12}  {gpu.get('name') or 'Unknown GPU'}"
        )
        details = [
            value
            for value in (
                _memory_label(gpu.get("memory_mib")),
                f"driver {gpu.get('driver_version')}"
                if gpu.get("driver_version")
                else "",
                f"compute {gpu.get('compute_capability')}"
                if gpu.get("compute_capability")
                else "",
            )
            if value
        ]
        if details:
            lines.append(f"                         {', '.join(details)}")
        identity = [
            value
            for value in (
                f"PCI {gpu.get('pci_bus_id')}" if gpu.get("pci_bus_id") else "",
                f"UUID {gpu.get('uuid')}" if gpu.get("uuid") else "",
            )
            if value
        ]
        if identity:
            lines.append(f"                         {', '.join(identity)}")
    selected_gpu = hardware["selected_gpu"]
    lines.extend([
        "", " GPU = nvidia-smi index; JAX* = expected index after selection.",
        section_divider, " PLANNED DEVICE", section_divider,
    ])
    if selected_gpu:
        lines.append(
            " Expected JAX GPU 0 --> "
            f"physical GPU {selected_gpu['index']} "
            f"({selected_gpu['name']})"
        )
    else:
        lines.append(" No physical GPU mapping is available.")
    if hardware["selection_note"]:
        lines.append(f" NOTE: {hardware['selection_note']}")
    return lines


def runtime_lines(info: dict[str, object]) -> list[str]:
    return [f" XLA memory preallocation: {info['runtime'].get('xla_preallocate', 'false')}"]
