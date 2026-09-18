#!/usr/bin/env python3
"""Read-only runtime inspection for the self-contained JAX launcher."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path


# Requirements for JAX's pip-managed CUDA 13 wheels:
# https://docs.jax.dev/en/latest/installation.html#nvidia-gpu
CUDA_MAJOR = 13
CUDA_MIN_DRIVER = 580
CUDA_MIN_COMPUTE_CAPABILITY = 7.5
METAL_MIN_MACOS = (14, 4)


def _first_error_line(value: str) -> str:
    return next((line.strip() for line in value.splitlines() if line.strip()), "")


def _cpu_info() -> dict[str, object]:
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


def _query_nvidia_gpus() -> tuple[list[dict[str, object]], str]:
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
        detail = _first_error_line(result.stderr or result.stdout)
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


def _query_metal_gpus() -> tuple[list[dict[str, object]], str]:
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
        detail = _first_error_line(result.stderr or result.stdout)
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


def resolve_visible_devices(visible_devices: str) -> str:
    """Use the same selection policy as inspection before CUDA initializes."""
    gpus, error = _query_nvidia_gpus()
    if error:
        raise ValueError(error)
    return select_gpus(gpus, visible_devices).resolved_devices


def _numeric_prefix(value: object) -> float | None:
    match = re.match(r"\d+(?:\.\d+)?", str(value).strip())
    return float(match.group()) if match else None


@dataclass(frozen=True)
class GpuSelection:
    """A preflight plan, not confirmation that JAX initialized these devices."""

    devices: tuple[dict[str, object], ...]
    resolved_devices: str | None
    order_known: bool = True
    note: str = ""


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


def _gpu_compatibility(
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


def _metal_compatibility(
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


def _installed_runtime(venv: Path) -> dict[str, object]:
    python = venv / "bin" / "python"
    if not python.is_file() or not os.access(python, os.X_OK):
        return {}
    script = """
import json
import platform
from importlib.metadata import PackageNotFoundError, version

packages = {}
for name in ("jax", "jaxlib", "jax-cuda13-plugin", "jax-cuda13-pjrt", "jax-metal"):
    try:
        packages[name] = version(name)
    except PackageNotFoundError:
        packages[name] = None
print(json.dumps({"python": platform.python_version(), "packages": packages}))
"""
    try:
        result = subprocess.run(
            [str(python), "-c", script], capture_output=True, text=True, timeout=3
        )
        if result.returncode == 0:
            return json.loads(result.stdout)
    except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError):
        pass
    return {}


def inspect_runtime(
    profile: str,
    accelerator: str,
    requirements: Path,
    venv: Path,
    required_jax: str,
    planned_python: str | None = None,
) -> dict[str, object]:
    cpu = _cpu_info()
    if accelerator == "cuda13":
        gpus, gpu_error = _query_nvidia_gpus()
    elif accelerator == "metal":
        gpus, gpu_error = _query_metal_gpus()
    else:
        gpus, gpu_error = [], ""
    requested = os.environ.get("CUDA_VISIBLE_DEVICES") if accelerator == "cuda13" else None
    device_order = os.environ.get("CUDA_DEVICE_ORDER", "FASTEST_FIRST")
    selection = GpuSelection((), requested, order_known=False)
    selectable, compatibility_reason = True, ""
    if accelerator == "cuda13":
        try:
            if gpu_error:
                raise ValueError(gpu_error)
            selection = select_gpus(gpus, requested, device_order)
            selectable, compatibility_reason = _gpu_compatibility(selection.devices, "")
        except ValueError as exc:
            selectable, compatibility_reason = False, str(exc)
    elif accelerator == "metal":
        selectable, compatibility_reason = _metal_compatibility(gpus, gpu_error)
        if selectable:
            selection = GpuSelection(tuple(gpus), None)
    for gpu in gpus:
        gpu["cuda_ordinal"] = None
        gpu["jax_device_index"] = None
        gpu["visible"] = gpu in selection.devices
    if selection.order_known:
        for index, gpu in enumerate(selection.devices):
            # These are expected *visible* ordinals, not nvidia-smi indices or
            # an attempt to reproduce CUDA's pre-filter hardware enumeration.
            gpu["cuda_ordinal"] = index if accelerator == "cuda13" else None
            gpu["jax_device_index"] = index
    selected_gpu = selection.devices[0] if selection.devices and selection.order_known else None
    installed = _installed_runtime(venv)
    packages = (
        installed.get("packages")
        if isinstance(installed.get("packages"), dict)
        else {}
    )
    required_packages = ["jax", "jaxlib"]
    expected_packages = {"jax": required_jax, "jaxlib": required_jax}
    if accelerator == "cuda13":
        required_packages.extend(("jax-cuda13-plugin", "jax-cuda13-pjrt"))
        expected_packages.update(
            {"jax-cuda13-plugin": required_jax, "jax-cuda13-pjrt": required_jax}
        )
    elif accelerator == "metal":
        required_packages.append("jax-metal")
        expected_packages["jax-metal"] = "0.1.1"

    requirements_hash = ""
    try:
        requirements_hash = hashlib.sha256(requirements.read_bytes()).hexdigest()
    except OSError:
        selectable = False
        compatibility_reason = f"Requirements file is missing: {requirements}"
    try:
        installed_hash = (venv / ".clubb-jax-requirements.sha256").read_text(
            encoding="utf-8"
        ).strip()
    except OSError:
        installed_hash = ""
    environment_ready = bool(
        installed
        and requirements_hash
        and requirements_hash == installed_hash
        and all(packages.get(name) == expected_packages[name] for name in required_packages)
    )

    if not selectable:
        status = "unavailable"
        reason = compatibility_reason
    elif environment_ready:
        status = "ready"
        reason = ""
    else:
        status = "setup_required"
        reason = "The managed environment will be prepared when the run starts."

    python_version = str(
        installed.get("python") or planned_python or platform.python_version()
    )
    return {
        "schema_version": 1,
        "profile": profile,
        "selectable": selectable,
        "status": status,
        "reason": reason,
        "hardware": {
            "cpu": cpu,
            "gpus": gpus,
            "cuda_visible_devices": selection.resolved_devices,
            "requested_cuda_visible_devices": requested,
            "cuda_device_order": device_order,
            "cuda_device_order_source": os.environ.get(
                "CLUBB_JAX_CUDA_DEVICE_ORDER_SOURCE", "user setting"
            ),
            "selection_note": selection.note,
            "selected_gpu": selected_gpu,
        },
        "runtime": {
            "xla_preallocate": (
                os.environ.get("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
                if accelerator == "cuda13"
                else None
            ),
            "backend_verified": False,
            "accelerator": accelerator,
            "cuda_major": CUDA_MAJOR if accelerator == "cuda13" else None,
            "minimum_driver": CUDA_MIN_DRIVER if accelerator == "cuda13" else None,
            "minimum_compute_capability": (
                CUDA_MIN_COMPUTE_CAPABILITY if accelerator == "cuda13" else None
            ),
            "venv": str(venv),
            "python": {
                "version": python_version,
                "installed": bool(installed),
            },
            "jax": {
                "required": required_jax,
                "installed": packages.get("jax"),
            },
            "jaxlib": {
                "required": required_jax,
                "installed": packages.get("jaxlib"),
            },
        },
    }


def _memory_label(memory_mib: object) -> str:
    if not isinstance(memory_mib, int):
        return ""
    return f"{memory_mib / 1024:.1f} GiB"


def format_human(info: dict[str, object]) -> str:
    runtime = info["runtime"]
    hardware = info["hardware"]
    profile = str(info["profile"]).upper()
    status = {
        "ready": "ENVIRONMENT READY",
        "setup_required": "SETUP REQUIRED",
        "unavailable": "UNAVAILABLE",
    }.get(info["status"], str(info["status"]).upper())
    divider = "=" * 78
    section_divider = "-" * 78
    lines = [
        divider,
        f" CLUBB JAX {profile} PRE-RUN CHECK",
        divider,
        f" STATUS: {status}",
    ]
    if info["profile"] == "cpu":
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
    elif runtime["accelerator"] == "metal":
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
    else:
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
    python = runtime["python"]
    jax = runtime["jax"]
    lines.extend(
        [
            section_divider,
            " RUNTIME",
            section_divider,
            " Backend: "
            + (
                "CUDA " + str(runtime["cuda_major"])
                if runtime["accelerator"] == "cuda13"
                else str(runtime["accelerator"]).upper()
            ),
            *([f" XLA memory preallocation: {runtime.get('xla_preallocate', 'false')}"]
              if runtime["accelerator"] == "cuda13" else []),
            f" Python: {python['version']} "
            f"({'installed' if python['installed'] else 'planned'})",
            f" JAX: {jax['installed'] or jax['required']} "
            f"({'installed' if jax['installed'] else 'planned'})",
            f" Environment: {runtime['venv']}",
            " Backend initialization: not checked (read-only preflight)",
        ]
    )
    if info["reason"]:
        lines.extend([section_divider, f" REASON: {info['reason']}"])
    lines.append(divider)
    return "\n".join(lines)


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--resolve-visible-devices":
        try:
            print(resolve_visible_devices(sys.argv[2]))
        except ValueError as exc:
            print(exc, file=sys.stderr)
            raise SystemExit(1) from exc
        return

    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("cpu", "gpu"), required=True)
    parser.add_argument("--accelerator", required=True)
    parser.add_argument("--requirements", type=Path, required=True)
    parser.add_argument("--venv", type=Path, required=True)
    parser.add_argument("--required-jax", required=True)
    parser.add_argument("--python-version")
    parser.add_argument("--format", choices=("human", "json"), default="human")
    parser.add_argument(
        "--require-selectable",
        action="store_true",
        help="exit unsuccessfully when the selected runtime cannot run on this host",
    )
    args = parser.parse_args()
    info = inspect_runtime(
        args.profile,
        args.accelerator,
        args.requirements,
        args.venv,
        args.required_jax,
        args.python_version,
    )
    if args.require_selectable and not info["selectable"]:
        print(info["reason"], file=sys.stderr)
        raise SystemExit(1)
    if args.format == "json":
        print(json.dumps(info, sort_keys=True))
    else:
        print(format_human(info))


if __name__ == "__main__":
    main()
