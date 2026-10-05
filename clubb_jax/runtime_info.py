#!/usr/bin/env python3
"""Read-only runtime inspection for the self-contained JAX launcher."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path


if __package__:
    from . import backends
else:
    import backends


def _installed_runtime(venv: Path) -> dict[str, object]:
    python = venv / "bin" / "python"
    if not python.is_file() or not os.access(python, os.X_OK):
        return {}
    script = """
import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version

packages = {}
for name in json.loads(sys.argv[1]):
    try:
        packages[name] = version(name)
    except PackageNotFoundError:
        packages[name] = None
print(json.dumps({"python": platform.python_version(), "packages": packages}))
"""
    try:
        result = subprocess.run(
            [str(python), "-c", script, json.dumps(backends.package_names())], capture_output=True, text=True, timeout=3
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
    backend = backends.get_backend(accelerator)
    cpu = backends.cpu.hardware_info()
    gpus, selection, selectable, compatibility_reason = backend.inspect_devices()
    visible_variable = getattr(backend, "VISIBLE_DEVICES_VARIABLE", None)
    requested = os.environ.get(visible_variable) if visible_variable else None
    device_order = os.environ.get("CUDA_DEVICE_ORDER", "FASTEST_FIRST")
    for gpu in gpus:
        gpu["cuda_ordinal"] = None
        gpu["jax_device_index"] = None
        gpu["visible"] = gpu in selection.devices
    if selection.order_known:
        for index, gpu in enumerate(selection.devices):
            # These are expected *visible* ordinals, not nvidia-smi indices or
            # an attempt to reproduce CUDA's pre-filter hardware enumeration.
            gpu["cuda_ordinal"] = index if getattr(backend, "CUDA_ORDINALS", False) else None
            gpu["jax_device_index"] = index
    selected_gpu = selection.devices[0] if selection.devices and selection.order_known else None
    installed = _installed_runtime(venv)
    packages = (
        installed.get("packages")
        if isinstance(installed.get("packages"), dict)
        else {}
    )
    expected_packages = backend.expected_packages(required_jax)
    required_packages = list(expected_packages)

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
    runtime_metadata = {
        "xla_preallocate": None, "cuda_major": None,
        "minimum_driver": None, "minimum_compute_capability": None,
    }
    if hasattr(backend, "runtime_metadata"):
        runtime_metadata.update(backend.runtime_metadata())
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
            "xla_flags": os.environ.get("XLA_FLAGS", ""),
            **runtime_metadata,
            "backend_verified": False,
            "accelerator": accelerator,
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


def format_human(info: dict[str, object]) -> str:
    runtime = info["runtime"]
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
    backend = backends.get_backend(str(runtime["accelerator"]))
    lines.extend(backend.format_devices(info, section_divider))
    python = runtime["python"]
    jax = runtime["jax"]
    lines.extend(
        [
            section_divider,
            " RUNTIME",
            section_divider,
            " Backend: " + backend.DISPLAY_NAME,
            *(backend.runtime_lines(info) if hasattr(backend, "runtime_lines") else []),
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
    if len(sys.argv) == 3 and sys.argv[1] == "-resolve_visible_devices":
        try:
            print(backends.cuda.resolve_visible_devices(sys.argv[2]))
        except ValueError as exc:
            print(exc, file=sys.stderr)
            raise SystemExit(1) from exc
        return

    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument('-profile', dest='profile', choices=("cpu", "gpu"), required=True)
    parser.add_argument('-accelerator', dest='accelerator', required=True)
    parser.add_argument('-requirements', dest='requirements', type=Path, required=True)
    parser.add_argument('-venv', dest='venv', type=Path, required=True)
    parser.add_argument('-required_jax', dest='required_jax', required=True)
    parser.add_argument('-python_version', dest='python_version')
    parser.add_argument('-format', dest='format', choices=("human", "json"), default="human")
    parser.add_argument(
        '-require_selectable', dest='require_selectable',
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
