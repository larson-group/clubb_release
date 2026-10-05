"""Static backend dispatch; importing backends never initializes JAX or installs packages."""

from __future__ import annotations

from . import cpu, cuda, metal, rocm

BACKENDS = {backend.NAME: backend for backend in (cpu, cuda, rocm, metal)}


def get_backend(accelerator: str):
    return BACKENDS[accelerator]


def native_gpu_accelerator() -> str:
    if metal.is_native_host():
        return metal.NAME
    if not cuda.has_gpu() and rocm.has_gpu():
        return rocm.NAME
    return cuda.NAME


def package_names() -> tuple[str, ...]:
    return tuple(dict.fromkeys(
        name for backend in BACKENDS.values() for name in backend.expected_packages("")
    ))
