"""Standard-library-only contracts shared by JAX backends."""

from __future__ import annotations

from dataclasses import dataclass

PYTHON_NAMES = ("python3.14", "python3.13", "python3.12", "python3.11", "python3", "python")
PYTHON_REQUIREMENT = "Python 3.11+"
MIN_PYTHON = (3, 11)


class LauncherError(RuntimeError):
    """A user-facing launcher configuration or setup error."""


def fail(message: str) -> None:
    raise LauncherError(message)


def jax_version(python_version: tuple[int, int]) -> str:
    return "0.11.0" if python_version >= (3, 12) else "0.10.0"


def expected_packages(required_jax: str) -> dict[str, str]:
    return {"jax": required_jax, "jaxlib": required_jax}


def first_error_line(value: str) -> str:
    return next((line.strip() for line in value.splitlines() if line.strip()), "")


@dataclass(frozen=True)
class GpuSelection:
    """A preflight plan, not confirmation that JAX initialized these devices."""

    devices: tuple[dict[str, object], ...]
    resolved_devices: str | None
    order_known: bool = True
    note: str = ""
