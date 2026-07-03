"""System-dependent defaults for tuner front ends."""

from __future__ import annotations

import os


def available_logical_cpu_count() -> int:
    """Return the logical CPUs available to this process."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def default_max_workers() -> int:
    """Return a conservative physical-core worker default."""
    return max(1, available_logical_cpu_count() // 2)
