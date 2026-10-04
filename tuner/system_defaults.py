"""CPU availability and worker defaults shared by scripts and tuner front ends."""

from __future__ import annotations

import os


def available_logical_cpu_count() -> int:
    """Count logical CPUs available to this process, respecting CPU affinity."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def default_max_workers() -> int:
    """Use half the available logical CPUs, rounded down, with at least one worker."""
    return max(1, available_logical_cpu_count() // 2)
