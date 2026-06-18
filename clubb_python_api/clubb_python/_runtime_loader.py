"""Runtime helpers for importing CLUBB's F2PY extension safely."""

from __future__ import annotations

import ctypes
import subprocess
import sys
from pathlib import Path


_NETCDF_PRELOADED = False


def _pkg_config_output(*args: str) -> list[str]:
    """Return non-empty pkg-config output lines, or an empty list on failure."""
    try:
        result = subprocess.run(
            ["pkg-config", *args],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []

    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _candidate_netcdf_paths() -> list[Path]:
    """Return likely libnetcdf paths from the active pkg-config environment."""
    candidates: list[Path] = []

    for libdir in _pkg_config_output("--variable=libdir", "netcdf"):
        libdir_path = Path(libdir)
        for libname in ("libnetcdf.so.19", "libnetcdf.so"):
            candidate = libdir_path / libname
            if candidate.exists():
                candidates.append(candidate)

    return candidates


def preload_matching_netcdf_c() -> None:
    """Preload the pkg-config-selected NetCDF-C library before importing clubb_f2py.

    Some module environments provide a valid Intel netcdf-fortran package whose
    embedded RPATH still points at a mismatched NetCDF-C stack.  Preloading the
    NetCDF-C library selected by `pkg-config netcdf` ensures the Python F2PY
    extension resolves against the matching dependency chain instead of pulling
    in an incompatible fallback provider at import time.
    """
    global _NETCDF_PRELOADED

    if _NETCDF_PRELOADED or not sys.platform.startswith("linux"):
        return

    rtld_global = getattr(ctypes, "RTLD_GLOBAL", None)
    if rtld_global is None:
        return

    for candidate in _candidate_netcdf_paths():
        try:
            ctypes.CDLL(str(candidate), mode=rtld_global)
        except OSError:
            continue
        _NETCDF_PRELOADED = True
        return
