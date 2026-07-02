"""Shared constants and mutable runtime state for the compile tab."""

from __future__ import annotations

import os
import threading

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TOOLCHAIN_DIR = os.path.join(REPO_ROOT, "cmake", "toolchains")
BUILD_DIR = os.path.join(REPO_ROOT, "build")
INSTALL_DIR = os.path.join(REPO_ROOT, "install")
MAX_UI_LOG_LINES = 3000

COMPILER_NAME_MAP = {
    "intel": "intel",
    "intel-oneapi": "intel",
    "intel-classic": "intel",
    "ifx": "intel",
    "ifort": "intel",
    "nvhpc": "nvhpc",
    "nvfortran": "nvhpc",
    "gcc": "gcc",
    "gnu": "gcc",
    "gfortran": "gcc",
    "ftn": "cce",
    "cce": "cce",
    "crayftn": "cce",
}

COMPILER_COMMANDS = ("gfortran", "ifx", "ifort", "nvfortran", "ftn", "crayftn")

COMPILE_LOCK = threading.Lock()
COMPILE_PROC = {"proc": None, "job": None}
COMPILE_STREAM_LOCK = threading.Lock()

