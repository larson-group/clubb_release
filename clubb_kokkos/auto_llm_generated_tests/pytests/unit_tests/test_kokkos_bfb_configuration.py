"""Exercise the local BFB policy without a Fortran compiler or model build.

This isolated CMake consumer checks that requesting the optional matrix backend cannot falsely enable it
when its compiler or either required library target is absent. Manual route:
human_reviewed_tests/run_pytests.sh -unit with this component's tests selected.
"""

import subprocess
from pathlib import Path

import pytest

POLICY = Path(__file__).resolve().parents[3] / "cmake" / "KokkosBfb.cmake"


@pytest.mark.parametrize("requested,fortran,bridge,lapack,effective", [
    (False, False, False, False, False),
    (True, False, False, False, False),
    (True, False, True, True, False),
    (True, True, False, True, False),
    (True, True, True, False, False),
    (True, True, True, True, True),
])
def test_kokkos_bfb_requires_available_fortran_and_bridge(
    tmp_path, requested, fortran, bridge, lapack, effective,
):
    """Configure the real policy; a missing dependency warns and disables the macro."""
    source = tmp_path / "CMakeLists.txt"
    source.write_text("\n".join([
        "cmake_minimum_required(VERSION 3.18)",
        "project(BfbPolicy LANGUAGES NONE)",
        f"set(KOKKOS_BFB {'ON' if requested else 'OFF'})",
        f"set(CMAKE_Fortran_COMPILER_LOADED {'ON' if fortran else 'OFF'})",
        "add_library(clubb_f2k_lib INTERFACE)" if bridge else "",
        "add_library(lapack_lib INTERFACE)" if lapack else "",
        f'include("{POLICY}")',
        'file(WRITE "${CMAKE_BINARY_DIR}/effective.txt" "${KOKKOS_BFB_ENABLED}")',
    ]) + "\n")
    result = subprocess.run(
        ["cmake", "-S", str(tmp_path), "-B", str(tmp_path / "build")],
        text=True, capture_output=True, timeout=30, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "build/effective.txt").read_text() == ("ON" if effective else "OFF")
    output = result.stdout + result.stderr
    assert ("backend is unavailable" in output) is (requested and not effective)
    if requested and not effective:
        assert "Using custom Cholesky" in output
        assert "effective KOKKOS_BFB is OFF" in output
