"""Unit tests for fixed-width string conversion helpers."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.string_conversion import (
    fortran_char_matrix_to_python_strings,
    python_strings_to_fortran_char_matrix,
)


def test_python_strings_to_fortran_char_matrix():
    """Python strings should encode to a padded rank-2 byte matrix."""
    result = python_strings_to_fortran_char_matrix(["C1", "C2rt"], width=6)

    assert result.shape == (2, 6)
    assert result.dtype == np.dtype("S1")
    assert result.tolist()[0] == [b"C", b"1", b" ", b" ", b" ", b" "]
    assert result.tolist()[1] == [b"C", b"2", b"r", b"t", b" ", b" "]


def test_fortran_char_matrix_to_python_strings():
    """Rank-2 byte matrices should decode to stripped Python strings."""
    raw = np.array(
        [
            [b"C", b"1", b" ", b" ", b" ", b" "],
            [b"C", b"2", b"r", b"t", b" ", b" "],
        ],
        dtype="S1",
        order="F",
    )

    result = fortran_char_matrix_to_python_strings(raw)

    assert result == ["C1", "C2rt"]


def test_python_strings_to_fortran_char_matrix_truncates():
    """Strings longer than the target width should be truncated."""
    result = python_strings_to_fortran_char_matrix(["abcdef"], width=4)

    assert result.tolist()[0] == [b"a", b"b", b"c", b"d"]


def test_fortran_char_matrix_to_python_strings_handles_empty():
    """All-space rows should decode to empty strings."""
    raw = np.array([[b" ", b" ", b" "]], dtype="S1", order="F")

    assert fortran_char_matrix_to_python_strings(raw) == [""]
