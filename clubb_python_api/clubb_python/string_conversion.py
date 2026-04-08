"""Helpers for converting fixed-width Fortran string arrays at the F2PY boundary."""

import numpy as np


def fortran_char_matrix_to_python_strings(raw_values) -> list[str]:
    """Decode F2PY-returned fixed-width Fortran char data into Python strings."""
    raw_values = np.asarray(raw_values)

    if raw_values.ndim == 2 and raw_values.dtype.kind == "S":
        return [
            b"".join(row.tolist()).decode("ascii", errors="replace").strip()
            for row in raw_values
        ]

    return [
        value.decode("ascii", errors="replace").strip()
        if isinstance(value, bytes) else str(value).strip()
        for value in raw_values
    ]


def python_strings_to_fortran_char_matrix(values, width: int) -> np.ndarray:
    """Encode Python strings as a rank-2 byte matrix for F2PY fixed-char arrays."""
    matrix = np.full((len(values), width), b" ", dtype="S1", order="F")
    for i, value in enumerate(values):
        encoded = str(value).encode("ascii", errors="replace")[:width]
        if encoded:
            matrix[i, :len(encoded)] = np.frombuffer(encoded, dtype="S1")
    return matrix
