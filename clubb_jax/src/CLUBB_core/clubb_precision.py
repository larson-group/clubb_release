"""Process-wide floating-point precision control for CLUBB-JAX.

The Fortran clubb_precision module selects the working real kind at build
time. JAX instead selects precision through the global jax_enable_x64 setting,
which must be configured consistently before arrays are created.

Set CLUBB_JAX_PRECISION before starting Python:

* double (default) enables float64 for Fortran-faithful validation.
* single or a supported alias disables x64 for lower-memory float32 runs.

Single precision is intended for performance and memory experiments. It is not
expected to match a double-precision Fortran build.
"""

import os

import jax


_SINGLE_ALIASES = ("single", "float32", "f32", "32", "real4", "sp")

CLUBB_JAX_PRECISION = os.environ.get(
    "CLUBB_JAX_PRECISION", "double"
).strip().lower()
USE_X64 = CLUBB_JAX_PRECISION not in _SINGLE_ALIASES


def configure_jax_precision() -> None:
    """Apply the process-wide precision selected by CLUBB_JAX_PRECISION."""
    jax.config.update("jax_enable_x64", USE_X64)
