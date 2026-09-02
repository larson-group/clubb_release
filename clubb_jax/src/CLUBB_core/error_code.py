"""JAX-side mirror of `src/CLUBB_core/error_code.F90`.

Fortran comments:
  Description:
    Since f90/95 lacks enumeration, we're stuck numbering each
    error code by hand like this.

    We are "enumerating" error codes to be used with CLUBB. Adding
    additional codes is as simple adding an additional integer
    parameter. The error codes are ranked by severity, the higher
    number being more servere. When two errors occur, assign the
    most servere to the output.

    This code also handles subroutines related to debug_level. See
    the 'set_clubb_debug_level_api' description for more detail.

  References:
    None

Porting deviations:
- The Fortran module stores ``err_header`` and initializes process/thread
  headers for formatted error output.  The JAX mirror only needs the error
  constants and debug-level predicate used by the Python/JAX code, so the
  formatted Fortran I/O header routine is intentionally omitted.
- Fortran API suffixes are dropped from Python function names.
"""

CLUBB_NO_ERROR = 0
CLUBB_GENERALIZED_GRD_TEST_ERR = 50
CLUBB_FATAL_ERROR = 99

_debug_level = 0


def set_debug_level(level: int):
    """Accessor for clubb_debug_level.

    Fortran comments:
      Description:
        Accessor for clubb_debug_level

        0 => Print no debug messages to the screen
        1 => Print lightweight debug messages, e.g. print statements
        2 => Print debug messages that require extra testing,
             e.g. checks for NaNs and spurious negative values.
      References:
        None
    """
    global _debug_level
    # ---- Begin Code ----
    _debug_level = int(level)


def clubb_at_least_debug_level(level: int) -> bool:
    """Checks to see if clubb has been set to a specified debug level."""
    # ---- Begin Code ----
    return int(level) <= _debug_level
