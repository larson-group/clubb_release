"""User-facing wrappers for routines from CLUBB_core/error_code.F90."""

import clubb_f2py

_debug_level = 0


def set_debug_level(level: int):
    """Set CLUBB debug verbosity level."""
    global _debug_level
    _debug_level = int(level)
    clubb_f2py.f2py_set_clubb_debug_level(_debug_level)


def clubb_at_least_debug_level(level: int):
    """Return whether CLUBB debug verbosity is at least level."""
    return _debug_level >= int(level)


def reset_err_code():
    """Reset error codes to no-error without re-allocating."""
    clubb_f2py.f2py_reset_err_code()


def initialize_error_headers():
    """Initialize the shared CLUBB error header string."""
    clubb_f2py.f2py_initialize_error_headers()
