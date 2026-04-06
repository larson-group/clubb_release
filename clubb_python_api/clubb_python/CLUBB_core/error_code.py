"""User-facing wrappers for routines from CLUBB_core/error_code.F90."""

import clubb_f2py


def set_debug_level(level: int):
    """Set CLUBB debug verbosity level."""
    clubb_f2py.f2py_set_clubb_debug_level(int(level))


def reset_err_code():
    """Reset error codes to no-error without re-allocating."""
    clubb_f2py.f2py_reset_err_code()


def initialize_error_headers():
    """Initialize the shared CLUBB error header string."""
    clubb_f2py.f2py_initialize_error_headers()
