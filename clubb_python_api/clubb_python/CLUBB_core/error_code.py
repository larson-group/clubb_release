"""User-facing wrappers for routines from CLUBB_core/error_code.F90."""

import clubb_f2py

CLUBB_NO_ERROR = 0
CLUBB_GENERALIZED_GRD_TEST_ERR = 50
CLUBB_FATAL_ERROR = 99

_debug_level = 0


def set_debug_level(level: int):
    """Set CLUBB debug verbosity level."""
    global _debug_level
    _debug_level = int(level)
    clubb_f2py.f2py_set_clubb_debug_level(_debug_level)


def clubb_at_least_debug_level(level: int):
    """Return whether CLUBB debug verbosity is at least level."""
    return _debug_level >= int(level)


def reset_err_code(err_info=None):
    """Reset error codes to no-error without re-allocating."""
    clubb_f2py.f2py_reset_err_code()
    from clubb_python.derived_types.err_info_converter import reset_cached_err_code

    reset_cached_err_code()
    if err_info is None:
        return None
    return err_info.reset_code()


def initialize_error_headers():
    """Initialize the shared CLUBB error header string."""
    clubb_f2py.f2py_initialize_error_headers()
