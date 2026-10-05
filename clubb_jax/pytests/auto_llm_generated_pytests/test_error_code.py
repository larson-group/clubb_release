"""Checks for the current Fortran debug-level semantics."""

from clubb_jax.src.CLUBB_core.error_code import (
    clubb_at_least_debug_level,
    set_debug_level,
)


def test_negative_debug_level_disables_level_zero_checks():
    try:
        set_debug_level(-1)
        assert not clubb_at_least_debug_level(0)
        assert clubb_at_least_debug_level(-1)
    finally:
        set_debug_level(0)
