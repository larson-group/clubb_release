"""Contracts for Fortran-bound-aware Tune default intervals."""

from tuner.parameter_ranges import (
    compiled_parameter_hard_bounds,
    default_range_for_parameter,
    resolve_parameter_specs,
)
from utilities.clubb_settings_validation import parameter_hard_bounds


def test_tune_hard_bounds_do_not_require_the_compiled_python_api():
    assert compiled_parameter_hard_bounds() == parameter_hard_bounds()


def test_fully_bounded_parameter_uses_exact_fortran_interval():
    assert default_range_for_parameter(0.4, {"min": 0.0, "max": 1.0}) == (0.0, 1.0)


def test_open_side_keeps_legacy_relative_envelope():
    assert default_range_for_parameter(2.0, {"min": 0.0, "max": None}) == (0.0, 8.0)
    assert default_range_for_parameter(2.0, {"min": None, "max": None}) == (0.5, 8.0)


def test_zero_default_with_open_upper_bound_gets_a_useful_suggestion():
    assert default_range_for_parameter(0.0, {"min": 0.0, "max": None}) == (0.0, 4.0)


def test_zero_default_still_uses_an_exact_finite_fortran_upper_bound():
    assert default_range_for_parameter(0.0, {"min": 0.0, "max": 1.0}) == (0.0, 1.0)


def test_linked_targets_intersect_their_derived_ranges():
    resolved = resolve_parameter_specs(
        [{"name": "C6rt", "targets": ["C6rt", "C6thl"]}],
        {
            "C6rt": {"default": 2.0, "min": 0.0, "max": 8.0},
            "C6thl": {"default": 2.0, "min": 0.0, "max": 6.0},
        },
    )
    assert resolved == [{"name": "C6rt", "targets": ["C6rt", "C6thl"], "min": 0.0, "max": 6.0}]
