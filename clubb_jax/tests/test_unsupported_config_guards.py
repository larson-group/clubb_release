"""Regression checks for the active split driver's fail-loud feature gates."""

from __future__ import annotations

from pathlib import Path

import pytest

from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
from clubb_jax.src.Radiation.parameters_radiation import initialize_radiation_parameters
from clubb_jax.src.clubb_case_initalization import _check_unsupported_features


_BASE = {
    "microphys_scheme": "none",
    "rad_scheme": "none",
    "l_calc_thlp2_rad": False,
}


def _check(cfg=None, **overrides):
    arguments = {**_BASE, **overrides}
    _check_unsupported_features(
        {} if cfg is None else cfg,
        get_default_config_flags(),
        **arguments,
    )


def test_baseline_and_supported_radiation_schemes_pass():
    _check()
    for scheme in ("simplified", "simplified_bomex", "lba"):
        _check(rad_scheme=scheme)


@pytest.mark.parametrize(
    ("cfg", "overrides", "message"),
    [
        ({}, {"microphys_scheme": "morrison"}, "microphys_scheme"),
        ({"l_cloud_sed": True}, {}, "l_cloud_sed"),
        ({}, {"rad_scheme": "bugsrad"}, "rad_scheme"),
        ({}, {"l_calc_thlp2_rad": True}, "l_calc_thlp2_rad"),
        ({"wp2_sponge_damp_settings%l_sponge_damping": True}, {}, "Sponge damping"),
        ({"lh_microphys_type": "cluster"}, {}, "lh_microphys_type"),
        ({"l_silhs_rad": True}, {}, "l_silhs_rad"),
        ({"l_restart": True}, {}, "l_restart"),
        ({"l_input_fields": True}, {}, "l_input_fields"),
        ({"l_test_grid_generalization": True}, {}, "l_test_grid_generalization"),
        ({"grid_adapt_in_time_method": 1}, {}, "grid_adapt_in_time_method"),
    ],
)
def test_unsupported_active_driver_features_fail_loud(cfg, overrides, message):
    with pytest.raises(ValueError, match=message):
        _check(cfg, **overrides)


def test_invalid_solar_date_uses_the_source_calendar_failure():
    with pytest.raises(ValueError, match="gregorian2julian_day"):
        initialize_radiation_parameters(
            {"l_sw_radiation": True, "day": 32, "month": 12, "year": 2023},
            "simplified",
            Path("."),
        )
