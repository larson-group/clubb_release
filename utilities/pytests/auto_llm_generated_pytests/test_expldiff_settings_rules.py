"""The no-build validator mirrors the l_do_expldiff_rtm_thlm eddy-scalar check."""

import pytest

from utilities.clubb_settings_validation import validate_clubb_settings

EXPLDIFF = {"l_do_expldiff_rtm_thlm": True}
RULE = "expldiff_eddy_scalar_indices"


def _codes(**kwargs):
    return {issue.code for issue in validate_clubb_settings(parameters={}, **kwargs)} & {RULE}


@pytest.mark.parametrize("indices", [(-1, -1), (-1, 2), (1, 0), (1, 5), (2, 2)])
def test_expldiff_requires_distinct_valid_indices(indices):
    thl, rt = indices
    assert _codes(flags=EXPLDIFF, edsclr_dim=4, iiedsclr_thl=thl, iiedsclr_rt=rt) == {RULE}


@pytest.mark.parametrize("kwargs", [
    dict(flags=EXPLDIFF, edsclr_dim=2, iiedsclr_thl=1, iiedsclr_rt=2),
    dict(flags=EXPLDIFF, edsclr_dim=4, iiedsclr_thl=4, iiedsclr_rt=2),
    dict(flags=EXPLDIFF, edsclr_dim=0),  # fewer than two eddy scalars: the flag has no effect
    dict(flags=EXPLDIFF, edsclr_dim=1, iiedsclr_thl=1, iiedsclr_rt=1),
    dict(flags={}, edsclr_dim=0),
    dict(flags=EXPLDIFF),  # edsclr_dim unknown: the rule is not evaluated
])
def test_expldiff_accepts_valid_disabled_or_unknown_eddy_scalars(kwargs):
    assert _codes(**kwargs) == set()
