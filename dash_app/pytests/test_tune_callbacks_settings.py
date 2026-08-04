"""Regression checks for Tune's pattern-matching callback return values."""

from dash import no_update

from dash_app.tune_tab.callbacks_settings import removal_click_is_real, wildcard_no_update


def test_wildcard_no_update_matches_every_dynamic_component():
    """Dash requires an element per ALL-pattern output, even for no_update."""
    values = wildcard_no_update(["arm", "bomex"])

    assert values == [no_update, no_update]


def test_newly_hydrated_remove_buttons_do_not_delete_their_rows():
    """Dynamic ALL inputs arrive at zero clicks during workspace hydration."""
    row_order = [0, 1]

    assert removal_click_is_real([0, 0], row_order, 0) is False
    assert removal_click_is_real([1, 0], row_order, 1) is False
    assert removal_click_is_real([1, 1], row_order, 1) is True
