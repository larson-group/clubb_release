"""Focused coverage for linked custom-hypergrid coordinates."""

import numpy as np

from input.scripts.create_multi_col_params import custom_hypergrid, parse_hypergrid_range_spec


def test_linked_hypergrid_coordinate_assigns_identical_columns():
    specs = parse_hypergrid_range_spec("C6rtb=C6thlb/0:4/3")
    generated, count = custom_hypergrid({"C6rtb": 2.0, "C6thlb": 2.0, "C1": 1.0}, specs)

    assert count == 3
    assert np.array_equal(generated["C6rtb"], generated["C6thlb"])
    assert np.array_equal(generated["C6rtb"], np.array([0.0, 2.0, 4.0]))


def test_linked_hypergrid_rejects_overlapping_targets():
    specs = parse_hypergrid_range_spec("C6rtb=C6thlb/0:4/3,C6thlb/0:4/3")

    try:
        custom_hypergrid({"C6rtb": 2.0, "C6thlb": 2.0}, specs)
    except ValueError as exc:
        assert "more than one" in str(exc)
    else:
        raise AssertionError("overlapping linked hypergrid targets were accepted")
