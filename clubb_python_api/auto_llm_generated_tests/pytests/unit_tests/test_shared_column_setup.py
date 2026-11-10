"""Check shared standalone grid setup and retain supplied host geometry through F2PY."""

import numpy as np
import pytest

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info


@pytest.mark.parametrize("grid_type", [1, 2, 3])
def test_standalone_grid_setup_with_common_geometry(grid_type):
    """All standalone grid types produce a full batch with shared level counts and heights."""
    ngrdcol = 3
    zm = np.tile([0., 100., 200., 300., 400.], (ngrdcol, 1))
    zt = 0.5 * (zm[:, 1:] + zm[:, :-1])
    clubb_api.init_err_info(ngrdcol)
    gr, err = clubb_api.setup_grid(
        5, ngrdcol, np.zeros(ngrdcol), False, True, grid_type,
        np.full(ngrdcol, 100.), np.zeros(ngrdcol), np.full(ngrdcol, 400.),
        zm, zt, err_info=ErrInfo(ngrdcol=ngrdcol),
    )
    assert not np.any(err.err_code)
    assert gr.zm.shape == (ngrdcol, 5)
    assert gr.zt.shape == (ngrdcol, 4)
    np.testing.assert_array_equal(gr.zm, zm)
    np.testing.assert_array_equal(gr.zt, zt)


def test_host_supplied_column_heights_are_retained():
    """Shared standalone crop decisions do not replace host-supplied column geometry."""
    zm = np.array([[0., 100., 200., 300., 400.], [0., 120., 250., 380., 500.]])
    zt = 0.5 * (zm[:, 1:] + zm[:, :-1])
    clubb_api.init_err_info(2)
    gr, err = clubb_api.setup_grid(
        5, 2, np.zeros(2), True, True, 3,
        np.zeros(2), np.zeros(2), zm[:, -1], zm, zt,
        err_info=ErrInfo(ngrdcol=2),
    )
    assert not np.any(err.err_code)
    assert gr.zm.shape == (2, 5)
    assert gr.zt.shape == (2, 4)
    np.testing.assert_array_equal(gr.zm, zm)
    np.testing.assert_array_equal(gr.zt, zt)


@pytest.mark.parametrize("ascending", [True, False])
def test_stretched_grid_uses_each_columns_crop_bounds(ascending):
    """Different crop offsets retain equal counts for ascending and descending input."""
    zm = np.array([[0., 100., 200., 300., 400.], [50., 150., 250., 350., 450.]])
    if not ascending:
        zm = zm[:, ::-1]
    zt = 0.5 * (zm[:, 1:] + zm[:, :-1])
    clubb_api.init_err_info(2)
    gr, err = clubb_api.setup_grid(
        5, 2, np.zeros(2), False, ascending, 3,
        np.zeros(2), np.array([100., 50.]), np.array([300., 250.]),
        zm, zt, err_info=ErrInfo(ngrdcol=2),
    )
    assert not np.any(err.err_code)
    expected = np.array([[100., 200., 300.], [50., 150., 250.]])
    if not ascending:
        expected = expected[:, ::-1]
    assert gr.zm.shape == (2, 3)
    np.testing.assert_array_equal(gr.zm, expected)


@pytest.mark.parametrize("grid_type", [1, 2, 3])
def test_setup_rejects_unequal_retained_counts(grid_type):
    """Unequal counts fail before allocating a misleading common-size grid."""
    zm = np.tile([0., 100., 200., 300., 400.], (2, 1))
    zt = np.array([[50., 150., 250., 350.], [50., 150., 250., 450.]])
    clubb_api.init_err_info(2)
    with pytest.raises(RuntimeError, match="not initialized"):
        clubb_api.setup_grid(
            5, 2, np.zeros(2), False, True, grid_type,
            np.array([100., 200.]), np.zeros(2), np.array([400., 300.]),
            zm, zt, err_info=ErrInfo(ngrdcol=2),
        )
    assert np.all(get_fortran_err_info().err_code != 0)


@pytest.mark.parametrize("grid_type", [1, 2])
def test_setup_allows_different_domain_bounds_with_equal_counts(grid_type):
    """Equal level counts allow different start heights, top heights and spacing."""
    zm = np.array([[0., 100., 200., 300., 400.], [50., 250., 450., 650., 850.]])
    zt = 0.5 * (zm[:, 1:] + zm[:, :-1])
    clubb_api.init_err_info(2)
    gr, err = clubb_api.setup_grid(
        5, 2, np.zeros(2), False, True, grid_type,
        np.array([100., 200.]), zm[:, 0], zm[:, -1],
        zm, zt, err_info=ErrInfo(ngrdcol=2),
    )
    assert not np.any(err.err_code)
    assert gr.zm.shape == (2, 5)
    assert gr.zt.shape == (2, 4)
    np.testing.assert_array_equal(gr.zm, zm)
    np.testing.assert_array_equal(gr.zt, zt)
