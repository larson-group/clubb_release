"""ASTEX tendencies must retain distinct column grids and profile shapes."""
from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np

from clubb_jax.src.Benchmark_cases.astex_a209 import astex_a209_tndcy
from clubb_jax.src.CLUBB_core.grid_class import GRID_TYPE_STRETCHED_ZM, setup_grid


def test_subsidence_uses_each_column_grid_and_keeps_scalar_axes():
    """Linear subsidence follows each grid, with zero momentum boundaries."""
    heights = np.array([[0., 1., 2., 3., 4.], [10., 12., 14., 16., 18.]])
    gr = setup_grid(
        2, np.zeros(2), heights[:, 0], heights[:, -1],
        grid_type=GRID_TYPE_STRETCHED_ZM, momentum_heights=heights,
    )
    idx = SimpleNamespace(iisclr_thl=1, iisclr_rt=2, iiedsclr_thl=1, iiedsclr_rt=2)
    wm_zt, wm_zm, thlm, rtm, sclrm, edsclrm = astex_a209_tndcy(2, 2, 2, idx, gr)
    assert wm_zt.shape == thlm.shape == rtm.shape == (2, 4)
    assert wm_zm.shape == (2, 5)
    assert sclrm.shape == edsclrm.shape == (2, 4, 2)
    np.testing.assert_array_equal(wm_zt, -5.e-6 * np.asarray(gr.zt))
    expected_zm = -5.e-6 * heights
    expected_zm[:, [0, -1]] = 0.
    np.testing.assert_allclose(wm_zm, expected_zm, rtol=2.e-16, atol=0.)
    for field in (thlm, rtm, sclrm, edsclrm):
        np.testing.assert_array_equal(field, jnp.zeros_like(field))
