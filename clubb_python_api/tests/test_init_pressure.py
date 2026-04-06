"""Test init_pressure: hydrostatic pressure calculation via F2PY.

Validates that the Fortran init_pressure routine produces physically
reasonable pressure and Exner function profiles when called through the
Python interface layer.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

# Add clubb_python_api/ to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.derived_types.grid_class import Grid
from clubb_python.clubb_api import init_pressure


def _make_simple_grid(ngrdcol=1, nzt=10):
    """Create a simple evenly-spaced ascending grid for testing.

    Returns a Grid with nzt thermodynamic levels and nzm = nzt + 1
    momentum levels, evenly spaced from 0 to ~5000m.
    """
    nzm = nzt + 1
    dz = 500.0  # 500m spacing

    # zm levels: 0, 500, 1000, ..., 5000
    zm_1d = np.arange(nzm, dtype=np.float64) * dz
    # zt levels: midpoints: 250, 750, 1250, ..., 4750
    zt_1d = 0.5 * (zm_1d[:-1] + zm_1d[1:])

    zm = np.tile(zm_1d, (ngrdcol, 1))
    zt = np.tile(zt_1d, (ngrdcol, 1))

    dzm = np.full((ngrdcol, nzm), dz, dtype=np.float64)
    dzt = np.full((ngrdcol, nzt), dz, dtype=np.float64)
    invrs_dzm = np.full((ngrdcol, nzm), 1.0 / dz, dtype=np.float64)
    invrs_dzt = np.full((ngrdcol, nzt), 1.0 / dz, dtype=np.float64)

    # Equal interpolation weights
    weights_zt2zm = np.full((ngrdcol, nzm, 2), 0.5, dtype=np.float64)
    weights_zm2zt = np.full((ngrdcol, nzt, 2), 0.5, dtype=np.float64)

    return Grid(
        nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
        zm=zm, zt=zt, dzm=dzm, dzt=dzt,
        invrs_dzm=invrs_dzm, invrs_dzt=invrs_dzt,
        weights_zt2zm=weights_zt2zm, weights_zm2zt=weights_zm2zt,
        k_lb_zm=0, k_ub_zm=nzm - 1, k_lb_zt=0, k_ub_zt=nzt - 1,
        grid_dir_indx=1, grid_dir=1.0,
    )


class TestInitPressure:
    """Tests for the init_pressure wrapper."""

    def test_pressure_decreases_with_height(self):
        """Pressure should monotonically decrease with altitude."""
        gr = _make_simple_grid(ngrdcol=1, nzt=10)
        thvm = np.full((1, gr.nzt), 300.0, dtype=np.float64)
        p_sfc = np.array([101325.0], dtype=np.float64)

        p_zt, exner_zt, p_zm, exner_zm = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc
        )

        # zt-level pressure should decrease (ascending grid)
        for col in range(gr.ngrdcol):
            assert np.all(np.diff(p_zt[col, :]) < 0), \
                "Pressure on zt levels should decrease with height"

    def test_exner_decreases_with_height(self):
        """Exner function should decrease with altitude."""
        gr = _make_simple_grid(ngrdcol=1, nzt=10)
        thvm = np.full((1, gr.nzt), 300.0, dtype=np.float64)
        p_sfc = np.array([101325.0], dtype=np.float64)

        p_zt, exner_zt, p_zm, exner_zm = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc
        )

        for col in range(gr.ngrdcol):
            assert np.all(np.diff(exner_zt[col, :]) < 0), \
                "Exner on zt levels should decrease with height"

    def test_physically_reasonable_values(self):
        """Pressure and Exner should be in physically reasonable ranges."""
        gr = _make_simple_grid(ngrdcol=1, nzt=10)
        thvm = np.full((1, gr.nzt), 300.0, dtype=np.float64)
        p_sfc = np.array([101325.0], dtype=np.float64)

        p_zt, exner_zt, p_zm, exner_zm = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc
        )

        # Pressure: should be between ~40000 Pa (top ~5km) and ~105000 Pa
        assert np.all(p_zt > 30000.0), "Pressure too low"
        assert np.all(p_zt < 110000.0), "Pressure too high"

        # Exner: (p/p0)^(Rd/Cp) — should be between ~0.6 and ~1.05
        assert np.all(exner_zt > 0.5), "Exner too low"
        assert np.all(exner_zt < 1.1), "Exner too high"

        # zm-level values should also be reasonable
        assert np.all(p_zm > 25000.0), "zm pressure too low"
        assert np.all(p_zm < 110000.0), "zm pressure too high"
        assert np.all(exner_zm > 0.5), "zm Exner too low"
        assert np.all(exner_zm < 1.1), "zm Exner too high"

    def test_output_shapes(self):
        """Output arrays should have correct shapes."""
        ngrdcol, nzt = 2, 8
        gr = _make_simple_grid(ngrdcol=ngrdcol, nzt=nzt)
        thvm = np.full((ngrdcol, nzt), 300.0, dtype=np.float64)
        p_sfc = np.full(ngrdcol, 101325.0, dtype=np.float64)

        p_zt, exner_zt, p_zm, exner_zm = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc
        )

        assert p_zt.shape == (ngrdcol, nzt)
        assert exner_zt.shape == (ngrdcol, nzt)
        assert p_zm.shape == (ngrdcol, gr.nzm)
        assert exner_zm.shape == (ngrdcol, gr.nzm)

    def test_higher_surface_pressure(self):
        """Higher surface pressure should give higher pressures throughout."""
        gr = _make_simple_grid(ngrdcol=1, nzt=10)
        thvm = np.full((1, gr.nzt), 300.0, dtype=np.float64)

        p_sfc_lo = np.array([95000.0], dtype=np.float64)
        p_sfc_hi = np.array([105000.0], dtype=np.float64)

        p_lo, _, _, _ = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc_lo
        )
        p_hi, _, _, _ = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc_hi
        )

        assert np.all(p_hi > p_lo), \
            "Higher surface pressure should give higher pressure at all levels"

    def test_multi_column(self):
        """Multiple columns with different surface pressures."""
        ngrdcol = 3
        gr = _make_simple_grid(ngrdcol=ngrdcol, nzt=10)
        thvm = np.full((ngrdcol, gr.nzt), 300.0, dtype=np.float64)
        p_sfc = np.array([101325.0, 95000.0, 105000.0], dtype=np.float64)

        p_zt, exner_zt, p_zm, exner_zm = init_pressure(
            gr=gr, ngrdcol=gr.ngrdcol, nzt=gr.nzt, nzm=gr.nzm, thvm=thvm, p_sfc=p_sfc
        )

        # Column with highest surface pressure should have highest pressure
        # at every level
        assert np.all(p_zt[2, :] > p_zt[0, :]), \
            "Column 2 (p_sfc=105000) should have higher p than column 0 (101325)"
        assert np.all(p_zt[0, :] > p_zt[1, :]), \
            "Column 0 (p_sfc=101325) should have higher p than column 1 (95000)"
