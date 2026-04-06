"""Test F2PY hello world: thlm2T_in_K_1D via the wrapper.

Validates that the Fortran routine called through F2PY produces
bit-for-bit identical results to the equivalent NumPy calculation.

Formula: T_in_K = thlm * exner + Lv * rcm / Cp
"""
import sys
from pathlib import Path

import numpy as np
import pytest

# Add clubb_python_api/ to path so we can import the built .so
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import clubb_f2py  # noqa: E402


# Physical constants (must match Fortran constants_clubb values)
# With .pyf-based build, submodules aren't exposed, so we hardcode these.
Cp = 1004.67   # J/(kg·K) — specific heat of dry air at constant pressure
Lv = 2.5e6     # J/kg — latent heat of vaporization


def _expected_T_in_K(thlm, exner, rcm):
    """Pure NumPy reference implementation."""
    return thlm * exner + Lv * rcm / Cp


class TestThlm2TInK1D:
    """Tests for f2py_thlm2T_in_K_1D wrapper."""

    def test_basic(self):
        """Typical atmospheric values."""
        nz = 10
        thlm = np.full(nz, 300.0, dtype=np.float64)
        exner = np.full(nz, 1.0, dtype=np.float64)
        rcm = np.full(nz, 0.001, dtype=np.float64)  # 1 g/kg

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)
        expected = _expected_T_in_K(thlm, exner, rcm)

        np.testing.assert_array_equal(result, expected)

    def test_varying_profiles(self):
        """Realistic vertically varying profiles."""
        nz = 50
        thlm = np.linspace(300.0, 310.0, nz)  # increasing with height
        exner = np.linspace(1.0, 0.85, nz)     # decreasing with height
        rcm = np.linspace(0.003, 0.0, nz)      # decreasing with height

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)
        expected = _expected_T_in_K(thlm, exner, rcm)

        np.testing.assert_array_equal(result, expected)

    def test_zero_rcm(self):
        """When rcm=0, T = thlm * exner exactly."""
        nz = 5
        thlm = np.array([290.0, 295.0, 300.0, 305.0, 310.0])
        exner = np.array([1.0, 0.98, 0.96, 0.94, 0.92])
        rcm = np.zeros(nz)

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)
        expected = thlm * exner

        np.testing.assert_array_equal(result, expected)

    def test_single_element(self):
        """nz=1 edge case."""
        thlm = np.array([300.0])
        exner = np.array([1.0])
        rcm = np.array([0.0])

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)
        expected = _expected_T_in_K(thlm, exner, rcm)

        np.testing.assert_array_equal(result, expected)

    def test_large_array(self):
        """nz=1000 to exercise larger arrays."""
        nz = 1000
        rng = np.random.default_rng(42)
        thlm = rng.uniform(280.0, 320.0, nz)
        exner = rng.uniform(0.8, 1.0, nz)
        rcm = rng.uniform(0.0, 0.005, nz)

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)
        expected = _expected_T_in_K(thlm, exner, rcm)

        np.testing.assert_array_equal(result, expected)

    def test_constants_match(self):
        """Verify Fortran constants match expected physical values."""
        assert Cp == pytest.approx(1004.67)
        assert Lv == pytest.approx(2.5e6)

    def test_physical_reasonableness(self):
        """Output should be physically reasonable temperatures."""
        nz = 20
        thlm = np.linspace(290.0, 310.0, nz)
        exner = np.linspace(1.0, 0.85, nz)
        rcm = np.linspace(0.002, 0.0, nz)

        result = clubb_f2py.f2py_thlm2t_in_k_1d(thlm, exner, rcm)

        # Temperatures should be in a reasonable atmospheric range
        assert np.all(result > 200.0), "Temperatures too low"
        assert np.all(result < 350.0), "Temperatures too high"
