"""Regression checks for the COBRA surface-flux Fortran mirror."""

import numpy as np

from clubb_jax.src.Benchmark_cases.cobra import cobra_sfclyr
from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx


def test_cobra_assigns_co2_heat_and_moisture_scalar_fluxes():
    time_dependent_input.time_sfc_given = np.linspace(0.0, 10.0, 49)
    time_dependent_input.sens_ht_given = np.linspace(100.0, 200.0, 49)
    time_dependent_input.latent_ht_given = np.linspace(50.0, 70.0, 49)
    time_dependent_input.CO2_sfc_given = np.linspace(-12.0, -10.0, 49)
    time_dependent_input.T_sfc_given = np.linspace(284.0, 286.0, 49)
    sclr_idx = SclrIdx(
        iisclr_rt=2, iisclr_thl=3, iisclr_CO2=1,
        iiedsclr_rt=3, iiedsclr_thl=2, iiedsclr_CO2=1,
    )
    rho = np.array([1.0, 1.25])

    heat, moisture, _ustar, passive, eddy, temperature = cobra_sfclyr(
        2,
        3,
        3,
        sclr_idx,
        5.0,
        np.array([20.0, 20.0]),
        rho,
        np.array([290.0, 291.0]),
        ubar=np.array([5.0, 6.0]),
    )

    expected_heat = 150.0 / (rho * Cp)
    expected_moisture = 60.0 / (rho * Lv)
    expected_co2 = -11.0 * 0.02897 / rho
    np.testing.assert_allclose(heat, expected_heat, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(moisture, expected_moisture, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(passive[:, 0], expected_co2, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(passive[:, 1], expected_moisture, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(passive[:, 2], expected_heat, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(eddy[:, 0], expected_co2, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(eddy[:, 1], expected_heat, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(eddy[:, 2], expected_moisture, rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(temperature, np.array([285.0, 285.0]))
