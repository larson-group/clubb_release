"""Checks for pressure-coordinate forcing interpolation and application."""

from types import SimpleNamespace

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.time_dependent_input import (
    apply_time_dependent_forcings_from_array,
    read_to_grid,
)
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.Input_fields.input_names import rt_f_name, temperature_f_name


def _forcing_arrays(coordinates, values):
    ntimes = values.shape[1]
    coordinate_values = np.broadcast_to(np.asarray(coordinates)[:, None], (len(coordinates), ntimes))
    return coordinate_values.copy(), np.asarray(values)


def test_pressure_coordinate_interp():
    source_pressure = np.array([-90000.0, -50000.0])
    model_pressure = -np.linspace(95000.0, 40000.0, 25)
    values = np.array([[1.0e-4, 2.0e-4], [3.0e-4, 4.0e-4]])
    variables = _forcing_arrays(source_pressure, values)

    result = read_to_grid(2, 2, 25, model_pressure, *variables)
    for time_index in range(2):
        expected = np.interp(
            model_pressure,
            source_pressure,
            values[:, time_index],
            left=0.0,
            right=0.0,
        )
        np.testing.assert_allclose(result[:, time_index], expected, rtol=3.0e-16, atol=0.0)


def test_height_coordinate_interp_is_independent_of_pressure():
    source_height = np.array([500.0, 8000.0])
    model_height = np.linspace(100.0, 12000.0, 25)
    values = np.array([[1.0e-4], [3.0e-4]])
    variables = _forcing_arrays(source_height, values)

    result = read_to_grid(2, 1, 25, model_height, *variables)
    expected = np.interp(model_height, source_height, values[:, 0], left=0.0, right=0.0)
    np.testing.assert_allclose(result[:, 0], expected, rtol=3.0e-16, atol=0.0)


def test_apply_temperature_forcing_converts_to_thlm(monkeypatch):
    nzt = 8
    ngrdcol = 1
    temperature_forcing = jnp.linspace(-2.0e-4, -1.0e-4, nzt)
    names = ("z[m]", temperature_f_name) + (rt_f_name,) * (time_dependent_input.nforcings - 2)
    monkeypatch.setattr(time_dependent_input, "forcing_names", names)

    forcings_array = jnp.full((ngrdcol, nzt, time_dependent_input.nforcings), -999.9)
    forcings_array = forcings_array.at[:, :, 1].set(temperature_forcing)
    zeros_zt = jnp.zeros((ngrdcol, nzt))
    zeros_zm = jnp.zeros((ngrdcol, nzt - 1))
    exner = jnp.linspace(1.0, 0.6, nzt)[None, :]
    gr = SimpleNamespace(nzt=nzt, nzm=nzt - 1, ngrdcol=ngrdcol)
    sclr_idx = SclrIdx(-1, -1, -1, -1, -1, -1)

    result = apply_time_dependent_forcings_from_array(
        ngrdcol,
        nzt - 1,
        nzt,
        0,
        0,
        sclr_idx,
        False,
        gr,
        zeros_zt,
        jnp.ones_like(zeros_zt),
        exner,
        forcings_array,
        zeros_zt,
        zeros_zt,
        zeros_zt,
        zeros_zt,
        zeros_zt,
        zeros_zt,
        zeros_zt,
        zeros_zm,
        zeros_zt,
        zeros_zt,
        jnp.zeros((ngrdcol, nzt, 0)),
        jnp.zeros((ngrdcol, nzt, 0)),
    )

    np.testing.assert_allclose(result[0], temperature_forcing[None, :] / exner)
    assert isinstance(result[0], jax.Array)


def test_unknown_forcing_name_is_rejected(monkeypatch):
    nzt = 2
    names = ("z[m]", "unknown[K/s]") + (rt_f_name,) * (time_dependent_input.nforcings - 2)
    monkeypatch.setattr(time_dependent_input, "forcing_names", names)

    zeros_zt = jnp.zeros((1, nzt))
    with pytest.raises(ValueError, match="Incompatible forcing type"):
        apply_time_dependent_forcings_from_array(
            1,
            1,
            nzt,
            0,
            0,
            SclrIdx(-1, -1, -1, -1, -1, -1),
            False,
            SimpleNamespace(nzt=nzt, nzm=1, ngrdcol=1),
            zeros_zt,
            jnp.ones_like(zeros_zt),
            jnp.ones_like(zeros_zt),
            jnp.ones((1, nzt, time_dependent_input.nforcings)),
            zeros_zt,
            zeros_zt,
            zeros_zt,
            zeros_zt,
            zeros_zt,
            zeros_zt,
            zeros_zt,
            jnp.zeros((1, 1)),
            zeros_zt,
            zeros_zt,
            jnp.zeros((1, nzt, 0)),
            jnp.zeros((1, nzt, 0)),
        )
