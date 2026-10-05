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
from clubb_jax.src.Input_fields.input_reader import two_dim_read_var
from clubb_jax.src.Input_fields.input_names import rt_f_name, temperature_f_name


def _forcing_vars(coordinate_name, coordinates, values):
    ntimes = values.shape[1]
    coordinate_values = np.broadcast_to(np.asarray(coordinates)[:, None], (len(coordinates), ntimes))
    return [
        two_dim_read_var(coordinate_name, "Time[s]", coordinate_name, coordinate_values.copy()),
        two_dim_read_var(temperature_f_name, "Time[s]", coordinate_name, np.asarray(values)),
    ]


def test_pressure_coordinate_interp():
    source_pressure = np.array([-90000.0, -50000.0])
    model_pressure = -np.linspace(95000.0, 40000.0, 25)
    values = np.array([[1.0e-4, 2.0e-4], [3.0e-4, 4.0e-4]])
    variables = _forcing_vars("Press[Pa]", source_pressure, values)

    result = read_to_grid(2, 2, 2, 25, model_pressure, variables, temperature_f_name)
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
    variables = _forcing_vars("z[m]", source_height, values)

    result = read_to_grid(2, 2, 1, 25, model_height, variables, temperature_f_name)
    expected = np.interp(model_height, source_height, values[:, 0], left=0.0, right=0.0)
    np.testing.assert_allclose(result[:, 0], expected, rtol=3.0e-16, atol=0.0)


def test_apply_temperature_forcing_converts_to_thlm(monkeypatch):
    nzt = 8
    ngrdcol = 1
    temperature_forcing = jnp.linspace(-2.0e-4, -1.0e-4, nzt)
    forcing_data = [
        two_dim_read_var("z[m]", "Time[s]", "z[m]", np.zeros((nzt, 1))),
        two_dim_read_var(temperature_f_name, "Time[s]", "z[m]", np.zeros((nzt, 1))),
    ]
    forcing_data.extend(
        two_dim_read_var(rt_f_name, "Time[s]", "z[m]", np.zeros((nzt, 1)))
        for _ in range(time_dependent_input.nforcings - len(forcing_data))
    )
    monkeypatch.setattr(time_dependent_input, "t_dependent_forcing_data", forcing_data)

    forcings_array = jnp.full((time_dependent_input.nforcings, nzt), -999.9)
    forcings_array = forcings_array.at[1].set(temperature_forcing)
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
    forcing_data = [
        two_dim_read_var("z[m]", "Time[s]", "z[m]", np.zeros((nzt, 1))),
        two_dim_read_var("unknown[K/s]", "Time[s]", "z[m]", np.zeros((nzt, 1))),
    ]
    forcing_data.extend(
        two_dim_read_var(rt_f_name, "Time[s]", "z[m]", np.zeros((nzt, 1)))
        for _ in range(time_dependent_input.nforcings - len(forcing_data))
    )
    monkeypatch.setattr(time_dependent_input, "t_dependent_forcing_data", forcing_data)

    zeros_zt = jnp.zeros((1, nzt))
    with pytest.raises(ValueError, match="Incompatible forcing type"):
        apply_time_dependent_forcings_from_array(
            1,
            1,
            nzt,
            0,
            0,
            SclrIdx(-1, -1, -1, -1, -1, -1),
            SimpleNamespace(nzt=nzt, nzm=1, ngrdcol=1),
            zeros_zt,
            jnp.ones_like(zeros_zt),
            jnp.ones_like(zeros_zt),
            jnp.ones((time_dependent_input.nforcings, nzt)),
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
