"""Tests for surface-boundary input extraction."""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.Benchmark_cases.prescribe_forcings import read_surface_var_for_bc
from clubb_jax.src.CLUBB_core.constants_clubb import kappa, p0
from clubb_jax.src.CLUBB_core.grid_class import setup_grid, zt2zm
from clubb_jax.src.CLUBB_core.interpolation import mono_cubic_interp


def test_interpolated_surface_values_match_column_reference():
    """The vectorized path must preserve the former per-column calculation."""
    ngrdcol = 3
    gr = setup_grid(
        ngrdcol,
        deltaz=np.array([10.0, 12.0, 15.0]),
        zm_init=0.0,
        zm_top=100.0,
    )
    offsets = jnp.arange(ngrdcol, dtype=jnp.float64)[:, None]
    um = 2.0 + 0.03 * gr.zt + offsets
    vm = -1.0 + 0.02 * gr.zt - offsets
    rtm = 0.01 - 1.0e-5 * gr.zt + offsets * 1.0e-4
    thlm = 285.0 + 0.01 * gr.zt + offsets
    exner = 0.99 - 1.0e-4 * gr.zt
    rho_zm = 1.2 - 1.0e-3 * gr.zm + offsets * 0.01
    p_sfc = jnp.array([100000.0, 99000.0, 98000.0])

    compiled_read_surface = jax.jit(
        read_surface_var_for_bc, static_argnums=(1, 9)
    )
    actual = compiled_read_surface(
        gr, ngrdcol, um, vm, rtm, thlm, rho_zm, exner, p_sfc, True
    )

    fields_zm = [
        zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, field)
        for field in (um, vm, rtm, thlm, exner)
    ]
    fields_zm[-1] = fields_zm[-1].at[:, 0].set((p_sfc / p0) ** kappa)
    k_min = jnp.argmin(jnp.abs(gr.zt - 25.0), axis=1)
    expected_fields = [[] for _ in fields_zm]
    expected_rho = []
    for column in range(ngrdcol):
        k00 = int(k_min[column])
        km1 = max(k00 - 1, 0)
        kp1 = k00 + 1
        kp2 = k00 + 2
        expected_rho.append(rho_zm[column, k00])
        for expected, values in zip(expected_fields, fields_zm):
            expected.append(
                mono_cubic_interp(
                    25.0,
                    km1,
                    k00,
                    kp1,
                    kp2,
                    gr.zm[column, km1],
                    gr.zm[column, k00],
                    gr.zm[column, kp1],
                    gr.zm[column, kp2],
                    values[column, km1],
                    values[column, k00],
                    values[column, kp1],
                    values[column, kp2],
                )
            )

    expected = (
        jnp.full(ngrdcol, 25.0),
        jnp.asarray(expected_fields[0]),
        jnp.asarray(expected_fields[1]),
        jnp.asarray(expected_fields[2]),
        jnp.asarray(expected_fields[3]),
        jnp.asarray(expected_rho),
        jnp.asarray(expected_fields[4]),
    )
    for actual_field, expected_field in zip(actual, expected):
        np.testing.assert_allclose(
            actual_field,
            expected_field,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
