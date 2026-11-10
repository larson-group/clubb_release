"""Check shared standalone setup with independently evolving column fields."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.sponge_layer_damping import (
    initialize_tau_sponge_damp_api, sponge_damp_settings, sponge_damp_xm,
)
from clubb_jax.src.Benchmark_cases import time_dependent_input as forcing
from clubb_jax.src.Input_fields.input_reader import two_dim_read_var, one_dim_read_var
from clubb_jax.src.Input_fields.input_names import z_name, thetal_f_name, time_name


@pytest.mark.parametrize("grid_type", [1, 2, 3])
def test_grid_setup_uses_common_bounds(grid_type):
    """Even and stretched input grids produce the same geometry in every column."""
    gr = setup_grid(3, 100., 0., 400., grid_type=grid_type,
                    momentum_heights=np.tile([0., 100., 200., 300., 400.], (3, 1)),
                    thermodynamic_heights=np.tile([50., 150., 250., 350.], (3, 1)))
    assert gr.zm.shape == (3, 5)
    assert gr.zt.shape == (3, 4)
    np.testing.assert_array_equal(gr.zm, np.tile([0., 100., 200., 300., 400.], (3, 1)))
    np.testing.assert_array_equal(gr.zt, np.tile([50., 150., 250., 350.], (3, 1)))


def test_sponge_uses_column_heights_with_shared_domain():
    """Equal endpoints and level counts allow different sponge levels and time scales."""
    heights = np.array([[0., 100., 600., 800., 1000.], [0., 200., 850., 950., 1000.]])
    gr = setup_grid(2, 0., 0., 1000., grid_type=3, momentum_heights=heights)
    settings = sponge_damp_settings(60., 600., 0.25)
    profile = initialize_tau_sponge_damp_api(2, gr.nzt, gr, 10., gr.zt, settings)
    assert profile.tau_sponge_damp.shape == (2, gr.nzt)
    np.testing.assert_array_equal(profile.sponge_layer_depth, [250., 250.])
    distance = 1000. - np.asarray(gr.zt)
    inside = distance < 250.
    expected_tau = np.where(inside, 60. * 10.**(distance / 250.), 600.)
    np.testing.assert_allclose(profile.tau_sponge_damp, expected_tau, rtol=1e-15)
    xm = jnp.broadcast_to(jnp.array([[300.], [320.]]), (2, gr.nzt))
    ref = jnp.full_like(xm, 280.)
    result = jax.jit(lambda: sponge_damp_xm(2, gr.nzm, gr.nzt, gr, 10., gr.zt, gr.zm, ref, xm, profile))()
    ratio = 10. / expected_tau
    expected = np.where(inside, (np.asarray(xm) + ratio * 280.) / (1. + ratio), np.asarray(xm))
    assert result.shape == xm.shape
    np.testing.assert_allclose(result, expected, rtol=1e-15)
    assert result[0, -2] == xm[0, -2]
    assert result[1, -2] != xm[1, -2]


def test_shared_forcing_file_interpolates_each_columns_grid(monkeypatch):
    """A shared file produces different profiles on grids with common endpoints."""
    z = np.array([[0., 0.], [1000., 1000.]])
    raw = [two_dim_read_var(z_name, time_name, z_name, z.copy())]
    raw += [two_dim_read_var(thetal_f_name, time_name, z_name, z.copy() / 1000.)
            for _ in range(forcing.nforcings - 1)]
    monkeypatch.setattr(forcing, "t_dependent_forcing_data_f_grid", raw)
    monkeypatch.setattr(forcing, "t_dependent_forcing_data", [])
    monkeypatch.setattr(forcing, "dimension_var", one_dim_read_var(time_name, time_name, np.array([0., 10.])))
    grids = np.array([[0., 100., 300., 1000.], [0., 200., 600., 1000.]])
    forcing.initialize_t_dependent_forcings(2, 0, "unused", 4, grids, np.ones((2, 4)), 0)
    values = np.asarray(forcing.t_dependent_forcing_data[1].values)
    assert values.shape == (2, 4, 2)
    np.testing.assert_allclose(values, np.repeat(grids[:, :, None] / 1000., 2, axis=2), rtol=1e-15)


def test_grid_crops_different_bounds_to_equal_level_counts():
    """Each column can retain a different slice of the same-sized input arrays."""
    heights = np.array([[0., 100., 200., 300., 400.], [50., 150., 250., 350., 450.]])
    gr = setup_grid(2, 0., [100., 50.], [300., 250.], grid_type=3, momentum_heights=heights)
    assert gr.zm.shape == (2, 3)
    np.testing.assert_array_equal(gr.zm, [[100., 200., 300.], [50., 150., 250.]])


@pytest.mark.parametrize("grid_type", [1, 2, 3])
def test_grid_rejects_unequal_retained_level_counts(grid_type):
    """Different domain bounds must still produce a rectangular column batch."""
    with pytest.raises(ValueError, match="same number of levels"):
        setup_grid(2, [100., 200.], [0., 0.], [400., 300.], grid_type=grid_type,
                   momentum_heights=np.tile([0., 100., 200., 300., 400.], (2, 1)),
                   thermodynamic_heights=np.array([[50., 150., 250., 350.], [50., 150., 250., 450.]]))


def test_sponge_depth_follows_each_columns_top():
    """Different model tops produce distinct layer depths and damping time scales."""
    heights = np.array([[0., 100., 600., 800., 1000.], [0., 200., 800., 1500., 2000.]])
    gr = setup_grid(2, 0., 0., [1000., 2000.], grid_type=3, momentum_heights=heights)
    profile = initialize_tau_sponge_damp_api(2, gr.nzt, gr, 10., gr.zt, sponge_damp_settings(60., 600., 0.25))
    np.testing.assert_array_equal(profile.sponge_layer_depth, [250., 500.])
    for i in range(2):
        one = setup_grid(1, 0., 0., heights[i,-1], grid_type=3, momentum_heights=heights[i:i+1])
        expected = initialize_tau_sponge_damp_api(1, one.nzt, one, 10., one.zt, sponge_damp_settings(60., 600., 0.25))
        np.testing.assert_array_equal(profile.tau_sponge_damp[i:i+1], expected.tau_sponge_damp)


@pytest.mark.parametrize("grid_type", [1, 2])
def test_setup_allows_different_domain_bounds_with_equal_counts(grid_type):
    """Equal level counts allow different start heights, top heights and spacing."""
    zm = np.array([[0., 100., 200., 300., 400.], [50., 250., 450., 650., 850.]])
    zt = 0.5 * (zm[:, 1:] + zm[:, :-1])
    gr = setup_grid(2, [100., 200.], zm[:, 0], zm[:, -1], grid_type=grid_type,
                    momentum_heights=zm, thermodynamic_heights=zt)
    assert gr.zm.shape == (2, 5)
    assert gr.zt.shape == (2, 4)
    np.testing.assert_array_equal(gr.zm, zm)
    np.testing.assert_array_equal(gr.zt, zt)
