"""Check replayed timesteps against the normal driver with a small numerical kernel.

The real StatsWriter supplies spin-up, sampling, resets and output cadence.
Only the physical timestep is replaced by a two-column recurrence, so these
checks exercise checkpoint scheduling and gradients without compiling CLUBB.
"""
import importlib

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.stats_netcdf import StatsWriter


@pytest.mark.parametrize('split_window', [False, True])
def test_checkpoint_replay_preserves_stats_gradients_and_single_host_output(tmp_path, monkeypatch, split_window):
    """Compare spin-up, split advancement and three output windows with normal advancement.

    Check two independent columns, finite differences, repeated pullbacks, and
    reuse across parameter updates. Host output must occur once per window,
    even when the same saved trajectory is differentiated more than once.
    """
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')
    from clubb_jax.src.CLUBB_core import code_timer_module
    monkeypatch.setattr(code_timer_module, 'l_timer_enabled', False)
    monkeypatch.setattr(driver, 'clubb_at_least_debug_level', lambda level: False)
    registry = tmp_path / 'stats.in'
    registry.write_text('&clubb_stats_nl\n entry(1) = "profile | zt | 1 | State"\n/\n')
    writer = StatsWriter(str(registry), '', 3, 4, 2, np.array([1., 2., 3.]), np.arange(4.),
                         2., 4., 1., 1, 1, 2000, 0., stats_tstart=2., stats_tend=14.)
    writes = []
    original_write = writer._write_timestep

    def write(time, stats):
        """Count actual host output calls, including any accidental replay calls."""
        writes.append(float(time))
        original_write(time, stats)

    monkeypatch.setattr(writer, '_write_timestep', write)

    @jax.jit
    def kernel(x, parameter, time_current, radiation, unchanged):
        """Update fields while returning an unchanged JIT output alongside them."""
        return (.9 + .01 * parameter) * x + parameter * jnp.cos(.1 * time_current) + .02 * radiation, unchanged

    def step(state, itime, time_current, l_sample, l_rad_itime):
        """Advance a stable parameter-dependent recurrence and the real stats bank."""
        state['l_sample'] = l_sample
        p = state['parameter'][:, None]
        state['x'], state['unchanged'] = kernel(state['x'], p, time_current, l_rad_itime, state['unchanged'])
        state['_jax_stats'] = state['_jax_stats'].update('profile', state['x'])

    monkeypatch.setattr(driver, '_advance_timestep', step)
    initial = dict(l_stats=True, l_sample=False, stats_output_path='', stats_writer=writer,
                   ngrdcol=2, dt_main=1., dt_rad=3., time_initial=0., time_final=14., ifinal=14,
                   # This read-only table has string leaves and must stay host-owned.
                   forcing_table={'time': jnp.arange(3.), 'units': 's'})
    cache = {}

    def score(parameter, checkpointed):
        """Return final-state/profile scores while resetting the ordinary host schedule."""
        writer.reset()
        state = dict(initial, parameter=parameter, unchanged=jnp.array([2., 3.]),
                     x=jnp.array([[1., 2., 3.], [2., 3., 4.]]))
        options = {'checkpoint_cache': cache if checkpointed else None}
        if split_window:
            driver.advance_clubb_to_end(state, False, itime_end=7, **options)
            driver.advance_clubb_to_end(state, False, itime_start=8, itime_end=14, **options)
        else:
            driver.advance_clubb_to_end(state, False, **options)
        assert state['l_sample']
        stats = state['_jax_stats']
        averages = stats.buffers[0][0] / jnp.maximum(stats.nsamples[0][0], 1)
        return jnp.sum(averages, axis=1) + jnp.sum(state['x'], axis=1)

    p = jnp.array([.4, .8])
    expected, tangent = jax.jvp(lambda x: score(x, False), (p,), (jnp.ones_like(p),))
    jax.block_until_ready((expected, tangent))
    jax.effects_barrier()
    writes.clear()
    value, pullback = jax.vjp(lambda x: score(x, True), p)
    gradient, = pullback(jnp.ones_like(value))
    second_gradient, = pullback(jnp.array([1., 0.]))
    jax.block_until_ready((value, gradient, second_gradient))
    jax.effects_barrier()
    assert writes == [6., 10., 14.]
    np.testing.assert_allclose(value, expected, rtol=1e-13)
    np.testing.assert_allclose(gradient, tangent, rtol=1e-13)
    np.testing.assert_allclose(second_gradient, [tangent[0], 0.], rtol=1e-13, atol=1e-14)
    h = 1e-5
    finite_difference = (score(p + h, False) - score(p - h, False)) / (2 * h)
    np.testing.assert_allclose(gradient, finite_difference, rtol=1e-8, atol=1e-9)
    count = len(cache)
    repeated = score(p + .01, True)
    jax.block_until_ready(repeated)
    assert len(cache) == count
    np.testing.assert_allclose(repeated, score(p + .01, False), rtol=1e-13)


@pytest.mark.parametrize('fault', ['debug', 'stats_disabled', 'stats_suppressed', 'file_output'])
def test_checkpointing_rejects_host_only_modes_before_scheduling(monkeypatch, fault):
    """Reverse replay cannot include debug checks, disabled stats or file output."""
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')
    monkeypatch.setattr(driver, 'clubb_at_least_debug_level', lambda level: fault == 'debug')
    state = {'l_stats': fault != 'stats_disabled',
             'stats_output_path': 'stats.nc' if fault == 'file_output' else ''}
    with pytest.raises(ValueError, match='enabled in-memory statistics'):
        driver.advance_clubb_to_end(state, False, l_suppress_stats=fault == 'stats_suppressed',
                                    checkpoint_cache={})


def test_delayed_stats_output_cannot_reset_a_later_timestep_schedule(tmp_path, monkeypatch):
    """A traced output's delayed host callback must not overwrite newer cadence flags."""
    registry = tmp_path / 'stats.in'
    registry.write_text('&clubb_stats_nl\n entry(1) = "profile | zt | 1 | State"\n/\n')
    writer = StatsWriter(str(registry), '', 3, 4, 2, np.array([1., 2., 3.]), np.arange(4.),
                         1., 4., 1., 1, 1, 2000, 0., stats_tstart=0., stats_tend=8.)
    from clubb_jax.src.CLUBB_core.stats_netcdf import JaxStats
    stats = JaxStats.from_layout(writer.get_jax_layout(), ncol=2).begin_timestep(
        l_sample=True, l_reset_accumulators=True)
    pending = []

    def enqueue(callback, *args, **kwargs):
        """Delay host reporting until the next output step has been scheduled."""
        pending.append(callback)

    monkeypatch.setattr(jax.debug, 'callback', enqueue)
    writer.begin_timestep(3)

    def finish(value):
        """Trace the real output boundary while deferring its host report."""
        writer.end_timestep(4., stats.update('profile', value))
        return value * value

    jax.jvp(finish, (jnp.ones((2, 3)),), (jnp.ones((2, 3)),))
    assert len(pending) == 1
    writer.begin_timestep(7)
    assert writer.l_sample and writer.l_last_sample
    pending[0](4., stats)
    assert writer.l_sample and writer.l_last_sample
