"""Check debug error gates at absolute driver iterations.

Adapted from the radiation/error branches of advance_clubb_to_end in
src/clubb_driver.F90. Mocked physics isolates whether debug levels -1 through 2
check radiation errors at the correct cadence for cold and resumed windows.
"""

import jax.numpy as jnp
import pytest

@pytest.mark.parametrize('debug_level', [-1, 0, 1, 2])
@pytest.mark.parametrize('iinit', [1, 2, 3])
def test_driver_error_checks_follow_debug_level(monkeypatch, debug_level, iinit):
    import importlib
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')
    monkeypatch.setattr(error_code, '_debug_level', debug_level)
    state = dict(dt_main=60., dt_rad=120., time_initial=0., iinit=iinit, ifinal=3, time_final=180.,
                 l_stats=False, l_calc_thlp2_rad=False, ngrdcol=1, nzm=2, nzt=1,
                 err_info=ErrInfo.initialized(1))
    for key in ('thlm', 'rtm', 'rcm', 'exner', 'thv_ds_zt', 'thlm_forcing', 'radht'):
        state[key] = jnp.ones((1, 1))
    # This fixture isolates radiation error checks; the normal driver also
    # carries previous microphysics tendencies even when microphysics is off.
    for key in ('rtm_forcing', 'rcm_mc', 'rvm_mc', 'thlm_mc'):
        state[key] = jnp.zeros((1, 1))
    for field in ('wprtp', 'wpthlp', 'rtp2', 'thlp2', 'rtpthlp'):
        state[field + '_forcing'] = jnp.zeros((1, 2))
        state[field + '_mc'] = jnp.zeros((1, 2))
    monkeypatch.setattr(driver, 'calculate_thvm', lambda **kwargs: kwargs['thlm'])
    monkeypatch.setattr(driver, '_prescribe_forcings', lambda state, time: None)
    monkeypatch.setattr(driver, '_advance_clubb_core', lambda state: None)
    monkeypatch.setattr(driver, '_advance_microphysics', lambda state, itime, time, l_rad: None)
    radiation_calls = []
    def fatal_radiation(state, time_current, l_rad_itime):
        radiation_calls.append(l_rad_itime)
        state['err_info'] = state['err_info'].set_fatal()
    monkeypatch.setattr(driver, '_advance_radiation', fatal_radiation)
    if debug_level >= 0:
        with pytest.raises(RuntimeError, match='Fatal error in radiation'):
            driver.advance_clubb_to_end(state, l_stdout=False)
        assert radiation_calls == [(iinit % 2 == 0) or (iinit == 1)]
    else:
        driver.advance_clubb_to_end(state, l_stdout=False)
        assert state['err_info'].is_fatal()  # Flags remain available to callers.
        assert radiation_calls == [
            (itime % 2 == 0) or (itime == 1) for itime in range(iinit, 4)
        ]
