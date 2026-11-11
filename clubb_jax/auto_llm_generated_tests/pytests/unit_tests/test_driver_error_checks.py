"""The driver honors configured error checks without running a model case."""
import importlib
import jax.numpy as jnp
import pytest

@pytest.mark.parametrize('debug_level', [-1, 0, 1, 2])
def test_driver_error_checks_follow_debug_level(monkeypatch, debug_level):
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    driver = importlib.import_module('clubb_jax.src.advance_clubb_to_end')
    monkeypatch.setattr(error_code, '_debug_level', debug_level)
    state = dict(dt_main=60., dt_rad=120., time_initial=0., ifinal=3, time_final=180.,
                 l_stats=False, l_calc_thlp2_rad=False, ngrdcol=1, nzm=2, nzt=1,
                 err_info=ErrInfo.initialized(1))
    for key in ('thlm', 'rtm', 'rcm', 'exner', 'thv_ds_zt', 'thlm_forcing', 'radht'):
        state[key] = jnp.ones((1, 1))
    # Unrelated tendencies stay inert so this check isolates the radiation error gate.
    for key in ('rtm_forcing', 'wprtp_forcing', 'wpthlp_forcing', 'rtp2_forcing',
                'thlp2_forcing', 'rtpthlp_forcing', 'rcm_mc', 'rvm_mc', 'thlm_mc',
                'wprtp_mc', 'wpthlp_mc', 'rtp2_mc', 'thlp2_mc', 'rtpthlp_mc'):
        state[key] = jnp.zeros((1, 1))
    monkeypatch.setattr(driver, '_advance_microphysics', lambda *args, **kwargs: None)
    monkeypatch.setattr(driver, 'calculate_thvm', lambda **kwargs: kwargs['thlm'])
    monkeypatch.setattr(driver, '_prescribe_forcings', lambda state, time: None)
    monkeypatch.setattr(driver, '_advance_clubb_core', lambda state: None)
    radiation_calls = []
    def fatal_radiation(state, time_current, l_rad_itime):
        radiation_calls.append(l_rad_itime)
        state['err_info'] = state['err_info'].set_fatal()
    monkeypatch.setattr(driver, '_advance_radiation', fatal_radiation)
    if debug_level >= 0:
        with pytest.raises(RuntimeError, match='Fatal error in radiation'):
            driver.advance_clubb_to_end(state, l_stdout=False)
        assert radiation_calls == [True]
    else:
        driver.advance_clubb_to_end(state, l_stdout=False)
        assert state['err_info'].is_fatal()  # Flags remain available to callers.
        assert radiation_calls == [True, True, False]
