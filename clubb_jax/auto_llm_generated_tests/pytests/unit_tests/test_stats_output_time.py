"""Regression checks for the inlined CLUBB stats lifecycle."""


from clubb_jax.src import advance_clubb_to_end as driver


def test_stats_timestamp_uses_end_of_current_model_step(monkeypatch):
    calls = []
    updates = []

    class Stats:
        def begin_timestep(self, *, l_sample, reset_accumulators):
            assert l_sample
            assert not reset_accumulators
            return self

        def update(self, name, values):
            updates.append((name, values))
            return self

    stats = Stats()

    class Writer:
        l_reset = False

        def begin_timestep(self, itime_idx):
            assert itime_idx == 0
            return True, True

        def end_timestep(self, stats_time, jax_stats):
            calls.append((stats_time, jax_stats))
            return jax_stats

    class ErrorInfo:
        @staticmethod
        def is_fatal():
            return False

    monkeypatch.setattr(driver, "calculate_thvm", lambda **kwargs: None)
    monkeypatch.setattr(driver, "_prescribe_forcings", lambda state, time_current: None)
    monkeypatch.setattr(driver, "_advance_clubb_core", lambda state: None)
    monkeypatch.setattr(driver, "_advance_radiation", lambda **kwargs: None)

    state = {
        "dt_main": 60.0,
        "dt_rad": 60.0,
        "time_initial": 120.0,
        "time_final": 180.0,
        "ifinal": 1,
        "l_stats": True,
        "stats_writer": Writer(),
        "_jax_stats": stats,
        "thlm_forcing": 0.0,
        "radht": 0.0,
        "l_calc_thlp2_rad": False,
        "err_info": ErrorInfo(),
        "thlm": None,
        "rtm": None,
        "rcm": None,
        "exner": None,
        "thv_ds_zt": None,
        "Nc_in_cloud": 2.0,
        "cloud_frac": 0.25,
        "nzm": 1,
        "nzt": 1,
        "ngrdcol": 1,
    }
    state.update(gr=None, microphys_scheme='none', l_cloud_sed=False)
    for name in ('rcm','rvm','thlm','wprtp','wpthlp','rtp2','thlp2','rtpthlp'):
        state[name+'_mc']=0.0
    for name in ('rtm','wprtp','wpthlp','rtp2','thlp2','rtpthlp'):
        state[name+'_forcing']=0.0
    driver.advance_clubb_to_end(state, l_stdout=False)

    assert calls == [(180.0, stats)]
    assert updates[:2] == [("Ncm", 0.5), ("Nc_in_cloud", 2.0)]
    assert dict(updates)["rtm_mc"] == 0.0
    assert state["Ncm"] == 0.5
    assert state["_jax_stats"] is stats
