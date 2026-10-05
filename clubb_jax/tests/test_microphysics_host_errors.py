"""Fortran debug guards at the host's compiled microphysics boundaries."""

from collections import defaultdict
from types import SimpleNamespace

import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src import advance_clubb_to_end as driver
from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.Microphys import advance_microphys_module


@pytest.mark.parametrize('boundary', ['prep', 'tendencies', 'transport'])
@pytest.mark.parametrize('debug_level', [-1, 0])
def test_fatal_column_respects_source_host_debug_guard(monkeypatch, boundary, debug_level):
    # Inject only kernel returns: the real host routine must propagate a
    # per-column error at debug=-1 and stop at the source boundary at debug=0.
    values = jnp.zeros((2, 3))
    state = defaultdict(lambda: values)
    stats = SimpleNamespace(update=lambda *args: stats)
    state.update(
        gr=SimpleNamespace(ngrdcol=2), microphys_scheme='khairoutdinov_kogan',
        err_info=ErrInfo.initialized(2), _jax_stats=stats,
        flags=SimpleNamespace(tridiag_solve_method=1, fill_holes_type=1, l_upwind_xm_ma=True),
        silhs_config_flags=SimpleNamespace(
            l_lh_importance_sampling=False, l_lh_instant_var_covar_src=False,
        ),
    )
    monkeypatch.setattr(error_code, '_debug_level', debug_level)
    calls = []

    def prep(*args):
        calls.append('prep')
        err = state['err_info']
        if boundary == 'prep':
            err = err.set_fatal(jnp.array([True, False]))
        return stats, err, *((values,) * 23)

    def tendencies(*args):
        calls.append('tendencies')
        error = jnp.array([boundary == 'tendencies', False])
        return stats, *((values,) * 15), error

    def transport(*args):
        calls.append('transport')
        err = state['err_info']
        if boundary == 'transport':
            err = err.set_fatal(jnp.array([True, False]))
        return stats, *((values,) * 8), err, values, values

    monkeypatch.setattr(driver, 'pdf_hydromet_microphys_prep', prep)
    monkeypatch.setattr(driver, 'calc_microphys_scheme_tendcies', tendencies)
    monkeypatch.setattr(driver, 'advance_microphys', transport)
    monkeypatch.setattr(advance_microphys_module, 'write_adv_micro_errors', lambda *args: None)
    if debug_level >= 0:
        message = {
            'prep': 'pdf_hydromet_microphys_prep',
            'tendencies': 'calc_microphys_scheme_tendcies',
            'transport': 'advance_microphys',
        }[boundary]
        with pytest.raises(RuntimeError, match=message):
            driver._advance_microphysics(state, 1, 0., True)
        expected = ['prep', 'tendencies', 'transport']
        assert calls == expected[:expected.index(boundary) + 1]
    else:
        driver._advance_microphysics(state, 1, 0., True)
        assert calls == ['prep', 'tendencies', 'transport']
    np.testing.assert_array_equal(state['err_info'].fatal_mask(), [True, False])
