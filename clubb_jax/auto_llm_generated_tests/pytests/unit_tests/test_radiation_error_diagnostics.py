"""Radiation failures print directly while retaining functional fatal status.

Synthetic inputs exercise the lookup failure in radiation_module.F90 and the
invalid-time stop in cos_solar_zen_module.F90. Cached calls check that messages
depend on runtime values, without initializing or advancing an SCM case.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Radiation.parameters_radiation import initialize_radiation_parameters
from clubb_jax.src.Radiation.radiation_module import radiation_driver


@pytest.mark.parametrize("fixed_solar_angle", [False, True])
def test_radiation_errors_print_on_each_failed_runtime_call(fixed_solar_angle, capsys, tmp_path):
    """Valid calls stay quiet; failed cached calls print and retain fatal flags."""
    gr = setup_grid(ngrdcol=2, deltaz=100., zm_init=0., zm_top=500., grid_type=1)
    values = jnp.ones((gr.ngrdcol, gr.nzt))
    params = initialize_radiation_parameters(
        dict(
            l_sw_radiation=True, l_fix_cos_solar_zen=fixed_solar_angle,
            cos_solar_zen_times=[0., 10.], cos_solar_zen_values=[0.2, 0.4],
            day=1, month=1, year=2020,
        ),
        "none", tmp_path,
    )
    info = ErrInfo.initialized(gr.ngrdcol)
    stats = JaxStats.empty(l_sample=False, names=(), ncol=gr.ngrdcol, max_nlev=gr.nzm)
    bad_time = 11. if fixed_solar_angle else jnp.nan
    message = (
        "Time not found in cos_solar_zen_times" if fixed_solar_angle
        else "cos_solar_zen received an invalid current_time."
    )
    for failed in (False, True, True, False):
        result = radiation_driver(
            gr, jnp.asarray(bad_time if failed else 5.), 0., 0,
            gr.ngrdcol, 1, 1, 2020, 0., 0.,
            values, jnp.ones((gr.ngrdcol, gr.nzm)), values, values, values, values,
            values, values, values,
            jnp.zeros((gr.ngrdcol, gr.nzt, 0)), None, stats, info, params,
        )
        jax.block_until_ready(result)
        jax.effects_barrier()
        text = capsys.readouterr().err
        np.testing.assert_array_equal(result[1].fatal_mask(), [failed, failed])
        if failed:
            assert message in text
        else:
            assert text == ""
